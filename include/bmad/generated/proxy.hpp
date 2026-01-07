#pragma once

#include "bmad/fortran_arrays.hpp"
#include "bmad/proxy_base.hpp"

#include <complex>
#include <memory>
#include <string>

extern "C" {
// Forward declarations for Fortran interface
void spline_struct_get_x0(const void* struct_obj, double* value_out);
void spline_struct_set_x0(void* struct_obj, double value_in);
void spline_struct_get_y0(const void* struct_obj, double* value_out);
void spline_struct_set_y0(void* struct_obj, double value_in);
void spline_struct_get_x1(const void* struct_obj, double* value_out);
void spline_struct_set_x1(void* struct_obj, double value_in);
void spline_struct_get_coef_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void spin_polar_struct_get_polarization(
    const void* struct_obj,
    double* value_out);
void spin_polar_struct_set_polarization(void* struct_obj, double value_in);
void spin_polar_struct_get_theta(const void* struct_obj, double* value_out);
void spin_polar_struct_set_theta(void* struct_obj, double value_in);
void spin_polar_struct_get_phi(const void* struct_obj, double* value_out);
void spin_polar_struct_set_phi(void* struct_obj, double value_in);
void spin_polar_struct_get_xi(const void* struct_obj, double* value_out);
void spin_polar_struct_set_xi(void* struct_obj, double value_in);
void ac_kicker_time_struct_get_amp(const void* struct_obj, double* value_out);
void ac_kicker_time_struct_set_amp(void* struct_obj, double value_in);
void ac_kicker_time_struct_get_time(const void* struct_obj, double* value_out);
void ac_kicker_time_struct_set_time(void* struct_obj, double value_in);
void ac_kicker_time_struct_get_spline(const void* struct_obj, void** ptr_out);
void ac_kicker_time_struct_set_spline(void* struct_obj, const void* src_ptr);
void ac_kicker_freq_struct_get_f(const void* struct_obj, double* value_out);
void ac_kicker_freq_struct_set_f(void* struct_obj, double value_in);
void ac_kicker_freq_struct_get_amp(const void* struct_obj, double* value_out);
void ac_kicker_freq_struct_set_amp(void* struct_obj, double value_in);
void ac_kicker_freq_struct_get_phi(const void* struct_obj, double* value_out);
void ac_kicker_freq_struct_set_phi(void* struct_obj, double value_in);

void ac_kicker_struct_get_amp_vs_time_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void ac_kicker_struct_get_frequency_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void interval1_coef_struct_get_c0(const void* struct_obj, double* value_out);
void interval1_coef_struct_set_c0(void* struct_obj, double value_in);
void interval1_coef_struct_get_c1(const void* struct_obj, double* value_out);
void interval1_coef_struct_set_c1(void* struct_obj, double value_in);
void interval1_coef_struct_get_n_exp(const void* struct_obj, double* value_out);
void interval1_coef_struct_set_n_exp(void* struct_obj, double value_in);
void photon_reflect_table_struct_get_angle_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void photon_reflect_table_struct_get_energy_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);

void photon_reflect_table_struct_get_int1_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void photon_reflect_table_struct_get_p_reflect_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void photon_reflect_table_struct_get_max_energy(
    const void* struct_obj,
    double* value_out);
void photon_reflect_table_struct_set_max_energy(
    void* struct_obj,
    double value_in);
void photon_reflect_table_struct_get_p_reflect_scratch_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void photon_reflect_table_struct_get_bragg_angle_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void photon_reflect_surface_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void photon_reflect_surface_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void photon_reflect_surface_struct_get_description_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void photon_reflect_surface_struct_set_description(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void photon_reflect_surface_struct_get_reflectivity_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void photon_reflect_surface_struct_set_reflectivity_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void photon_reflect_surface_struct_get_table_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void photon_reflect_surface_struct_get_surface_roughness_rms(
    const void* struct_obj,
    double* value_out);
void photon_reflect_surface_struct_set_surface_roughness_rms(
    void* struct_obj,
    double value_in);
void photon_reflect_surface_struct_get_roughness_correlation_len(
    const void* struct_obj,
    double* value_out);
void photon_reflect_surface_struct_set_roughness_correlation_len(
    void* struct_obj,
    double value_in);
void photon_reflect_surface_struct_get_ix_surface(
    const void* struct_obj,
    int* value_out);
void photon_reflect_surface_struct_set_ix_surface(
    void* struct_obj,
    int value_in);
void coord_struct_get_vec_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void coord_struct_get_s(const void* struct_obj, double* value_out);
void coord_struct_set_s(void* struct_obj, double value_in);
void coord_struct_get_t(const void* struct_obj, long double* value_out);
void coord_struct_set_t(void* struct_obj, long double value_in);
void coord_struct_get_spin_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void coord_struct_get_field_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void coord_struct_get_phase_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void coord_struct_get_charge(const void* struct_obj, double* value_out);
void coord_struct_set_charge(void* struct_obj, double value_in);
void coord_struct_get_dt_ref(const void* struct_obj, double* value_out);
void coord_struct_set_dt_ref(void* struct_obj, double value_in);
void coord_struct_get_r(const void* struct_obj, double* value_out);
void coord_struct_set_r(void* struct_obj, double value_in);
void coord_struct_get_p0c(const void* struct_obj, double* value_out);
void coord_struct_set_p0c(void* struct_obj, double value_in);
void coord_struct_get_E_potential(const void* struct_obj, double* value_out);
void coord_struct_set_E_potential(void* struct_obj, double value_in);
void coord_struct_get_beta(const void* struct_obj, double* value_out);
void coord_struct_set_beta(void* struct_obj, double value_in);
void coord_struct_get_ix_ele(const void* struct_obj, int* value_out);
void coord_struct_set_ix_ele(void* struct_obj, int value_in);
void coord_struct_get_ix_branch(const void* struct_obj, int* value_out);
void coord_struct_set_ix_branch(void* struct_obj, int value_in);
void coord_struct_get_ix_turn(const void* struct_obj, int* value_out);
void coord_struct_set_ix_turn(void* struct_obj, int value_in);
void coord_struct_get_ix_user(const void* struct_obj, int* value_out);
void coord_struct_set_ix_user(void* struct_obj, int value_in);
void coord_struct_get_state(const void* struct_obj, int* value_out);
void coord_struct_set_state(void* struct_obj, int value_in);
void coord_struct_get_direction(const void* struct_obj, int* value_out);
void coord_struct_set_direction(void* struct_obj, int value_in);
void coord_struct_get_time_dir(const void* struct_obj, int* value_out);
void coord_struct_set_time_dir(void* struct_obj, int value_in);
void coord_struct_get_species(const void* struct_obj, int* value_out);
void coord_struct_set_species(void* struct_obj, int value_in);
void coord_struct_get_location(const void* struct_obj, int* value_out);
void coord_struct_set_location(void* struct_obj, int value_in);

void coord_array_struct_get_orbit_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bpm_phase_coupling_struct_get_K_22a(
    const void* struct_obj,
    double* value_out);
void bpm_phase_coupling_struct_set_K_22a(void* struct_obj, double value_in);
void bpm_phase_coupling_struct_get_K_12a(
    const void* struct_obj,
    double* value_out);
void bpm_phase_coupling_struct_set_K_12a(void* struct_obj, double value_in);
void bpm_phase_coupling_struct_get_K_11b(
    const void* struct_obj,
    double* value_out);
void bpm_phase_coupling_struct_set_K_11b(void* struct_obj, double value_in);
void bpm_phase_coupling_struct_get_K_12b(
    const void* struct_obj,
    double* value_out);
void bpm_phase_coupling_struct_set_K_12b(void* struct_obj, double value_in);
void bpm_phase_coupling_struct_get_Cbar22_a(
    const void* struct_obj,
    double* value_out);
void bpm_phase_coupling_struct_set_Cbar22_a(void* struct_obj, double value_in);
void bpm_phase_coupling_struct_get_Cbar12_a(
    const void* struct_obj,
    double* value_out);
void bpm_phase_coupling_struct_set_Cbar12_a(void* struct_obj, double value_in);
void bpm_phase_coupling_struct_get_Cbar11_b(
    const void* struct_obj,
    double* value_out);
void bpm_phase_coupling_struct_set_Cbar11_b(void* struct_obj, double value_in);
void bpm_phase_coupling_struct_get_Cbar12_b(
    const void* struct_obj,
    double* value_out);
void bpm_phase_coupling_struct_set_Cbar12_b(void* struct_obj, double value_in);
void bpm_phase_coupling_struct_get_phi_a(
    const void* struct_obj,
    double* value_out);
void bpm_phase_coupling_struct_set_phi_a(void* struct_obj, double value_in);
void bpm_phase_coupling_struct_get_phi_b(
    const void* struct_obj,
    double* value_out);
void bpm_phase_coupling_struct_set_phi_b(void* struct_obj, double value_in);
void expression_atom_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void expression_atom_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void expression_atom_struct_get_type(const void* struct_obj, int* value_out);
void expression_atom_struct_set_type(void* struct_obj, int value_in);
void expression_atom_struct_get_value(
    const void* struct_obj,
    double* value_out);
void expression_atom_struct_set_value(void* struct_obj, double value_in);
void wake_sr_z_long_struct_get_w_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void wake_sr_z_long_struct_get_fw_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    bool* is_alloc);
void wake_sr_z_long_struct_get_fbunch_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    bool* is_alloc);
void wake_sr_z_long_struct_get_w_out_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    bool* is_alloc);
void wake_sr_z_long_struct_get_dz(const void* struct_obj, double* value_out);
void wake_sr_z_long_struct_set_dz(void* struct_obj, double value_in);
void wake_sr_z_long_struct_get_z0(const void* struct_obj, double* value_out);
void wake_sr_z_long_struct_set_z0(void* struct_obj, double value_in);
void wake_sr_z_long_struct_get_smoothing_sigma(
    const void* struct_obj,
    double* value_out);
void wake_sr_z_long_struct_set_smoothing_sigma(
    void* struct_obj,
    double value_in);
void wake_sr_z_long_struct_get_position_dependence(
    const void* struct_obj,
    int* value_out);
void wake_sr_z_long_struct_set_position_dependence(
    void* struct_obj,
    int value_in);
void wake_sr_z_long_struct_get_time_based(
    const void* struct_obj,
    bool* value_out);
void wake_sr_z_long_struct_set_time_based(void* struct_obj, bool value_in);
void wake_sr_mode_struct_get_amp(const void* struct_obj, double* value_out);
void wake_sr_mode_struct_set_amp(void* struct_obj, double value_in);
void wake_sr_mode_struct_get_damp(const void* struct_obj, double* value_out);
void wake_sr_mode_struct_set_damp(void* struct_obj, double value_in);
void wake_sr_mode_struct_get_k(const void* struct_obj, double* value_out);
void wake_sr_mode_struct_set_k(void* struct_obj, double value_in);
void wake_sr_mode_struct_get_phi(const void* struct_obj, double* value_out);
void wake_sr_mode_struct_set_phi(void* struct_obj, double value_in);
void wake_sr_mode_struct_get_b_sin(const void* struct_obj, double* value_out);
void wake_sr_mode_struct_set_b_sin(void* struct_obj, double value_in);
void wake_sr_mode_struct_get_b_cos(const void* struct_obj, double* value_out);
void wake_sr_mode_struct_set_b_cos(void* struct_obj, double value_in);
void wake_sr_mode_struct_get_a_sin(const void* struct_obj, double* value_out);
void wake_sr_mode_struct_set_a_sin(void* struct_obj, double value_in);
void wake_sr_mode_struct_get_a_cos(const void* struct_obj, double* value_out);
void wake_sr_mode_struct_set_a_cos(void* struct_obj, double value_in);
void wake_sr_mode_struct_get_polarization(
    const void* struct_obj,
    int* value_out);
void wake_sr_mode_struct_set_polarization(void* struct_obj, int value_in);
void wake_sr_mode_struct_get_position_dependence(
    const void* struct_obj,
    int* value_out);
void wake_sr_mode_struct_set_position_dependence(
    void* struct_obj,
    int value_in);
void wake_sr_struct_get_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void wake_sr_struct_set_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void wake_sr_struct_get_z_long(const void* struct_obj, void** ptr_out);
void wake_sr_struct_set_z_long(void* struct_obj, const void* src_ptr);

void wake_sr_struct_get_long_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void wake_sr_struct_get_trans_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void wake_sr_struct_get_z_ref_long(const void* struct_obj, double* value_out);
void wake_sr_struct_set_z_ref_long(void* struct_obj, double value_in);
void wake_sr_struct_get_z_ref_trans(const void* struct_obj, double* value_out);
void wake_sr_struct_set_z_ref_trans(void* struct_obj, double value_in);
void wake_sr_struct_get_z_max(const void* struct_obj, double* value_out);
void wake_sr_struct_set_z_max(void* struct_obj, double value_in);
void wake_sr_struct_get_amp_scale(const void* struct_obj, double* value_out);
void wake_sr_struct_set_amp_scale(void* struct_obj, double value_in);
void wake_sr_struct_get_z_scale(const void* struct_obj, double* value_out);
void wake_sr_struct_set_z_scale(void* struct_obj, double value_in);
void wake_sr_struct_get_scale_with_length(
    const void* struct_obj,
    bool* value_out);
void wake_sr_struct_set_scale_with_length(void* struct_obj, bool value_in);
void wake_lr_mode_struct_get_freq(const void* struct_obj, double* value_out);
void wake_lr_mode_struct_set_freq(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_freq_in(const void* struct_obj, double* value_out);
void wake_lr_mode_struct_set_freq_in(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_R_over_Q(
    const void* struct_obj,
    double* value_out);
void wake_lr_mode_struct_set_R_over_Q(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_Q(const void* struct_obj, double* value_out);
void wake_lr_mode_struct_set_Q(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_damp(const void* struct_obj, double* value_out);
void wake_lr_mode_struct_set_damp(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_phi(const void* struct_obj, double* value_out);
void wake_lr_mode_struct_set_phi(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_angle(const void* struct_obj, double* value_out);
void wake_lr_mode_struct_set_angle(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_b_sin(const void* struct_obj, double* value_out);
void wake_lr_mode_struct_set_b_sin(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_b_cos(const void* struct_obj, double* value_out);
void wake_lr_mode_struct_set_b_cos(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_a_sin(const void* struct_obj, double* value_out);
void wake_lr_mode_struct_set_a_sin(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_a_cos(const void* struct_obj, double* value_out);
void wake_lr_mode_struct_set_a_cos(void* struct_obj, double value_in);
void wake_lr_mode_struct_get_m(const void* struct_obj, int* value_out);
void wake_lr_mode_struct_set_m(void* struct_obj, int value_in);
void wake_lr_mode_struct_get_polarized(const void* struct_obj, bool* value_out);
void wake_lr_mode_struct_set_polarized(void* struct_obj, bool value_in);
void wake_lr_struct_get_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void wake_lr_struct_set_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void wake_lr_struct_get_mode_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void wake_lr_struct_get_t_ref(const void* struct_obj, double* value_out);
void wake_lr_struct_set_t_ref(void* struct_obj, double value_in);
void wake_lr_struct_get_freq_spread(const void* struct_obj, double* value_out);
void wake_lr_struct_set_freq_spread(void* struct_obj, double value_in);
void wake_lr_struct_get_amp_scale(const void* struct_obj, double* value_out);
void wake_lr_struct_set_amp_scale(void* struct_obj, double value_in);
void wake_lr_struct_get_time_scale(const void* struct_obj, double* value_out);
void wake_lr_struct_set_time_scale(void* struct_obj, double value_in);
void wake_lr_struct_get_self_wake_on(const void* struct_obj, bool* value_out);
void wake_lr_struct_set_self_wake_on(void* struct_obj, bool value_in);
void lat_ele_loc_struct_get_ix_ele(const void* struct_obj, int* value_out);
void lat_ele_loc_struct_set_ix_ele(void* struct_obj, int value_in);
void lat_ele_loc_struct_get_ix_branch(const void* struct_obj, int* value_out);
void lat_ele_loc_struct_set_ix_branch(void* struct_obj, int value_in);
void wake_struct_get_sr(const void* struct_obj, void** ptr_out);
void wake_struct_set_sr(void* struct_obj, const void* src_ptr);
void wake_struct_get_lr(const void* struct_obj, void** ptr_out);
void wake_struct_set_lr(void* struct_obj, const void* src_ptr);
void taylor_term_struct_get_coef(const void* struct_obj, double* value_out);
void taylor_term_struct_set_coef(void* struct_obj, double value_in);
void taylor_term_struct_get_expn_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void taylor_struct_get_ref(const void* struct_obj, double* value_out);
void taylor_struct_set_ref(void* struct_obj, double value_in);

void taylor_struct_get_term_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void em_taylor_term_struct_get_coef(const void* struct_obj, double* value_out);
void em_taylor_term_struct_set_coef(void* struct_obj, double value_in);
void em_taylor_term_struct_get_expn_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void em_taylor_struct_get_ref(const void* struct_obj, double* value_out);
void em_taylor_struct_set_ref(void* struct_obj, double value_in);

void em_taylor_struct_get_term_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void cartesian_map_term1_struct_get_coef(
    const void* struct_obj,
    double* value_out);
void cartesian_map_term1_struct_set_coef(void* struct_obj, double value_in);
void cartesian_map_term1_struct_get_kx(
    const void* struct_obj,
    double* value_out);
void cartesian_map_term1_struct_set_kx(void* struct_obj, double value_in);
void cartesian_map_term1_struct_get_ky(
    const void* struct_obj,
    double* value_out);
void cartesian_map_term1_struct_set_ky(void* struct_obj, double value_in);
void cartesian_map_term1_struct_get_kz(
    const void* struct_obj,
    double* value_out);
void cartesian_map_term1_struct_set_kz(void* struct_obj, double value_in);
void cartesian_map_term1_struct_get_x0(
    const void* struct_obj,
    double* value_out);
void cartesian_map_term1_struct_set_x0(void* struct_obj, double value_in);
void cartesian_map_term1_struct_get_y0(
    const void* struct_obj,
    double* value_out);
void cartesian_map_term1_struct_set_y0(void* struct_obj, double value_in);
void cartesian_map_term1_struct_get_phi_z(
    const void* struct_obj,
    double* value_out);
void cartesian_map_term1_struct_set_phi_z(void* struct_obj, double value_in);
void cartesian_map_term1_struct_get_family(
    const void* struct_obj,
    int* value_out);
void cartesian_map_term1_struct_set_family(void* struct_obj, int value_in);
void cartesian_map_term1_struct_get_form(
    const void* struct_obj,
    int* value_out);
void cartesian_map_term1_struct_set_form(void* struct_obj, int value_in);
void cartesian_map_term_struct_get_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void cartesian_map_term_struct_set_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void cartesian_map_term_struct_get_n_link(
    const void* struct_obj,
    int* value_out);
void cartesian_map_term_struct_set_n_link(void* struct_obj, int value_in);

void cartesian_map_term_struct_get_term_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void cartesian_map_struct_get_field_scale(
    const void* struct_obj,
    double* value_out);
void cartesian_map_struct_set_field_scale(void* struct_obj, double value_in);
void cartesian_map_struct_get_r0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void cartesian_map_struct_get_master_parameter(
    const void* struct_obj,
    int* value_out);
void cartesian_map_struct_set_master_parameter(void* struct_obj, int value_in);
void cartesian_map_struct_get_ele_anchor_pt(
    const void* struct_obj,
    int* value_out);
void cartesian_map_struct_set_ele_anchor_pt(void* struct_obj, int value_in);
void cartesian_map_struct_get_field_type(
    const void* struct_obj,
    int* value_out);
void cartesian_map_struct_set_field_type(void* struct_obj, int value_in);
void cartesian_map_struct_get_ptr(const void* struct_obj, void** ptr_out);
void cartesian_map_struct_set_ptr(void* struct_obj, const void* src_ptr);
void cylindrical_map_term1_struct_get_e_coef(
    const void* struct_obj,
    std::complex<double>* value_out);
void cylindrical_map_term1_struct_set_e_coef(
    void* struct_obj,
    std::complex<double> value_in);
void cylindrical_map_term1_struct_get_b_coef(
    const void* struct_obj,
    std::complex<double>* value_out);
void cylindrical_map_term1_struct_set_b_coef(
    void* struct_obj,
    std::complex<double> value_in);
void cylindrical_map_term_struct_get_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void cylindrical_map_term_struct_set_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void cylindrical_map_term_struct_get_n_link(
    const void* struct_obj,
    int* value_out);
void cylindrical_map_term_struct_set_n_link(void* struct_obj, int value_in);

void cylindrical_map_term_struct_get_term_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void cylindrical_map_struct_get_m(const void* struct_obj, int* value_out);
void cylindrical_map_struct_set_m(void* struct_obj, int value_in);
void cylindrical_map_struct_get_harmonic(
    const void* struct_obj,
    int* value_out);
void cylindrical_map_struct_set_harmonic(void* struct_obj, int value_in);
void cylindrical_map_struct_get_phi0_fieldmap(
    const void* struct_obj,
    double* value_out);
void cylindrical_map_struct_set_phi0_fieldmap(
    void* struct_obj,
    double value_in);
void cylindrical_map_struct_get_theta0_azimuth(
    const void* struct_obj,
    double* value_out);
void cylindrical_map_struct_set_theta0_azimuth(
    void* struct_obj,
    double value_in);
void cylindrical_map_struct_get_field_scale(
    const void* struct_obj,
    double* value_out);
void cylindrical_map_struct_set_field_scale(void* struct_obj, double value_in);
void cylindrical_map_struct_get_master_parameter(
    const void* struct_obj,
    int* value_out);
void cylindrical_map_struct_set_master_parameter(
    void* struct_obj,
    int value_in);
void cylindrical_map_struct_get_ele_anchor_pt(
    const void* struct_obj,
    int* value_out);
void cylindrical_map_struct_set_ele_anchor_pt(void* struct_obj, int value_in);
void cylindrical_map_struct_get_dz(const void* struct_obj, double* value_out);
void cylindrical_map_struct_set_dz(void* struct_obj, double value_in);
void cylindrical_map_struct_get_r0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void cylindrical_map_struct_get_ptr(const void* struct_obj, void** ptr_out);
void cylindrical_map_struct_set_ptr(void* struct_obj, const void* src_ptr);
void bicubic_cmplx_coef_struct_get_coef_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void bicubic_cmplx_coef_struct_get_i_box_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void tricubic_cmplx_coef_struct_get_coef_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void tricubic_cmplx_coef_struct_get_i_box_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void grid_field_pt1_struct_get_E_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    bool* is_alloc);
void grid_field_pt1_struct_get_B_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    bool* is_alloc);
void grid_field_pt_struct_get_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void grid_field_pt_struct_set_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void grid_field_pt_struct_get_n_link(const void* struct_obj, int* value_out);
void grid_field_pt_struct_set_n_link(void* struct_obj, int value_in);

void grid_field_pt_struct_get_pt_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void grid_field_struct_get_geometry(const void* struct_obj, int* value_out);
void grid_field_struct_set_geometry(void* struct_obj, int value_in);
void grid_field_struct_get_harmonic(const void* struct_obj, int* value_out);
void grid_field_struct_set_harmonic(void* struct_obj, int value_in);
void grid_field_struct_get_phi0_fieldmap(
    const void* struct_obj,
    double* value_out);
void grid_field_struct_set_phi0_fieldmap(void* struct_obj, double value_in);
void grid_field_struct_get_field_scale(
    const void* struct_obj,
    double* value_out);
void grid_field_struct_set_field_scale(void* struct_obj, double value_in);
void grid_field_struct_get_field_type(const void* struct_obj, int* value_out);
void grid_field_struct_set_field_type(void* struct_obj, int value_in);
void grid_field_struct_get_master_parameter(
    const void* struct_obj,
    int* value_out);
void grid_field_struct_set_master_parameter(void* struct_obj, int value_in);
void grid_field_struct_get_ele_anchor_pt(
    const void* struct_obj,
    int* value_out);
void grid_field_struct_set_ele_anchor_pt(void* struct_obj, int value_in);
void grid_field_struct_get_interpolation_order(
    const void* struct_obj,
    int* value_out);
void grid_field_struct_set_interpolation_order(void* struct_obj, int value_in);
void grid_field_struct_get_dr_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void grid_field_struct_get_r0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void grid_field_struct_get_curved_ref_frame(
    const void* struct_obj,
    bool* value_out);
void grid_field_struct_set_curved_ref_frame(void* struct_obj, bool value_in);
void grid_field_struct_get_ptr(const void* struct_obj, void** ptr_out);
void grid_field_struct_set_ptr(void* struct_obj, const void* src_ptr);

void grid_field_struct_get_bi_coef_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void grid_field_struct_get_tri_coef_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void floor_position_struct_get_r_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void floor_position_struct_get_w_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void floor_position_struct_get_theta(const void* struct_obj, double* value_out);
void floor_position_struct_set_theta(void* struct_obj, double value_in);
void floor_position_struct_get_phi(const void* struct_obj, double* value_out);
void floor_position_struct_set_phi(void* struct_obj, double value_in);
void floor_position_struct_get_psi(const void* struct_obj, double* value_out);
void floor_position_struct_set_psi(void* struct_obj, double value_in);
void high_energy_space_charge_struct_get_closed_orb(
    const void* struct_obj,
    void** ptr_out);
void high_energy_space_charge_struct_set_closed_orb(
    void* struct_obj,
    const void* src_ptr);
void high_energy_space_charge_struct_get_kick_const(
    const void* struct_obj,
    double* value_out);
void high_energy_space_charge_struct_set_kick_const(
    void* struct_obj,
    double value_in);
void high_energy_space_charge_struct_get_sig_x(
    const void* struct_obj,
    double* value_out);
void high_energy_space_charge_struct_set_sig_x(
    void* struct_obj,
    double value_in);
void high_energy_space_charge_struct_get_sig_y(
    const void* struct_obj,
    double* value_out);
void high_energy_space_charge_struct_set_sig_y(
    void* struct_obj,
    double value_in);
void high_energy_space_charge_struct_get_phi(
    const void* struct_obj,
    double* value_out);
void high_energy_space_charge_struct_set_phi(void* struct_obj, double value_in);
void high_energy_space_charge_struct_get_sin_phi(
    const void* struct_obj,
    double* value_out);
void high_energy_space_charge_struct_set_sin_phi(
    void* struct_obj,
    double value_in);
void high_energy_space_charge_struct_get_cos_phi(
    const void* struct_obj,
    double* value_out);
void high_energy_space_charge_struct_set_cos_phi(
    void* struct_obj,
    double value_in);
void high_energy_space_charge_struct_get_sig_z(
    const void* struct_obj,
    double* value_out);
void high_energy_space_charge_struct_set_sig_z(
    void* struct_obj,
    double value_in);
void xy_disp_struct_get_eta(const void* struct_obj, double* value_out);
void xy_disp_struct_set_eta(void* struct_obj, double value_in);
void xy_disp_struct_get_etap(const void* struct_obj, double* value_out);
void xy_disp_struct_set_etap(void* struct_obj, double value_in);
void xy_disp_struct_get_deta_ds(const void* struct_obj, double* value_out);
void xy_disp_struct_set_deta_ds(void* struct_obj, double value_in);
void xy_disp_struct_get_sigma(const void* struct_obj, double* value_out);
void xy_disp_struct_set_sigma(void* struct_obj, double value_in);
void xy_disp_struct_get_deta_dpz(const void* struct_obj, double* value_out);
void xy_disp_struct_set_deta_dpz(void* struct_obj, double value_in);
void xy_disp_struct_get_detap_dpz(const void* struct_obj, double* value_out);
void xy_disp_struct_set_detap_dpz(void* struct_obj, double value_in);
void twiss_struct_get_beta(const void* struct_obj, double* value_out);
void twiss_struct_set_beta(void* struct_obj, double value_in);
void twiss_struct_get_alpha(const void* struct_obj, double* value_out);
void twiss_struct_set_alpha(void* struct_obj, double value_in);
void twiss_struct_get_gamma(const void* struct_obj, double* value_out);
void twiss_struct_set_gamma(void* struct_obj, double value_in);
void twiss_struct_get_phi(const void* struct_obj, double* value_out);
void twiss_struct_set_phi(void* struct_obj, double value_in);
void twiss_struct_get_eta(const void* struct_obj, double* value_out);
void twiss_struct_set_eta(void* struct_obj, double value_in);
void twiss_struct_get_etap(const void* struct_obj, double* value_out);
void twiss_struct_set_etap(void* struct_obj, double value_in);
void twiss_struct_get_deta_ds(const void* struct_obj, double* value_out);
void twiss_struct_set_deta_ds(void* struct_obj, double value_in);
void twiss_struct_get_sigma(const void* struct_obj, double* value_out);
void twiss_struct_set_sigma(void* struct_obj, double value_in);
void twiss_struct_get_sigma_p(const void* struct_obj, double* value_out);
void twiss_struct_set_sigma_p(void* struct_obj, double value_in);
void twiss_struct_get_emit(const void* struct_obj, double* value_out);
void twiss_struct_set_emit(void* struct_obj, double value_in);
void twiss_struct_get_norm_emit(const void* struct_obj, double* value_out);
void twiss_struct_set_norm_emit(void* struct_obj, double value_in);
void twiss_struct_get_chrom(const void* struct_obj, double* value_out);
void twiss_struct_set_chrom(void* struct_obj, double value_in);
void twiss_struct_get_dbeta_dpz(const void* struct_obj, double* value_out);
void twiss_struct_set_dbeta_dpz(void* struct_obj, double value_in);
void twiss_struct_get_dalpha_dpz(const void* struct_obj, double* value_out);
void twiss_struct_set_dalpha_dpz(void* struct_obj, double value_in);
void twiss_struct_get_deta_dpz(const void* struct_obj, double* value_out);
void twiss_struct_set_deta_dpz(void* struct_obj, double value_in);
void twiss_struct_get_detap_dpz(const void* struct_obj, double* value_out);
void twiss_struct_set_detap_dpz(void* struct_obj, double value_in);
void mode3_struct_get_v_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void mode3_struct_get_a(const void* struct_obj, void** ptr_out);
void mode3_struct_set_a(void* struct_obj, const void* src_ptr);
void mode3_struct_get_b(const void* struct_obj, void** ptr_out);
void mode3_struct_set_b(void* struct_obj, const void* src_ptr);
void mode3_struct_get_c(const void* struct_obj, void** ptr_out);
void mode3_struct_set_c(void* struct_obj, const void* src_ptr);
void mode3_struct_get_x(const void* struct_obj, void** ptr_out);
void mode3_struct_set_x(void* struct_obj, const void* src_ptr);
void mode3_struct_get_y(const void* struct_obj, void** ptr_out);
void mode3_struct_set_y(void* struct_obj, const void* src_ptr);
void bookkeeping_state_struct_get_attributes(
    const void* struct_obj,
    int* value_out);
void bookkeeping_state_struct_set_attributes(void* struct_obj, int value_in);
void bookkeeping_state_struct_get_control(
    const void* struct_obj,
    int* value_out);
void bookkeeping_state_struct_set_control(void* struct_obj, int value_in);
void bookkeeping_state_struct_get_floor_position(
    const void* struct_obj,
    int* value_out);
void bookkeeping_state_struct_set_floor_position(
    void* struct_obj,
    int value_in);
void bookkeeping_state_struct_get_s_position(
    const void* struct_obj,
    int* value_out);
void bookkeeping_state_struct_set_s_position(void* struct_obj, int value_in);
void bookkeeping_state_struct_get_ref_energy(
    const void* struct_obj,
    int* value_out);
void bookkeeping_state_struct_set_ref_energy(void* struct_obj, int value_in);
void bookkeeping_state_struct_get_mat6(const void* struct_obj, int* value_out);
void bookkeeping_state_struct_set_mat6(void* struct_obj, int value_in);
void bookkeeping_state_struct_get_rad_int(
    const void* struct_obj,
    int* value_out);
void bookkeeping_state_struct_set_rad_int(void* struct_obj, int value_in);
void bookkeeping_state_struct_get_ptc(const void* struct_obj, int* value_out);
void bookkeeping_state_struct_set_ptc(void* struct_obj, int value_in);
void bookkeeping_state_struct_get_has_misalign(
    const void* struct_obj,
    bool* value_out);
void bookkeeping_state_struct_set_has_misalign(void* struct_obj, bool value_in);
void rad_map_struct_get_ref_orb_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void rad_map_struct_get_damp_dmat_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void rad_map_struct_get_xfer_damp_vec_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void rad_map_struct_get_xfer_damp_mat_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void rad_map_struct_get_stoc_mat_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void rad_map_ele_struct_get_rm0(const void* struct_obj, void** ptr_out);
void rad_map_ele_struct_set_rm0(void* struct_obj, const void* src_ptr);
void rad_map_ele_struct_get_rm1(const void* struct_obj, void** ptr_out);
void rad_map_ele_struct_set_rm1(void* struct_obj, const void* src_ptr);
void rad_map_ele_struct_get_stale(const void* struct_obj, bool* value_out);
void rad_map_ele_struct_set_stale(void* struct_obj, bool value_in);
void gen_grad1_struct_get_m(const void* struct_obj, int* value_out);
void gen_grad1_struct_set_m(void* struct_obj, int value_in);
void gen_grad1_struct_get_sincos(const void* struct_obj, int* value_out);
void gen_grad1_struct_set_sincos(void* struct_obj, int value_in);
void gen_grad1_struct_get_n_deriv_max(const void* struct_obj, int* value_out);
void gen_grad1_struct_set_n_deriv_max(void* struct_obj, int value_in);
void gen_grad1_struct_get_deriv_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void gen_grad_map_struct_get_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void gen_grad_map_struct_set_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void gen_grad_map_struct_get_gg_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void gen_grad_map_struct_get_ele_anchor_pt(
    const void* struct_obj,
    int* value_out);
void gen_grad_map_struct_set_ele_anchor_pt(void* struct_obj, int value_in);
void gen_grad_map_struct_get_field_type(const void* struct_obj, int* value_out);
void gen_grad_map_struct_set_field_type(void* struct_obj, int value_in);
void gen_grad_map_struct_get_iz0(const void* struct_obj, int* value_out);
void gen_grad_map_struct_set_iz0(void* struct_obj, int value_in);
void gen_grad_map_struct_get_iz1(const void* struct_obj, int* value_out);
void gen_grad_map_struct_set_iz1(void* struct_obj, int value_in);
void gen_grad_map_struct_get_dz(const void* struct_obj, double* value_out);
void gen_grad_map_struct_set_dz(void* struct_obj, double value_in);
void gen_grad_map_struct_get_r0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void gen_grad_map_struct_get_field_scale(
    const void* struct_obj,
    double* value_out);
void gen_grad_map_struct_set_field_scale(void* struct_obj, double value_in);
void gen_grad_map_struct_get_master_parameter(
    const void* struct_obj,
    int* value_out);
void gen_grad_map_struct_set_master_parameter(void* struct_obj, int value_in);
void gen_grad_map_struct_get_curved_ref_frame(
    const void* struct_obj,
    bool* value_out);
void gen_grad_map_struct_set_curved_ref_frame(void* struct_obj, bool value_in);
void surface_segmented_pt_struct_get_x0(
    const void* struct_obj,
    double* value_out);
void surface_segmented_pt_struct_set_x0(void* struct_obj, double value_in);
void surface_segmented_pt_struct_get_y0(
    const void* struct_obj,
    double* value_out);
void surface_segmented_pt_struct_set_y0(void* struct_obj, double value_in);
void surface_segmented_pt_struct_get_z0(
    const void* struct_obj,
    double* value_out);
void surface_segmented_pt_struct_set_z0(void* struct_obj, double value_in);
void surface_segmented_pt_struct_get_dz_dx(
    const void* struct_obj,
    double* value_out);
void surface_segmented_pt_struct_set_dz_dx(void* struct_obj, double value_in);
void surface_segmented_pt_struct_get_dz_dy(
    const void* struct_obj,
    double* value_out);
void surface_segmented_pt_struct_set_dz_dy(void* struct_obj, double value_in);
void surface_segmented_struct_get_active(
    const void* struct_obj,
    bool* value_out);
void surface_segmented_struct_set_active(void* struct_obj, bool value_in);
void surface_segmented_struct_get_dr_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void surface_segmented_struct_get_r0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);

void surface_segmented_struct_get_pt_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void surface_h_misalign_pt_struct_get_x0(
    const void* struct_obj,
    double* value_out);
void surface_h_misalign_pt_struct_set_x0(void* struct_obj, double value_in);
void surface_h_misalign_pt_struct_get_y0(
    const void* struct_obj,
    double* value_out);
void surface_h_misalign_pt_struct_set_y0(void* struct_obj, double value_in);
void surface_h_misalign_pt_struct_get_rot_y(
    const void* struct_obj,
    double* value_out);
void surface_h_misalign_pt_struct_set_rot_y(void* struct_obj, double value_in);
void surface_h_misalign_pt_struct_get_rot_t(
    const void* struct_obj,
    double* value_out);
void surface_h_misalign_pt_struct_set_rot_t(void* struct_obj, double value_in);
void surface_h_misalign_pt_struct_get_rot_y_rms(
    const void* struct_obj,
    double* value_out);
void surface_h_misalign_pt_struct_set_rot_y_rms(
    void* struct_obj,
    double value_in);
void surface_h_misalign_pt_struct_get_rot_t_rms(
    const void* struct_obj,
    double* value_out);
void surface_h_misalign_pt_struct_set_rot_t_rms(
    void* struct_obj,
    double value_in);
void surface_h_misalign_struct_get_active(
    const void* struct_obj,
    bool* value_out);
void surface_h_misalign_struct_set_active(void* struct_obj, bool value_in);
void surface_h_misalign_struct_get_dr_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void surface_h_misalign_struct_get_r0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);

void surface_h_misalign_struct_get_pt_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void surface_displacement_pt_struct_get_x0(
    const void* struct_obj,
    double* value_out);
void surface_displacement_pt_struct_set_x0(void* struct_obj, double value_in);
void surface_displacement_pt_struct_get_y0(
    const void* struct_obj,
    double* value_out);
void surface_displacement_pt_struct_set_y0(void* struct_obj, double value_in);
void surface_displacement_pt_struct_get_z0(
    const void* struct_obj,
    double* value_out);
void surface_displacement_pt_struct_set_z0(void* struct_obj, double value_in);
void surface_displacement_pt_struct_get_dz_dx(
    const void* struct_obj,
    double* value_out);
void surface_displacement_pt_struct_set_dz_dx(
    void* struct_obj,
    double value_in);
void surface_displacement_pt_struct_get_dz_dy(
    const void* struct_obj,
    double* value_out);
void surface_displacement_pt_struct_set_dz_dy(
    void* struct_obj,
    double value_in);
void surface_displacement_pt_struct_get_d2z_dxdy(
    const void* struct_obj,
    double* value_out);
void surface_displacement_pt_struct_set_d2z_dxdy(
    void* struct_obj,
    double value_in);
void surface_displacement_struct_get_active(
    const void* struct_obj,
    bool* value_out);
void surface_displacement_struct_set_active(void* struct_obj, bool value_in);
void surface_displacement_struct_get_dr_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void surface_displacement_struct_get_r0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);

void surface_displacement_struct_get_pt_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void target_point_struct_get_r_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void surface_curvature_struct_get_xy_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void surface_curvature_struct_get_spherical(
    const void* struct_obj,
    double* value_out);
void surface_curvature_struct_set_spherical(void* struct_obj, double value_in);
void surface_curvature_struct_get_elliptical_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void surface_curvature_struct_get_has_curvature(
    const void* struct_obj,
    bool* value_out);
void surface_curvature_struct_set_has_curvature(
    void* struct_obj,
    bool value_in);
void photon_target_struct_get_type(const void* struct_obj, int* value_out);
void photon_target_struct_set_type(void* struct_obj, int value_in);
void photon_target_struct_get_n_corner(const void* struct_obj, int* value_out);
void photon_target_struct_set_n_corner(void* struct_obj, int value_in);
void photon_target_struct_get_ele_loc(const void* struct_obj, void** ptr_out);
void photon_target_struct_set_ele_loc(void* struct_obj, const void* src_ptr);

void photon_target_struct_get_corner_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void photon_target_struct_get_center(const void* struct_obj, void** ptr_out);
void photon_target_struct_set_center(void* struct_obj, const void* src_ptr);
void photon_material_struct_get_f0_m1(
    const void* struct_obj,
    std::complex<double>* value_out);
void photon_material_struct_set_f0_m1(
    void* struct_obj,
    std::complex<double> value_in);
void photon_material_struct_get_f0_m2(
    const void* struct_obj,
    std::complex<double>* value_out);
void photon_material_struct_set_f0_m2(
    void* struct_obj,
    std::complex<double> value_in);
void photon_material_struct_get_f_0(
    const void* struct_obj,
    std::complex<double>* value_out);
void photon_material_struct_set_f_0(
    void* struct_obj,
    std::complex<double> value_in);
void photon_material_struct_get_f_h(
    const void* struct_obj,
    std::complex<double>* value_out);
void photon_material_struct_set_f_h(
    void* struct_obj,
    std::complex<double> value_in);
void photon_material_struct_get_f_hbar(
    const void* struct_obj,
    std::complex<double>* value_out);
void photon_material_struct_set_f_hbar(
    void* struct_obj,
    std::complex<double> value_in);
void photon_material_struct_get_f_hkl(
    const void* struct_obj,
    std::complex<double>* value_out);
void photon_material_struct_set_f_hkl(
    void* struct_obj,
    std::complex<double> value_in);
void photon_material_struct_get_h_norm_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void photon_material_struct_get_l_ref_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void pixel_pt_struct_get_n_photon(const void* struct_obj, int64_t* value_out);
void pixel_pt_struct_set_n_photon(void* struct_obj, int64_t value_in);
void pixel_pt_struct_get_E_x(
    const void* struct_obj,
    std::complex<double>* value_out);
void pixel_pt_struct_set_E_x(void* struct_obj, std::complex<double> value_in);
void pixel_pt_struct_get_E_y(
    const void* struct_obj,
    std::complex<double>* value_out);
void pixel_pt_struct_set_E_y(void* struct_obj, std::complex<double> value_in);
void pixel_pt_struct_get_intensity_x(const void* struct_obj, double* value_out);
void pixel_pt_struct_set_intensity_x(void* struct_obj, double value_in);
void pixel_pt_struct_get_intensity_y(const void* struct_obj, double* value_out);
void pixel_pt_struct_set_intensity_y(void* struct_obj, double value_in);
void pixel_pt_struct_get_intensity(const void* struct_obj, double* value_out);
void pixel_pt_struct_set_intensity(void* struct_obj, double value_in);
void pixel_pt_struct_get_orbit_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void pixel_pt_struct_get_orbit_rms_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void pixel_pt_struct_get_init_orbit_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void pixel_pt_struct_get_init_orbit_rms_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void pixel_detec_struct_get_dr_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void pixel_detec_struct_get_r0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void pixel_detec_struct_get_n_track_tot(
    const void* struct_obj,
    int64_t* value_out);
void pixel_detec_struct_set_n_track_tot(void* struct_obj, int64_t value_in);
void pixel_detec_struct_get_n_hit_detec(
    const void* struct_obj,
    int64_t* value_out);
void pixel_detec_struct_set_n_hit_detec(void* struct_obj, int64_t value_in);
void pixel_detec_struct_get_n_hit_pixel(
    const void* struct_obj,
    int64_t* value_out);
void pixel_detec_struct_set_n_hit_pixel(void* struct_obj, int64_t value_in);

void pixel_detec_struct_get_pt_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void photon_element_struct_get_curvature(
    const void* struct_obj,
    void** ptr_out);
void photon_element_struct_set_curvature(void* struct_obj, const void* src_ptr);
void photon_element_struct_get_target(const void* struct_obj, void** ptr_out);
void photon_element_struct_set_target(void* struct_obj, const void* src_ptr);
void photon_element_struct_get_material(const void* struct_obj, void** ptr_out);
void photon_element_struct_set_material(void* struct_obj, const void* src_ptr);
void photon_element_struct_get_segmented(
    const void* struct_obj,
    void** ptr_out);
void photon_element_struct_set_segmented(void* struct_obj, const void* src_ptr);
void photon_element_struct_get_h_misalign(
    const void* struct_obj,
    void** ptr_out);
void photon_element_struct_set_h_misalign(
    void* struct_obj,
    const void* src_ptr);
void photon_element_struct_get_displacement(
    const void* struct_obj,
    void** ptr_out);
void photon_element_struct_set_displacement(
    void* struct_obj,
    const void* src_ptr);
void photon_element_struct_get_pixel(const void* struct_obj, void** ptr_out);
void photon_element_struct_set_pixel(void* struct_obj, const void* src_ptr);
void photon_element_struct_get_reflectivity_table_type(
    const void* struct_obj,
    int* value_out);
void photon_element_struct_set_reflectivity_table_type(
    void* struct_obj,
    int value_in);
void photon_element_struct_get_reflectivity_table_sigma(
    const void* struct_obj,
    void** ptr_out);
void photon_element_struct_set_reflectivity_table_sigma(
    void* struct_obj,
    const void* src_ptr);
void photon_element_struct_get_reflectivity_table_pi(
    const void* struct_obj,
    void** ptr_out);
void photon_element_struct_set_reflectivity_table_pi(
    void* struct_obj,
    const void* src_ptr);

void photon_element_struct_get_init_energy_prob_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void photon_element_struct_get_integrated_init_energy_prob_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void wall3d_vertex_struct_get_x(const void* struct_obj, double* value_out);
void wall3d_vertex_struct_set_x(void* struct_obj, double value_in);
void wall3d_vertex_struct_get_y(const void* struct_obj, double* value_out);
void wall3d_vertex_struct_set_y(void* struct_obj, double value_in);
void wall3d_vertex_struct_get_radius_x(
    const void* struct_obj,
    double* value_out);
void wall3d_vertex_struct_set_radius_x(void* struct_obj, double value_in);
void wall3d_vertex_struct_get_radius_y(
    const void* struct_obj,
    double* value_out);
void wall3d_vertex_struct_set_radius_y(void* struct_obj, double value_in);
void wall3d_vertex_struct_get_tilt(const void* struct_obj, double* value_out);
void wall3d_vertex_struct_set_tilt(void* struct_obj, double value_in);
void wall3d_vertex_struct_get_angle(const void* struct_obj, double* value_out);
void wall3d_vertex_struct_set_angle(void* struct_obj, double value_in);
void wall3d_vertex_struct_get_x0(const void* struct_obj, double* value_out);
void wall3d_vertex_struct_set_x0(void* struct_obj, double value_in);
void wall3d_vertex_struct_get_y0(const void* struct_obj, double* value_out);
void wall3d_vertex_struct_set_y0(void* struct_obj, double value_in);
void wall3d_vertex_struct_get_type(const void* struct_obj, int* value_out);
void wall3d_vertex_struct_set_type(void* struct_obj, int value_in);
void wall3d_section_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void wall3d_section_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void wall3d_section_struct_get_material_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void wall3d_section_struct_set_material(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void wall3d_section_struct_get_v_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void wall3d_section_struct_get_surface(const void* struct_obj, void** ptr_out);
void wall3d_section_struct_set_surface(void* struct_obj, const void* src_ptr);
void wall3d_section_struct_get_type(const void* struct_obj, int* value_out);
void wall3d_section_struct_set_type(void* struct_obj, int value_in);
void wall3d_section_struct_get_n_vertex_input(
    const void* struct_obj,
    int* value_out);
void wall3d_section_struct_set_n_vertex_input(void* struct_obj, int value_in);
void wall3d_section_struct_get_ix_ele(const void* struct_obj, int* value_out);
void wall3d_section_struct_set_ix_ele(void* struct_obj, int value_in);
void wall3d_section_struct_get_ix_branch(
    const void* struct_obj,
    int* value_out);
void wall3d_section_struct_set_ix_branch(void* struct_obj, int value_in);
void wall3d_section_struct_get_vertices_state(
    const void* struct_obj,
    int* value_out);
void wall3d_section_struct_set_vertices_state(void* struct_obj, int value_in);
void wall3d_section_struct_get_patch_in_region(
    const void* struct_obj,
    bool* value_out);
void wall3d_section_struct_set_patch_in_region(void* struct_obj, bool value_in);
void wall3d_section_struct_get_thickness(
    const void* struct_obj,
    double* value_out);
void wall3d_section_struct_set_thickness(void* struct_obj, double value_in);
void wall3d_section_struct_get_s(const void* struct_obj, double* value_out);
void wall3d_section_struct_set_s(void* struct_obj, double value_in);
void wall3d_section_struct_get_r0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void wall3d_section_struct_get_dx0_ds(
    const void* struct_obj,
    double* value_out);
void wall3d_section_struct_set_dx0_ds(void* struct_obj, double value_in);
void wall3d_section_struct_get_dy0_ds(
    const void* struct_obj,
    double* value_out);
void wall3d_section_struct_set_dy0_ds(void* struct_obj, double value_in);
void wall3d_section_struct_get_x0_coef_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void wall3d_section_struct_get_y0_coef_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void wall3d_section_struct_get_dr_ds(const void* struct_obj, double* value_out);
void wall3d_section_struct_set_dr_ds(void* struct_obj, double value_in);
void wall3d_section_struct_get_p1_coef_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void wall3d_section_struct_get_p2_coef_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void wall3d_struct_get_name_info(const void* s, char** d, int* bounds, bool* a);
void wall3d_struct_set_name(void* struct_obj, const char* str_ptr, int str_len);
void wall3d_struct_get_type(const void* struct_obj, int* value_out);
void wall3d_struct_set_type(void* struct_obj, int value_in);
void wall3d_struct_get_ix_wall3d(const void* struct_obj, int* value_out);
void wall3d_struct_set_ix_wall3d(void* struct_obj, int value_in);
void wall3d_struct_get_n_link(const void* struct_obj, int* value_out);
void wall3d_struct_set_n_link(void* struct_obj, int value_in);
void wall3d_struct_get_thickness(const void* struct_obj, double* value_out);
void wall3d_struct_set_thickness(void* struct_obj, double value_in);
void wall3d_struct_get_clear_material_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void wall3d_struct_set_clear_material(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void wall3d_struct_get_opaque_material_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void wall3d_struct_set_opaque_material(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void wall3d_struct_get_superimpose(const void* struct_obj, bool* value_out);
void wall3d_struct_set_superimpose(void* struct_obj, bool value_in);
void wall3d_struct_get_ele_anchor_pt(const void* struct_obj, int* value_out);
void wall3d_struct_set_ele_anchor_pt(void* struct_obj, int value_in);

void wall3d_struct_get_section_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void ramper_lord_struct_get_ix_ele(const void* struct_obj, int* value_out);
void ramper_lord_struct_set_ix_ele(void* struct_obj, int value_in);
void ramper_lord_struct_get_ix_con(const void* struct_obj, int* value_out);
void ramper_lord_struct_set_ix_con(void* struct_obj, int value_in);
void ramper_lord_struct_get_attrib_ptr(
    const void* struct_obj,
    double** ptr_out);
void ramper_lord_struct_set_attrib_ptr(void* struct_obj, double value_in);
void control_struct_get_value(const void* struct_obj, double* value_out);
void control_struct_set_value(void* struct_obj, double value_in);
void control_struct_get_y_knot_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);

void control_struct_get_stack_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void control_struct_get_slave(const void* struct_obj, void** ptr_out);
void control_struct_set_slave(void* struct_obj, const void* src_ptr);
void control_struct_get_lord(const void* struct_obj, void** ptr_out);
void control_struct_set_lord(void* struct_obj, const void* src_ptr);
void control_struct_get_slave_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void control_struct_set_slave_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void control_struct_get_attribute_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void control_struct_set_attribute(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void control_struct_get_ix_attrib(const void* struct_obj, int* value_out);
void control_struct_set_ix_attrib(void* struct_obj, int value_in);
void control_var1_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void control_var1_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void control_var1_struct_get_value(const void* struct_obj, double* value_out);
void control_var1_struct_set_value(void* struct_obj, double value_in);
void control_var1_struct_get_old_value(
    const void* struct_obj,
    double* value_out);
void control_var1_struct_set_old_value(void* struct_obj, double value_in);
void control_ramp1_struct_get_y_knot_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);

void control_ramp1_struct_get_stack_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void control_ramp1_struct_get_attribute_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void control_ramp1_struct_set_attribute(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void control_ramp1_struct_get_slave_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void control_ramp1_struct_set_slave_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void control_ramp1_struct_get_is_controller(
    const void* struct_obj,
    bool* value_out);
void control_ramp1_struct_set_is_controller(void* struct_obj, bool value_in);

void controller_struct_get_var_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void controller_struct_get_ramp_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void controller_struct_get_ramper_lord_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void controller_struct_get_x_knot_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ellipse_beam_init_struct_get_part_per_ellipse(
    const void* struct_obj,
    int* value_out);
void ellipse_beam_init_struct_set_part_per_ellipse(
    void* struct_obj,
    int value_in);
void ellipse_beam_init_struct_get_n_ellipse(
    const void* struct_obj,
    int* value_out);
void ellipse_beam_init_struct_set_n_ellipse(void* struct_obj, int value_in);
void ellipse_beam_init_struct_get_sigma_cutoff(
    const void* struct_obj,
    double* value_out);
void ellipse_beam_init_struct_set_sigma_cutoff(
    void* struct_obj,
    double value_in);
void kv_beam_init_struct_get_part_per_phi_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void kv_beam_init_struct_get_n_I2(const void* struct_obj, int* value_out);
void kv_beam_init_struct_set_n_I2(void* struct_obj, int value_in);
void kv_beam_init_struct_get_A(const void* struct_obj, double* value_out);
void kv_beam_init_struct_set_A(void* struct_obj, double value_in);
void grid_beam_init_struct_get_n_x(const void* struct_obj, int* value_out);
void grid_beam_init_struct_set_n_x(void* struct_obj, int value_in);
void grid_beam_init_struct_get_n_px(const void* struct_obj, int* value_out);
void grid_beam_init_struct_set_n_px(void* struct_obj, int value_in);
void grid_beam_init_struct_get_x_min(const void* struct_obj, double* value_out);
void grid_beam_init_struct_set_x_min(void* struct_obj, double value_in);
void grid_beam_init_struct_get_x_max(const void* struct_obj, double* value_out);
void grid_beam_init_struct_set_x_max(void* struct_obj, double value_in);
void grid_beam_init_struct_get_px_min(
    const void* struct_obj,
    double* value_out);
void grid_beam_init_struct_set_px_min(void* struct_obj, double value_in);
void grid_beam_init_struct_get_px_max(
    const void* struct_obj,
    double* value_out);
void grid_beam_init_struct_set_px_max(void* struct_obj, double value_in);
void beam_init_struct_get_position_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void beam_init_struct_set_position_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void beam_init_struct_get_distribution_type_info(
    const void* s,
    char** d,
    int* bounds, // [lower, upper]
    int* str_len,
    bool* is_alloc);

void beam_init_struct_get_spin_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);

void beam_init_struct_get_ellipse_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void beam_init_struct_get_KV(const void* struct_obj, void** ptr_out);
void beam_init_struct_set_KV(void* struct_obj, const void* src_ptr);

void beam_init_struct_get_grid_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void beam_init_struct_get_center_jitter_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void beam_init_struct_get_emit_jitter_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void beam_init_struct_get_sig_z_jitter(
    const void* struct_obj,
    double* value_out);
void beam_init_struct_set_sig_z_jitter(void* struct_obj, double value_in);
void beam_init_struct_get_sig_pz_jitter(
    const void* struct_obj,
    double* value_out);
void beam_init_struct_set_sig_pz_jitter(void* struct_obj, double value_in);
void beam_init_struct_get_n_particle(const void* struct_obj, int* value_out);
void beam_init_struct_set_n_particle(void* struct_obj, int value_in);
void beam_init_struct_get_renorm_center(
    const void* struct_obj,
    bool* value_out);
void beam_init_struct_set_renorm_center(void* struct_obj, bool value_in);
void beam_init_struct_get_renorm_sigma(const void* struct_obj, bool* value_out);
void beam_init_struct_set_renorm_sigma(void* struct_obj, bool value_in);
void beam_init_struct_get_random_engine_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void beam_init_struct_set_random_engine(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void beam_init_struct_get_random_gauss_converter_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void beam_init_struct_set_random_gauss_converter(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void beam_init_struct_get_random_sigma_cutoff(
    const void* struct_obj,
    double* value_out);
void beam_init_struct_set_random_sigma_cutoff(
    void* struct_obj,
    double value_in);
void beam_init_struct_get_a_norm_emit(
    const void* struct_obj,
    double* value_out);
void beam_init_struct_set_a_norm_emit(void* struct_obj, double value_in);
void beam_init_struct_get_b_norm_emit(
    const void* struct_obj,
    double* value_out);
void beam_init_struct_set_b_norm_emit(void* struct_obj, double value_in);
void beam_init_struct_get_a_emit(const void* struct_obj, double* value_out);
void beam_init_struct_set_a_emit(void* struct_obj, double value_in);
void beam_init_struct_get_b_emit(const void* struct_obj, double* value_out);
void beam_init_struct_set_b_emit(void* struct_obj, double value_in);
void beam_init_struct_get_dPz_dz(const void* struct_obj, double* value_out);
void beam_init_struct_set_dPz_dz(void* struct_obj, double value_in);
void beam_init_struct_get_center_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void beam_init_struct_get_t_offset(const void* struct_obj, double* value_out);
void beam_init_struct_set_t_offset(void* struct_obj, double value_in);
void beam_init_struct_get_dt_bunch(const void* struct_obj, double* value_out);
void beam_init_struct_set_dt_bunch(void* struct_obj, double value_in);
void beam_init_struct_get_sig_z(const void* struct_obj, double* value_out);
void beam_init_struct_set_sig_z(void* struct_obj, double value_in);
void beam_init_struct_get_sig_pz(const void* struct_obj, double* value_out);
void beam_init_struct_set_sig_pz(void* struct_obj, double value_in);
void beam_init_struct_get_bunch_charge(
    const void* struct_obj,
    double* value_out);
void beam_init_struct_set_bunch_charge(void* struct_obj, double value_in);
void beam_init_struct_get_n_bunch(const void* struct_obj, int* value_out);
void beam_init_struct_set_n_bunch(void* struct_obj, int value_in);
void beam_init_struct_get_ix_turn(const void* struct_obj, int* value_out);
void beam_init_struct_set_ix_turn(void* struct_obj, int value_in);
void beam_init_struct_get_species_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void beam_init_struct_set_species(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void beam_init_struct_get_full_6D_coupling_calc(
    const void* struct_obj,
    bool* value_out);
void beam_init_struct_set_full_6D_coupling_calc(
    void* struct_obj,
    bool value_in);
void beam_init_struct_get_use_particle_start(
    const void* struct_obj,
    bool* value_out);
void beam_init_struct_set_use_particle_start(void* struct_obj, bool value_in);
void beam_init_struct_get_use_t_coords(const void* struct_obj, bool* value_out);
void beam_init_struct_set_use_t_coords(void* struct_obj, bool value_in);
void beam_init_struct_get_use_z_as_t(const void* struct_obj, bool* value_out);
void beam_init_struct_set_use_z_as_t(void* struct_obj, bool value_in);
void beam_init_struct_get_file_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void beam_init_struct_set_file_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void lat_param_struct_get_n_part(const void* struct_obj, double* value_out);
void lat_param_struct_set_n_part(void* struct_obj, double value_in);
void lat_param_struct_get_total_length(
    const void* struct_obj,
    double* value_out);
void lat_param_struct_set_total_length(void* struct_obj, double value_in);
void lat_param_struct_get_unstable_factor(
    const void* struct_obj,
    double* value_out);
void lat_param_struct_set_unstable_factor(void* struct_obj, double value_in);
void lat_param_struct_get_t1_with_RF_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void lat_param_struct_get_t1_no_RF_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void lat_param_struct_get_spin_tune(const void* struct_obj, double* value_out);
void lat_param_struct_set_spin_tune(void* struct_obj, double value_in);
void lat_param_struct_get_particle(const void* struct_obj, int* value_out);
void lat_param_struct_set_particle(void* struct_obj, int value_in);
void lat_param_struct_get_default_tracking_species(
    const void* struct_obj,
    int* value_out);
void lat_param_struct_set_default_tracking_species(
    void* struct_obj,
    int value_in);
void lat_param_struct_get_geometry(const void* struct_obj, int* value_out);
void lat_param_struct_set_geometry(void* struct_obj, int value_in);
void lat_param_struct_get_ixx(const void* struct_obj, int* value_out);
void lat_param_struct_set_ixx(void* struct_obj, int value_in);
void lat_param_struct_get_stable(const void* struct_obj, bool* value_out);
void lat_param_struct_set_stable(void* struct_obj, bool value_in);
void lat_param_struct_get_live_branch(const void* struct_obj, bool* value_out);
void lat_param_struct_set_live_branch(void* struct_obj, bool value_in);
void lat_param_struct_get_g1_integral(
    const void* struct_obj,
    double* value_out);
void lat_param_struct_set_g1_integral(void* struct_obj, double value_in);
void lat_param_struct_get_g2_integral(
    const void* struct_obj,
    double* value_out);
void lat_param_struct_set_g2_integral(void* struct_obj, double value_in);
void lat_param_struct_get_g3_integral(
    const void* struct_obj,
    double* value_out);
void lat_param_struct_set_g3_integral(void* struct_obj, double value_in);
void lat_param_struct_get_bookkeeping_state(
    const void* struct_obj,
    void** ptr_out);
void lat_param_struct_set_bookkeeping_state(
    void* struct_obj,
    const void* src_ptr);
void lat_param_struct_get_beam_init(const void* struct_obj, void** ptr_out);
void lat_param_struct_set_beam_init(void* struct_obj, const void* src_ptr);
void mode_info_struct_get_stable(const void* struct_obj, bool* value_out);
void mode_info_struct_set_stable(void* struct_obj, bool value_in);
void mode_info_struct_get_tune(const void* struct_obj, double* value_out);
void mode_info_struct_set_tune(void* struct_obj, double value_in);
void mode_info_struct_get_emit(const void* struct_obj, double* value_out);
void mode_info_struct_set_emit(void* struct_obj, double value_in);
void mode_info_struct_get_chrom(const void* struct_obj, double* value_out);
void mode_info_struct_set_chrom(void* struct_obj, double value_in);
void mode_info_struct_get_sigma(const void* struct_obj, double* value_out);
void mode_info_struct_set_sigma(void* struct_obj, double value_in);
void mode_info_struct_get_sigmap(const void* struct_obj, double* value_out);
void mode_info_struct_set_sigmap(void* struct_obj, double value_in);
void pre_tracker_struct_get_who(const void* struct_obj, int* value_out);
void pre_tracker_struct_set_who(void* struct_obj, int value_in);
void pre_tracker_struct_get_ix_ele_start(
    const void* struct_obj,
    int* value_out);
void pre_tracker_struct_set_ix_ele_start(void* struct_obj, int value_in);
void pre_tracker_struct_get_ix_ele_end(const void* struct_obj, int* value_out);
void pre_tracker_struct_set_ix_ele_end(void* struct_obj, int value_in);
void pre_tracker_struct_get_input_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void pre_tracker_struct_set_input_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void anormal_mode_struct_get_emittance(
    const void* struct_obj,
    double* value_out);
void anormal_mode_struct_set_emittance(void* struct_obj, double value_in);
void anormal_mode_struct_get_emittance_no_vert(
    const void* struct_obj,
    double* value_out);
void anormal_mode_struct_set_emittance_no_vert(
    void* struct_obj,
    double value_in);
void anormal_mode_struct_get_synch_int_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void anormal_mode_struct_get_j_damp(const void* struct_obj, double* value_out);
void anormal_mode_struct_set_j_damp(void* struct_obj, double value_in);
void anormal_mode_struct_get_alpha_damp(
    const void* struct_obj,
    double* value_out);
void anormal_mode_struct_set_alpha_damp(void* struct_obj, double value_in);
void anormal_mode_struct_get_chrom(const void* struct_obj, double* value_out);
void anormal_mode_struct_set_chrom(void* struct_obj, double value_in);
void anormal_mode_struct_get_tune(const void* struct_obj, double* value_out);
void anormal_mode_struct_set_tune(void* struct_obj, double value_in);
void linac_normal_mode_struct_get_i2_E4(
    const void* struct_obj,
    double* value_out);
void linac_normal_mode_struct_set_i2_E4(void* struct_obj, double value_in);
void linac_normal_mode_struct_get_i3_E7(
    const void* struct_obj,
    double* value_out);
void linac_normal_mode_struct_set_i3_E7(void* struct_obj, double value_in);
void linac_normal_mode_struct_get_i5a_E6(
    const void* struct_obj,
    double* value_out);
void linac_normal_mode_struct_set_i5a_E6(void* struct_obj, double value_in);
void linac_normal_mode_struct_get_i5b_E6(
    const void* struct_obj,
    double* value_out);
void linac_normal_mode_struct_set_i5b_E6(void* struct_obj, double value_in);
void linac_normal_mode_struct_get_sig_E1(
    const void* struct_obj,
    double* value_out);
void linac_normal_mode_struct_set_sig_E1(void* struct_obj, double value_in);
void linac_normal_mode_struct_get_a_emittance_end(
    const void* struct_obj,
    double* value_out);
void linac_normal_mode_struct_set_a_emittance_end(
    void* struct_obj,
    double value_in);
void linac_normal_mode_struct_get_b_emittance_end(
    const void* struct_obj,
    double* value_out);
void linac_normal_mode_struct_set_b_emittance_end(
    void* struct_obj,
    double value_in);
void normal_modes_struct_get_synch_int_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void normal_modes_struct_get_sigE_E(const void* struct_obj, double* value_out);
void normal_modes_struct_set_sigE_E(void* struct_obj, double value_in);
void normal_modes_struct_get_sig_z(const void* struct_obj, double* value_out);
void normal_modes_struct_set_sig_z(void* struct_obj, double value_in);
void normal_modes_struct_get_e_loss(const void* struct_obj, double* value_out);
void normal_modes_struct_set_e_loss(void* struct_obj, double value_in);
void normal_modes_struct_get_rf_voltage(
    const void* struct_obj,
    double* value_out);
void normal_modes_struct_set_rf_voltage(void* struct_obj, double value_in);
void normal_modes_struct_get_pz_aperture(
    const void* struct_obj,
    double* value_out);
void normal_modes_struct_set_pz_aperture(void* struct_obj, double value_in);
void normal_modes_struct_get_pz_average(
    const void* struct_obj,
    double* value_out);
void normal_modes_struct_set_pz_average(void* struct_obj, double value_in);
void normal_modes_struct_get_momentum_compaction(
    const void* struct_obj,
    double* value_out);
void normal_modes_struct_set_momentum_compaction(
    void* struct_obj,
    double value_in);
void normal_modes_struct_get_dpz_damp(
    const void* struct_obj,
    double* value_out);
void normal_modes_struct_set_dpz_damp(void* struct_obj, double value_in);
void normal_modes_struct_get_a(const void* struct_obj, void** ptr_out);
void normal_modes_struct_set_a(void* struct_obj, const void* src_ptr);
void normal_modes_struct_get_b(const void* struct_obj, void** ptr_out);
void normal_modes_struct_set_b(void* struct_obj, const void* src_ptr);
void normal_modes_struct_get_z(const void* struct_obj, void** ptr_out);
void normal_modes_struct_set_z(void* struct_obj, const void* src_ptr);
void normal_modes_struct_get_lin(const void* struct_obj, void** ptr_out);
void normal_modes_struct_set_lin(void* struct_obj, const void* src_ptr);
void em_field_struct_get_E_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void em_field_struct_get_B_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void em_field_struct_get_dE_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void em_field_struct_get_dB_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void em_field_struct_get_phi(const void* struct_obj, double* value_out);
void em_field_struct_set_phi(void* struct_obj, double value_in);
void em_field_struct_get_phi_B(const void* struct_obj, double* value_out);
void em_field_struct_set_phi_B(void* struct_obj, double value_in);
void em_field_struct_get_A_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void strong_beam_struct_get_ix_slice(const void* struct_obj, int* value_out);
void strong_beam_struct_set_ix_slice(void* struct_obj, int value_in);
void strong_beam_struct_get_x_center(const void* struct_obj, double* value_out);
void strong_beam_struct_set_x_center(void* struct_obj, double value_in);
void strong_beam_struct_get_y_center(const void* struct_obj, double* value_out);
void strong_beam_struct_set_y_center(void* struct_obj, double value_in);
void strong_beam_struct_get_x_sigma(const void* struct_obj, double* value_out);
void strong_beam_struct_set_x_sigma(void* struct_obj, double value_in);
void strong_beam_struct_get_y_sigma(const void* struct_obj, double* value_out);
void strong_beam_struct_set_y_sigma(void* struct_obj, double value_in);
void strong_beam_struct_get_dx(const void* struct_obj, double* value_out);
void strong_beam_struct_set_dx(void* struct_obj, double value_in);
void strong_beam_struct_get_dy(const void* struct_obj, double* value_out);
void strong_beam_struct_set_dy(void* struct_obj, double value_in);
void track_point_struct_get_s_lab(const void* struct_obj, double* value_out);
void track_point_struct_set_s_lab(void* struct_obj, double value_in);
void track_point_struct_get_s_body(const void* struct_obj, double* value_out);
void track_point_struct_set_s_body(void* struct_obj, double value_in);
void track_point_struct_get_orb(const void* struct_obj, void** ptr_out);
void track_point_struct_set_orb(void* struct_obj, const void* src_ptr);
void track_point_struct_get_field(const void* struct_obj, void** ptr_out);
void track_point_struct_set_field(void* struct_obj, const void* src_ptr);
void track_point_struct_get_strong_beam(const void* struct_obj, void** ptr_out);
void track_point_struct_set_strong_beam(void* struct_obj, const void* src_ptr);
void track_point_struct_get_vec0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void track_point_struct_get_mat6_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);

void track_struct_get_pt_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void track_struct_get_ds_save(const void* struct_obj, double* value_out);
void track_struct_set_ds_save(void* struct_obj, double value_in);
void track_struct_get_n_pt(const void* struct_obj, int* value_out);
void track_struct_set_n_pt(void* struct_obj, int value_in);
void track_struct_get_n_bad(const void* struct_obj, int* value_out);
void track_struct_set_n_bad(void* struct_obj, int value_in);
void track_struct_get_n_ok(const void* struct_obj, int* value_out);
void track_struct_set_n_ok(void* struct_obj, int value_in);
void space_charge_common_struct_get_ds_track_step(
    const void* struct_obj,
    double* value_out);
void space_charge_common_struct_set_ds_track_step(
    void* struct_obj,
    double value_in);
void space_charge_common_struct_get_dt_track_step(
    const void* struct_obj,
    double* value_out);
void space_charge_common_struct_set_dt_track_step(
    void* struct_obj,
    double value_in);
void space_charge_common_struct_get_cathode_strength_cutoff(
    const void* struct_obj,
    double* value_out);
void space_charge_common_struct_set_cathode_strength_cutoff(
    void* struct_obj,
    double value_in);
void space_charge_common_struct_get_rel_tol_tracking(
    const void* struct_obj,
    double* value_out);
void space_charge_common_struct_set_rel_tol_tracking(
    void* struct_obj,
    double value_in);
void space_charge_common_struct_get_abs_tol_tracking(
    const void* struct_obj,
    double* value_out);
void space_charge_common_struct_set_abs_tol_tracking(
    void* struct_obj,
    double value_in);
void space_charge_common_struct_get_beam_chamber_height(
    const void* struct_obj,
    double* value_out);
void space_charge_common_struct_set_beam_chamber_height(
    void* struct_obj,
    double value_in);
void space_charge_common_struct_get_lsc_sigma_cutoff(
    const void* struct_obj,
    double* value_out);
void space_charge_common_struct_set_lsc_sigma_cutoff(
    void* struct_obj,
    double value_in);
void space_charge_common_struct_get_particle_sigma_cutoff(
    const void* struct_obj,
    double* value_out);
void space_charge_common_struct_set_particle_sigma_cutoff(
    void* struct_obj,
    double value_in);
void space_charge_common_struct_get_space_charge_mesh_size_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void space_charge_common_struct_get_csr3d_mesh_size_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void space_charge_common_struct_get_n_bin(
    const void* struct_obj,
    int* value_out);
void space_charge_common_struct_set_n_bin(void* struct_obj, int value_in);
void space_charge_common_struct_get_particle_bin_span(
    const void* struct_obj,
    int* value_out);
void space_charge_common_struct_set_particle_bin_span(
    void* struct_obj,
    int value_in);
void space_charge_common_struct_get_n_shield_images(
    const void* struct_obj,
    int* value_out);
void space_charge_common_struct_set_n_shield_images(
    void* struct_obj,
    int value_in);
void space_charge_common_struct_get_sc_min_in_bin(
    const void* struct_obj,
    int* value_out);
void space_charge_common_struct_set_sc_min_in_bin(
    void* struct_obj,
    int value_in);
void space_charge_common_struct_get_lsc_kick_transverse_dependence(
    const void* struct_obj,
    bool* value_out);
void space_charge_common_struct_set_lsc_kick_transverse_dependence(
    void* struct_obj,
    bool value_in);
void space_charge_common_struct_get_debug(
    const void* struct_obj,
    bool* value_out);
void space_charge_common_struct_set_debug(void* struct_obj, bool value_in);
void space_charge_common_struct_get_diagnostic_output_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void space_charge_common_struct_set_diagnostic_output_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void bmad_common_struct_get_max_aperture_limit(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_max_aperture_limit(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_d_orb_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void bmad_common_struct_get_default_ds_step(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_default_ds_step(void* struct_obj, double value_in);
void bmad_common_struct_get_significant_length(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_significant_length(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_rel_tol_tracking(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_rel_tol_tracking(void* struct_obj, double value_in);
void bmad_common_struct_get_abs_tol_tracking(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_abs_tol_tracking(void* struct_obj, double value_in);
void bmad_common_struct_get_rel_tol_adaptive_tracking(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_rel_tol_adaptive_tracking(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_abs_tol_adaptive_tracking(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_abs_tol_adaptive_tracking(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_init_ds_adaptive_tracking(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_init_ds_adaptive_tracking(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_min_ds_adaptive_tracking(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_min_ds_adaptive_tracking(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_fatal_ds_adaptive_tracking(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_fatal_ds_adaptive_tracking(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_autoscale_amp_abs_tol(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_autoscale_amp_abs_tol(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_autoscale_amp_rel_tol(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_autoscale_amp_rel_tol(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_autoscale_phase_tol(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_autoscale_phase_tol(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_electric_dipole_moment(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_electric_dipole_moment(
    void* struct_obj,
    double value_in);
void bmad_common_struct_get_synch_rad_scale(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_synch_rad_scale(void* struct_obj, double value_in);
void bmad_common_struct_get_sad_eps_scale(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_sad_eps_scale(void* struct_obj, double value_in);
void bmad_common_struct_get_sad_amp_max(
    const void* struct_obj,
    double* value_out);
void bmad_common_struct_set_sad_amp_max(void* struct_obj, double value_in);
void bmad_common_struct_get_sad_n_div_max(
    const void* struct_obj,
    int* value_out);
void bmad_common_struct_set_sad_n_div_max(void* struct_obj, int value_in);
void bmad_common_struct_get_taylor_order(
    const void* struct_obj,
    int* value_out);
void bmad_common_struct_set_taylor_order(void* struct_obj, int value_in);
void bmad_common_struct_get_runge_kutta_order(
    const void* struct_obj,
    int* value_out);
void bmad_common_struct_set_runge_kutta_order(void* struct_obj, int value_in);
void bmad_common_struct_get_default_integ_order(
    const void* struct_obj,
    int* value_out);
void bmad_common_struct_set_default_integ_order(void* struct_obj, int value_in);
void bmad_common_struct_get_max_num_runge_kutta_step(
    const void* struct_obj,
    int* value_out);
void bmad_common_struct_set_max_num_runge_kutta_step(
    void* struct_obj,
    int value_in);
void bmad_common_struct_get_rf_phase_below_transition_ref(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_rf_phase_below_transition_ref(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_sr_wakes_on(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_sr_wakes_on(void* struct_obj, bool value_in);
void bmad_common_struct_get_lr_wakes_on(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_lr_wakes_on(void* struct_obj, bool value_in);
void bmad_common_struct_get_auto_bookkeeper(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_auto_bookkeeper(void* struct_obj, bool value_in);
void bmad_common_struct_get_high_energy_space_charge_on(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_high_energy_space_charge_on(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_csr_and_space_charge_on(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_csr_and_space_charge_on(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_spin_tracking_on(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_spin_tracking_on(void* struct_obj, bool value_in);
void bmad_common_struct_get_spin_sokolov_ternov_flipping_on(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_spin_sokolov_ternov_flipping_on(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_radiation_damping_on(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_radiation_damping_on(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_radiation_zero_average(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_radiation_zero_average(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_radiation_fluctuations_on(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_radiation_fluctuations_on(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_conserve_taylor_maps(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_conserve_taylor_maps(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_absolute_time_tracking(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_absolute_time_tracking(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_absolute_time_ref_shift(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_absolute_time_ref_shift(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_convert_to_kinetic_momentum(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_convert_to_kinetic_momentum(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_normalize_twiss(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_normalize_twiss(void* struct_obj, bool value_in);
void bmad_common_struct_get_aperture_limit_on(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_aperture_limit_on(void* struct_obj, bool value_in);
void bmad_common_struct_get_spin_n0_direction_user_set(
    const void* struct_obj,
    bool* value_out);
void bmad_common_struct_set_spin_n0_direction_user_set(
    void* struct_obj,
    bool value_in);
void bmad_common_struct_get_debug(const void* struct_obj, bool* value_out);
void bmad_common_struct_set_debug(void* struct_obj, bool value_in);
void rad_int1_struct_get_i0(const void* struct_obj, double* value_out);
void rad_int1_struct_set_i0(void* struct_obj, double value_in);
void rad_int1_struct_get_i1(const void* struct_obj, double* value_out);
void rad_int1_struct_set_i1(void* struct_obj, double value_in);
void rad_int1_struct_get_i2(const void* struct_obj, double* value_out);
void rad_int1_struct_set_i2(void* struct_obj, double value_in);
void rad_int1_struct_get_i3(const void* struct_obj, double* value_out);
void rad_int1_struct_set_i3(void* struct_obj, double value_in);
void rad_int1_struct_get_i4a(const void* struct_obj, double* value_out);
void rad_int1_struct_set_i4a(void* struct_obj, double value_in);
void rad_int1_struct_get_i4b(const void* struct_obj, double* value_out);
void rad_int1_struct_set_i4b(void* struct_obj, double value_in);
void rad_int1_struct_get_i4z(const void* struct_obj, double* value_out);
void rad_int1_struct_set_i4z(void* struct_obj, double value_in);
void rad_int1_struct_get_i5a(const void* struct_obj, double* value_out);
void rad_int1_struct_set_i5a(void* struct_obj, double value_in);
void rad_int1_struct_get_i5b(const void* struct_obj, double* value_out);
void rad_int1_struct_set_i5b(void* struct_obj, double value_in);
void rad_int1_struct_get_i6b(const void* struct_obj, double* value_out);
void rad_int1_struct_set_i6b(void* struct_obj, double value_in);
void rad_int1_struct_get_lin_i2_E4(const void* struct_obj, double* value_out);
void rad_int1_struct_set_lin_i2_E4(void* struct_obj, double value_in);
void rad_int1_struct_get_lin_i3_E7(const void* struct_obj, double* value_out);
void rad_int1_struct_set_lin_i3_E7(void* struct_obj, double value_in);
void rad_int1_struct_get_lin_i5a_E6(const void* struct_obj, double* value_out);
void rad_int1_struct_set_lin_i5a_E6(void* struct_obj, double value_in);
void rad_int1_struct_get_lin_i5b_E6(const void* struct_obj, double* value_out);
void rad_int1_struct_set_lin_i5b_E6(void* struct_obj, double value_in);
void rad_int1_struct_get_lin_norm_emit_a(
    const void* struct_obj,
    double* value_out);
void rad_int1_struct_set_lin_norm_emit_a(void* struct_obj, double value_in);
void rad_int1_struct_get_lin_norm_emit_b(
    const void* struct_obj,
    double* value_out);
void rad_int1_struct_set_lin_norm_emit_b(void* struct_obj, double value_in);
void rad_int1_struct_get_lin_sig_E(const void* struct_obj, double* value_out);
void rad_int1_struct_set_lin_sig_E(void* struct_obj, double value_in);
void rad_int1_struct_get_n_steps(const void* struct_obj, double* value_out);
void rad_int1_struct_set_n_steps(void* struct_obj, double value_in);

void rad_int_branch_struct_get_ele_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void rad_int_all_ele_struct_get_branch_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void rf_stair_step_struct_get_E_tot0(const void* struct_obj, double* value_out);
void rf_stair_step_struct_set_E_tot0(void* struct_obj, double value_in);
void rf_stair_step_struct_get_E_tot1(const void* struct_obj, double* value_out);
void rf_stair_step_struct_set_E_tot1(void* struct_obj, double value_in);
void rf_stair_step_struct_get_p0c(const void* struct_obj, double* value_out);
void rf_stair_step_struct_set_p0c(void* struct_obj, double value_in);
void rf_stair_step_struct_get_p1c(const void* struct_obj, double* value_out);
void rf_stair_step_struct_set_p1c(void* struct_obj, double value_in);
void rf_stair_step_struct_get_scale(const void* struct_obj, double* value_out);
void rf_stair_step_struct_set_scale(void* struct_obj, double value_in);
void rf_stair_step_struct_get_time(const void* struct_obj, double* value_out);
void rf_stair_step_struct_set_time(void* struct_obj, double value_in);
void rf_stair_step_struct_get_s0(const void* struct_obj, double* value_out);
void rf_stair_step_struct_set_s0(void* struct_obj, double value_in);
void rf_stair_step_struct_get_s(const void* struct_obj, double* value_out);
void rf_stair_step_struct_set_s(void* struct_obj, double value_in);
void rf_stair_step_struct_get_ix_step(const void* struct_obj, int* value_out);
void rf_stair_step_struct_set_ix_step(void* struct_obj, int value_in);

void rf_ele_struct_get_steps_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void rf_ele_struct_get_ds_step(const void* struct_obj, double* value_out);
void rf_ele_struct_set_ds_step(void* struct_obj, double value_in);
void ele_struct_get_name_info(const void* s, char** d, int* bounds, bool* a);
void ele_struct_set_name(void* struct_obj, const char* str_ptr, int str_len);
void ele_struct_get_type_info(const void* s, char** d, int* bounds, bool* a);
void ele_struct_set_type(void* struct_obj, const char* str_ptr, int str_len);
void ele_struct_get_alias_info(const void* s, char** d, int* bounds, bool* a);
void ele_struct_set_alias(void* struct_obj, const char* str_ptr, int str_len);
void ele_struct_get_component_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void ele_struct_set_component_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void ele_struct_get_descrip_info(
    const void* s,
    char** d,
    int* len,
    bool* is_alloc);

void ele_struct_set_descrip(void* struct_obj, const char* str_ptr, int str_len);
void ele_struct_get_a(const void* struct_obj, void** ptr_out);
void ele_struct_set_a(void* struct_obj, const void* src_ptr);
void ele_struct_get_b(const void* struct_obj, void** ptr_out);
void ele_struct_set_b(void* struct_obj, const void* src_ptr);
void ele_struct_get_z(const void* struct_obj, void** ptr_out);
void ele_struct_set_z(void* struct_obj, const void* src_ptr);
void ele_struct_get_x(const void* struct_obj, void** ptr_out);
void ele_struct_set_x(void* struct_obj, const void* src_ptr);
void ele_struct_get_y(const void* struct_obj, void** ptr_out);
void ele_struct_set_y(void* struct_obj, const void* src_ptr);
void ele_struct_get_ac_kick(const void* struct_obj, void** ptr_out);
void ele_struct_set_ac_kick(void* struct_obj, const void* src_ptr);
void ele_struct_get_bookkeeping_state(const void* struct_obj, void** ptr_out);
void ele_struct_set_bookkeeping_state(void* struct_obj, const void* src_ptr);
void ele_struct_get_branch(const void* struct_obj, void** ptr_out);
void ele_struct_set_branch(void* struct_obj, const void* src_ptr);
void ele_struct_get_control(const void* struct_obj, void** ptr_out);
void ele_struct_set_control(void* struct_obj, const void* src_ptr);
void ele_struct_get_rf(const void* struct_obj, void** ptr_out);
void ele_struct_set_rf(void* struct_obj, const void* src_ptr);
void ele_struct_get_lord(const void* struct_obj, void** ptr_out);
void ele_struct_set_lord(void* struct_obj, const void* src_ptr);
void ele_struct_get_floor(const void* struct_obj, void** ptr_out);
void ele_struct_set_floor(void* struct_obj, const void* src_ptr);
void ele_struct_get_high_energy_space_charge(
    const void* struct_obj,
    void** ptr_out);
void ele_struct_set_high_energy_space_charge(
    void* struct_obj,
    const void* src_ptr);
void ele_struct_get_mode3(const void* struct_obj, void** ptr_out);
void ele_struct_set_mode3(void* struct_obj, const void* src_ptr);
void ele_struct_get_photon(const void* struct_obj, void** ptr_out);
void ele_struct_set_photon(void* struct_obj, const void* src_ptr);
void ele_struct_get_rad_map(const void* struct_obj, void** ptr_out);
void ele_struct_set_rad_map(void* struct_obj, const void* src_ptr);

void ele_struct_get_taylor_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void ele_struct_get_spin_taylor_ref_orb_in_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);

void ele_struct_get_spin_taylor_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void ele_struct_get_wake(const void* struct_obj, void** ptr_out);
void ele_struct_set_wake(void* struct_obj, const void* src_ptr);

void ele_struct_get_wall3d_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void ele_struct_get_cartesian_map_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void ele_struct_get_cylindrical_map_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void ele_struct_get_gen_grad_map_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void ele_struct_get_grid_field_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void ele_struct_get_map_ref_orb_in(const void* struct_obj, void** ptr_out);
void ele_struct_set_map_ref_orb_in(void* struct_obj, const void* src_ptr);
void ele_struct_get_map_ref_orb_out(const void* struct_obj, void** ptr_out);
void ele_struct_set_map_ref_orb_out(void* struct_obj, const void* src_ptr);
void ele_struct_get_time_ref_orb_in(const void* struct_obj, void** ptr_out);
void ele_struct_set_time_ref_orb_in(void* struct_obj, const void* src_ptr);
void ele_struct_get_time_ref_orb_out(const void* struct_obj, void** ptr_out);
void ele_struct_set_time_ref_orb_out(void* struct_obj, const void* src_ptr);
void ele_struct_get_value_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ele_struct_get_old_value_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ele_struct_get_spin_q_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void ele_struct_get_vec0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ele_struct_get_mat6_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void ele_struct_get_c_mat_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void ele_struct_get_dc_mat_dpz_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void ele_struct_get_gamma_c(const void* struct_obj, double* value_out);
void ele_struct_set_gamma_c(void* struct_obj, double value_in);
void ele_struct_get_s_start(const void* struct_obj, double* value_out);
void ele_struct_set_s_start(void* struct_obj, double value_in);
void ele_struct_get_s(const void* struct_obj, double* value_out);
void ele_struct_set_s(void* struct_obj, double value_in);
void ele_struct_get_ref_time(const void* struct_obj, double* value_out);
void ele_struct_set_ref_time(void* struct_obj, double value_in);
void ele_struct_get_a_pole_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ele_struct_get_b_pole_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ele_struct_get_a_pole_elec_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ele_struct_get_b_pole_elec_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ele_struct_get_custom_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ele_struct_get_r_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void ele_struct_get_key(const void* struct_obj, int* value_out);
void ele_struct_set_key(void* struct_obj, int value_in);
void ele_struct_get_sub_key(const void* struct_obj, int* value_out);
void ele_struct_set_sub_key(void* struct_obj, int value_in);
void ele_struct_get_ix_ele(const void* struct_obj, int* value_out);
void ele_struct_set_ix_ele(void* struct_obj, int value_in);
void ele_struct_get_ix_branch(const void* struct_obj, int* value_out);
void ele_struct_set_ix_branch(void* struct_obj, int value_in);
void ele_struct_get_lord_status(const void* struct_obj, int* value_out);
void ele_struct_set_lord_status(void* struct_obj, int value_in);
void ele_struct_get_n_slave(const void* struct_obj, int* value_out);
void ele_struct_set_n_slave(void* struct_obj, int value_in);
void ele_struct_get_n_slave_field(const void* struct_obj, int* value_out);
void ele_struct_set_n_slave_field(void* struct_obj, int value_in);
void ele_struct_get_ix1_slave(const void* struct_obj, int* value_out);
void ele_struct_set_ix1_slave(void* struct_obj, int value_in);
void ele_struct_get_slave_status(const void* struct_obj, int* value_out);
void ele_struct_set_slave_status(void* struct_obj, int value_in);
void ele_struct_get_n_lord(const void* struct_obj, int* value_out);
void ele_struct_set_n_lord(void* struct_obj, int value_in);
void ele_struct_get_n_lord_field(const void* struct_obj, int* value_out);
void ele_struct_set_n_lord_field(void* struct_obj, int value_in);
void ele_struct_get_n_lord_ramper(const void* struct_obj, int* value_out);
void ele_struct_set_n_lord_ramper(void* struct_obj, int value_in);
void ele_struct_get_ic1_lord(const void* struct_obj, int* value_out);
void ele_struct_set_ic1_lord(void* struct_obj, int value_in);
void ele_struct_get_ix_pointer(const void* struct_obj, int* value_out);
void ele_struct_set_ix_pointer(void* struct_obj, int value_in);
void ele_struct_get_ixx(const void* struct_obj, int* value_out);
void ele_struct_set_ixx(void* struct_obj, int value_in);
void ele_struct_get_iyy(const void* struct_obj, int* value_out);
void ele_struct_set_iyy(void* struct_obj, int value_in);
void ele_struct_get_izz(const void* struct_obj, int* value_out);
void ele_struct_set_izz(void* struct_obj, int value_in);
void ele_struct_get_mat6_calc_method(const void* struct_obj, int* value_out);
void ele_struct_set_mat6_calc_method(void* struct_obj, int value_in);
void ele_struct_get_tracking_method(const void* struct_obj, int* value_out);
void ele_struct_set_tracking_method(void* struct_obj, int value_in);
void ele_struct_get_spin_tracking_method(
    const void* struct_obj,
    int* value_out);
void ele_struct_set_spin_tracking_method(void* struct_obj, int value_in);
void ele_struct_get_csr_method(const void* struct_obj, int* value_out);
void ele_struct_set_csr_method(void* struct_obj, int value_in);
void ele_struct_get_space_charge_method(const void* struct_obj, int* value_out);
void ele_struct_set_space_charge_method(void* struct_obj, int value_in);
void ele_struct_get_ptc_integration_type(
    const void* struct_obj,
    int* value_out);
void ele_struct_set_ptc_integration_type(void* struct_obj, int value_in);
void ele_struct_get_field_calc(const void* struct_obj, int* value_out);
void ele_struct_set_field_calc(void* struct_obj, int value_in);
void ele_struct_get_aperture_at(const void* struct_obj, int* value_out);
void ele_struct_set_aperture_at(void* struct_obj, int value_in);
void ele_struct_get_aperture_type(const void* struct_obj, int* value_out);
void ele_struct_set_aperture_type(void* struct_obj, int value_in);
void ele_struct_get_ref_species(const void* struct_obj, int* value_out);
void ele_struct_set_ref_species(void* struct_obj, int value_in);
void ele_struct_get_orientation(const void* struct_obj, int* value_out);
void ele_struct_set_orientation(void* struct_obj, int value_in);
void ele_struct_get_symplectify(const void* struct_obj, bool* value_out);
void ele_struct_set_symplectify(void* struct_obj, bool value_in);
void ele_struct_get_mode_flip(const void* struct_obj, bool* value_out);
void ele_struct_set_mode_flip(void* struct_obj, bool value_in);
void ele_struct_get_multipoles_on(const void* struct_obj, bool* value_out);
void ele_struct_set_multipoles_on(void* struct_obj, bool value_in);
void ele_struct_get_scale_multipoles(const void* struct_obj, bool* value_out);
void ele_struct_set_scale_multipoles(void* struct_obj, bool value_in);
void ele_struct_get_taylor_map_includes_offsets(
    const void* struct_obj,
    bool* value_out);
void ele_struct_set_taylor_map_includes_offsets(
    void* struct_obj,
    bool value_in);
void ele_struct_get_field_master(const void* struct_obj, bool* value_out);
void ele_struct_set_field_master(void* struct_obj, bool value_in);
void ele_struct_get_is_on(const void* struct_obj, bool* value_out);
void ele_struct_set_is_on(void* struct_obj, bool value_in);
void ele_struct_get_logic(const void* struct_obj, bool* value_out);
void ele_struct_set_logic(void* struct_obj, bool value_in);
void ele_struct_get_bmad_logic(const void* struct_obj, bool* value_out);
void ele_struct_set_bmad_logic(void* struct_obj, bool value_in);
void ele_struct_get_select(const void* struct_obj, bool* value_out);
void ele_struct_set_select(void* struct_obj, bool value_in);
void ele_struct_get_offset_moves_aperture(
    const void* struct_obj,
    bool* value_out);
void ele_struct_set_offset_moves_aperture(void* struct_obj, bool value_in);
void complex_taylor_term_struct_get_coef(
    const void* struct_obj,
    std::complex<double>* value_out);
void complex_taylor_term_struct_set_coef(
    void* struct_obj,
    std::complex<double> value_in);
void complex_taylor_term_struct_get_expn_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void complex_taylor_struct_get_ref(
    const void* struct_obj,
    std::complex<double>* value_out);
void complex_taylor_struct_set_ref(
    void* struct_obj,
    std::complex<double> value_in);

void complex_taylor_struct_get_term_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void branch_struct_get_name_info(const void* s, char** d, int* bounds, bool* a);
void branch_struct_set_name(void* struct_obj, const char* str_ptr, int str_len);
void branch_struct_get_ix_branch(const void* struct_obj, int* value_out);
void branch_struct_set_ix_branch(void* struct_obj, int value_in);
void branch_struct_get_ix_from_branch(const void* struct_obj, int* value_out);
void branch_struct_set_ix_from_branch(void* struct_obj, int value_in);
void branch_struct_get_ix_from_ele(const void* struct_obj, int* value_out);
void branch_struct_set_ix_from_ele(void* struct_obj, int value_in);
void branch_struct_get_ix_to_ele(const void* struct_obj, int* value_out);
void branch_struct_set_ix_to_ele(void* struct_obj, int value_in);
void branch_struct_get_ix_fixer(const void* struct_obj, int* value_out);
void branch_struct_set_ix_fixer(void* struct_obj, int value_in);
void branch_struct_get_n_ele_track(const void* struct_obj, int* value_out);
void branch_struct_set_n_ele_track(void* struct_obj, int value_in);
void branch_struct_get_n_ele_max(const void* struct_obj, int* value_out);
void branch_struct_set_n_ele_max(void* struct_obj, int value_in);
void branch_struct_get_lat(const void* struct_obj, void** ptr_out);
void branch_struct_set_lat(void* struct_obj, const void* src_ptr);
void branch_struct_get_a(const void* struct_obj, void** ptr_out);
void branch_struct_set_a(void* struct_obj, const void* src_ptr);
void branch_struct_get_b(const void* struct_obj, void** ptr_out);
void branch_struct_set_b(void* struct_obj, const void* src_ptr);
void branch_struct_get_z(const void* struct_obj, void** ptr_out);
void branch_struct_set_z(void* struct_obj, const void* src_ptr);

void branch_struct_get_ele_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void branch_struct_get_param(const void* struct_obj, void** ptr_out);
void branch_struct_set_param(void* struct_obj, const void* src_ptr);
void branch_struct_get_particle_start(const void* struct_obj, void** ptr_out);
void branch_struct_set_particle_start(void* struct_obj, const void* src_ptr);

void branch_struct_get_wall3d_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void lat_struct_get_use_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void lat_struct_set_use_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void lat_struct_get_lattice_info(const void* s, char** d, int* bounds, bool* a);
void lat_struct_set_lattice(void* struct_obj, const char* str_ptr, int str_len);
void lat_struct_get_machine_info(const void* s, char** d, int* bounds, bool* a);
void lat_struct_set_machine(void* struct_obj, const char* str_ptr, int str_len);
void lat_struct_get_input_file_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void lat_struct_set_input_file_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void lat_struct_get_title_info(const void* s, char** d, int* bounds, bool* a);
void lat_struct_set_title(void* struct_obj, const char* str_ptr, int str_len);

void lat_struct_get_print_str_info(
    const void* s,
    char** d,
    int* bounds, // [lower, upper]
    int* str_len,
    bool* is_alloc);

void lat_struct_get_constant_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void lat_struct_get_a(const void* struct_obj, void** ptr_out);
void lat_struct_set_a(void* struct_obj, const void* src_ptr);
void lat_struct_get_b(const void* struct_obj, void** ptr_out);
void lat_struct_set_b(void* struct_obj, const void* src_ptr);
void lat_struct_get_z(const void* struct_obj, void** ptr_out);
void lat_struct_set_z(void* struct_obj, const void* src_ptr);
void lat_struct_get_param(const void* struct_obj, void** ptr_out);
void lat_struct_set_param(void* struct_obj, const void* src_ptr);
void lat_struct_get_lord_state(const void* struct_obj, void** ptr_out);
void lat_struct_set_lord_state(void* struct_obj, const void* src_ptr);
void lat_struct_get_ele_init(const void* struct_obj, void** ptr_out);
void lat_struct_set_ele_init(void* struct_obj, const void* src_ptr);

void lat_struct_get_ele_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void lat_struct_get_branch_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void lat_struct_get_control_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void lat_struct_get_particle_start(const void* struct_obj, void** ptr_out);
void lat_struct_set_particle_start(void* struct_obj, const void* src_ptr);
void lat_struct_get_beam_init(const void* struct_obj, void** ptr_out);
void lat_struct_set_beam_init(void* struct_obj, const void* src_ptr);
void lat_struct_get_pre_tracker(const void* struct_obj, void** ptr_out);
void lat_struct_set_pre_tracker(void* struct_obj, const void* src_ptr);
void lat_struct_get_custom_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void lat_struct_get_version(const void* struct_obj, int* value_out);
void lat_struct_set_version(void* struct_obj, int value_in);
void lat_struct_get_n_ele_track(const void* struct_obj, int** ptr_out);
void lat_struct_set_n_ele_track(void* struct_obj, int value_in);
void lat_struct_get_n_ele_max(const void* struct_obj, int** ptr_out);
void lat_struct_set_n_ele_max(void* struct_obj, int value_in);
void lat_struct_get_n_control_max(const void* struct_obj, int* value_out);
void lat_struct_set_n_control_max(void* struct_obj, int value_in);
void lat_struct_get_n_ic_max(const void* struct_obj, int* value_out);
void lat_struct_set_n_ic_max(void* struct_obj, int value_in);
void lat_struct_get_input_taylor_order(const void* struct_obj, int* value_out);
void lat_struct_set_input_taylor_order(void* struct_obj, int value_in);
void lat_struct_get_ic_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void lat_struct_get_photon_type(const void* struct_obj, int* value_out);
void lat_struct_set_photon_type(void* struct_obj, int value_in);
void lat_struct_get_creation_hash(const void* struct_obj, int* value_out);
void lat_struct_set_creation_hash(void* struct_obj, int value_in);
void lat_struct_get_ramper_slave_bookkeeping(
    const void* struct_obj,
    int* value_out);
void lat_struct_set_ramper_slave_bookkeeping(void* struct_obj, int value_in);

void bunch_struct_get_particle_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bunch_struct_get_ix_z_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void bunch_struct_get_charge_tot(const void* struct_obj, double* value_out);
void bunch_struct_set_charge_tot(void* struct_obj, double value_in);
void bunch_struct_get_charge_live(const void* struct_obj, double* value_out);
void bunch_struct_set_charge_live(void* struct_obj, double value_in);
void bunch_struct_get_z_center(const void* struct_obj, double* value_out);
void bunch_struct_set_z_center(void* struct_obj, double value_in);
void bunch_struct_get_t_center(const void* struct_obj, double* value_out);
void bunch_struct_set_t_center(void* struct_obj, double value_in);
void bunch_struct_get_t0(const void* struct_obj, double* value_out);
void bunch_struct_set_t0(void* struct_obj, double value_in);
void bunch_struct_get_drift_between_t_and_s(
    const void* struct_obj,
    bool* value_out);
void bunch_struct_set_drift_between_t_and_s(void* struct_obj, bool value_in);
void bunch_struct_get_ix_ele(const void* struct_obj, int* value_out);
void bunch_struct_set_ix_ele(void* struct_obj, int value_in);
void bunch_struct_get_ix_bunch(const void* struct_obj, int* value_out);
void bunch_struct_set_ix_bunch(void* struct_obj, int value_in);
void bunch_struct_get_ix_turn(const void* struct_obj, int* value_out);
void bunch_struct_set_ix_turn(void* struct_obj, int value_in);
void bunch_struct_get_n_live(const void* struct_obj, int* value_out);
void bunch_struct_set_n_live(void* struct_obj, int value_in);
void bunch_struct_get_n_good(const void* struct_obj, int* value_out);
void bunch_struct_set_n_good(void* struct_obj, int value_in);
void bunch_struct_get_n_bad(const void* struct_obj, int* value_out);
void bunch_struct_set_n_bad(void* struct_obj, int value_in);
void bunch_params_struct_get_centroid(const void* struct_obj, void** ptr_out);
void bunch_params_struct_set_centroid(void* struct_obj, const void* src_ptr);
void bunch_params_struct_get_x(const void* struct_obj, void** ptr_out);
void bunch_params_struct_set_x(void* struct_obj, const void* src_ptr);
void bunch_params_struct_get_y(const void* struct_obj, void** ptr_out);
void bunch_params_struct_set_y(void* struct_obj, const void* src_ptr);
void bunch_params_struct_get_z(const void* struct_obj, void** ptr_out);
void bunch_params_struct_set_z(void* struct_obj, const void* src_ptr);
void bunch_params_struct_get_a(const void* struct_obj, void** ptr_out);
void bunch_params_struct_set_a(void* struct_obj, const void* src_ptr);
void bunch_params_struct_get_b(const void* struct_obj, void** ptr_out);
void bunch_params_struct_set_b(void* struct_obj, const void* src_ptr);
void bunch_params_struct_get_c(const void* struct_obj, void** ptr_out);
void bunch_params_struct_set_c(void* struct_obj, const void* src_ptr);
void bunch_params_struct_get_sigma_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void bunch_params_struct_get_rel_max_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void bunch_params_struct_get_rel_min_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void bunch_params_struct_get_s(const void* struct_obj, double* value_out);
void bunch_params_struct_set_s(void* struct_obj, double value_in);
void bunch_params_struct_get_t(const void* struct_obj, double* value_out);
void bunch_params_struct_set_t(void* struct_obj, double value_in);
void bunch_params_struct_get_sigma_t(const void* struct_obj, double* value_out);
void bunch_params_struct_set_sigma_t(void* struct_obj, double value_in);
void bunch_params_struct_get_charge_live(
    const void* struct_obj,
    double* value_out);
void bunch_params_struct_set_charge_live(void* struct_obj, double value_in);
void bunch_params_struct_get_charge_tot(
    const void* struct_obj,
    double* value_out);
void bunch_params_struct_set_charge_tot(void* struct_obj, double value_in);
void bunch_params_struct_get_n_particle_tot(
    const void* struct_obj,
    int* value_out);
void bunch_params_struct_set_n_particle_tot(void* struct_obj, int value_in);
void bunch_params_struct_get_n_particle_live(
    const void* struct_obj,
    int* value_out);
void bunch_params_struct_set_n_particle_live(void* struct_obj, int value_in);
void bunch_params_struct_get_n_particle_lost_in_ele(
    const void* struct_obj,
    int* value_out);
void bunch_params_struct_set_n_particle_lost_in_ele(
    void* struct_obj,
    int value_in);
void bunch_params_struct_get_n_good_steps(
    const void* struct_obj,
    int* value_out);
void bunch_params_struct_set_n_good_steps(void* struct_obj, int value_in);
void bunch_params_struct_get_n_bad_steps(
    const void* struct_obj,
    int* value_out);
void bunch_params_struct_set_n_bad_steps(void* struct_obj, int value_in);
void bunch_params_struct_get_ix_ele(const void* struct_obj, int* value_out);
void bunch_params_struct_set_ix_ele(void* struct_obj, int value_in);
void bunch_params_struct_get_location(const void* struct_obj, int* value_out);
void bunch_params_struct_set_location(void* struct_obj, int value_in);
void bunch_params_struct_get_twiss_valid(
    const void* struct_obj,
    bool* value_out);
void bunch_params_struct_set_twiss_valid(void* struct_obj, bool value_in);

void beam_struct_get_bunch_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void aperture_point_struct_get_x(const void* struct_obj, double* value_out);
void aperture_point_struct_set_x(void* struct_obj, double value_in);
void aperture_point_struct_get_y(const void* struct_obj, double* value_out);
void aperture_point_struct_set_y(void* struct_obj, double value_in);
void aperture_point_struct_get_plane(const void* struct_obj, int* value_out);
void aperture_point_struct_set_plane(void* struct_obj, int value_in);
void aperture_point_struct_get_ix_ele(const void* struct_obj, int* value_out);
void aperture_point_struct_set_ix_ele(void* struct_obj, int value_in);
void aperture_point_struct_get_i_turn(const void* struct_obj, int* value_out);
void aperture_point_struct_set_i_turn(void* struct_obj, int value_in);
void aperture_param_struct_get_min_angle(
    const void* struct_obj,
    double* value_out);
void aperture_param_struct_set_min_angle(void* struct_obj, double value_in);
void aperture_param_struct_get_max_angle(
    const void* struct_obj,
    double* value_out);
void aperture_param_struct_set_max_angle(void* struct_obj, double value_in);
void aperture_param_struct_get_n_angle(const void* struct_obj, int* value_out);
void aperture_param_struct_set_n_angle(void* struct_obj, int value_in);
void aperture_param_struct_get_n_turn(const void* struct_obj, int* value_out);
void aperture_param_struct_set_n_turn(void* struct_obj, int value_in);
void aperture_param_struct_get_x_init(
    const void* struct_obj,
    double* value_out);
void aperture_param_struct_set_x_init(void* struct_obj, double value_in);
void aperture_param_struct_get_y_init(
    const void* struct_obj,
    double* value_out);
void aperture_param_struct_set_y_init(void* struct_obj, double value_in);
void aperture_param_struct_get_rel_accuracy(
    const void* struct_obj,
    double* value_out);
void aperture_param_struct_set_rel_accuracy(void* struct_obj, double value_in);
void aperture_param_struct_get_abs_accuracy(
    const void* struct_obj,
    double* value_out);
void aperture_param_struct_set_abs_accuracy(void* struct_obj, double value_in);
void aperture_param_struct_get_start_ele_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void aperture_param_struct_set_start_ele(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void aperture_scan_struct_get_point_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void aperture_scan_struct_get_ref_orb(const void* struct_obj, void** ptr_out);
void aperture_scan_struct_set_ref_orb(void* struct_obj, const void* src_ptr);
void aperture_scan_struct_get_pz_start(
    const void* struct_obj,
    double* value_out);
void aperture_scan_struct_set_pz_start(void* struct_obj, double value_in);
void ele_pointer_struct_get_ele(const void* struct_obj, void** ptr_out);
void ele_pointer_struct_set_ele(void* struct_obj, const void* src_ptr);
void ele_pointer_struct_get_loc(const void* struct_obj, void** ptr_out);
void ele_pointer_struct_set_loc(void* struct_obj, const void* src_ptr);
void ele_pointer_struct_get_id(const void* struct_obj, int* value_out);
void ele_pointer_struct_set_id(void* struct_obj, int value_in);
void expression_tree_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void expression_tree_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void expression_tree_struct_get_type(const void* struct_obj, int* value_out);
void expression_tree_struct_set_type(void* struct_obj, int value_in);
void expression_tree_struct_get_value(
    const void* struct_obj,
    double* value_out);
void expression_tree_struct_set_value(void* struct_obj, double value_in);

void expression_tree_struct_get_node_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void nametable_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds, // [lower, upper]
    int* str_len,
    bool* is_alloc);

void nametable_struct_get_index_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void nametable_struct_get_n_min(const void* struct_obj, int* value_out);
void nametable_struct_set_n_min(void* struct_obj, int value_in);
void nametable_struct_get_n_max(const void* struct_obj, int* value_out);
void nametable_struct_set_n_max(void* struct_obj, int value_in);
void tao_spin_dn_dpz_struct_get_vec_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_spin_dn_dpz_struct_get_partial_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void tao_spin_dn_dpz_struct_get_partial2_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void resonance_h_struct_get_id_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void resonance_h_struct_set_id(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void resonance_h_struct_get_c_val(
    const void* struct_obj,
    std::complex<double>* value_out);
void resonance_h_struct_set_c_val(
    void* struct_obj,
    std::complex<double> value_in);
void spin_orbit_map1_struct_get_orb_mat_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void spin_orbit_map1_struct_get_vec0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void spin_orbit_map1_struct_get_spin_q_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void spin_axis_struct_get_l_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void spin_axis_struct_get_n0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void spin_axis_struct_get_m_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ptc_normal_form_struct_get_ele_origin(
    const void* struct_obj,
    void** ptr_out);
void ptc_normal_form_struct_set_ele_origin(
    void* struct_obj,
    const void* src_ptr);
void ptc_normal_form_struct_get_orb0_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void ptc_normal_form_struct_get_valid_map(
    const void* struct_obj,
    bool* value_out);
void ptc_normal_form_struct_set_valid_map(void* struct_obj, bool value_in);
void bmad_normal_form_struct_get_ele_origin(
    const void* struct_obj,
    void** ptr_out);
void bmad_normal_form_struct_set_ele_origin(
    void* struct_obj,
    const void* src_ptr);

void bmad_normal_form_struct_get_M_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bmad_normal_form_struct_get_A_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bmad_normal_form_struct_get_A_inv_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bmad_normal_form_struct_get_dhdj_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bmad_normal_form_struct_get_F_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bmad_normal_form_struct_get_L_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bmad_normal_form_struct_get_h_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bunch_track_struct_get_pt_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bunch_track_struct_get_ds_save(const void* struct_obj, double* value_out);
void bunch_track_struct_set_ds_save(void* struct_obj, double value_in);
void bunch_track_struct_get_n_pt(const void* struct_obj, int* value_out);
void bunch_track_struct_set_n_pt(void* struct_obj, int value_in);
void summation_rdt_struct_get_h11001(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h11001(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h00111(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h00111(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h20001(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h20001(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h00201(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h00201(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h10002(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h10002(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h21000(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h21000(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h30000(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h30000(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h10110(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h10110(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h10020(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h10020(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h10200(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h10200(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h31000(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h31000(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h40000(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h40000(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h20110(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h20110(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h11200(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h11200(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h20020(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h20020(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h20200(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h20200(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h00310(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h00310(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h00400(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h00400(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h22000(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h22000(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h00220(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h00220(
    void* struct_obj,
    std::complex<double> value_in);
void summation_rdt_struct_get_h11110(
    const void* struct_obj,
    std::complex<double>* value_out);
void summation_rdt_struct_set_h11110(
    void* struct_obj,
    std::complex<double> value_in);
void tao_ele_shape_struct_get_ele_id_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_ele_shape_struct_set_ele_id(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_ele_shape_struct_get_shape_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_ele_shape_struct_set_shape(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_ele_shape_struct_get_color_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_ele_shape_struct_set_color(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_ele_shape_struct_get_size(const void* struct_obj, double* value_out);
void tao_ele_shape_struct_set_size(void* struct_obj, double value_in);
void tao_ele_shape_struct_get_label_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_ele_shape_struct_set_label(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_ele_shape_struct_get_draw(const void* struct_obj, bool* value_out);
void tao_ele_shape_struct_set_draw(void* struct_obj, bool value_in);
void tao_ele_shape_struct_get_multi(const void* struct_obj, bool* value_out);
void tao_ele_shape_struct_set_multi(void* struct_obj, bool value_in);
void tao_ele_shape_struct_get_line_width(
    const void* struct_obj,
    int* value_out);
void tao_ele_shape_struct_set_line_width(void* struct_obj, int value_in);
void tao_ele_shape_struct_get_offset(const void* struct_obj, double* value_out);
void tao_ele_shape_struct_set_offset(void* struct_obj, double value_in);
void tao_ele_shape_struct_get_ix_key(const void* struct_obj, int* value_out);
void tao_ele_shape_struct_set_ix_key(void* struct_obj, int value_in);
void tao_ele_shape_struct_get_name_ele_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_ele_shape_struct_set_name_ele(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void tao_ele_shape_struct_get_uni_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_ele_pointer_struct_get_eles_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_ele_pointer_struct_get_n_loc(const void* struct_obj, int* value_out);
void tao_ele_pointer_struct_set_n_loc(void* struct_obj, int value_in);
void tao_curve_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_curve_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_curve_struct_get_data_source_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_curve_struct_set_data_source(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_curve_struct_get_data_index_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_curve_struct_set_data_index(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_curve_struct_get_data_type_x_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_curve_struct_set_data_type_x(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void tao_curve_struct_get_data_type_info(
    const void* s,
    char** d,
    int* len,
    bool* is_alloc);

void tao_curve_struct_set_data_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_curve_struct_get_ele_ref_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_curve_struct_set_ele_ref_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_curve_struct_get_legend_text_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_curve_struct_set_legend_text(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_curve_struct_get_message_text_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_curve_struct_set_message_text(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_curve_struct_get_component_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_curve_struct_set_component(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_curve_struct_get_why_invalid_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_curve_struct_set_why_invalid(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_curve_struct_get_g(const void* struct_obj, void** ptr_out);
void tao_curve_struct_set_g(void* struct_obj, const void* src_ptr);
void tao_curve_struct_get_hist(const void* struct_obj, void** ptr_out);
void tao_curve_struct_set_hist(void* struct_obj, const void* src_ptr);
void tao_curve_struct_get_z_color(const void* struct_obj, void** ptr_out);
void tao_curve_struct_set_z_color(void* struct_obj, const void* src_ptr);
void tao_curve_struct_get_x_line_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_curve_struct_get_y_line_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_curve_struct_get_y2_line_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_curve_struct_get_ix_line_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void tao_curve_struct_get_x_symb_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_curve_struct_get_y_symb_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_curve_struct_get_z_symb_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_curve_struct_get_err_symb_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_curve_struct_get_symb_size_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_curve_struct_get_ix_symb_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void tao_curve_struct_get_y_axis_scale_factor(
    const void* struct_obj,
    double* value_out);
void tao_curve_struct_set_y_axis_scale_factor(
    void* struct_obj,
    double value_in);
void tao_curve_struct_get_line(const void* struct_obj, void** ptr_out);
void tao_curve_struct_set_line(void* struct_obj, const void* src_ptr);
void tao_curve_struct_get_symbol(const void* struct_obj, void** ptr_out);
void tao_curve_struct_set_symbol(void* struct_obj, const void* src_ptr);
void tao_curve_struct_get_orbit(const void* struct_obj, void** ptr_out);
void tao_curve_struct_set_orbit(void* struct_obj, const void* src_ptr);
void tao_curve_struct_get_ix_universe(const void* struct_obj, int* value_out);
void tao_curve_struct_set_ix_universe(void* struct_obj, int value_in);
void tao_curve_struct_get_symbol_every(const void* struct_obj, int* value_out);
void tao_curve_struct_set_symbol_every(void* struct_obj, int value_in);
void tao_curve_struct_get_ix_branch(const void* struct_obj, int* value_out);
void tao_curve_struct_set_ix_branch(void* struct_obj, int value_in);
void tao_curve_struct_get_ix_bunch(const void* struct_obj, int* value_out);
void tao_curve_struct_set_ix_bunch(void* struct_obj, int value_in);
void tao_curve_struct_get_n_turn(const void* struct_obj, int* value_out);
void tao_curve_struct_set_n_turn(void* struct_obj, int value_in);
void tao_curve_struct_get_use_y2(const void* struct_obj, bool* value_out);
void tao_curve_struct_set_use_y2(void* struct_obj, bool value_in);
void tao_curve_struct_get_draw_line(const void* struct_obj, bool* value_out);
void tao_curve_struct_set_draw_line(void* struct_obj, bool value_in);
void tao_curve_struct_get_draw_symbols(const void* struct_obj, bool* value_out);
void tao_curve_struct_set_draw_symbols(void* struct_obj, bool value_in);
void tao_curve_struct_get_draw_symbol_index(
    const void* struct_obj,
    bool* value_out);
void tao_curve_struct_set_draw_symbol_index(void* struct_obj, bool value_in);
void tao_curve_struct_get_draw_error_bars(
    const void* struct_obj,
    bool* value_out);
void tao_curve_struct_set_draw_error_bars(void* struct_obj, bool value_in);
void tao_curve_struct_get_smooth_line_calc(
    const void* struct_obj,
    bool* value_out);
void tao_curve_struct_set_smooth_line_calc(void* struct_obj, bool value_in);
void tao_curve_struct_get_valid(const void* struct_obj, bool* value_out);
void tao_curve_struct_set_valid(void* struct_obj, bool value_in);
void tao_curve_color_struct_get_data_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_curve_color_struct_set_data_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_curve_color_struct_get_is_on(const void* struct_obj, bool* value_out);
void tao_curve_color_struct_set_is_on(void* struct_obj, bool value_in);
void tao_curve_color_struct_get_min(const void* struct_obj, double* value_out);
void tao_curve_color_struct_set_min(void* struct_obj, double value_in);
void tao_curve_color_struct_get_max(const void* struct_obj, double* value_out);
void tao_curve_color_struct_set_max(void* struct_obj, double value_in);
void tao_curve_color_struct_get_autoscale(
    const void* struct_obj,
    bool* value_out);
void tao_curve_color_struct_set_autoscale(void* struct_obj, bool value_in);
void tao_curve_orbit_struct_get_x(const void* struct_obj, double* value_out);
void tao_curve_orbit_struct_set_x(void* struct_obj, double value_in);
void tao_curve_orbit_struct_get_y(const void* struct_obj, double* value_out);
void tao_curve_orbit_struct_set_y(void* struct_obj, double value_in);
void tao_curve_orbit_struct_get_t(const void* struct_obj, double* value_out);
void tao_curve_orbit_struct_set_t(void* struct_obj, double value_in);
void tao_histogram_struct_get_density_normalized(
    const void* struct_obj,
    bool* value_out);
void tao_histogram_struct_set_density_normalized(
    void* struct_obj,
    bool value_in);
void tao_histogram_struct_get_weight_by_charge(
    const void* struct_obj,
    bool* value_out);
void tao_histogram_struct_set_weight_by_charge(void* struct_obj, bool value_in);
void tao_histogram_struct_get_minimum(
    const void* struct_obj,
    double* value_out);
void tao_histogram_struct_set_minimum(void* struct_obj, double value_in);
void tao_histogram_struct_get_maximum(
    const void* struct_obj,
    double* value_out);
void tao_histogram_struct_set_maximum(void* struct_obj, double value_in);
void tao_histogram_struct_get_width(const void* struct_obj, double* value_out);
void tao_histogram_struct_set_width(void* struct_obj, double value_in);
void tao_histogram_struct_get_center(const void* struct_obj, double* value_out);
void tao_histogram_struct_set_center(void* struct_obj, double value_in);
void tao_histogram_struct_get_number(const void* struct_obj, int* value_out);
void tao_histogram_struct_set_number(void* struct_obj, int value_in);
void lat_ele_order1_struct_get_ix_branch(
    const void* struct_obj,
    int* value_out);
void lat_ele_order1_struct_set_ix_branch(void* struct_obj, int value_in);
void lat_ele_order1_struct_get_ix_order(const void* struct_obj, int* value_out);
void lat_ele_order1_struct_set_ix_order(void* struct_obj, int value_in);

void lat_ele_order_array_struct_get_ele_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lat_sigma_struct_get_mat_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void tao_spin_ele_struct_get_dn_dpz(const void* struct_obj, void** ptr_out);
void tao_spin_ele_struct_set_dn_dpz(void* struct_obj, const void* src_ptr);
void tao_spin_ele_struct_get_orb_eigen_val_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_spin_ele_struct_get_orb_eigen_vec_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void tao_spin_ele_struct_get_spin_eigen_vec_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void tao_spin_ele_struct_get_valid(const void* struct_obj, bool* value_out);
void tao_spin_ele_struct_set_valid(void* struct_obj, bool value_in);
void tao_plot_cache_struct_get_ele_to_s(const void* struct_obj, void** ptr_out);
void tao_plot_cache_struct_set_ele_to_s(void* struct_obj, const void* src_ptr);
void tao_plot_cache_struct_get_orbit(const void* struct_obj, void** ptr_out);
void tao_plot_cache_struct_set_orbit(void* struct_obj, const void* src_ptr);
void tao_plot_cache_struct_get_err(const void* struct_obj, bool* value_out);
void tao_plot_cache_struct_set_err(void* struct_obj, bool value_in);
void tao_spin_polarization_struct_get_tune(
    const void* struct_obj,
    double* value_out);
void tao_spin_polarization_struct_set_tune(void* struct_obj, double value_in);
void tao_spin_polarization_struct_get_pol_limit_st(
    const void* struct_obj,
    double* value_out);
void tao_spin_polarization_struct_set_pol_limit_st(
    void* struct_obj,
    double value_in);
void tao_spin_polarization_struct_get_pol_limit_dk(
    const void* struct_obj,
    double* value_out);
void tao_spin_polarization_struct_set_pol_limit_dk(
    void* struct_obj,
    double value_in);
void tao_spin_polarization_struct_get_pol_limit_dk_partial_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_spin_polarization_struct_get_pol_limit_dk_partial2_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_spin_polarization_struct_get_pol_rate_bks(
    const void* struct_obj,
    double* value_out);
void tao_spin_polarization_struct_set_pol_rate_bks(
    void* struct_obj,
    double value_in);
void tao_spin_polarization_struct_get_depol_rate(
    const void* struct_obj,
    double* value_out);
void tao_spin_polarization_struct_set_depol_rate(
    void* struct_obj,
    double value_in);
void tao_spin_polarization_struct_get_depol_rate_partial_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_spin_polarization_struct_get_depol_rate_partial2_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_spin_polarization_struct_get_integral_bn(
    const void* struct_obj,
    double* value_out);
void tao_spin_polarization_struct_set_integral_bn(
    void* struct_obj,
    double value_in);
void tao_spin_polarization_struct_get_integral_bdn(
    const void* struct_obj,
    double* value_out);
void tao_spin_polarization_struct_set_integral_bdn(
    void* struct_obj,
    double value_in);
void tao_spin_polarization_struct_get_integral_1ns(
    const void* struct_obj,
    double* value_out);
void tao_spin_polarization_struct_set_integral_1ns(
    void* struct_obj,
    double value_in);
void tao_spin_polarization_struct_get_integral_dn2(
    const void* struct_obj,
    double* value_out);
void tao_spin_polarization_struct_set_integral_dn2(
    void* struct_obj,
    double value_in);
void tao_spin_polarization_struct_get_valid(
    const void* struct_obj,
    bool* value_out);
void tao_spin_polarization_struct_set_valid(void* struct_obj, bool value_in);
void tao_spin_polarization_struct_get_q_1turn(
    const void* struct_obj,
    void** ptr_out);
void tao_spin_polarization_struct_set_q_1turn(
    void* struct_obj,
    const void* src_ptr);

void tao_spin_polarization_struct_get_q_ele_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lattice_branch_struct_get_tao_lat(
    const void* struct_obj,
    void** ptr_out);
void tao_lattice_branch_struct_set_tao_lat(
    void* struct_obj,
    const void* src_ptr);

void tao_lattice_branch_struct_get_lat_sigma_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lattice_branch_struct_get_spin_ele_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lattice_branch_struct_get_bunch_params_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lattice_branch_struct_get_bunch_params_comb_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lattice_branch_struct_get_orbit_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lattice_branch_struct_get_plot_cache_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lattice_branch_struct_get_spin(const void* struct_obj, void** ptr_out);
void tao_lattice_branch_struct_set_spin(void* struct_obj, const void* src_ptr);
void tao_lattice_branch_struct_get_srdt(const void* struct_obj, void** ptr_out);
void tao_lattice_branch_struct_set_srdt(void* struct_obj, const void* src_ptr);
void tao_lattice_branch_struct_get_orb0(const void* struct_obj, void** ptr_out);
void tao_lattice_branch_struct_set_orb0(void* struct_obj, const void* src_ptr);
void tao_lattice_branch_struct_get_modes_ri(
    const void* struct_obj,
    void** ptr_out);
void tao_lattice_branch_struct_set_modes_ri(
    void* struct_obj,
    const void* src_ptr);
void tao_lattice_branch_struct_get_modes_6d(
    const void* struct_obj,
    void** ptr_out);
void tao_lattice_branch_struct_set_modes_6d(
    void* struct_obj,
    const void* src_ptr);
void tao_lattice_branch_struct_get_ptc_normal_form(
    const void* struct_obj,
    void** ptr_out);
void tao_lattice_branch_struct_set_ptc_normal_form(
    void* struct_obj,
    const void* src_ptr);
void tao_lattice_branch_struct_get_bmad_normal_form(
    const void* struct_obj,
    void** ptr_out);
void tao_lattice_branch_struct_set_bmad_normal_form(
    void* struct_obj,
    const void* src_ptr);

void tao_lattice_branch_struct_get_high_E_orb_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lattice_branch_struct_get_low_E_orb_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lattice_branch_struct_get_taylor_save_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_lattice_branch_struct_get_cache_x_min(
    const void* struct_obj,
    double* value_out);
void tao_lattice_branch_struct_set_cache_x_min(
    void* struct_obj,
    double value_in);
void tao_lattice_branch_struct_get_cache_x_max(
    const void* struct_obj,
    double* value_out);
void tao_lattice_branch_struct_set_cache_x_max(
    void* struct_obj,
    double value_in);
void tao_lattice_branch_struct_get_comb_ds_save(
    const void* struct_obj,
    double* value_out);
void tao_lattice_branch_struct_set_comb_ds_save(
    void* struct_obj,
    double value_in);
void tao_lattice_branch_struct_get_ix_ref_taylor(
    const void* struct_obj,
    int* value_out);
void tao_lattice_branch_struct_set_ix_ref_taylor(
    void* struct_obj,
    int value_in);
void tao_lattice_branch_struct_get_ix_ele_taylor(
    const void* struct_obj,
    int* value_out);
void tao_lattice_branch_struct_set_ix_ele_taylor(
    void* struct_obj,
    int value_in);
void tao_lattice_branch_struct_get_track_state(
    const void* struct_obj,
    int* value_out);
void tao_lattice_branch_struct_set_track_state(void* struct_obj, int value_in);
void tao_lattice_branch_struct_get_cache_n_pts(
    const void* struct_obj,
    int* value_out);
void tao_lattice_branch_struct_set_cache_n_pts(void* struct_obj, int value_in);
void tao_lattice_branch_struct_get_ix_rad_int_cache(
    const void* struct_obj,
    int* value_out);
void tao_lattice_branch_struct_set_ix_rad_int_cache(
    void* struct_obj,
    int value_in);
void tao_lattice_branch_struct_get_has_open_match_element(
    const void* struct_obj,
    bool* value_out);
void tao_lattice_branch_struct_set_has_open_match_element(
    void* struct_obj,
    bool value_in);
void tao_lattice_branch_struct_get_plot_cache_valid(
    const void* struct_obj,
    bool* value_out);
void tao_lattice_branch_struct_set_plot_cache_valid(
    void* struct_obj,
    bool value_in);
void tao_lattice_branch_struct_get_spin_map_valid(
    const void* struct_obj,
    bool* value_out);
void tao_lattice_branch_struct_set_spin_map_valid(
    void* struct_obj,
    bool value_in);
void tao_lattice_branch_struct_get_twiss_valid(
    const void* struct_obj,
    bool* value_out);
void tao_lattice_branch_struct_set_twiss_valid(void* struct_obj, bool value_in);
void tao_lattice_branch_struct_get_mode_flip_here(
    const void* struct_obj,
    bool* value_out);
void tao_lattice_branch_struct_set_mode_flip_here(
    void* struct_obj,
    bool value_in);
void tao_lattice_branch_struct_get_chrom_calc_ok(
    const void* struct_obj,
    bool* value_out);
void tao_lattice_branch_struct_set_chrom_calc_ok(
    void* struct_obj,
    bool value_in);
void tao_lattice_branch_struct_get_rad_int_calc_ok(
    const void* struct_obj,
    bool* value_out);
void tao_lattice_branch_struct_set_rad_int_calc_ok(
    void* struct_obj,
    bool value_in);
void tao_lattice_branch_struct_get_emit_6d_calc_ok(
    const void* struct_obj,
    bool* value_out);
void tao_lattice_branch_struct_set_emit_6d_calc_ok(
    void* struct_obj,
    bool value_in);
void tao_lattice_branch_struct_get_sigma_track_ok(
    const void* struct_obj,
    bool* value_out);
void tao_lattice_branch_struct_set_sigma_track_ok(
    void* struct_obj,
    bool value_in);
void tao_model_element_struct_get_beam(const void* struct_obj, void** ptr_out);
void tao_model_element_struct_set_beam(void* struct_obj, const void* src_ptr);
void tao_model_element_struct_get_save_beam_internally(
    const void* struct_obj,
    bool* value_out);
void tao_model_element_struct_set_save_beam_internally(
    void* struct_obj,
    bool value_in);
void tao_model_element_struct_get_save_beam_to_file(
    const void* struct_obj,
    bool* value_out);
void tao_model_element_struct_set_save_beam_to_file(
    void* struct_obj,
    bool value_in);
void tao_beam_branch_struct_get_beam_at_start(
    const void* struct_obj,
    void** ptr_out);
void tao_beam_branch_struct_set_beam_at_start(
    void* struct_obj,
    const void* src_ptr);
void tao_beam_branch_struct_get_beam_init(
    const void* struct_obj,
    void** ptr_out);
void tao_beam_branch_struct_set_beam_init(
    void* struct_obj,
    const void* src_ptr);
void tao_beam_branch_struct_get_beam_init_used(
    const void* struct_obj,
    void** ptr_out);
void tao_beam_branch_struct_set_beam_init_used(
    void* struct_obj,
    const void* src_ptr);
void tao_beam_branch_struct_get_init_starting_distribution(
    const void* struct_obj,
    bool* value_out);
void tao_beam_branch_struct_set_init_starting_distribution(
    void* struct_obj,
    bool value_in);
void tao_beam_branch_struct_get_track_start_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_beam_branch_struct_set_track_start(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_beam_branch_struct_get_track_end_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_beam_branch_struct_set_track_end(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_beam_branch_struct_get_ix_branch(
    const void* struct_obj,
    int* value_out);
void tao_beam_branch_struct_set_ix_branch(void* struct_obj, int value_in);
void tao_beam_branch_struct_get_ix_track_start(
    const void* struct_obj,
    int* value_out);
void tao_beam_branch_struct_set_ix_track_start(void* struct_obj, int value_in);
void tao_beam_branch_struct_get_ix_track_end(
    const void* struct_obj,
    int* value_out);
void tao_beam_branch_struct_set_ix_track_end(void* struct_obj, int value_in);
void tao_d1_data_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_d1_data_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_d1_data_struct_get_d2(const void* struct_obj, void** ptr_out);
void tao_d1_data_struct_set_d2(void* struct_obj, const void* src_ptr);

void tao_d1_data_struct_get_d_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_d2_data_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_d2_data_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_d2_data_struct_get_data_file_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_d2_data_struct_set_data_file_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_d2_data_struct_get_ref_file_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_d2_data_struct_set_ref_file_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_d2_data_struct_get_data_date_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_d2_data_struct_set_data_date(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_d2_data_struct_get_ref_date_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_d2_data_struct_set_ref_date(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void tao_d2_data_struct_get_descrip_info(
    const void* s,
    char** d,
    int* bounds, // [lower, upper]
    int* str_len,
    bool* is_alloc);

void tao_d2_data_struct_get_d1_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_d2_data_struct_get_ix_universe(const void* struct_obj, int* value_out);
void tao_d2_data_struct_set_ix_universe(void* struct_obj, int value_in);
void tao_d2_data_struct_get_ix_d2_data(const void* struct_obj, int* value_out);
void tao_d2_data_struct_set_ix_d2_data(void* struct_obj, int value_in);
void tao_d2_data_struct_get_ix_ref(const void* struct_obj, int* value_out);
void tao_d2_data_struct_set_ix_ref(void* struct_obj, int value_in);
void tao_d2_data_struct_get_data_read_in(
    const void* struct_obj,
    bool* value_out);
void tao_d2_data_struct_set_data_read_in(void* struct_obj, bool value_in);
void tao_d2_data_struct_get_ref_read_in(
    const void* struct_obj,
    bool* value_out);
void tao_d2_data_struct_set_ref_read_in(void* struct_obj, bool value_in);
void tao_data_var_component_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_data_var_component_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_data_var_component_struct_get_sign(
    const void* struct_obj,
    double* value_out);
void tao_data_var_component_struct_set_sign(void* struct_obj, double value_in);
void tao_graph_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_graph_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_graph_struct_get_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_graph_struct_set_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_graph_struct_get_title_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_graph_struct_set_title(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_graph_struct_get_title_suffix_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_graph_struct_set_title_suffix(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void tao_graph_struct_get_text_legend_info(
    const void* s,
    char** d,
    int* bounds, // [lower, upper]
    int* str_len,
    bool* is_alloc);

void tao_graph_struct_get_text_legend_out_info(
    const void* s,
    char** d,
    int* bounds, // [lower, upper]
    int* str_len,
    bool* is_alloc);

void tao_graph_struct_get_why_invalid_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_graph_struct_set_why_invalid(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void tao_graph_struct_get_curve_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_graph_struct_get_p(const void* struct_obj, void** ptr_out);
void tao_graph_struct_set_p(void* struct_obj, const void* src_ptr);
void tao_graph_struct_get_floor_plan(const void* struct_obj, void** ptr_out);
void tao_graph_struct_set_floor_plan(void* struct_obj, const void* src_ptr);
void tao_graph_struct_get_text_legend_origin(
    const void* struct_obj,
    void** ptr_out);
void tao_graph_struct_set_text_legend_origin(
    void* struct_obj,
    const void* src_ptr);
void tao_graph_struct_get_curve_legend_origin(
    const void* struct_obj,
    void** ptr_out);
void tao_graph_struct_set_curve_legend_origin(
    void* struct_obj,
    const void* src_ptr);
void tao_graph_struct_get_curve_legend(const void* struct_obj, void** ptr_out);
void tao_graph_struct_set_curve_legend(void* struct_obj, const void* src_ptr);
void tao_graph_struct_get_x(const void* struct_obj, void** ptr_out);
void tao_graph_struct_set_x(void* struct_obj, const void* src_ptr);
void tao_graph_struct_get_y(const void* struct_obj, void** ptr_out);
void tao_graph_struct_set_y(void* struct_obj, const void* src_ptr);
void tao_graph_struct_get_x2(const void* struct_obj, void** ptr_out);
void tao_graph_struct_set_x2(void* struct_obj, const void* src_ptr);
void tao_graph_struct_get_y2(const void* struct_obj, void** ptr_out);
void tao_graph_struct_set_y2(void* struct_obj, const void* src_ptr);
void tao_graph_struct_get_margin(const void* struct_obj, void** ptr_out);
void tao_graph_struct_set_margin(void* struct_obj, const void* src_ptr);
void tao_graph_struct_get_scale_margin(const void* struct_obj, void** ptr_out);
void tao_graph_struct_set_scale_margin(void* struct_obj, const void* src_ptr);
void tao_graph_struct_get_x_axis_scale_factor(
    const void* struct_obj,
    double* value_out);
void tao_graph_struct_set_x_axis_scale_factor(
    void* struct_obj,
    double value_in);
void tao_graph_struct_get_symbol_size_scale(
    const void* struct_obj,
    double* value_out);
void tao_graph_struct_set_symbol_size_scale(void* struct_obj, double value_in);
void tao_graph_struct_get_box_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void tao_graph_struct_get_ix_branch(const void* struct_obj, int* value_out);
void tao_graph_struct_set_ix_branch(void* struct_obj, int value_in);
void tao_graph_struct_get_ix_universe(const void* struct_obj, int* value_out);
void tao_graph_struct_set_ix_universe(void* struct_obj, int value_in);
void tao_graph_struct_get_clip(const void* struct_obj, bool* value_out);
void tao_graph_struct_set_clip(void* struct_obj, bool value_in);
void tao_graph_struct_get_y2_mirrors_y(const void* struct_obj, bool* value_out);
void tao_graph_struct_set_y2_mirrors_y(void* struct_obj, bool value_in);
void tao_graph_struct_get_limited(const void* struct_obj, bool* value_out);
void tao_graph_struct_set_limited(void* struct_obj, bool value_in);
void tao_graph_struct_get_draw_axes(const void* struct_obj, bool* value_out);
void tao_graph_struct_set_draw_axes(void* struct_obj, bool value_in);
void tao_graph_struct_get_draw_curve_legend(
    const void* struct_obj,
    bool* value_out);
void tao_graph_struct_set_draw_curve_legend(void* struct_obj, bool value_in);
void tao_graph_struct_get_draw_grid(const void* struct_obj, bool* value_out);
void tao_graph_struct_set_draw_grid(void* struct_obj, bool value_in);
void tao_graph_struct_get_draw_title(const void* struct_obj, bool* value_out);
void tao_graph_struct_set_draw_title(void* struct_obj, bool value_in);
void tao_graph_struct_get_draw_only_good_user_data_or_vars(
    const void* struct_obj,
    bool* value_out);
void tao_graph_struct_set_draw_only_good_user_data_or_vars(
    void* struct_obj,
    bool value_in);
void tao_graph_struct_get_allow_wrap_around(
    const void* struct_obj,
    bool* value_out);
void tao_graph_struct_set_allow_wrap_around(void* struct_obj, bool value_in);
void tao_graph_struct_get_is_valid(const void* struct_obj, bool* value_out);
void tao_graph_struct_set_is_valid(void* struct_obj, bool value_in);
void tao_plot_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_plot_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_plot_struct_get_description_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_plot_struct_set_description(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void tao_plot_struct_get_graph_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_plot_struct_get_r(const void* struct_obj, void** ptr_out);
void tao_plot_struct_set_r(void* struct_obj, const void* src_ptr);
void tao_plot_struct_get_ix_plot(const void* struct_obj, int* value_out);
void tao_plot_struct_set_ix_plot(void* struct_obj, int value_in);
void tao_plot_struct_get_n_curve_pts(const void* struct_obj, int* value_out);
void tao_plot_struct_set_n_curve_pts(void* struct_obj, int value_in);
void tao_plot_struct_get_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_plot_struct_set_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_plot_struct_get_x_axis_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_plot_struct_set_x_axis_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_plot_struct_get_autoscale_x(const void* struct_obj, bool* value_out);
void tao_plot_struct_set_autoscale_x(void* struct_obj, bool value_in);
void tao_plot_struct_get_autoscale_y(const void* struct_obj, bool* value_out);
void tao_plot_struct_set_autoscale_y(void* struct_obj, bool value_in);
void tao_plot_struct_get_autoscale_gang_x(
    const void* struct_obj,
    bool* value_out);
void tao_plot_struct_set_autoscale_gang_x(void* struct_obj, bool value_in);
void tao_plot_struct_get_autoscale_gang_y(
    const void* struct_obj,
    bool* value_out);
void tao_plot_struct_set_autoscale_gang_y(void* struct_obj, bool value_in);
void tao_plot_struct_get_list_with_show_plot_command(
    const void* struct_obj,
    bool* value_out);
void tao_plot_struct_set_list_with_show_plot_command(
    void* struct_obj,
    bool value_in);
void tao_plot_struct_get_phantom(const void* struct_obj, bool* value_out);
void tao_plot_struct_set_phantom(void* struct_obj, bool value_in);
void tao_plot_struct_get_default_plot(const void* struct_obj, bool* value_out);
void tao_plot_struct_set_default_plot(void* struct_obj, bool value_in);
void tao_plot_region_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_plot_region_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_plot_region_struct_get_plot(const void* struct_obj, void** ptr_out);
void tao_plot_region_struct_set_plot(void* struct_obj, const void* src_ptr);
void tao_plot_region_struct_get_location_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_plot_region_struct_get_visible(
    const void* struct_obj,
    bool* value_out);
void tao_plot_region_struct_set_visible(void* struct_obj, bool value_in);
void tao_plot_region_struct_get_list_with_show_plot_command(
    const void* struct_obj,
    bool* value_out);
void tao_plot_region_struct_set_list_with_show_plot_command(
    void* struct_obj,
    bool value_in);
void tao_plot_region_struct_get_setup_done(
    const void* struct_obj,
    bool* value_out);
void tao_plot_region_struct_set_setup_done(void* struct_obj, bool value_in);
void tao_universe_pointer_struct_get_u(const void* struct_obj, void** ptr_out);
void tao_universe_pointer_struct_set_u(void* struct_obj, const void* src_ptr);
void tao_super_universe_struct_get_global(
    const void* struct_obj,
    void** ptr_out);
void tao_super_universe_struct_set_global(
    void* struct_obj,
    const void* src_ptr);
void tao_super_universe_struct_get_init(const void* struct_obj, void** ptr_out);
void tao_super_universe_struct_set_init(void* struct_obj, const void* src_ptr);
void tao_super_universe_struct_get_com(const void* struct_obj, void** ptr_out);
void tao_super_universe_struct_set_com(void* struct_obj, const void* src_ptr);
void tao_super_universe_struct_get_plot_page(
    const void* struct_obj,
    void** ptr_out);
void tao_super_universe_struct_set_plot_page(
    void* struct_obj,
    const void* src_ptr);

void tao_super_universe_struct_get_v1_var_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_super_universe_struct_get_var_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_super_universe_struct_get_u_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_super_universe_struct_get_key_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void tao_super_universe_struct_get_building_wall(
    const void* struct_obj,
    void** ptr_out);
void tao_super_universe_struct_set_building_wall(
    void* struct_obj,
    const void* src_ptr);
void tao_super_universe_struct_get_wave(const void* struct_obj, void** ptr_out);
void tao_super_universe_struct_set_wave(void* struct_obj, const void* src_ptr);
void tao_super_universe_struct_get_n_var_used(
    const void* struct_obj,
    int* value_out);
void tao_super_universe_struct_set_n_var_used(void* struct_obj, int value_in);
void tao_super_universe_struct_get_n_v1_var_used(
    const void* struct_obj,
    int* value_out);
void tao_super_universe_struct_set_n_v1_var_used(
    void* struct_obj,
    int value_in);

void tao_super_universe_struct_get_history_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_super_universe_struct_get_initialized(
    const void* struct_obj,
    bool* value_out);
void tao_super_universe_struct_set_initialized(void* struct_obj, bool value_in);
void tao_var_struct_get_ele_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_var_struct_set_ele_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_var_struct_get_attrib_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_var_struct_set_attrib_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_var_struct_get_id_info(const void* s, char** d, int* bounds, bool* a);
void tao_var_struct_set_id(void* struct_obj, const char* str_ptr, int str_len);

void tao_var_struct_get_slave_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_var_struct_get_ix_v1(const void* struct_obj, int* value_out);
void tao_var_struct_set_ix_v1(void* struct_obj, int value_in);
void tao_var_struct_get_ix_var(const void* struct_obj, int* value_out);
void tao_var_struct_set_ix_var(void* struct_obj, int value_in);
void tao_var_struct_get_ix_dvar(const void* struct_obj, int* value_out);
void tao_var_struct_set_ix_dvar(void* struct_obj, int value_in);
void tao_var_struct_get_ix_attrib(const void* struct_obj, int* value_out);
void tao_var_struct_set_ix_attrib(void* struct_obj, int value_in);
void tao_var_struct_get_ix_key_table(const void* struct_obj, int* value_out);
void tao_var_struct_set_ix_key_table(void* struct_obj, int value_in);
void tao_var_struct_get_model_value(const void* struct_obj, double** ptr_out);
void tao_var_struct_set_model_value(void* struct_obj, double value_in);
void tao_var_struct_get_base_value(const void* struct_obj, double** ptr_out);
void tao_var_struct_set_base_value(void* struct_obj, double value_in);
void tao_var_struct_get_design_value(const void* struct_obj, double* value_out);
void tao_var_struct_set_design_value(void* struct_obj, double value_in);
void tao_var_struct_get_scratch_value(
    const void* struct_obj,
    double* value_out);
void tao_var_struct_set_scratch_value(void* struct_obj, double value_in);
void tao_var_struct_get_old_value(const void* struct_obj, double* value_out);
void tao_var_struct_set_old_value(void* struct_obj, double value_in);
void tao_var_struct_get_meas_value(const void* struct_obj, double* value_out);
void tao_var_struct_set_meas_value(void* struct_obj, double value_in);
void tao_var_struct_get_ref_value(const void* struct_obj, double* value_out);
void tao_var_struct_set_ref_value(void* struct_obj, double value_in);
void tao_var_struct_get_correction_value(
    const void* struct_obj,
    double* value_out);
void tao_var_struct_set_correction_value(void* struct_obj, double value_in);
void tao_var_struct_get_high_lim(const void* struct_obj, double* value_out);
void tao_var_struct_set_high_lim(void* struct_obj, double value_in);
void tao_var_struct_get_low_lim(const void* struct_obj, double* value_out);
void tao_var_struct_set_low_lim(void* struct_obj, double value_in);
void tao_var_struct_get_step(const void* struct_obj, double* value_out);
void tao_var_struct_set_step(void* struct_obj, double value_in);
void tao_var_struct_get_weight(const void* struct_obj, double* value_out);
void tao_var_struct_set_weight(void* struct_obj, double value_in);
void tao_var_struct_get_delta_merit(const void* struct_obj, double* value_out);
void tao_var_struct_set_delta_merit(void* struct_obj, double value_in);
void tao_var_struct_get_merit(const void* struct_obj, double* value_out);
void tao_var_struct_set_merit(void* struct_obj, double value_in);
void tao_var_struct_get_dMerit_dVar(const void* struct_obj, double* value_out);
void tao_var_struct_set_dMerit_dVar(void* struct_obj, double value_in);
void tao_var_struct_get_key_val0(const void* struct_obj, double* value_out);
void tao_var_struct_set_key_val0(void* struct_obj, double value_in);
void tao_var_struct_get_key_delta(const void* struct_obj, double* value_out);
void tao_var_struct_set_key_delta(void* struct_obj, double value_in);
void tao_var_struct_get_s(const void* struct_obj, double* value_out);
void tao_var_struct_set_s(void* struct_obj, double value_in);
void tao_var_struct_get_extend_val(const void* struct_obj, double* value_out);
void tao_var_struct_set_extend_val(void* struct_obj, double value_in);
void tao_var_struct_get_merit_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_var_struct_set_merit_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_var_struct_get_exists(const void* struct_obj, bool* value_out);
void tao_var_struct_set_exists(void* struct_obj, bool value_in);
void tao_var_struct_get_good_var(const void* struct_obj, bool* value_out);
void tao_var_struct_set_good_var(void* struct_obj, bool value_in);
void tao_var_struct_get_good_user(const void* struct_obj, bool* value_out);
void tao_var_struct_set_good_user(void* struct_obj, bool value_in);
void tao_var_struct_get_good_opt(const void* struct_obj, bool* value_out);
void tao_var_struct_set_good_opt(void* struct_obj, bool value_in);
void tao_var_struct_get_good_plot(const void* struct_obj, bool* value_out);
void tao_var_struct_set_good_plot(void* struct_obj, bool value_in);
void tao_var_struct_get_useit_opt(const void* struct_obj, bool* value_out);
void tao_var_struct_set_useit_opt(void* struct_obj, bool value_in);
void tao_var_struct_get_useit_plot(const void* struct_obj, bool* value_out);
void tao_var_struct_set_useit_plot(void* struct_obj, bool value_in);
void tao_var_struct_get_key_bound(const void* struct_obj, bool* value_out);
void tao_var_struct_set_key_bound(void* struct_obj, bool value_in);
void tao_var_struct_get_v1(const void* struct_obj, void** ptr_out);
void tao_var_struct_set_v1(void* struct_obj, const void* src_ptr);
void tao_var_slave_struct_get_ix_uni(const void* struct_obj, int* value_out);
void tao_var_slave_struct_set_ix_uni(void* struct_obj, int value_in);
void tao_var_slave_struct_get_ix_branch(const void* struct_obj, int* value_out);
void tao_var_slave_struct_set_ix_branch(void* struct_obj, int value_in);
void tao_var_slave_struct_get_ix_ele(const void* struct_obj, int* value_out);
void tao_var_slave_struct_set_ix_ele(void* struct_obj, int value_in);
void tao_var_slave_struct_get_model_value(
    const void* struct_obj,
    double** ptr_out);
void tao_var_slave_struct_set_model_value(void* struct_obj, double value_in);
void tao_var_slave_struct_get_base_value(
    const void* struct_obj,
    double** ptr_out);
void tao_var_slave_struct_set_base_value(void* struct_obj, double value_in);
void tao_lattice_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_lattice_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_lattice_struct_get_lat(const void* struct_obj, void** ptr_out);
void tao_lattice_struct_set_lat(void* struct_obj, const void* src_ptr);
void tao_lattice_struct_get_high_E_lat(const void* struct_obj, void** ptr_out);
void tao_lattice_struct_set_high_E_lat(void* struct_obj, const void* src_ptr);
void tao_lattice_struct_get_low_E_lat(const void* struct_obj, void** ptr_out);
void tao_lattice_struct_set_low_E_lat(void* struct_obj, const void* src_ptr);
void tao_lattice_struct_get_rad_int_by_ele_ri(
    const void* struct_obj,
    void** ptr_out);
void tao_lattice_struct_set_rad_int_by_ele_ri(
    void* struct_obj,
    const void* src_ptr);
void tao_lattice_struct_get_rad_int_by_ele_6d(
    const void* struct_obj,
    void** ptr_out);
void tao_lattice_struct_set_rad_int_by_ele_6d(
    void* struct_obj,
    const void* src_ptr);

void tao_lattice_struct_get_tao_branch_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_beam_uni_struct_get_saved_at_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_beam_uni_struct_set_saved_at(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_beam_uni_struct_get_dump_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_beam_uni_struct_set_dump_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_beam_uni_struct_get_dump_at_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_beam_uni_struct_set_dump_at(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_beam_uni_struct_get_track_beam_in_universe(
    const void* struct_obj,
    bool* value_out);
void tao_beam_uni_struct_set_track_beam_in_universe(
    void* struct_obj,
    bool value_in);
void tao_beam_uni_struct_get_always_reinit(
    const void* struct_obj,
    bool* value_out);
void tao_beam_uni_struct_set_always_reinit(void* struct_obj, bool value_in);
void tao_dynamic_aperture_struct_get_param(
    const void* struct_obj,
    void** ptr_out);
void tao_dynamic_aperture_struct_set_param(
    void* struct_obj,
    const void* src_ptr);

void tao_dynamic_aperture_struct_get_scan_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_dynamic_aperture_struct_get_pz_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_dynamic_aperture_struct_get_ellipse_scale(
    const void* struct_obj,
    double* value_out);
void tao_dynamic_aperture_struct_set_ellipse_scale(
    void* struct_obj,
    double value_in);
void tao_dynamic_aperture_struct_get_a_emit(
    const void* struct_obj,
    double* value_out);
void tao_dynamic_aperture_struct_set_a_emit(void* struct_obj, double value_in);
void tao_dynamic_aperture_struct_get_b_emit(
    const void* struct_obj,
    double* value_out);
void tao_dynamic_aperture_struct_set_b_emit(void* struct_obj, double value_in);

void tao_model_branch_struct_get_ele_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_model_branch_struct_get_beam(const void* struct_obj, void** ptr_out);
void tao_model_branch_struct_set_beam(void* struct_obj, const void* src_ptr);
void tao_spin_map_struct_get_valid(const void* struct_obj, bool* value_out);
void tao_spin_map_struct_set_valid(void* struct_obj, bool value_in);
void tao_spin_map_struct_get_map1(const void* struct_obj, void** ptr_out);
void tao_spin_map_struct_set_map1(void* struct_obj, const void* src_ptr);
void tao_spin_map_struct_get_axis_input(const void* struct_obj, void** ptr_out);
void tao_spin_map_struct_set_axis_input(void* struct_obj, const void* src_ptr);
void tao_spin_map_struct_get_axis0(const void* struct_obj, void** ptr_out);
void tao_spin_map_struct_set_axis0(void* struct_obj, const void* src_ptr);
void tao_spin_map_struct_get_axis1(const void* struct_obj, void** ptr_out);
void tao_spin_map_struct_set_axis1(void* struct_obj, const void* src_ptr);
void tao_spin_map_struct_get_ix_ele(const void* struct_obj, int* value_out);
void tao_spin_map_struct_set_ix_ele(void* struct_obj, int value_in);
void tao_spin_map_struct_get_ix_ref(const void* struct_obj, int* value_out);
void tao_spin_map_struct_set_ix_ref(void* struct_obj, int value_in);
void tao_spin_map_struct_get_ix_uni(const void* struct_obj, int* value_out);
void tao_spin_map_struct_set_ix_uni(void* struct_obj, int value_in);
void tao_spin_map_struct_get_ix_branch(const void* struct_obj, int* value_out);
void tao_spin_map_struct_set_ix_branch(void* struct_obj, int value_in);
void tao_spin_map_struct_get_mat8_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void tao_data_struct_get_ele_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_data_struct_set_ele_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_data_struct_get_ele_start_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_data_struct_set_ele_start_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_data_struct_get_ele_ref_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_data_struct_set_ele_ref_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void tao_data_struct_get_data_type_info(
    const void* s,
    char** d,
    int* len,
    bool* is_alloc);

void tao_data_struct_set_data_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_data_struct_get_merit_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_data_struct_set_merit_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_data_struct_get_id_info(const void* s, char** d, int* bounds, bool* a);
void tao_data_struct_set_id(void* struct_obj, const char* str_ptr, int str_len);
void tao_data_struct_get_data_source_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_data_struct_set_data_source(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_data_struct_get_why_invalid_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_data_struct_set_why_invalid(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_data_struct_get_ix_uni(const void* struct_obj, int* value_out);
void tao_data_struct_set_ix_uni(void* struct_obj, int value_in);
void tao_data_struct_get_ix_bunch(const void* struct_obj, int* value_out);
void tao_data_struct_set_ix_bunch(void* struct_obj, int value_in);
void tao_data_struct_get_ix_branch(const void* struct_obj, int* value_out);
void tao_data_struct_set_ix_branch(void* struct_obj, int value_in);
void tao_data_struct_get_ix_ele(const void* struct_obj, int* value_out);
void tao_data_struct_set_ix_ele(void* struct_obj, int value_in);
void tao_data_struct_get_ix_ele_start(const void* struct_obj, int* value_out);
void tao_data_struct_set_ix_ele_start(void* struct_obj, int value_in);
void tao_data_struct_get_ix_ele_ref(const void* struct_obj, int* value_out);
void tao_data_struct_set_ix_ele_ref(void* struct_obj, int value_in);
void tao_data_struct_get_ix_ele_merit(const void* struct_obj, int* value_out);
void tao_data_struct_set_ix_ele_merit(void* struct_obj, int value_in);
void tao_data_struct_get_ix_d1(const void* struct_obj, int* value_out);
void tao_data_struct_set_ix_d1(void* struct_obj, int value_in);
void tao_data_struct_get_ix_data(const void* struct_obj, int* value_out);
void tao_data_struct_set_ix_data(void* struct_obj, int value_in);
void tao_data_struct_get_ix_dModel(const void* struct_obj, int* value_out);
void tao_data_struct_set_ix_dModel(void* struct_obj, int value_in);
void tao_data_struct_get_eval_point(const void* struct_obj, int* value_out);
void tao_data_struct_set_eval_point(void* struct_obj, int value_in);
void tao_data_struct_get_meas_value(const void* struct_obj, double* value_out);
void tao_data_struct_set_meas_value(void* struct_obj, double value_in);
void tao_data_struct_get_ref_value(const void* struct_obj, double* value_out);
void tao_data_struct_set_ref_value(void* struct_obj, double value_in);
void tao_data_struct_get_model_value(const void* struct_obj, double* value_out);
void tao_data_struct_set_model_value(void* struct_obj, double value_in);
void tao_data_struct_get_design_value(
    const void* struct_obj,
    double* value_out);
void tao_data_struct_set_design_value(void* struct_obj, double value_in);
void tao_data_struct_get_old_value(const void* struct_obj, double* value_out);
void tao_data_struct_set_old_value(void* struct_obj, double value_in);
void tao_data_struct_get_base_value(const void* struct_obj, double* value_out);
void tao_data_struct_set_base_value(void* struct_obj, double value_in);
void tao_data_struct_get_error_rms(const void* struct_obj, double* value_out);
void tao_data_struct_set_error_rms(void* struct_obj, double value_in);
void tao_data_struct_get_delta_merit(const void* struct_obj, double* value_out);
void tao_data_struct_set_delta_merit(void* struct_obj, double value_in);
void tao_data_struct_get_weight(const void* struct_obj, double* value_out);
void tao_data_struct_set_weight(void* struct_obj, double value_in);
void tao_data_struct_get_invalid_value(
    const void* struct_obj,
    double* value_out);
void tao_data_struct_set_invalid_value(void* struct_obj, double value_in);
void tao_data_struct_get_merit(const void* struct_obj, double* value_out);
void tao_data_struct_set_merit(void* struct_obj, double value_in);
void tao_data_struct_get_s(const void* struct_obj, double* value_out);
void tao_data_struct_set_s(void* struct_obj, double value_in);
void tao_data_struct_get_s_offset(const void* struct_obj, double* value_out);
void tao_data_struct_set_s_offset(void* struct_obj, double value_in);
void tao_data_struct_get_ref_s_offset(
    const void* struct_obj,
    double* value_out);
void tao_data_struct_set_ref_s_offset(void* struct_obj, double value_in);
void tao_data_struct_get_err_message_printed(
    const void* struct_obj,
    bool* value_out);
void tao_data_struct_set_err_message_printed(void* struct_obj, bool value_in);
void tao_data_struct_get_exists(const void* struct_obj, bool* value_out);
void tao_data_struct_set_exists(void* struct_obj, bool value_in);
void tao_data_struct_get_good_model(const void* struct_obj, bool* value_out);
void tao_data_struct_set_good_model(void* struct_obj, bool value_in);
void tao_data_struct_get_good_base(const void* struct_obj, bool* value_out);
void tao_data_struct_set_good_base(void* struct_obj, bool value_in);
void tao_data_struct_get_good_design(const void* struct_obj, bool* value_out);
void tao_data_struct_set_good_design(void* struct_obj, bool value_in);
void tao_data_struct_get_good_meas(const void* struct_obj, bool* value_out);
void tao_data_struct_set_good_meas(void* struct_obj, bool value_in);
void tao_data_struct_get_good_ref(const void* struct_obj, bool* value_out);
void tao_data_struct_set_good_ref(void* struct_obj, bool value_in);
void tao_data_struct_get_good_user(const void* struct_obj, bool* value_out);
void tao_data_struct_set_good_user(void* struct_obj, bool value_in);
void tao_data_struct_get_good_opt(const void* struct_obj, bool* value_out);
void tao_data_struct_set_good_opt(void* struct_obj, bool value_in);
void tao_data_struct_get_good_plot(const void* struct_obj, bool* value_out);
void tao_data_struct_set_good_plot(void* struct_obj, bool value_in);
void tao_data_struct_get_useit_plot(const void* struct_obj, bool* value_out);
void tao_data_struct_set_useit_plot(void* struct_obj, bool value_in);
void tao_data_struct_get_useit_opt(const void* struct_obj, bool* value_out);
void tao_data_struct_set_useit_opt(void* struct_obj, bool value_in);
void tao_data_struct_get_spin_map(const void* struct_obj, void** ptr_out);
void tao_data_struct_set_spin_map(void* struct_obj, const void* src_ptr);
void tao_data_struct_get_d1(const void* struct_obj, void** ptr_out);
void tao_data_struct_set_d1(void* struct_obj, const void* src_ptr);
void tao_ping_scale_struct_get_a_mode_meas(
    const void* struct_obj,
    double* value_out);
void tao_ping_scale_struct_set_a_mode_meas(void* struct_obj, double value_in);
void tao_ping_scale_struct_get_a_mode_ref(
    const void* struct_obj,
    double* value_out);
void tao_ping_scale_struct_set_a_mode_ref(void* struct_obj, double value_in);
void tao_ping_scale_struct_get_b_mode_meas(
    const void* struct_obj,
    double* value_out);
void tao_ping_scale_struct_set_b_mode_meas(void* struct_obj, double value_in);
void tao_ping_scale_struct_get_b_mode_ref(
    const void* struct_obj,
    double* value_out);
void tao_ping_scale_struct_set_b_mode_ref(void* struct_obj, double value_in);
void tao_universe_calc_struct_get_srdt_for_data(
    const void* struct_obj,
    int* value_out);
void tao_universe_calc_struct_set_srdt_for_data(void* struct_obj, int value_in);
void tao_universe_calc_struct_get_rad_int_for_data(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_rad_int_for_data(
    void* struct_obj,
    bool value_in);
void tao_universe_calc_struct_get_rad_int_for_plotting(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_rad_int_for_plotting(
    void* struct_obj,
    bool value_in);
void tao_universe_calc_struct_get_chrom_for_data(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_chrom_for_data(
    void* struct_obj,
    bool value_in);
void tao_universe_calc_struct_get_chrom_for_plotting(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_chrom_for_plotting(
    void* struct_obj,
    bool value_in);
void tao_universe_calc_struct_get_lat_sigma_for_data(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_lat_sigma_for_data(
    void* struct_obj,
    bool value_in);
void tao_universe_calc_struct_get_lat_sigma_for_plotting(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_lat_sigma_for_plotting(
    void* struct_obj,
    bool value_in);
void tao_universe_calc_struct_get_dynamic_aperture(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_dynamic_aperture(
    void* struct_obj,
    bool value_in);
void tao_universe_calc_struct_get_one_turn_map(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_one_turn_map(void* struct_obj, bool value_in);
void tao_universe_calc_struct_get_lattice(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_lattice(void* struct_obj, bool value_in);
void tao_universe_calc_struct_get_twiss(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_twiss(void* struct_obj, bool value_in);
void tao_universe_calc_struct_get_track(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_track(void* struct_obj, bool value_in);
void tao_universe_calc_struct_get_spin_matrices(
    const void* struct_obj,
    bool* value_out);
void tao_universe_calc_struct_set_spin_matrices(
    void* struct_obj,
    bool value_in);

void lat_ele_order_struct_get_branch_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_title_struct_get_string_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_title_struct_set_string(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_title_struct_get_x(const void* struct_obj, double* value_out);
void tao_title_struct_set_x(void* struct_obj, double value_in);
void tao_title_struct_get_y(const void* struct_obj, double* value_out);
void tao_title_struct_set_y(void* struct_obj, double value_in);
void tao_title_struct_get_units_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_title_struct_set_units(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_title_struct_get_justify_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_title_struct_set_justify(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_title_struct_get_draw_it(const void* struct_obj, bool* value_out);
void tao_title_struct_set_draw_it(void* struct_obj, bool value_in);
void qp_rect_struct_get_x1(const void* struct_obj, double* value_out);
void qp_rect_struct_set_x1(void* struct_obj, double value_in);
void qp_rect_struct_get_x2(const void* struct_obj, double* value_out);
void qp_rect_struct_set_x2(void* struct_obj, double value_in);
void qp_rect_struct_get_y1(const void* struct_obj, double* value_out);
void qp_rect_struct_set_y1(void* struct_obj, double value_in);
void qp_rect_struct_get_y2(const void* struct_obj, double* value_out);
void qp_rect_struct_set_y2(void* struct_obj, double value_in);
void qp_rect_struct_get_units_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_rect_struct_set_units(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void tao_drawing_struct_get_ele_shape_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_shape_pattern_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_shape_pattern_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_shape_pattern_struct_get_line(const void* struct_obj, void** ptr_out);
void tao_shape_pattern_struct_set_line(void* struct_obj, const void* src_ptr);

void tao_shape_pattern_struct_get_pt_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_shape_pattern_point_struct_get_s(
    const void* struct_obj,
    double* value_out);
void tao_shape_pattern_point_struct_set_s(void* struct_obj, double value_in);
void tao_shape_pattern_point_struct_get_y(
    const void* struct_obj,
    double* value_out);
void tao_shape_pattern_point_struct_set_y(void* struct_obj, double value_in);
void tao_shape_pattern_point_struct_get_radius(
    const void* struct_obj,
    double* value_out);
void tao_shape_pattern_point_struct_set_radius(
    void* struct_obj,
    double value_in);
void qp_axis_struct_get_label_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_axis_struct_set_label(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void qp_axis_struct_get_min(const void* struct_obj, double* value_out);
void qp_axis_struct_set_min(void* struct_obj, double value_in);
void qp_axis_struct_get_max(const void* struct_obj, double* value_out);
void qp_axis_struct_set_max(void* struct_obj, double value_in);
void qp_axis_struct_get_tick_min(const void* struct_obj, double* value_out);
void qp_axis_struct_set_tick_min(void* struct_obj, double value_in);
void qp_axis_struct_get_tick_max(const void* struct_obj, double* value_out);
void qp_axis_struct_set_tick_max(void* struct_obj, double value_in);
void qp_axis_struct_get_eval_min(const void* struct_obj, double* value_out);
void qp_axis_struct_set_eval_min(void* struct_obj, double value_in);
void qp_axis_struct_get_eval_max(const void* struct_obj, double* value_out);
void qp_axis_struct_set_eval_max(void* struct_obj, double value_in);
void qp_axis_struct_get_dtick(const void* struct_obj, double* value_out);
void qp_axis_struct_set_dtick(void* struct_obj, double value_in);
void qp_axis_struct_get_number_offset(
    const void* struct_obj,
    double* value_out);
void qp_axis_struct_set_number_offset(void* struct_obj, double value_in);
void qp_axis_struct_get_label_offset(const void* struct_obj, double* value_out);
void qp_axis_struct_set_label_offset(void* struct_obj, double value_in);
void qp_axis_struct_get_major_tick_len(
    const void* struct_obj,
    double* value_out);
void qp_axis_struct_set_major_tick_len(void* struct_obj, double value_in);
void qp_axis_struct_get_minor_tick_len(
    const void* struct_obj,
    double* value_out);
void qp_axis_struct_set_minor_tick_len(void* struct_obj, double value_in);
void qp_axis_struct_get_label_color_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_axis_struct_set_label_color(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void qp_axis_struct_get_major_div(const void* struct_obj, int* value_out);
void qp_axis_struct_set_major_div(void* struct_obj, int value_in);
void qp_axis_struct_get_major_div_nominal(
    const void* struct_obj,
    int* value_out);
void qp_axis_struct_set_major_div_nominal(void* struct_obj, int value_in);
void qp_axis_struct_get_minor_div(const void* struct_obj, int* value_out);
void qp_axis_struct_set_minor_div(void* struct_obj, int value_in);
void qp_axis_struct_get_minor_div_max(const void* struct_obj, int* value_out);
void qp_axis_struct_set_minor_div_max(void* struct_obj, int value_in);
void qp_axis_struct_get_places(const void* struct_obj, int* value_out);
void qp_axis_struct_set_places(void* struct_obj, int value_in);
void qp_axis_struct_get_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_axis_struct_set_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void qp_axis_struct_get_bounds_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_axis_struct_set_bounds(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void qp_axis_struct_get_tick_side(const void* struct_obj, int* value_out);
void qp_axis_struct_set_tick_side(void* struct_obj, int value_in);
void qp_axis_struct_get_number_side(const void* struct_obj, int* value_out);
void qp_axis_struct_set_number_side(void* struct_obj, int value_in);
void qp_axis_struct_get_draw_label(const void* struct_obj, bool* value_out);
void qp_axis_struct_set_draw_label(void* struct_obj, bool value_in);
void qp_axis_struct_get_draw_numbers(const void* struct_obj, bool* value_out);
void qp_axis_struct_set_draw_numbers(void* struct_obj, bool value_in);
void qp_legend_struct_get_row_spacing(
    const void* struct_obj,
    double* value_out);
void qp_legend_struct_set_row_spacing(void* struct_obj, double value_in);
void qp_legend_struct_get_line_length(
    const void* struct_obj,
    double* value_out);
void qp_legend_struct_set_line_length(void* struct_obj, double value_in);
void qp_legend_struct_get_text_offset(
    const void* struct_obj,
    double* value_out);
void qp_legend_struct_set_text_offset(void* struct_obj, double value_in);
void qp_legend_struct_get_draw_line(const void* struct_obj, bool* value_out);
void qp_legend_struct_set_draw_line(void* struct_obj, bool value_in);
void qp_legend_struct_get_draw_symbol(const void* struct_obj, bool* value_out);
void qp_legend_struct_set_draw_symbol(void* struct_obj, bool value_in);
void qp_legend_struct_get_draw_text(const void* struct_obj, bool* value_out);
void qp_legend_struct_set_draw_text(void* struct_obj, bool value_in);
void qp_point_struct_get_x(const void* struct_obj, double* value_out);
void qp_point_struct_set_x(void* struct_obj, double value_in);
void qp_point_struct_get_y(const void* struct_obj, double* value_out);
void qp_point_struct_set_y(void* struct_obj, double value_in);
void qp_point_struct_get_units_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_point_struct_set_units(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void qp_line_struct_get_width(const void* struct_obj, int* value_out);
void qp_line_struct_set_width(void* struct_obj, int value_in);
void qp_line_struct_get_color_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_line_struct_set_color(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void qp_line_struct_get_pattern_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_line_struct_set_pattern(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void qp_symbol_struct_get_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_symbol_struct_set_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void qp_symbol_struct_get_height(const void* struct_obj, double* value_out);
void qp_symbol_struct_set_height(void* struct_obj, double value_in);
void qp_symbol_struct_get_color_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_symbol_struct_set_color(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void qp_symbol_struct_get_fill_pattern_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void qp_symbol_struct_set_fill_pattern(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void qp_symbol_struct_get_line_width(const void* struct_obj, int* value_out);
void qp_symbol_struct_set_line_width(void* struct_obj, int value_in);
void tao_floor_plan_struct_get_view_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_floor_plan_struct_set_view(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_floor_plan_struct_get_rotation(
    const void* struct_obj,
    double* value_out);
void tao_floor_plan_struct_set_rotation(void* struct_obj, double value_in);
void tao_floor_plan_struct_get_correct_distortion(
    const void* struct_obj,
    bool* value_out);
void tao_floor_plan_struct_set_correct_distortion(
    void* struct_obj,
    bool value_in);
void tao_floor_plan_struct_get_flip_label_side(
    const void* struct_obj,
    bool* value_out);
void tao_floor_plan_struct_set_flip_label_side(void* struct_obj, bool value_in);
void tao_floor_plan_struct_get_size_is_absolute(
    const void* struct_obj,
    bool* value_out);
void tao_floor_plan_struct_set_size_is_absolute(
    void* struct_obj,
    bool value_in);
void tao_floor_plan_struct_get_draw_only_first_pass(
    const void* struct_obj,
    bool* value_out);
void tao_floor_plan_struct_set_draw_only_first_pass(
    void* struct_obj,
    bool value_in);
void tao_floor_plan_struct_get_draw_building_wall(
    const void* struct_obj,
    bool* value_out);
void tao_floor_plan_struct_set_draw_building_wall(
    void* struct_obj,
    bool value_in);
void tao_floor_plan_struct_get_orbit_scale(
    const void* struct_obj,
    double* value_out);
void tao_floor_plan_struct_set_orbit_scale(void* struct_obj, double value_in);
void tao_floor_plan_struct_get_orbit_color_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_floor_plan_struct_set_orbit_color(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_floor_plan_struct_get_orbit_pattern_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_floor_plan_struct_set_orbit_pattern(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_floor_plan_struct_get_orbit_lattice_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_floor_plan_struct_set_orbit_lattice(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_floor_plan_struct_get_orbit_width(
    const void* struct_obj,
    int* value_out);
void tao_floor_plan_struct_set_orbit_width(void* struct_obj, int value_in);
void tao_v1_var_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_v1_var_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_v1_var_struct_get_ix_v1_var(const void* struct_obj, int* value_out);
void tao_v1_var_struct_set_ix_v1_var(void* struct_obj, int value_in);

void tao_v1_var_struct_get_v_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_global_struct_get_beam_dead_cutoff(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_beam_dead_cutoff(void* struct_obj, double value_in);
void tao_global_struct_get_lm_opt_deriv_reinit(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_lm_opt_deriv_reinit(
    void* struct_obj,
    double value_in);
void tao_global_struct_get_de_lm_step_ratio(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_de_lm_step_ratio(void* struct_obj, double value_in);
void tao_global_struct_get_de_var_to_population_factor(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_de_var_to_population_factor(
    void* struct_obj,
    double value_in);
void tao_global_struct_get_lmdif_eps(const void* struct_obj, double* value_out);
void tao_global_struct_set_lmdif_eps(void* struct_obj, double value_in);
void tao_global_struct_get_lmdif_negligible_merit(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_lmdif_negligible_merit(
    void* struct_obj,
    double value_in);
void tao_global_struct_get_svd_cutoff(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_svd_cutoff(void* struct_obj, double value_in);
void tao_global_struct_get_unstable_penalty(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_unstable_penalty(void* struct_obj, double value_in);
void tao_global_struct_get_merit_stop_value(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_merit_stop_value(void* struct_obj, double value_in);
void tao_global_struct_get_dmerit_stop_value(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_dmerit_stop_value(void* struct_obj, double value_in);
void tao_global_struct_get_random_sigma_cutoff(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_random_sigma_cutoff(
    void* struct_obj,
    double value_in);
void tao_global_struct_get_delta_e_chrom(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_delta_e_chrom(void* struct_obj, double value_in);
void tao_global_struct_get_max_plot_time(
    const void* struct_obj,
    double* value_out);
void tao_global_struct_set_max_plot_time(void* struct_obj, double value_in);
void tao_global_struct_get_default_universe(
    const void* struct_obj,
    int* value_out);
void tao_global_struct_set_default_universe(void* struct_obj, int value_in);
void tao_global_struct_get_default_branch(
    const void* struct_obj,
    int* value_out);
void tao_global_struct_set_default_branch(void* struct_obj, int value_in);
void tao_global_struct_get_n_opti_cycles(
    const void* struct_obj,
    int* value_out);
void tao_global_struct_set_n_opti_cycles(void* struct_obj, int value_in);
void tao_global_struct_get_n_opti_loops(const void* struct_obj, int* value_out);
void tao_global_struct_set_n_opti_loops(void* struct_obj, int value_in);
void tao_global_struct_get_n_threads(const void* struct_obj, int* value_out);
void tao_global_struct_set_n_threads(void* struct_obj, int value_in);
void tao_global_struct_get_phase_units(const void* struct_obj, int* value_out);
void tao_global_struct_set_phase_units(void* struct_obj, int value_in);
void tao_global_struct_get_bunch_to_plot(
    const void* struct_obj,
    int* value_out);
void tao_global_struct_set_bunch_to_plot(void* struct_obj, int value_in);
void tao_global_struct_get_random_seed(const void* struct_obj, int* value_out);
void tao_global_struct_set_random_seed(void* struct_obj, int value_in);
void tao_global_struct_get_n_top10_merit(
    const void* struct_obj,
    int* value_out);
void tao_global_struct_set_n_top10_merit(void* struct_obj, int value_in);
void tao_global_struct_get_srdt_gen_n_slices(
    const void* struct_obj,
    int* value_out);
void tao_global_struct_set_srdt_gen_n_slices(void* struct_obj, int value_in);
void tao_global_struct_get_datum_err_messages_max(
    const void* struct_obj,
    int* value_out);
void tao_global_struct_set_datum_err_messages_max(
    void* struct_obj,
    int value_in);
void tao_global_struct_get_srdt_sxt_n_slices(
    const void* struct_obj,
    int* value_out);
void tao_global_struct_set_srdt_sxt_n_slices(void* struct_obj, int value_in);
void tao_global_struct_get_srdt_use_cache(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_srdt_use_cache(void* struct_obj, bool value_in);
void tao_global_struct_get_quiet_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_quiet(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_random_engine_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_random_engine(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_random_gauss_converter_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_random_gauss_converter(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_track_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_track_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_lat_sigma_calc_uses_emit_from_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_lat_sigma_calc_uses_emit_from(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_prompt_string_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_prompt_string(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_prompt_color_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_prompt_color(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_optimizer_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_optimizer(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_print_command_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_print_command(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_var_out_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_var_out_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_history_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_global_struct_set_history_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_global_struct_get_beam_timer_on(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_beam_timer_on(void* struct_obj, bool value_in);
void tao_global_struct_get_box_plots(const void* struct_obj, bool* value_out);
void tao_global_struct_set_box_plots(void* struct_obj, bool value_in);
void tao_global_struct_get_blank_line_between_commands(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_blank_line_between_commands(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_cmd_file_abort_on_error(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_cmd_file_abort_on_error(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_concatenate_maps(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_concatenate_maps(void* struct_obj, bool value_in);
void tao_global_struct_get_derivative_recalc(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_derivative_recalc(void* struct_obj, bool value_in);
void tao_global_struct_get_derivative_uses_design(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_derivative_uses_design(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_disable_smooth_line_calc(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_disable_smooth_line_calc(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_draw_curve_off_scale_warn(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_draw_curve_off_scale_warn(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_external_plotting(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_external_plotting(void* struct_obj, bool value_in);
void tao_global_struct_get_label_lattice_elements(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_label_lattice_elements(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_label_keys(const void* struct_obj, bool* value_out);
void tao_global_struct_set_label_keys(void* struct_obj, bool value_in);
void tao_global_struct_get_lattice_calc_on(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_lattice_calc_on(void* struct_obj, bool value_in);
void tao_global_struct_get_only_limit_opt_vars(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_only_limit_opt_vars(void* struct_obj, bool value_in);
void tao_global_struct_get_opt_with_ref(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_opt_with_ref(void* struct_obj, bool value_in);
void tao_global_struct_get_opt_with_base(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_opt_with_base(void* struct_obj, bool value_in);
void tao_global_struct_get_opt_match_auto_recalc(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_opt_match_auto_recalc(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_opti_write_var_file(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_opti_write_var_file(void* struct_obj, bool value_in);
void tao_global_struct_get_optimizer_allow_user_abort(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_optimizer_allow_user_abort(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_optimizer_var_limit_warn(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_optimizer_var_limit_warn(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_plot_on(const void* struct_obj, bool* value_out);
void tao_global_struct_set_plot_on(void* struct_obj, bool value_in);
void tao_global_struct_get_rad_int_user_calc_on(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_rad_int_user_calc_on(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_rf_on(const void* struct_obj, bool* value_out);
void tao_global_struct_set_rf_on(void* struct_obj, bool value_in);
void tao_global_struct_get_single_step(const void* struct_obj, bool* value_out);
void tao_global_struct_set_single_step(void* struct_obj, bool value_in);
void tao_global_struct_get_stop_on_error(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_stop_on_error(void* struct_obj, bool value_in);
void tao_global_struct_get_svd_retreat_on_merit_increase(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_svd_retreat_on_merit_increase(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_var_limits_on(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_var_limits_on(void* struct_obj, bool value_in);
void tao_global_struct_get_wait_for_CR_in_single_mode(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_wait_for_CR_in_single_mode(
    void* struct_obj,
    bool value_in);
void tao_global_struct_get_symbol_import(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_symbol_import(void* struct_obj, bool value_in);
void tao_global_struct_get_debug_on(const void* struct_obj, bool* value_out);
void tao_global_struct_set_debug_on(void* struct_obj, bool value_in);
void tao_global_struct_get_expression_tree_on(
    const void* struct_obj,
    bool* value_out);
void tao_global_struct_set_expression_tree_on(void* struct_obj, bool value_in);
void tao_global_struct_get_verbose_on(const void* struct_obj, bool* value_out);
void tao_global_struct_set_verbose_on(void* struct_obj, bool value_in);
void tao_init_struct_get_parse_cmd_args(
    const void* struct_obj,
    bool* value_out);
void tao_init_struct_set_parse_cmd_args(void* struct_obj, bool value_in);
void tao_init_struct_get_debug_switch(const void* struct_obj, bool* value_out);
void tao_init_struct_set_debug_switch(void* struct_obj, bool value_in);
void tao_init_struct_get_external_plotting_switch(
    const void* struct_obj,
    bool* value_out);
void tao_init_struct_set_external_plotting_switch(
    void* struct_obj,
    bool value_in);
void tao_init_struct_get_init_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_init_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_hook_init_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_hook_init_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_hook_lat_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_hook_lat_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_hook_beam_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_hook_beam_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_hook_data_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_hook_data_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_hook_plot_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_hook_plot_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_hook_startup_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_hook_startup_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_hook_var_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_hook_var_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_hook_building_wall_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_hook_building_wall_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_init_file_arg_path_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_init_file_arg_path(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_lattice_file_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_lattice_file_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_hook_init_file_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_hook_init_file_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_init_file_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_init_file_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_beam_file_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_beam_file_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_beam_init_position_file_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_beam_init_position_file_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_command_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_command_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_data_file_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_data_file_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_plot_file_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_plot_file_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_startup_file_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_startup_file_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_var_file_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_var_file_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_building_wall_file_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_building_wall_file_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_geometry_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_geometry_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_slice_lattice_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_slice_lattice_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_start_branch_at_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_start_branch_at_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_log_startup_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_log_startup_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_no_stopping_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_no_stopping_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_noplot_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_noplot_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_no_rad_int_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_no_rad_int_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_reverse_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_reverse_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_debug_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_debug_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_disable_smooth_line_calc_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_disable_smooth_line_calc_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_rf_on_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_rf_on_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_prompt_color_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_prompt_color_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_quiet_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_quiet_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_noinit_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_noinit_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_nostartup_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_nostartup_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_symbol_import_arg_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_symbol_import_arg(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_init_struct_get_unique_name_suffix_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_init_struct_set_unique_name_suffix(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void tao_common_struct_get_plot_place_buffer_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_common_struct_get_covar_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void tao_common_struct_get_alpha_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void tao_common_struct_get_dummy_target(
    const void* struct_obj,
    double* value_out);
void tao_common_struct_set_dummy_target(void* struct_obj, double value_in);
void tao_common_struct_get_n_alias(const void* struct_obj, int* value_out);
void tao_common_struct_set_n_alias(void* struct_obj, int value_in);
void tao_common_struct_get_cmd_file_level(
    const void* struct_obj,
    int* value_out);
void tao_common_struct_set_cmd_file_level(void* struct_obj, int value_in);
void tao_common_struct_get_ix_key_bank(const void* struct_obj, int* value_out);
void tao_common_struct_set_ix_key_bank(void* struct_obj, int value_in);
void tao_common_struct_get_ix_history(const void* struct_obj, int* value_out);
void tao_common_struct_set_ix_history(void* struct_obj, int value_in);
void tao_common_struct_get_n_history(const void* struct_obj, int* value_out);
void tao_common_struct_set_n_history(void* struct_obj, int value_in);
void tao_common_struct_get_lev_loop(const void* struct_obj, int* value_out);
void tao_common_struct_set_lev_loop(void* struct_obj, int value_in);
void tao_common_struct_get_n_err_messages_printed(
    const void* struct_obj,
    int* value_out);
void tao_common_struct_set_n_err_messages_printed(
    void* struct_obj,
    int value_in);
void tao_common_struct_get_n_universes(const void* struct_obj, int* value_out);
void tao_common_struct_set_n_universes(void* struct_obj, int value_in);
void tao_common_struct_get_ix_beam_track_active_element(
    const void* struct_obj,
    int* value_out);
void tao_common_struct_set_ix_beam_track_active_element(
    void* struct_obj,
    int value_in);
void tao_common_struct_get_cmd_file_paused(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_cmd_file_paused(void* struct_obj, bool value_in);
void tao_common_struct_get_use_cmd_here(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_use_cmd_here(void* struct_obj, bool value_in);
void tao_common_struct_get_cmd_from_cmd_file(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_cmd_from_cmd_file(void* struct_obj, bool value_in);
void tao_common_struct_get_use_saved_beam_in_tracking(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_use_saved_beam_in_tracking(
    void* struct_obj,
    bool value_in);
void tao_common_struct_get_single_mode(const void* struct_obj, bool* value_out);
void tao_common_struct_set_single_mode(void* struct_obj, bool value_in);
void tao_common_struct_get_combine_consecutive_elements_of_like_name(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_combine_consecutive_elements_of_like_name(
    void* struct_obj,
    bool value_in);
void tao_common_struct_get_have_tracked_beam(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_have_tracked_beam(void* struct_obj, bool value_in);
void tao_common_struct_get_init_plot_needed(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_init_plot_needed(void* struct_obj, bool value_in);
void tao_common_struct_get_init_beam(const void* struct_obj, bool* value_out);
void tao_common_struct_set_init_beam(void* struct_obj, bool value_in);
void tao_common_struct_get_init_var(const void* struct_obj, bool* value_out);
void tao_common_struct_set_init_var(void* struct_obj, bool value_in);
void tao_common_struct_get_init_read_lat_info(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_init_read_lat_info(void* struct_obj, bool value_in);
void tao_common_struct_get_optimizer_running(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_optimizer_running(void* struct_obj, bool value_in);
void tao_common_struct_get_have_datums_using_expressions(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_have_datums_using_expressions(
    void* struct_obj,
    bool value_in);
void tao_common_struct_get_print_to_terminal(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_print_to_terminal(void* struct_obj, bool value_in);
void tao_common_struct_get_lattice_calc_done(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_lattice_calc_done(void* struct_obj, bool value_in);
void tao_common_struct_get_add_measurement_noise(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_add_measurement_noise(
    void* struct_obj,
    bool value_in);
void tao_common_struct_get_command_arg_has_been_executed(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_command_arg_has_been_executed(
    void* struct_obj,
    bool value_in);
void tao_common_struct_get_all_merit_weights_positive(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_all_merit_weights_positive(
    void* struct_obj,
    bool value_in);
void tao_common_struct_get_multi_turn_orbit_is_plotted(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_multi_turn_orbit_is_plotted(
    void* struct_obj,
    bool value_in);
void tao_common_struct_get_force_chrom_calc(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_force_chrom_calc(void* struct_obj, bool value_in);
void tao_common_struct_get_force_rad_int_calc(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_force_rad_int_calc(void* struct_obj, bool value_in);
void tao_common_struct_get_rad_int_ri_calc_on(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_rad_int_ri_calc_on(void* struct_obj, bool value_in);
void tao_common_struct_get_rad_int_6d_calc_on(
    const void* struct_obj,
    bool* value_out);
void tao_common_struct_set_rad_int_6d_calc_on(void* struct_obj, bool value_in);

void tao_common_struct_get_valid_plot_who_info(
    const void* s,
    char** d,
    int* bounds, // [lower, upper]
    int* str_len,
    bool* is_alloc);

void tao_common_struct_get_single_mode_buffer_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_common_struct_set_single_mode_buffer(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_common_struct_get_cmd_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_common_struct_set_cmd(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_common_struct_get_saved_cmd_line_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_common_struct_set_saved_cmd_line(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_plot_page_struct_get_title(const void* struct_obj, void** ptr_out);
void tao_plot_page_struct_set_title(void* struct_obj, const void* src_ptr);
void tao_plot_page_struct_get_subtitle(const void* struct_obj, void** ptr_out);
void tao_plot_page_struct_set_subtitle(void* struct_obj, const void* src_ptr);
void tao_plot_page_struct_get_border(const void* struct_obj, void** ptr_out);
void tao_plot_page_struct_set_border(void* struct_obj, const void* src_ptr);
void tao_plot_page_struct_get_floor_plan(
    const void* struct_obj,
    void** ptr_out);
void tao_plot_page_struct_set_floor_plan(void* struct_obj, const void* src_ptr);
void tao_plot_page_struct_get_lat_layout(
    const void* struct_obj,
    void** ptr_out);
void tao_plot_page_struct_set_lat_layout(void* struct_obj, const void* src_ptr);

void tao_plot_page_struct_get_pattern_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_plot_page_struct_get_template_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_plot_page_struct_get_region_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_plot_page_struct_get_plot_display_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_plot_page_struct_set_plot_display_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_plot_page_struct_get_size_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_plot_page_struct_get_text_height(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_text_height(void* struct_obj, double value_in);
void tao_plot_page_struct_get_main_title_text_scale(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_main_title_text_scale(
    void* struct_obj,
    double value_in);
void tao_plot_page_struct_get_graph_title_text_scale(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_graph_title_text_scale(
    void* struct_obj,
    double value_in);
void tao_plot_page_struct_get_axis_number_text_scale(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_axis_number_text_scale(
    void* struct_obj,
    double value_in);
void tao_plot_page_struct_get_axis_label_text_scale(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_axis_label_text_scale(
    void* struct_obj,
    double value_in);
void tao_plot_page_struct_get_legend_text_scale(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_legend_text_scale(
    void* struct_obj,
    double value_in);
void tao_plot_page_struct_get_key_table_text_scale(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_key_table_text_scale(
    void* struct_obj,
    double value_in);
void tao_plot_page_struct_get_floor_plan_shape_scale(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_floor_plan_shape_scale(
    void* struct_obj,
    double value_in);
void tao_plot_page_struct_get_floor_plan_text_scale(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_floor_plan_text_scale(
    void* struct_obj,
    double value_in);
void tao_plot_page_struct_get_lat_layout_shape_scale(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_lat_layout_shape_scale(
    void* struct_obj,
    double value_in);
void tao_plot_page_struct_get_lat_layout_text_scale(
    const void* struct_obj,
    double* value_out);
void tao_plot_page_struct_set_lat_layout_text_scale(
    void* struct_obj,
    double value_in);
void tao_plot_page_struct_get_n_curve_pts(
    const void* struct_obj,
    int* value_out);
void tao_plot_page_struct_set_n_curve_pts(void* struct_obj, int value_in);
void tao_plot_page_struct_get_id_window(const void* struct_obj, int* value_out);
void tao_plot_page_struct_set_id_window(void* struct_obj, int value_in);
void tao_plot_page_struct_get_delete_overlapping_plots(
    const void* struct_obj,
    bool* value_out);
void tao_plot_page_struct_set_delete_overlapping_plots(
    void* struct_obj,
    bool value_in);
void tao_plot_page_struct_get_draw_graph_title_suffix(
    const void* struct_obj,
    bool* value_out);
void tao_plot_page_struct_set_draw_graph_title_suffix(
    void* struct_obj,
    bool value_in);
void tao_building_wall_struct_get_orientation(
    const void* struct_obj,
    void** ptr_out);
void tao_building_wall_struct_set_orientation(
    void* struct_obj,
    const void* src_ptr);

void tao_building_wall_struct_get_section_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_building_wall_orientation_struct_get_theta(
    const void* struct_obj,
    double* value_out);
void tao_building_wall_orientation_struct_set_theta(
    void* struct_obj,
    double value_in);
void tao_building_wall_orientation_struct_get_x_offset(
    const void* struct_obj,
    double* value_out);
void tao_building_wall_orientation_struct_set_x_offset(
    void* struct_obj,
    double value_in);
void tao_building_wall_orientation_struct_get_z_offset(
    const void* struct_obj,
    double* value_out);
void tao_building_wall_orientation_struct_set_z_offset(
    void* struct_obj,
    double value_in);
void tao_building_wall_section_struct_get_name_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_building_wall_section_struct_set_name(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_building_wall_section_struct_get_constraint_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_building_wall_section_struct_set_constraint(
    void* struct_obj,
    const char* str_ptr,
    int str_len);

void tao_building_wall_section_struct_get_point_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_building_wall_point_struct_get_z(
    const void* struct_obj,
    double* value_out);
void tao_building_wall_point_struct_set_z(void* struct_obj, double value_in);
void tao_building_wall_point_struct_get_x(
    const void* struct_obj,
    double* value_out);
void tao_building_wall_point_struct_set_x(void* struct_obj, double value_in);
void tao_building_wall_point_struct_get_radius(
    const void* struct_obj,
    double* value_out);
void tao_building_wall_point_struct_set_radius(
    void* struct_obj,
    double value_in);
void tao_building_wall_point_struct_get_z_center(
    const void* struct_obj,
    double* value_out);
void tao_building_wall_point_struct_set_z_center(
    void* struct_obj,
    double value_in);
void tao_building_wall_point_struct_get_x_center(
    const void* struct_obj,
    double* value_out);
void tao_building_wall_point_struct_set_x_center(
    void* struct_obj,
    double value_in);
void tao_wave_struct_get_data_type_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void tao_wave_struct_set_data_type(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_wave_struct_get_rms_rel_a(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_rel_a(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_b(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_rel_b(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_as(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_rel_as(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_bs(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_rel_bs(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_ar(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_rel_ar(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_br(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_rel_br(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_k(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_rel_k(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_ks(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_rel_ks(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_kr(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_rel_kr(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_phi(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_phi(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_phi_s(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_phi_s(void* struct_obj, double value_in);
void tao_wave_struct_get_rms_phi_r(const void* struct_obj, double* value_out);
void tao_wave_struct_set_rms_phi_r(void* struct_obj, double value_in);
void tao_wave_struct_get_amp_ba_s(const void* struct_obj, double* value_out);
void tao_wave_struct_set_amp_ba_s(void* struct_obj, double value_in);
void tao_wave_struct_get_amp_ba_r(const void* struct_obj, double* value_out);
void tao_wave_struct_set_amp_ba_r(void* struct_obj, double value_in);
void tao_wave_struct_get_chi_a(const void* struct_obj, double* value_out);
void tao_wave_struct_set_chi_a(void* struct_obj, double value_in);
void tao_wave_struct_get_chi_c(const void* struct_obj, double* value_out);
void tao_wave_struct_set_chi_c(void* struct_obj, double value_in);
void tao_wave_struct_get_chi_ba(const void* struct_obj, double* value_out);
void tao_wave_struct_set_chi_ba(void* struct_obj, double value_in);
void tao_wave_struct_get_amp_a_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_wave_struct_get_amp_b_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_wave_struct_get_amp_ba_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_wave_struct_get_coef_a_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_wave_struct_get_coef_b_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_wave_struct_get_coef_ba_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void tao_wave_struct_get_n_func(const void* struct_obj, int* value_out);
void tao_wave_struct_set_n_func(void* struct_obj, int value_in);
void tao_wave_struct_get_ix_a1(const void* struct_obj, int* value_out);
void tao_wave_struct_set_ix_a1(void* struct_obj, int value_in);
void tao_wave_struct_get_ix_a2(const void* struct_obj, int* value_out);
void tao_wave_struct_set_ix_a2(void* struct_obj, int value_in);
void tao_wave_struct_get_ix_b1(const void* struct_obj, int* value_out);
void tao_wave_struct_set_ix_b1(void* struct_obj, int value_in);
void tao_wave_struct_get_ix_b2(const void* struct_obj, int* value_out);
void tao_wave_struct_set_ix_b2(void* struct_obj, int value_in);
void tao_wave_struct_get_i_a1(const void* struct_obj, int* value_out);
void tao_wave_struct_set_i_a1(void* struct_obj, int value_in);
void tao_wave_struct_get_i_a2(const void* struct_obj, int* value_out);
void tao_wave_struct_set_i_a2(void* struct_obj, int value_in);
void tao_wave_struct_get_i_b1(const void* struct_obj, int* value_out);
void tao_wave_struct_set_i_b1(void* struct_obj, int value_in);
void tao_wave_struct_get_i_b2(const void* struct_obj, int* value_out);
void tao_wave_struct_set_i_b2(void* struct_obj, int value_in);
void tao_wave_struct_get_n_a(const void* struct_obj, int* value_out);
void tao_wave_struct_set_n_a(void* struct_obj, int value_in);
void tao_wave_struct_get_n_b(const void* struct_obj, int* value_out);
void tao_wave_struct_set_n_b(void* struct_obj, int value_in);
void tao_wave_struct_get_i_curve_wrap_pt(
    const void* struct_obj,
    int* value_out);
void tao_wave_struct_set_i_curve_wrap_pt(void* struct_obj, int value_in);
void tao_wave_struct_get_ix_data_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void tao_wave_struct_get_n_kick(const void* struct_obj, int* value_out);
void tao_wave_struct_set_n_kick(void* struct_obj, int value_in);

void tao_wave_struct_get_kick_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_wave_struct_get_base_graph(const void* struct_obj, void** ptr_out);
void tao_wave_struct_set_base_graph(void* struct_obj, const void* src_ptr);
void tao_wave_struct_get_region(const void* struct_obj, void** ptr_out);
void tao_wave_struct_set_region(void* struct_obj, const void* src_ptr);
void tao_wave_struct_get_d1_dat(const void* struct_obj, void** ptr_out);
void tao_wave_struct_set_d1_dat(void* struct_obj, const void* src_ptr);
void tao_wave_kick_pt_struct_get_phi_s(
    const void* struct_obj,
    double* value_out);
void tao_wave_kick_pt_struct_set_phi_s(void* struct_obj, double value_in);
void tao_wave_kick_pt_struct_get_phi_r(
    const void* struct_obj,
    double* value_out);
void tao_wave_kick_pt_struct_set_phi_r(void* struct_obj, double value_in);
void tao_wave_kick_pt_struct_get_phi(const void* struct_obj, double* value_out);
void tao_wave_kick_pt_struct_set_phi(void* struct_obj, double value_in);
void tao_wave_kick_pt_struct_get_amp(const void* struct_obj, double* value_out);
void tao_wave_kick_pt_struct_set_amp(void* struct_obj, double value_in);
void tao_wave_kick_pt_struct_get_s(const void* struct_obj, double* value_out);
void tao_wave_kick_pt_struct_set_s(void* struct_obj, double value_in);
void tao_wave_kick_pt_struct_get_ix_dat_before_kick(
    const void* struct_obj,
    int* value_out);
void tao_wave_kick_pt_struct_set_ix_dat_before_kick(
    void* struct_obj,
    int value_in);
void tao_wave_kick_pt_struct_get_ele(const void* struct_obj, void** ptr_out);
void tao_wave_kick_pt_struct_set_ele(void* struct_obj, const void* src_ptr);

void tao_cmd_history_struct_get_cmd_info(
    const void* s,
    char** d,
    int* len,
    bool* is_alloc);

void tao_cmd_history_struct_set_cmd(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void tao_cmd_history_struct_get_ix(const void* struct_obj, int* value_out);
void tao_cmd_history_struct_set_ix(void* struct_obj, int value_in);
void tao_universe_struct_get_model(const void* struct_obj, void** ptr_out);
void tao_universe_struct_set_model(void* struct_obj, const void* src_ptr);
void tao_universe_struct_get_design(const void* struct_obj, void** ptr_out);
void tao_universe_struct_set_design(void* struct_obj, const void* src_ptr);
void tao_universe_struct_get_base(const void* struct_obj, void** ptr_out);
void tao_universe_struct_set_base(void* struct_obj, const void* src_ptr);
void tao_universe_struct_get_beam(const void* struct_obj, void** ptr_out);
void tao_universe_struct_set_beam(void* struct_obj, const void* src_ptr);
void tao_universe_struct_get_dynamic_aperture(
    const void* struct_obj,
    void** ptr_out);
void tao_universe_struct_set_dynamic_aperture(
    void* struct_obj,
    const void* src_ptr);

void tao_universe_struct_get_model_branch_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_universe_struct_get_d2_data_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_universe_struct_get_data_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void tao_universe_struct_get_ping_scale(const void* struct_obj, void** ptr_out);
void tao_universe_struct_set_ping_scale(void* struct_obj, const void* src_ptr);
void tao_universe_struct_get_scratch_lat(
    const void* struct_obj,
    void** ptr_out);
void tao_universe_struct_set_scratch_lat(void* struct_obj, const void* src_ptr);
void tao_universe_struct_get_calc(const void* struct_obj, void** ptr_out);
void tao_universe_struct_set_calc(void* struct_obj, const void* src_ptr);
void tao_universe_struct_get_ele_order(const void* struct_obj, void** ptr_out);
void tao_universe_struct_set_ele_order(void* struct_obj, const void* src_ptr);
void tao_universe_struct_get_spin_map(const void* struct_obj, void** ptr_out);
void tao_universe_struct_set_spin_map(void* struct_obj, const void* src_ptr);
void tao_universe_struct_get_dModel_dVar_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void tao_universe_struct_get_ix_uni(const void* struct_obj, int* value_out);
void tao_universe_struct_set_ix_uni(void* struct_obj, int value_in);
void tao_universe_struct_get_n_d2_data_used(
    const void* struct_obj,
    int* value_out);
void tao_universe_struct_set_n_d2_data_used(void* struct_obj, int value_in);
void tao_universe_struct_get_n_data_used(
    const void* struct_obj,
    int* value_out);
void tao_universe_struct_set_n_data_used(void* struct_obj, int value_in);
void tao_universe_struct_get_is_on(const void* struct_obj, bool* value_out);
void tao_universe_struct_set_is_on(void* struct_obj, bool value_in);
void tao_universe_struct_get_design_same_as_previous(
    const void* struct_obj,
    bool* value_out);
void tao_universe_struct_set_design_same_as_previous(
    void* struct_obj,
    bool value_in);
void tao_universe_struct_get_picked_uni(
    const void* struct_obj,
    bool* value_out);
void tao_universe_struct_set_picked_uni(void* struct_obj, bool value_in);
void mad_energy_struct_get_total(const void* struct_obj, double* value_out);
void mad_energy_struct_set_total(void* struct_obj, double value_in);
void mad_energy_struct_get_beta(const void* struct_obj, double* value_out);
void mad_energy_struct_set_beta(void* struct_obj, double value_in);
void mad_energy_struct_get_gamma(const void* struct_obj, double* value_out);
void mad_energy_struct_set_gamma(void* struct_obj, double value_in);
void mad_energy_struct_get_kinetic(const void* struct_obj, double* value_out);
void mad_energy_struct_set_kinetic(void* struct_obj, double value_in);
void mad_energy_struct_get_p0c(const void* struct_obj, double* value_out);
void mad_energy_struct_set_p0c(void* struct_obj, double value_in);
void mad_energy_struct_get_particle(const void* struct_obj, int* value_out);
void mad_energy_struct_set_particle(void* struct_obj, int value_in);
void mad_map_struct_get_k_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void mad_map_struct_get_r_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void mad_map_struct_get_t_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void random_state_struct_get_ix(const void* struct_obj, int64_t* value_out);
void random_state_struct_set_ix(void* struct_obj, int64_t value_in);
void random_state_struct_get_iy(const void* struct_obj, int64_t* value_out);
void random_state_struct_set_iy(void* struct_obj, int64_t value_in);
void random_state_struct_get_number_stored(
    const void* struct_obj,
    bool* value_out);
void random_state_struct_set_number_stored(void* struct_obj, bool value_in);
void random_state_struct_get_h_saved(const void* struct_obj, double* value_out);
void random_state_struct_set_h_saved(void* struct_obj, double value_in);
void random_state_struct_get_engine(const void* struct_obj, int* value_out);
void random_state_struct_set_engine(void* struct_obj, int value_in);
void random_state_struct_get_seed(const void* struct_obj, int* value_out);
void random_state_struct_set_seed(void* struct_obj, int value_in);
void random_state_struct_get_am(const void* struct_obj, double* value_out);
void random_state_struct_set_am(void* struct_obj, double value_in);
void random_state_struct_get_gauss_converter(
    const void* struct_obj,
    int* value_out);
void random_state_struct_set_gauss_converter(void* struct_obj, int value_in);
void random_state_struct_get_gauss_sigma_cut(
    const void* struct_obj,
    double* value_out);
void random_state_struct_set_gauss_sigma_cut(void* struct_obj, double value_in);
void random_state_struct_get_in_sobseq(
    const void* struct_obj,
    int64_t* value_out);
void random_state_struct_set_in_sobseq(void* struct_obj, int64_t value_in);
void random_state_struct_get_x_sobseq_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void bbu_stage_struct_get_ix_ele_lr_wake(
    const void* struct_obj,
    int* value_out);
void bbu_stage_struct_set_ix_ele_lr_wake(void* struct_obj, int value_in);
void bbu_stage_struct_get_ix_ele_stage_end(
    const void* struct_obj,
    int* value_out);
void bbu_stage_struct_set_ix_ele_stage_end(void* struct_obj, int value_in);
void bbu_stage_struct_get_ix_pass(const void* struct_obj, int* value_out);
void bbu_stage_struct_set_ix_pass(void* struct_obj, int value_in);
void bbu_stage_struct_get_ix_stage_pass1(
    const void* struct_obj,
    int* value_out);
void bbu_stage_struct_set_ix_stage_pass1(void* struct_obj, int value_in);
void bbu_stage_struct_get_ix_head_bunch(const void* struct_obj, int* value_out);
void bbu_stage_struct_set_ix_head_bunch(void* struct_obj, int value_in);
void bbu_stage_struct_get_ix_hom_max(const void* struct_obj, int* value_out);
void bbu_stage_struct_set_ix_hom_max(void* struct_obj, int value_in);
void bbu_stage_struct_get_hom_voltage_max(
    const void* struct_obj,
    double* value_out);
void bbu_stage_struct_set_hom_voltage_max(void* struct_obj, double value_in);
void bbu_stage_struct_get_time_at_wake_ele(
    const void* struct_obj,
    double* value_out);
void bbu_stage_struct_set_time_at_wake_ele(void* struct_obj, double value_in);
void bbu_stage_struct_get_ave_orb_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void bbu_stage_struct_get_rms_orb_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void bbu_stage_struct_get_min_orb_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void bbu_stage_struct_get_max_orb_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void bbu_stage_struct_get_n_orb(const void* struct_obj, int* value_out);
void bbu_stage_struct_set_n_orb(void* struct_obj, int value_in);

void bbu_beam_struct_get_bunch_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bbu_beam_struct_get_stage_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void bbu_beam_struct_get_ix_ele_bunch_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void bbu_beam_struct_get_ix_bunch_head(const void* struct_obj, int* value_out);
void bbu_beam_struct_set_ix_bunch_head(void* struct_obj, int value_in);
void bbu_beam_struct_get_ix_bunch_end(const void* struct_obj, int* value_out);
void bbu_beam_struct_set_ix_bunch_end(void* struct_obj, int value_in);
void bbu_beam_struct_get_n_bunch_in_lat(const void* struct_obj, int* value_out);
void bbu_beam_struct_set_n_bunch_in_lat(void* struct_obj, int value_in);
void bbu_beam_struct_get_ix_stage_voltage_max(
    const void* struct_obj,
    int* value_out);
void bbu_beam_struct_set_ix_stage_voltage_max(void* struct_obj, int value_in);
void bbu_beam_struct_get_hom_voltage_max(
    const void* struct_obj,
    double* value_out);
void bbu_beam_struct_set_hom_voltage_max(void* struct_obj, double value_in);
void bbu_beam_struct_get_time_now(const void* struct_obj, double* value_out);
void bbu_beam_struct_set_time_now(void* struct_obj, double value_in);
void bbu_beam_struct_get_one_turn_time(
    const void* struct_obj,
    double* value_out);
void bbu_beam_struct_set_one_turn_time(void* struct_obj, double value_in);
void bbu_beam_struct_get_rf_wavelength_max(
    const void* struct_obj,
    double* value_out);
void bbu_beam_struct_set_rf_wavelength_max(void* struct_obj, double value_in);
void bbu_param_struct_get_lat_filename_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void bbu_param_struct_set_lat_filename(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void bbu_param_struct_get_lat2_filename_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void bbu_param_struct_set_lat2_filename(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void bbu_param_struct_get_bunch_by_bunch_info_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void bbu_param_struct_set_bunch_by_bunch_info_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void bbu_param_struct_get_hybridize(const void* struct_obj, bool* value_out);
void bbu_param_struct_set_hybridize(void* struct_obj, bool value_in);
void bbu_param_struct_get_write_digested_hybrid_lat(
    const void* struct_obj,
    bool* value_out);
void bbu_param_struct_set_write_digested_hybrid_lat(
    void* struct_obj,
    bool value_in);
void bbu_param_struct_get_write_voltage_vs_time_dat(
    const void* struct_obj,
    bool* value_out);
void bbu_param_struct_set_write_voltage_vs_time_dat(
    void* struct_obj,
    bool value_in);
void bbu_param_struct_get_keep_overlays_and_groups(
    const void* struct_obj,
    bool* value_out);
void bbu_param_struct_set_keep_overlays_and_groups(
    void* struct_obj,
    bool value_in);
void bbu_param_struct_get_keep_all_lcavities(
    const void* struct_obj,
    bool* value_out);
void bbu_param_struct_set_keep_all_lcavities(void* struct_obj, bool value_in);
void bbu_param_struct_get_use_taylor_for_hybrids(
    const void* struct_obj,
    bool* value_out);
void bbu_param_struct_set_use_taylor_for_hybrids(
    void* struct_obj,
    bool value_in);
void bbu_param_struct_get_stable_orbit_anal(
    const void* struct_obj,
    bool* value_out);
void bbu_param_struct_set_stable_orbit_anal(void* struct_obj, bool value_in);
void bbu_param_struct_get_limit_factor(
    const void* struct_obj,
    double* value_out);
void bbu_param_struct_set_limit_factor(void* struct_obj, double value_in);
void bbu_param_struct_get_simulation_turns_max(
    const void* struct_obj,
    double* value_out);
void bbu_param_struct_set_simulation_turns_max(
    void* struct_obj,
    double value_in);
void bbu_param_struct_get_bunch_freq(const void* struct_obj, double* value_out);
void bbu_param_struct_set_bunch_freq(void* struct_obj, double value_in);
void bbu_param_struct_get_init_particle_offset(
    const void* struct_obj,
    double* value_out);
void bbu_param_struct_set_init_particle_offset(
    void* struct_obj,
    double value_in);
void bbu_param_struct_get_current(const void* struct_obj, double* value_out);
void bbu_param_struct_set_current(void* struct_obj, double value_in);
void bbu_param_struct_get_rel_tol(const void* struct_obj, double* value_out);
void bbu_param_struct_set_rel_tol(void* struct_obj, double value_in);
void bbu_param_struct_get_drscan(const void* struct_obj, bool* value_out);
void bbu_param_struct_set_drscan(void* struct_obj, bool value_in);
void bbu_param_struct_get_use_interpolated_threshold(
    const void* struct_obj,
    bool* value_out);
void bbu_param_struct_set_use_interpolated_threshold(
    void* struct_obj,
    bool value_in);
void bbu_param_struct_get_write_hom_info(
    const void* struct_obj,
    bool* value_out);
void bbu_param_struct_set_write_hom_info(void* struct_obj, bool value_in);
void bbu_param_struct_get_elindex(const void* struct_obj, int* value_out);
void bbu_param_struct_set_elindex(void* struct_obj, int value_in);
void bbu_param_struct_get_elname_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void bbu_param_struct_set_elname(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void bbu_param_struct_get_nstep(const void* struct_obj, int* value_out);
void bbu_param_struct_set_nstep(void* struct_obj, int value_in);
void bbu_param_struct_get_begdr(const void* struct_obj, double* value_out);
void bbu_param_struct_set_begdr(void* struct_obj, double value_in);
void bbu_param_struct_get_enddr(const void* struct_obj, double* value_out);
void bbu_param_struct_set_enddr(void* struct_obj, double value_in);
void bbu_param_struct_get_nrep(const void* struct_obj, int* value_out);
void bbu_param_struct_set_nrep(void* struct_obj, int value_in);
void bbu_param_struct_get_ran_seed(const void* struct_obj, int* value_out);
void bbu_param_struct_set_ran_seed(void* struct_obj, int value_in);
void bbu_param_struct_get_hom_order_cutoff(
    const void* struct_obj,
    int* value_out);
void bbu_param_struct_set_hom_order_cutoff(void* struct_obj, int value_in);
void bbu_param_struct_get_ran_gauss_sigma_cut(
    const void* struct_obj,
    double* value_out);
void bbu_param_struct_set_ran_gauss_sigma_cut(
    void* struct_obj,
    double value_in);
void bbu_param_struct_get_ele_track_end_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void bbu_param_struct_set_ele_track_end(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void bbu_param_struct_get_ix_ele_track_end(
    const void* struct_obj,
    int* value_out);
void bbu_param_struct_set_ix_ele_track_end(void* struct_obj, int value_in);
void bbu_param_struct_get_regression(const void* struct_obj, bool* value_out);
void bbu_param_struct_set_regression(void* struct_obj, bool value_in);
void bbu_param_struct_get_normalize_z_to_rf(
    const void* struct_obj,
    bool* value_out);
void bbu_param_struct_set_normalize_z_to_rf(void* struct_obj, bool value_in);
void bbu_param_struct_get_ramp_on(const void* struct_obj, bool* value_out);
void bbu_param_struct_set_ramp_on(void* struct_obj, bool value_in);
void bbu_param_struct_get_ramp_pattern_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void bbu_param_struct_get_ramp_n_start(const void* struct_obj, int* value_out);
void bbu_param_struct_set_ramp_n_start(void* struct_obj, int value_in);
void bbu_param_struct_get_n_ramp_pattern(
    const void* struct_obj,
    int* value_out);
void bbu_param_struct_set_n_ramp_pattern(void* struct_obj, int value_in);
void all_encompassing_struct_get_real_rp_0d(
    const void* struct_obj,
    double* value_out);
void all_encompassing_struct_set_real_rp_0d(void* struct_obj, double value_in);
void all_encompassing_struct_get_real_rp_1d_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_real_rp_2d_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_rp_3d_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_rp_0d_ptr(
    const void* struct_obj,
    double** ptr_out);
void all_encompassing_struct_set_real_rp_0d_ptr(
    void* struct_obj,
    double value_in);
void all_encompassing_struct_get_real_rp_1d_ptr_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_real_rp_2d_ptr_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_rp_3d_ptr_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_rp_1d_alloc_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_real_rp_2d_alloc_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_rp_3d_alloc_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_dp_0d(
    const void* struct_obj,
    double* value_out);
void all_encompassing_struct_set_real_dp_0d(void* struct_obj, double value_in);
void all_encompassing_struct_get_real_dp_1d_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_real_dp_2d_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_dp_3d_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_dp_0d_ptr(
    const void* struct_obj,
    double** ptr_out);
void all_encompassing_struct_set_real_dp_0d_ptr(
    void* struct_obj,
    double value_in);
void all_encompassing_struct_get_real_dp_1d_ptr_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_real_dp_2d_ptr_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_dp_3d_ptr_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_dp_1d_alloc_info(
    const void* s,
    double** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_real_dp_2d_alloc_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_real_dp_3d_alloc_info(
    const void* s,
    double** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_complex_dp_0d(
    const void* struct_obj,
    std::complex<double>* value_out);
void all_encompassing_struct_set_complex_dp_0d(
    void* struct_obj,
    std::complex<double> value_in);
void all_encompassing_struct_get_complex_dp_1d_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_complex_dp_2d_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_complex_dp_3d_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_complex_dp_1d_ptr_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_complex_dp_2d_ptr_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_complex_dp_3d_ptr_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_complex_dp_1d_alloc_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_complex_dp_2d_alloc_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_complex_dp_3d_alloc_info(
    const void* s,
    std::complex<double>** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_int_0d(const void* struct_obj, int* value_out);
void all_encompassing_struct_set_int_0d(void* struct_obj, int value_in);
void all_encompassing_struct_get_int_1d_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_int_2d_info(
    const void* s,
    int** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_int_3d_info(
    const void* s,
    int** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_int_0d_ptr(
    const void* struct_obj,
    int** ptr_out);
void all_encompassing_struct_set_int_0d_ptr(void* struct_obj, int value_in);
void all_encompassing_struct_get_int_1d_ptr_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_int_2d_ptr_info(
    const void* s,
    int** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_int_3d_ptr_info(
    const void* s,
    int** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_int_1d_alloc_info(
    const void* s,
    int** d,
    int* bounds,
    bool* is_alloc);
void all_encompassing_struct_get_int_2d_alloc_info(
    const void* s,
    int** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_int_3d_alloc_info(
    const void* s,
    int** d,
    int* bounds,
    int* strides,
    bool* is_alloc);
void all_encompassing_struct_get_int8_0d(
    const void* struct_obj,
    int64_t* value_out);
void all_encompassing_struct_set_int8_0d(void* struct_obj, int64_t value_in);
void all_encompassing_struct_get_int8_0d_ptr(
    const void* struct_obj,
    int64_t** ptr_out);
void all_encompassing_struct_set_int8_0d_ptr(
    void* struct_obj,
    int64_t value_in);
void all_encompassing_struct_get_logical_0d(
    const void* struct_obj,
    bool* value_out);
void all_encompassing_struct_set_logical_0d(void* struct_obj, bool value_in);
void all_encompassing_struct_get_logical_0d_ptr(
    const void* struct_obj,
    bool** ptr_out);
void all_encompassing_struct_set_logical_0d_ptr(
    void* struct_obj,
    bool value_in);
void all_encompassing_struct_get_type_0d(
    const void* struct_obj,
    void** ptr_out);
void all_encompassing_struct_set_type_0d(void* struct_obj, const void* src_ptr);

void all_encompassing_struct_get_type_1d_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void all_encompassing_struct_get_type_2d_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void all_encompassing_struct_get_type_3d_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void all_encompassing_struct_get_type_0d_ptr(
    const void* struct_obj,
    void** ptr_out);
void all_encompassing_struct_set_type_0d_ptr(
    void* struct_obj,
    const void* src_ptr);

void all_encompassing_struct_get_type_1d_ptr_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void all_encompassing_struct_get_type_2d_ptr_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void all_encompassing_struct_get_type_3d_ptr_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void all_encompassing_struct_get_type_1d_alloc_info(
    const void* s,
    void** d,
    int* bounds,
    bool* is_alloc,
    size_t* el_size);

void all_encompassing_struct_get_type_2d_alloc_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void all_encompassing_struct_get_type_3d_alloc_info(
    const void* s,
    void** d,
    int* bounds,
    int* strides,
    bool* a,
    size_t* es);

void test_sub_struct_get_sr(const void* struct_obj, void** ptr_out);
void test_sub_struct_set_sr(void* struct_obj, const void* src_ptr);
void test_sub_sub_struct_get_aaa(const void* struct_obj, int64_t* value_out);
void test_sub_sub_struct_set_aaa(void* struct_obj, int64_t value_in);
void test_sub_sub_struct_get_bbb(const void* struct_obj, int* value_out);
void test_sub_sub_struct_set_bbb(void* struct_obj, int value_in);
void test_sub_sub_struct_get_file_info(
    const void* s,
    char** d,
    int* bounds,
    bool* a);
void test_sub_sub_struct_set_file(
    void* struct_obj,
    const char* str_ptr,
    int str_len);
void test_sub_sub_struct_get_t_ref(const void* struct_obj, double* value_out);
void test_sub_sub_struct_set_t_ref(void* struct_obj, double value_in);
void test_sub_sub_struct_get_freq_spread(
    const void* struct_obj,
    double* value_out);
void test_sub_sub_struct_set_freq_spread(void* struct_obj, double value_in);
}

namespace Bmad {

extern "C" {

void* allocate_real_container();
void reallocate_real_container_data(void*, int, size_t) noexcept;
void deallocate_real_container(void*) noexcept;
void access_real_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_real16_container();
void reallocate_real16_container_data(void*, int, size_t) noexcept;
void deallocate_real16_container(void*) noexcept;
void access_real16_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_integer_container();
void reallocate_integer_container_data(void*, int, size_t) noexcept;
void deallocate_integer_container(void*) noexcept;
void access_integer_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_integer8_container();
void reallocate_integer8_container_data(void*, int, size_t) noexcept;
void deallocate_integer8_container(void*) noexcept;
void access_integer8_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_logical_container();
void reallocate_logical_container_data(void*, int, size_t) noexcept;
void deallocate_logical_container(void*) noexcept;
void access_logical_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_complex_container();
void reallocate_complex_container_data(void*, int, size_t) noexcept;
void deallocate_complex_container(void*) noexcept;
void access_complex_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_spline_struct(int n, size_t* element_size);
void deallocate_fortran_spline_struct(void* ptr, int n) noexcept;
void copy_fortran_spline_struct(const void* src, void* dst);

void* allocate_spline_struct_container();
void reallocate_spline_struct_container_data(void*, int, size_t) noexcept;
void deallocate_spline_struct_container(void*) noexcept;
void access_spline_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_spin_polar_struct(int n, size_t* element_size);
void deallocate_fortran_spin_polar_struct(void* ptr, int n) noexcept;
void copy_fortran_spin_polar_struct(const void* src, void* dst);

void* allocate_spin_polar_struct_container();
void reallocate_spin_polar_struct_container_data(void*, int, size_t) noexcept;
void deallocate_spin_polar_struct_container(void*) noexcept;
void access_spin_polar_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_ac_kicker_time_struct(int n, size_t* element_size);
void deallocate_fortran_ac_kicker_time_struct(void* ptr, int n) noexcept;
void copy_fortran_ac_kicker_time_struct(const void* src, void* dst);

void* allocate_ac_kicker_time_struct_container();
void reallocate_ac_kicker_time_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_ac_kicker_time_struct_container(void*) noexcept;
void access_ac_kicker_time_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_ac_kicker_freq_struct(int n, size_t* element_size);
void deallocate_fortran_ac_kicker_freq_struct(void* ptr, int n) noexcept;
void copy_fortran_ac_kicker_freq_struct(const void* src, void* dst);

void* allocate_ac_kicker_freq_struct_container();
void reallocate_ac_kicker_freq_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_ac_kicker_freq_struct_container(void*) noexcept;
void access_ac_kicker_freq_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_ac_kicker_struct(int n, size_t* element_size);
void deallocate_fortran_ac_kicker_struct(void* ptr, int n) noexcept;
void copy_fortran_ac_kicker_struct(const void* src, void* dst);

void* allocate_ac_kicker_struct_container();
void reallocate_ac_kicker_struct_container_data(void*, int, size_t) noexcept;
void deallocate_ac_kicker_struct_container(void*) noexcept;
void access_ac_kicker_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_interval1_coef_struct(int n, size_t* element_size);
void deallocate_fortran_interval1_coef_struct(void* ptr, int n) noexcept;
void copy_fortran_interval1_coef_struct(const void* src, void* dst);

void* allocate_interval1_coef_struct_container();
void reallocate_interval1_coef_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_interval1_coef_struct_container(void*) noexcept;
void access_interval1_coef_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_photon_reflect_table_struct(int n, size_t* element_size);
void deallocate_fortran_photon_reflect_table_struct(void* ptr, int n) noexcept;
void copy_fortran_photon_reflect_table_struct(const void* src, void* dst);

void* allocate_photon_reflect_table_struct_container();
void reallocate_photon_reflect_table_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_photon_reflect_table_struct_container(void*) noexcept;
void access_photon_reflect_table_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_photon_reflect_surface_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_photon_reflect_surface_struct(
    void* ptr,
    int n) noexcept;
void copy_fortran_photon_reflect_surface_struct(const void* src, void* dst);

void* allocate_photon_reflect_surface_struct_container();
void reallocate_photon_reflect_surface_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_photon_reflect_surface_struct_container(void*) noexcept;
void access_photon_reflect_surface_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_coord_struct(int n, size_t* element_size);
void deallocate_fortran_coord_struct(void* ptr, int n) noexcept;
void copy_fortran_coord_struct(const void* src, void* dst);

void* allocate_coord_struct_container();
void reallocate_coord_struct_container_data(void*, int, size_t) noexcept;
void deallocate_coord_struct_container(void*) noexcept;
void access_coord_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_coord_array_struct(int n, size_t* element_size);
void deallocate_fortran_coord_array_struct(void* ptr, int n) noexcept;
void copy_fortran_coord_array_struct(const void* src, void* dst);

void* allocate_coord_array_struct_container();
void reallocate_coord_array_struct_container_data(void*, int, size_t) noexcept;
void deallocate_coord_array_struct_container(void*) noexcept;
void access_coord_array_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bpm_phase_coupling_struct(int n, size_t* element_size);
void deallocate_fortran_bpm_phase_coupling_struct(void* ptr, int n) noexcept;
void copy_fortran_bpm_phase_coupling_struct(const void* src, void* dst);

void* allocate_bpm_phase_coupling_struct_container();
void reallocate_bpm_phase_coupling_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_bpm_phase_coupling_struct_container(void*) noexcept;
void access_bpm_phase_coupling_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_expression_atom_struct(int n, size_t* element_size);
void deallocate_fortran_expression_atom_struct(void* ptr, int n) noexcept;
void copy_fortran_expression_atom_struct(const void* src, void* dst);

void* allocate_expression_atom_struct_container();
void reallocate_expression_atom_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_expression_atom_struct_container(void*) noexcept;
void access_expression_atom_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_wake_sr_z_long_struct(int n, size_t* element_size);
void deallocate_fortran_wake_sr_z_long_struct(void* ptr, int n) noexcept;
void copy_fortran_wake_sr_z_long_struct(const void* src, void* dst);

void* allocate_wake_sr_z_long_struct_container();
void reallocate_wake_sr_z_long_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_wake_sr_z_long_struct_container(void*) noexcept;
void access_wake_sr_z_long_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_wake_sr_mode_struct(int n, size_t* element_size);
void deallocate_fortran_wake_sr_mode_struct(void* ptr, int n) noexcept;
void copy_fortran_wake_sr_mode_struct(const void* src, void* dst);

void* allocate_wake_sr_mode_struct_container();
void reallocate_wake_sr_mode_struct_container_data(void*, int, size_t) noexcept;
void deallocate_wake_sr_mode_struct_container(void*) noexcept;
void access_wake_sr_mode_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_wake_sr_struct(int n, size_t* element_size);
void deallocate_fortran_wake_sr_struct(void* ptr, int n) noexcept;
void copy_fortran_wake_sr_struct(const void* src, void* dst);

void* allocate_wake_sr_struct_container();
void reallocate_wake_sr_struct_container_data(void*, int, size_t) noexcept;
void deallocate_wake_sr_struct_container(void*) noexcept;
void access_wake_sr_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_wake_lr_mode_struct(int n, size_t* element_size);
void deallocate_fortran_wake_lr_mode_struct(void* ptr, int n) noexcept;
void copy_fortran_wake_lr_mode_struct(const void* src, void* dst);

void* allocate_wake_lr_mode_struct_container();
void reallocate_wake_lr_mode_struct_container_data(void*, int, size_t) noexcept;
void deallocate_wake_lr_mode_struct_container(void*) noexcept;
void access_wake_lr_mode_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_wake_lr_struct(int n, size_t* element_size);
void deallocate_fortran_wake_lr_struct(void* ptr, int n) noexcept;
void copy_fortran_wake_lr_struct(const void* src, void* dst);

void* allocate_wake_lr_struct_container();
void reallocate_wake_lr_struct_container_data(void*, int, size_t) noexcept;
void deallocate_wake_lr_struct_container(void*) noexcept;
void access_wake_lr_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_lat_ele_loc_struct(int n, size_t* element_size);
void deallocate_fortran_lat_ele_loc_struct(void* ptr, int n) noexcept;
void copy_fortran_lat_ele_loc_struct(const void* src, void* dst);

void* allocate_lat_ele_loc_struct_container();
void reallocate_lat_ele_loc_struct_container_data(void*, int, size_t) noexcept;
void deallocate_lat_ele_loc_struct_container(void*) noexcept;
void access_lat_ele_loc_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_wake_struct(int n, size_t* element_size);
void deallocate_fortran_wake_struct(void* ptr, int n) noexcept;
void copy_fortran_wake_struct(const void* src, void* dst);

void* allocate_wake_struct_container();
void reallocate_wake_struct_container_data(void*, int, size_t) noexcept;
void deallocate_wake_struct_container(void*) noexcept;
void access_wake_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_taylor_term_struct(int n, size_t* element_size);
void deallocate_fortran_taylor_term_struct(void* ptr, int n) noexcept;
void copy_fortran_taylor_term_struct(const void* src, void* dst);

void* allocate_taylor_term_struct_container();
void reallocate_taylor_term_struct_container_data(void*, int, size_t) noexcept;
void deallocate_taylor_term_struct_container(void*) noexcept;
void access_taylor_term_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_taylor_struct(int n, size_t* element_size);
void deallocate_fortran_taylor_struct(void* ptr, int n) noexcept;
void copy_fortran_taylor_struct(const void* src, void* dst);

void* allocate_taylor_struct_container();
void reallocate_taylor_struct_container_data(void*, int, size_t) noexcept;
void deallocate_taylor_struct_container(void*) noexcept;
void access_taylor_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_em_taylor_term_struct(int n, size_t* element_size);
void deallocate_fortran_em_taylor_term_struct(void* ptr, int n) noexcept;
void copy_fortran_em_taylor_term_struct(const void* src, void* dst);

void* allocate_em_taylor_term_struct_container();
void reallocate_em_taylor_term_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_em_taylor_term_struct_container(void*) noexcept;
void access_em_taylor_term_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_em_taylor_struct(int n, size_t* element_size);
void deallocate_fortran_em_taylor_struct(void* ptr, int n) noexcept;
void copy_fortran_em_taylor_struct(const void* src, void* dst);

void* allocate_em_taylor_struct_container();
void reallocate_em_taylor_struct_container_data(void*, int, size_t) noexcept;
void deallocate_em_taylor_struct_container(void*) noexcept;
void access_em_taylor_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_cartesian_map_term1_struct(int n, size_t* element_size);
void deallocate_fortran_cartesian_map_term1_struct(void* ptr, int n) noexcept;
void copy_fortran_cartesian_map_term1_struct(const void* src, void* dst);

void* allocate_cartesian_map_term1_struct_container();
void reallocate_cartesian_map_term1_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_cartesian_map_term1_struct_container(void*) noexcept;
void access_cartesian_map_term1_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_cartesian_map_term_struct(int n, size_t* element_size);
void deallocate_fortran_cartesian_map_term_struct(void* ptr, int n) noexcept;
void copy_fortran_cartesian_map_term_struct(const void* src, void* dst);

void* allocate_cartesian_map_term_struct_container();
void reallocate_cartesian_map_term_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_cartesian_map_term_struct_container(void*) noexcept;
void access_cartesian_map_term_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_cartesian_map_struct(int n, size_t* element_size);
void deallocate_fortran_cartesian_map_struct(void* ptr, int n) noexcept;
void copy_fortran_cartesian_map_struct(const void* src, void* dst);

void* allocate_cartesian_map_struct_container();
void reallocate_cartesian_map_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_cartesian_map_struct_container(void*) noexcept;
void access_cartesian_map_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_cylindrical_map_term1_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_cylindrical_map_term1_struct(void* ptr, int n) noexcept;
void copy_fortran_cylindrical_map_term1_struct(const void* src, void* dst);

void* allocate_cylindrical_map_term1_struct_container();
void reallocate_cylindrical_map_term1_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_cylindrical_map_term1_struct_container(void*) noexcept;
void access_cylindrical_map_term1_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_cylindrical_map_term_struct(int n, size_t* element_size);
void deallocate_fortran_cylindrical_map_term_struct(void* ptr, int n) noexcept;
void copy_fortran_cylindrical_map_term_struct(const void* src, void* dst);

void* allocate_cylindrical_map_term_struct_container();
void reallocate_cylindrical_map_term_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_cylindrical_map_term_struct_container(void*) noexcept;
void access_cylindrical_map_term_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_cylindrical_map_struct(int n, size_t* element_size);
void deallocate_fortran_cylindrical_map_struct(void* ptr, int n) noexcept;
void copy_fortran_cylindrical_map_struct(const void* src, void* dst);

void* allocate_cylindrical_map_struct_container();
void reallocate_cylindrical_map_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_cylindrical_map_struct_container(void*) noexcept;
void access_cylindrical_map_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bicubic_cmplx_coef_struct(int n, size_t* element_size);
void deallocate_fortran_bicubic_cmplx_coef_struct(void* ptr, int n) noexcept;
void copy_fortran_bicubic_cmplx_coef_struct(const void* src, void* dst);

void* allocate_bicubic_cmplx_coef_struct_container();
void reallocate_bicubic_cmplx_coef_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_bicubic_cmplx_coef_struct_container(void*) noexcept;
void access_bicubic_cmplx_coef_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tricubic_cmplx_coef_struct(int n, size_t* element_size);
void deallocate_fortran_tricubic_cmplx_coef_struct(void* ptr, int n) noexcept;
void copy_fortran_tricubic_cmplx_coef_struct(const void* src, void* dst);

void* allocate_tricubic_cmplx_coef_struct_container();
void reallocate_tricubic_cmplx_coef_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tricubic_cmplx_coef_struct_container(void*) noexcept;
void access_tricubic_cmplx_coef_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_grid_field_pt1_struct(int n, size_t* element_size);
void deallocate_fortran_grid_field_pt1_struct(void* ptr, int n) noexcept;
void copy_fortran_grid_field_pt1_struct(const void* src, void* dst);

void* allocate_grid_field_pt1_struct_container();
void reallocate_grid_field_pt1_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_grid_field_pt1_struct_container(void*) noexcept;
void access_grid_field_pt1_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_grid_field_pt_struct(int n, size_t* element_size);
void deallocate_fortran_grid_field_pt_struct(void* ptr, int n) noexcept;
void copy_fortran_grid_field_pt_struct(const void* src, void* dst);

void* allocate_grid_field_pt_struct_container();
void reallocate_grid_field_pt_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_grid_field_pt_struct_container(void*) noexcept;
void access_grid_field_pt_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_grid_field_struct(int n, size_t* element_size);
void deallocate_fortran_grid_field_struct(void* ptr, int n) noexcept;
void copy_fortran_grid_field_struct(const void* src, void* dst);

void* allocate_grid_field_struct_container();
void reallocate_grid_field_struct_container_data(void*, int, size_t) noexcept;
void deallocate_grid_field_struct_container(void*) noexcept;
void access_grid_field_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_floor_position_struct(int n, size_t* element_size);
void deallocate_fortran_floor_position_struct(void* ptr, int n) noexcept;
void copy_fortran_floor_position_struct(const void* src, void* dst);

void* allocate_floor_position_struct_container();
void reallocate_floor_position_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_floor_position_struct_container(void*) noexcept;
void access_floor_position_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_high_energy_space_charge_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_high_energy_space_charge_struct(
    void* ptr,
    int n) noexcept;
void copy_fortran_high_energy_space_charge_struct(const void* src, void* dst);

void* allocate_high_energy_space_charge_struct_container();
void reallocate_high_energy_space_charge_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_high_energy_space_charge_struct_container(void*) noexcept;
void access_high_energy_space_charge_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_xy_disp_struct(int n, size_t* element_size);
void deallocate_fortran_xy_disp_struct(void* ptr, int n) noexcept;
void copy_fortran_xy_disp_struct(const void* src, void* dst);

void* allocate_xy_disp_struct_container();
void reallocate_xy_disp_struct_container_data(void*, int, size_t) noexcept;
void deallocate_xy_disp_struct_container(void*) noexcept;
void access_xy_disp_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_twiss_struct(int n, size_t* element_size);
void deallocate_fortran_twiss_struct(void* ptr, int n) noexcept;
void copy_fortran_twiss_struct(const void* src, void* dst);

void* allocate_twiss_struct_container();
void reallocate_twiss_struct_container_data(void*, int, size_t) noexcept;
void deallocate_twiss_struct_container(void*) noexcept;
void access_twiss_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_mode3_struct(int n, size_t* element_size);
void deallocate_fortran_mode3_struct(void* ptr, int n) noexcept;
void copy_fortran_mode3_struct(const void* src, void* dst);

void* allocate_mode3_struct_container();
void reallocate_mode3_struct_container_data(void*, int, size_t) noexcept;
void deallocate_mode3_struct_container(void*) noexcept;
void access_mode3_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bookkeeping_state_struct(int n, size_t* element_size);
void deallocate_fortran_bookkeeping_state_struct(void* ptr, int n) noexcept;
void copy_fortran_bookkeeping_state_struct(const void* src, void* dst);

void* allocate_bookkeeping_state_struct_container();
void reallocate_bookkeeping_state_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_bookkeeping_state_struct_container(void*) noexcept;
void access_bookkeeping_state_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_rad_map_struct(int n, size_t* element_size);
void deallocate_fortran_rad_map_struct(void* ptr, int n) noexcept;
void copy_fortran_rad_map_struct(const void* src, void* dst);

void* allocate_rad_map_struct_container();
void reallocate_rad_map_struct_container_data(void*, int, size_t) noexcept;
void deallocate_rad_map_struct_container(void*) noexcept;
void access_rad_map_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_rad_map_ele_struct(int n, size_t* element_size);
void deallocate_fortran_rad_map_ele_struct(void* ptr, int n) noexcept;
void copy_fortran_rad_map_ele_struct(const void* src, void* dst);

void* allocate_rad_map_ele_struct_container();
void reallocate_rad_map_ele_struct_container_data(void*, int, size_t) noexcept;
void deallocate_rad_map_ele_struct_container(void*) noexcept;
void access_rad_map_ele_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_gen_grad1_struct(int n, size_t* element_size);
void deallocate_fortran_gen_grad1_struct(void* ptr, int n) noexcept;
void copy_fortran_gen_grad1_struct(const void* src, void* dst);

void* allocate_gen_grad1_struct_container();
void reallocate_gen_grad1_struct_container_data(void*, int, size_t) noexcept;
void deallocate_gen_grad1_struct_container(void*) noexcept;
void access_gen_grad1_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_gen_grad_map_struct(int n, size_t* element_size);
void deallocate_fortran_gen_grad_map_struct(void* ptr, int n) noexcept;
void copy_fortran_gen_grad_map_struct(const void* src, void* dst);

void* allocate_gen_grad_map_struct_container();
void reallocate_gen_grad_map_struct_container_data(void*, int, size_t) noexcept;
void deallocate_gen_grad_map_struct_container(void*) noexcept;
void access_gen_grad_map_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_surface_segmented_pt_struct(int n, size_t* element_size);
void deallocate_fortran_surface_segmented_pt_struct(void* ptr, int n) noexcept;
void copy_fortran_surface_segmented_pt_struct(const void* src, void* dst);

void* allocate_surface_segmented_pt_struct_container();
void reallocate_surface_segmented_pt_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_surface_segmented_pt_struct_container(void*) noexcept;
void access_surface_segmented_pt_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_surface_segmented_struct(int n, size_t* element_size);
void deallocate_fortran_surface_segmented_struct(void* ptr, int n) noexcept;
void copy_fortran_surface_segmented_struct(const void* src, void* dst);

void* allocate_surface_segmented_struct_container();
void reallocate_surface_segmented_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_surface_segmented_struct_container(void*) noexcept;
void access_surface_segmented_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_surface_h_misalign_pt_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_surface_h_misalign_pt_struct(void* ptr, int n) noexcept;
void copy_fortran_surface_h_misalign_pt_struct(const void* src, void* dst);

void* allocate_surface_h_misalign_pt_struct_container();
void reallocate_surface_h_misalign_pt_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_surface_h_misalign_pt_struct_container(void*) noexcept;
void access_surface_h_misalign_pt_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_surface_h_misalign_struct(int n, size_t* element_size);
void deallocate_fortran_surface_h_misalign_struct(void* ptr, int n) noexcept;
void copy_fortran_surface_h_misalign_struct(const void* src, void* dst);

void* allocate_surface_h_misalign_struct_container();
void reallocate_surface_h_misalign_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_surface_h_misalign_struct_container(void*) noexcept;
void access_surface_h_misalign_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_surface_displacement_pt_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_surface_displacement_pt_struct(
    void* ptr,
    int n) noexcept;
void copy_fortran_surface_displacement_pt_struct(const void* src, void* dst);

void* allocate_surface_displacement_pt_struct_container();
void reallocate_surface_displacement_pt_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_surface_displacement_pt_struct_container(void*) noexcept;
void access_surface_displacement_pt_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_surface_displacement_struct(int n, size_t* element_size);
void deallocate_fortran_surface_displacement_struct(void* ptr, int n) noexcept;
void copy_fortran_surface_displacement_struct(const void* src, void* dst);

void* allocate_surface_displacement_struct_container();
void reallocate_surface_displacement_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_surface_displacement_struct_container(void*) noexcept;
void access_surface_displacement_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_target_point_struct(int n, size_t* element_size);
void deallocate_fortran_target_point_struct(void* ptr, int n) noexcept;
void copy_fortran_target_point_struct(const void* src, void* dst);

void* allocate_target_point_struct_container();
void reallocate_target_point_struct_container_data(void*, int, size_t) noexcept;
void deallocate_target_point_struct_container(void*) noexcept;
void access_target_point_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_surface_curvature_struct(int n, size_t* element_size);
void deallocate_fortran_surface_curvature_struct(void* ptr, int n) noexcept;
void copy_fortran_surface_curvature_struct(const void* src, void* dst);

void* allocate_surface_curvature_struct_container();
void reallocate_surface_curvature_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_surface_curvature_struct_container(void*) noexcept;
void access_surface_curvature_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_photon_target_struct(int n, size_t* element_size);
void deallocate_fortran_photon_target_struct(void* ptr, int n) noexcept;
void copy_fortran_photon_target_struct(const void* src, void* dst);

void* allocate_photon_target_struct_container();
void reallocate_photon_target_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_photon_target_struct_container(void*) noexcept;
void access_photon_target_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_photon_material_struct(int n, size_t* element_size);
void deallocate_fortran_photon_material_struct(void* ptr, int n) noexcept;
void copy_fortran_photon_material_struct(const void* src, void* dst);

void* allocate_photon_material_struct_container();
void reallocate_photon_material_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_photon_material_struct_container(void*) noexcept;
void access_photon_material_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_pixel_pt_struct(int n, size_t* element_size);
void deallocate_fortran_pixel_pt_struct(void* ptr, int n) noexcept;
void copy_fortran_pixel_pt_struct(const void* src, void* dst);

void* allocate_pixel_pt_struct_container();
void reallocate_pixel_pt_struct_container_data(void*, int, size_t) noexcept;
void deallocate_pixel_pt_struct_container(void*) noexcept;
void access_pixel_pt_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_pixel_detec_struct(int n, size_t* element_size);
void deallocate_fortran_pixel_detec_struct(void* ptr, int n) noexcept;
void copy_fortran_pixel_detec_struct(const void* src, void* dst);

void* allocate_pixel_detec_struct_container();
void reallocate_pixel_detec_struct_container_data(void*, int, size_t) noexcept;
void deallocate_pixel_detec_struct_container(void*) noexcept;
void access_pixel_detec_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_photon_element_struct(int n, size_t* element_size);
void deallocate_fortran_photon_element_struct(void* ptr, int n) noexcept;
void copy_fortran_photon_element_struct(const void* src, void* dst);

void* allocate_photon_element_struct_container();
void reallocate_photon_element_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_photon_element_struct_container(void*) noexcept;
void access_photon_element_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_wall3d_vertex_struct(int n, size_t* element_size);
void deallocate_fortran_wall3d_vertex_struct(void* ptr, int n) noexcept;
void copy_fortran_wall3d_vertex_struct(const void* src, void* dst);

void* allocate_wall3d_vertex_struct_container();
void reallocate_wall3d_vertex_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_wall3d_vertex_struct_container(void*) noexcept;
void access_wall3d_vertex_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_wall3d_section_struct(int n, size_t* element_size);
void deallocate_fortran_wall3d_section_struct(void* ptr, int n) noexcept;
void copy_fortran_wall3d_section_struct(const void* src, void* dst);

void* allocate_wall3d_section_struct_container();
void reallocate_wall3d_section_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_wall3d_section_struct_container(void*) noexcept;
void access_wall3d_section_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_wall3d_struct(int n, size_t* element_size);
void deallocate_fortran_wall3d_struct(void* ptr, int n) noexcept;
void copy_fortran_wall3d_struct(const void* src, void* dst);

void* allocate_wall3d_struct_container();
void reallocate_wall3d_struct_container_data(void*, int, size_t) noexcept;
void deallocate_wall3d_struct_container(void*) noexcept;
void access_wall3d_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_ramper_lord_struct(int n, size_t* element_size);
void deallocate_fortran_ramper_lord_struct(void* ptr, int n) noexcept;
void copy_fortran_ramper_lord_struct(const void* src, void* dst);

void* allocate_ramper_lord_struct_container();
void reallocate_ramper_lord_struct_container_data(void*, int, size_t) noexcept;
void deallocate_ramper_lord_struct_container(void*) noexcept;
void access_ramper_lord_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_control_struct(int n, size_t* element_size);
void deallocate_fortran_control_struct(void* ptr, int n) noexcept;
void copy_fortran_control_struct(const void* src, void* dst);

void* allocate_control_struct_container();
void reallocate_control_struct_container_data(void*, int, size_t) noexcept;
void deallocate_control_struct_container(void*) noexcept;
void access_control_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_control_var1_struct(int n, size_t* element_size);
void deallocate_fortran_control_var1_struct(void* ptr, int n) noexcept;
void copy_fortran_control_var1_struct(const void* src, void* dst);

void* allocate_control_var1_struct_container();
void reallocate_control_var1_struct_container_data(void*, int, size_t) noexcept;
void deallocate_control_var1_struct_container(void*) noexcept;
void access_control_var1_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_control_ramp1_struct(int n, size_t* element_size);
void deallocate_fortran_control_ramp1_struct(void* ptr, int n) noexcept;
void copy_fortran_control_ramp1_struct(const void* src, void* dst);

void* allocate_control_ramp1_struct_container();
void reallocate_control_ramp1_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_control_ramp1_struct_container(void*) noexcept;
void access_control_ramp1_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_controller_struct(int n, size_t* element_size);
void deallocate_fortran_controller_struct(void* ptr, int n) noexcept;
void copy_fortran_controller_struct(const void* src, void* dst);

void* allocate_controller_struct_container();
void reallocate_controller_struct_container_data(void*, int, size_t) noexcept;
void deallocate_controller_struct_container(void*) noexcept;
void access_controller_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_ellipse_beam_init_struct(int n, size_t* element_size);
void deallocate_fortran_ellipse_beam_init_struct(void* ptr, int n) noexcept;
void copy_fortran_ellipse_beam_init_struct(const void* src, void* dst);

void* allocate_ellipse_beam_init_struct_container();
void reallocate_ellipse_beam_init_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_ellipse_beam_init_struct_container(void*) noexcept;
void access_ellipse_beam_init_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_kv_beam_init_struct(int n, size_t* element_size);
void deallocate_fortran_kv_beam_init_struct(void* ptr, int n) noexcept;
void copy_fortran_kv_beam_init_struct(const void* src, void* dst);

void* allocate_kv_beam_init_struct_container();
void reallocate_kv_beam_init_struct_container_data(void*, int, size_t) noexcept;
void deallocate_kv_beam_init_struct_container(void*) noexcept;
void access_kv_beam_init_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_grid_beam_init_struct(int n, size_t* element_size);
void deallocate_fortran_grid_beam_init_struct(void* ptr, int n) noexcept;
void copy_fortran_grid_beam_init_struct(const void* src, void* dst);

void* allocate_grid_beam_init_struct_container();
void reallocate_grid_beam_init_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_grid_beam_init_struct_container(void*) noexcept;
void access_grid_beam_init_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_beam_init_struct(int n, size_t* element_size);
void deallocate_fortran_beam_init_struct(void* ptr, int n) noexcept;
void copy_fortran_beam_init_struct(const void* src, void* dst);

void* allocate_beam_init_struct_container();
void reallocate_beam_init_struct_container_data(void*, int, size_t) noexcept;
void deallocate_beam_init_struct_container(void*) noexcept;
void access_beam_init_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_lat_param_struct(int n, size_t* element_size);
void deallocate_fortran_lat_param_struct(void* ptr, int n) noexcept;
void copy_fortran_lat_param_struct(const void* src, void* dst);

void* allocate_lat_param_struct_container();
void reallocate_lat_param_struct_container_data(void*, int, size_t) noexcept;
void deallocate_lat_param_struct_container(void*) noexcept;
void access_lat_param_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_mode_info_struct(int n, size_t* element_size);
void deallocate_fortran_mode_info_struct(void* ptr, int n) noexcept;
void copy_fortran_mode_info_struct(const void* src, void* dst);

void* allocate_mode_info_struct_container();
void reallocate_mode_info_struct_container_data(void*, int, size_t) noexcept;
void deallocate_mode_info_struct_container(void*) noexcept;
void access_mode_info_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_pre_tracker_struct(int n, size_t* element_size);
void deallocate_fortran_pre_tracker_struct(void* ptr, int n) noexcept;
void copy_fortran_pre_tracker_struct(const void* src, void* dst);

void* allocate_pre_tracker_struct_container();
void reallocate_pre_tracker_struct_container_data(void*, int, size_t) noexcept;
void deallocate_pre_tracker_struct_container(void*) noexcept;
void access_pre_tracker_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_anormal_mode_struct(int n, size_t* element_size);
void deallocate_fortran_anormal_mode_struct(void* ptr, int n) noexcept;
void copy_fortran_anormal_mode_struct(const void* src, void* dst);

void* allocate_anormal_mode_struct_container();
void reallocate_anormal_mode_struct_container_data(void*, int, size_t) noexcept;
void deallocate_anormal_mode_struct_container(void*) noexcept;
void access_anormal_mode_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_linac_normal_mode_struct(int n, size_t* element_size);
void deallocate_fortran_linac_normal_mode_struct(void* ptr, int n) noexcept;
void copy_fortran_linac_normal_mode_struct(const void* src, void* dst);

void* allocate_linac_normal_mode_struct_container();
void reallocate_linac_normal_mode_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_linac_normal_mode_struct_container(void*) noexcept;
void access_linac_normal_mode_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_normal_modes_struct(int n, size_t* element_size);
void deallocate_fortran_normal_modes_struct(void* ptr, int n) noexcept;
void copy_fortran_normal_modes_struct(const void* src, void* dst);

void* allocate_normal_modes_struct_container();
void reallocate_normal_modes_struct_container_data(void*, int, size_t) noexcept;
void deallocate_normal_modes_struct_container(void*) noexcept;
void access_normal_modes_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_em_field_struct(int n, size_t* element_size);
void deallocate_fortran_em_field_struct(void* ptr, int n) noexcept;
void copy_fortran_em_field_struct(const void* src, void* dst);

void* allocate_em_field_struct_container();
void reallocate_em_field_struct_container_data(void*, int, size_t) noexcept;
void deallocate_em_field_struct_container(void*) noexcept;
void access_em_field_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_strong_beam_struct(int n, size_t* element_size);
void deallocate_fortran_strong_beam_struct(void* ptr, int n) noexcept;
void copy_fortran_strong_beam_struct(const void* src, void* dst);

void* allocate_strong_beam_struct_container();
void reallocate_strong_beam_struct_container_data(void*, int, size_t) noexcept;
void deallocate_strong_beam_struct_container(void*) noexcept;
void access_strong_beam_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_track_point_struct(int n, size_t* element_size);
void deallocate_fortran_track_point_struct(void* ptr, int n) noexcept;
void copy_fortran_track_point_struct(const void* src, void* dst);

void* allocate_track_point_struct_container();
void reallocate_track_point_struct_container_data(void*, int, size_t) noexcept;
void deallocate_track_point_struct_container(void*) noexcept;
void access_track_point_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_track_struct(int n, size_t* element_size);
void deallocate_fortran_track_struct(void* ptr, int n) noexcept;
void copy_fortran_track_struct(const void* src, void* dst);

void* allocate_track_struct_container();
void reallocate_track_struct_container_data(void*, int, size_t) noexcept;
void deallocate_track_struct_container(void*) noexcept;
void access_track_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_space_charge_common_struct(int n, size_t* element_size);
void deallocate_fortran_space_charge_common_struct(void* ptr, int n) noexcept;
void copy_fortran_space_charge_common_struct(const void* src, void* dst);

void* allocate_space_charge_common_struct_container();
void reallocate_space_charge_common_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_space_charge_common_struct_container(void*) noexcept;
void access_space_charge_common_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bmad_common_struct(int n, size_t* element_size);
void deallocate_fortran_bmad_common_struct(void* ptr, int n) noexcept;
void copy_fortran_bmad_common_struct(const void* src, void* dst);

void* allocate_bmad_common_struct_container();
void reallocate_bmad_common_struct_container_data(void*, int, size_t) noexcept;
void deallocate_bmad_common_struct_container(void*) noexcept;
void access_bmad_common_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_rad_int1_struct(int n, size_t* element_size);
void deallocate_fortran_rad_int1_struct(void* ptr, int n) noexcept;
void copy_fortran_rad_int1_struct(const void* src, void* dst);

void* allocate_rad_int1_struct_container();
void reallocate_rad_int1_struct_container_data(void*, int, size_t) noexcept;
void deallocate_rad_int1_struct_container(void*) noexcept;
void access_rad_int1_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_rad_int_branch_struct(int n, size_t* element_size);
void deallocate_fortran_rad_int_branch_struct(void* ptr, int n) noexcept;
void copy_fortran_rad_int_branch_struct(const void* src, void* dst);

void* allocate_rad_int_branch_struct_container();
void reallocate_rad_int_branch_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_rad_int_branch_struct_container(void*) noexcept;
void access_rad_int_branch_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_rad_int_all_ele_struct(int n, size_t* element_size);
void deallocate_fortran_rad_int_all_ele_struct(void* ptr, int n) noexcept;
void copy_fortran_rad_int_all_ele_struct(const void* src, void* dst);

void* allocate_rad_int_all_ele_struct_container();
void reallocate_rad_int_all_ele_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_rad_int_all_ele_struct_container(void*) noexcept;
void access_rad_int_all_ele_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_rf_stair_step_struct(int n, size_t* element_size);
void deallocate_fortran_rf_stair_step_struct(void* ptr, int n) noexcept;
void copy_fortran_rf_stair_step_struct(const void* src, void* dst);

void* allocate_rf_stair_step_struct_container();
void reallocate_rf_stair_step_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_rf_stair_step_struct_container(void*) noexcept;
void access_rf_stair_step_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_rf_ele_struct(int n, size_t* element_size);
void deallocate_fortran_rf_ele_struct(void* ptr, int n) noexcept;
void copy_fortran_rf_ele_struct(const void* src, void* dst);

void* allocate_rf_ele_struct_container();
void reallocate_rf_ele_struct_container_data(void*, int, size_t) noexcept;
void deallocate_rf_ele_struct_container(void*) noexcept;
void access_rf_ele_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_ele_struct(int n, size_t* element_size);
void deallocate_fortran_ele_struct(void* ptr, int n) noexcept;
void copy_fortran_ele_struct(const void* src, void* dst);

void* allocate_ele_struct_container();
void reallocate_ele_struct_container_data(void*, int, size_t) noexcept;
void deallocate_ele_struct_container(void*) noexcept;
void access_ele_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_complex_taylor_term_struct(int n, size_t* element_size);
void deallocate_fortran_complex_taylor_term_struct(void* ptr, int n) noexcept;
void copy_fortran_complex_taylor_term_struct(const void* src, void* dst);

void* allocate_complex_taylor_term_struct_container();
void reallocate_complex_taylor_term_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_complex_taylor_term_struct_container(void*) noexcept;
void access_complex_taylor_term_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_complex_taylor_struct(int n, size_t* element_size);
void deallocate_fortran_complex_taylor_struct(void* ptr, int n) noexcept;
void copy_fortran_complex_taylor_struct(const void* src, void* dst);

void* allocate_complex_taylor_struct_container();
void reallocate_complex_taylor_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_complex_taylor_struct_container(void*) noexcept;
void access_complex_taylor_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_branch_struct(int n, size_t* element_size);
void deallocate_fortran_branch_struct(void* ptr, int n) noexcept;
void copy_fortran_branch_struct(const void* src, void* dst);

void* allocate_branch_struct_container();
void reallocate_branch_struct_container_data(void*, int, size_t) noexcept;
void deallocate_branch_struct_container(void*) noexcept;
void access_branch_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_lat_struct(int n, size_t* element_size);
void deallocate_fortran_lat_struct(void* ptr, int n) noexcept;
void copy_fortran_lat_struct(const void* src, void* dst);

void* allocate_lat_struct_container();
void reallocate_lat_struct_container_data(void*, int, size_t) noexcept;
void deallocate_lat_struct_container(void*) noexcept;
void access_lat_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bunch_struct(int n, size_t* element_size);
void deallocate_fortran_bunch_struct(void* ptr, int n) noexcept;
void copy_fortran_bunch_struct(const void* src, void* dst);

void* allocate_bunch_struct_container();
void reallocate_bunch_struct_container_data(void*, int, size_t) noexcept;
void deallocate_bunch_struct_container(void*) noexcept;
void access_bunch_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bunch_params_struct(int n, size_t* element_size);
void deallocate_fortran_bunch_params_struct(void* ptr, int n) noexcept;
void copy_fortran_bunch_params_struct(const void* src, void* dst);

void* allocate_bunch_params_struct_container();
void reallocate_bunch_params_struct_container_data(void*, int, size_t) noexcept;
void deallocate_bunch_params_struct_container(void*) noexcept;
void access_bunch_params_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_beam_struct(int n, size_t* element_size);
void deallocate_fortran_beam_struct(void* ptr, int n) noexcept;
void copy_fortran_beam_struct(const void* src, void* dst);

void* allocate_beam_struct_container();
void reallocate_beam_struct_container_data(void*, int, size_t) noexcept;
void deallocate_beam_struct_container(void*) noexcept;
void access_beam_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_aperture_point_struct(int n, size_t* element_size);
void deallocate_fortran_aperture_point_struct(void* ptr, int n) noexcept;
void copy_fortran_aperture_point_struct(const void* src, void* dst);

void* allocate_aperture_point_struct_container();
void reallocate_aperture_point_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_aperture_point_struct_container(void*) noexcept;
void access_aperture_point_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_aperture_param_struct(int n, size_t* element_size);
void deallocate_fortran_aperture_param_struct(void* ptr, int n) noexcept;
void copy_fortran_aperture_param_struct(const void* src, void* dst);

void* allocate_aperture_param_struct_container();
void reallocate_aperture_param_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_aperture_param_struct_container(void*) noexcept;
void access_aperture_param_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_aperture_scan_struct(int n, size_t* element_size);
void deallocate_fortran_aperture_scan_struct(void* ptr, int n) noexcept;
void copy_fortran_aperture_scan_struct(const void* src, void* dst);

void* allocate_aperture_scan_struct_container();
void reallocate_aperture_scan_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_aperture_scan_struct_container(void*) noexcept;
void access_aperture_scan_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_ele_pointer_struct(int n, size_t* element_size);
void deallocate_fortran_ele_pointer_struct(void* ptr, int n) noexcept;
void copy_fortran_ele_pointer_struct(const void* src, void* dst);

void* allocate_ele_pointer_struct_container();
void reallocate_ele_pointer_struct_container_data(void*, int, size_t) noexcept;
void deallocate_ele_pointer_struct_container(void*) noexcept;
void access_ele_pointer_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_expression_tree_struct(int n, size_t* element_size);
void deallocate_fortran_expression_tree_struct(void* ptr, int n) noexcept;
void copy_fortran_expression_tree_struct(const void* src, void* dst);

void* allocate_expression_tree_struct_container();
void reallocate_expression_tree_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_expression_tree_struct_container(void*) noexcept;
void access_expression_tree_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_nametable_struct(int n, size_t* element_size);
void deallocate_fortran_nametable_struct(void* ptr, int n) noexcept;
void copy_fortran_nametable_struct(const void* src, void* dst);

void* allocate_nametable_struct_container();
void reallocate_nametable_struct_container_data(void*, int, size_t) noexcept;
void deallocate_nametable_struct_container(void*) noexcept;
void access_nametable_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_spin_dn_dpz_struct(int n, size_t* element_size);
void deallocate_fortran_tao_spin_dn_dpz_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_spin_dn_dpz_struct(const void* src, void* dst);

void* allocate_tao_spin_dn_dpz_struct_container();
void reallocate_tao_spin_dn_dpz_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_spin_dn_dpz_struct_container(void*) noexcept;
void access_tao_spin_dn_dpz_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_resonance_h_struct(int n, size_t* element_size);
void deallocate_fortran_resonance_h_struct(void* ptr, int n) noexcept;
void copy_fortran_resonance_h_struct(const void* src, void* dst);

void* allocate_resonance_h_struct_container();
void reallocate_resonance_h_struct_container_data(void*, int, size_t) noexcept;
void deallocate_resonance_h_struct_container(void*) noexcept;
void access_resonance_h_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_spin_orbit_map1_struct(int n, size_t* element_size);
void deallocate_fortran_spin_orbit_map1_struct(void* ptr, int n) noexcept;
void copy_fortran_spin_orbit_map1_struct(const void* src, void* dst);

void* allocate_spin_orbit_map1_struct_container();
void reallocate_spin_orbit_map1_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_spin_orbit_map1_struct_container(void*) noexcept;
void access_spin_orbit_map1_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_spin_axis_struct(int n, size_t* element_size);
void deallocate_fortran_spin_axis_struct(void* ptr, int n) noexcept;
void copy_fortran_spin_axis_struct(const void* src, void* dst);

void* allocate_spin_axis_struct_container();
void reallocate_spin_axis_struct_container_data(void*, int, size_t) noexcept;
void deallocate_spin_axis_struct_container(void*) noexcept;
void access_spin_axis_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_ptc_normal_form_struct(int n, size_t* element_size);
void deallocate_fortran_ptc_normal_form_struct(void* ptr, int n) noexcept;
void copy_fortran_ptc_normal_form_struct(const void* src, void* dst);

void* allocate_ptc_normal_form_struct_container();
void reallocate_ptc_normal_form_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_ptc_normal_form_struct_container(void*) noexcept;
void access_ptc_normal_form_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bmad_normal_form_struct(int n, size_t* element_size);
void deallocate_fortran_bmad_normal_form_struct(void* ptr, int n) noexcept;
void copy_fortran_bmad_normal_form_struct(const void* src, void* dst);

void* allocate_bmad_normal_form_struct_container();
void reallocate_bmad_normal_form_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_bmad_normal_form_struct_container(void*) noexcept;
void access_bmad_normal_form_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bunch_track_struct(int n, size_t* element_size);
void deallocate_fortran_bunch_track_struct(void* ptr, int n) noexcept;
void copy_fortran_bunch_track_struct(const void* src, void* dst);

void* allocate_bunch_track_struct_container();
void reallocate_bunch_track_struct_container_data(void*, int, size_t) noexcept;
void deallocate_bunch_track_struct_container(void*) noexcept;
void access_bunch_track_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_summation_rdt_struct(int n, size_t* element_size);
void deallocate_fortran_summation_rdt_struct(void* ptr, int n) noexcept;
void copy_fortran_summation_rdt_struct(const void* src, void* dst);

void* allocate_summation_rdt_struct_container();
void reallocate_summation_rdt_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_summation_rdt_struct_container(void*) noexcept;
void access_summation_rdt_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_ele_shape_struct(int n, size_t* element_size);
void deallocate_fortran_tao_ele_shape_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_ele_shape_struct(const void* src, void* dst);

void* allocate_tao_ele_shape_struct_container();
void reallocate_tao_ele_shape_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_ele_shape_struct_container(void*) noexcept;
void access_tao_ele_shape_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_ele_pointer_struct(int n, size_t* element_size);
void deallocate_fortran_tao_ele_pointer_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_ele_pointer_struct(const void* src, void* dst);

void* allocate_tao_ele_pointer_struct_container();
void reallocate_tao_ele_pointer_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_ele_pointer_struct_container(void*) noexcept;
void access_tao_ele_pointer_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_curve_struct(int n, size_t* element_size);
void deallocate_fortran_tao_curve_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_curve_struct(const void* src, void* dst);

void* allocate_tao_curve_struct_container();
void reallocate_tao_curve_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_curve_struct_container(void*) noexcept;
void access_tao_curve_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_curve_color_struct(int n, size_t* element_size);
void deallocate_fortran_tao_curve_color_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_curve_color_struct(const void* src, void* dst);

void* allocate_tao_curve_color_struct_container();
void reallocate_tao_curve_color_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_curve_color_struct_container(void*) noexcept;
void access_tao_curve_color_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_curve_orbit_struct(int n, size_t* element_size);
void deallocate_fortran_tao_curve_orbit_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_curve_orbit_struct(const void* src, void* dst);

void* allocate_tao_curve_orbit_struct_container();
void reallocate_tao_curve_orbit_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_curve_orbit_struct_container(void*) noexcept;
void access_tao_curve_orbit_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_histogram_struct(int n, size_t* element_size);
void deallocate_fortran_tao_histogram_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_histogram_struct(const void* src, void* dst);

void* allocate_tao_histogram_struct_container();
void reallocate_tao_histogram_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_histogram_struct_container(void*) noexcept;
void access_tao_histogram_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_lat_ele_order1_struct(int n, size_t* element_size);
void deallocate_fortran_lat_ele_order1_struct(void* ptr, int n) noexcept;
void copy_fortran_lat_ele_order1_struct(const void* src, void* dst);

void* allocate_lat_ele_order1_struct_container();
void reallocate_lat_ele_order1_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_lat_ele_order1_struct_container(void*) noexcept;
void access_lat_ele_order1_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_lat_ele_order_array_struct(int n, size_t* element_size);
void deallocate_fortran_lat_ele_order_array_struct(void* ptr, int n) noexcept;
void copy_fortran_lat_ele_order_array_struct(const void* src, void* dst);

void* allocate_lat_ele_order_array_struct_container();
void reallocate_lat_ele_order_array_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_lat_ele_order_array_struct_container(void*) noexcept;
void access_lat_ele_order_array_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_lat_sigma_struct(int n, size_t* element_size);
void deallocate_fortran_tao_lat_sigma_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_lat_sigma_struct(const void* src, void* dst);

void* allocate_tao_lat_sigma_struct_container();
void reallocate_tao_lat_sigma_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_lat_sigma_struct_container(void*) noexcept;
void access_tao_lat_sigma_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_spin_ele_struct(int n, size_t* element_size);
void deallocate_fortran_tao_spin_ele_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_spin_ele_struct(const void* src, void* dst);

void* allocate_tao_spin_ele_struct_container();
void reallocate_tao_spin_ele_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_spin_ele_struct_container(void*) noexcept;
void access_tao_spin_ele_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_plot_cache_struct(int n, size_t* element_size);
void deallocate_fortran_tao_plot_cache_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_plot_cache_struct(const void* src, void* dst);

void* allocate_tao_plot_cache_struct_container();
void reallocate_tao_plot_cache_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_plot_cache_struct_container(void*) noexcept;
void access_tao_plot_cache_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_spin_polarization_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_tao_spin_polarization_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_spin_polarization_struct(const void* src, void* dst);

void* allocate_tao_spin_polarization_struct_container();
void reallocate_tao_spin_polarization_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_spin_polarization_struct_container(void*) noexcept;
void access_tao_spin_polarization_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_lattice_branch_struct(int n, size_t* element_size);
void deallocate_fortran_tao_lattice_branch_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_lattice_branch_struct(const void* src, void* dst);

void* allocate_tao_lattice_branch_struct_container();
void reallocate_tao_lattice_branch_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_lattice_branch_struct_container(void*) noexcept;
void access_tao_lattice_branch_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_model_element_struct(int n, size_t* element_size);
void deallocate_fortran_tao_model_element_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_model_element_struct(const void* src, void* dst);

void* allocate_tao_model_element_struct_container();
void reallocate_tao_model_element_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_model_element_struct_container(void*) noexcept;
void access_tao_model_element_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_beam_branch_struct(int n, size_t* element_size);
void deallocate_fortran_tao_beam_branch_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_beam_branch_struct(const void* src, void* dst);

void* allocate_tao_beam_branch_struct_container();
void reallocate_tao_beam_branch_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_beam_branch_struct_container(void*) noexcept;
void access_tao_beam_branch_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_d1_data_struct(int n, size_t* element_size);
void deallocate_fortran_tao_d1_data_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_d1_data_struct(const void* src, void* dst);

void* allocate_tao_d1_data_struct_container();
void reallocate_tao_d1_data_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_d1_data_struct_container(void*) noexcept;
void access_tao_d1_data_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_d2_data_struct(int n, size_t* element_size);
void deallocate_fortran_tao_d2_data_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_d2_data_struct(const void* src, void* dst);

void* allocate_tao_d2_data_struct_container();
void reallocate_tao_d2_data_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_d2_data_struct_container(void*) noexcept;
void access_tao_d2_data_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_data_var_component_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_tao_data_var_component_struct(
    void* ptr,
    int n) noexcept;
void copy_fortran_tao_data_var_component_struct(const void* src, void* dst);

void* allocate_tao_data_var_component_struct_container();
void reallocate_tao_data_var_component_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_data_var_component_struct_container(void*) noexcept;
void access_tao_data_var_component_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_graph_struct(int n, size_t* element_size);
void deallocate_fortran_tao_graph_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_graph_struct(const void* src, void* dst);

void* allocate_tao_graph_struct_container();
void reallocate_tao_graph_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_graph_struct_container(void*) noexcept;
void access_tao_graph_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_plot_struct(int n, size_t* element_size);
void deallocate_fortran_tao_plot_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_plot_struct(const void* src, void* dst);

void* allocate_tao_plot_struct_container();
void reallocate_tao_plot_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_plot_struct_container(void*) noexcept;
void access_tao_plot_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_plot_region_struct(int n, size_t* element_size);
void deallocate_fortran_tao_plot_region_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_plot_region_struct(const void* src, void* dst);

void* allocate_tao_plot_region_struct_container();
void reallocate_tao_plot_region_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_plot_region_struct_container(void*) noexcept;
void access_tao_plot_region_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_universe_pointer_struct(int n, size_t* element_size);
void deallocate_fortran_tao_universe_pointer_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_universe_pointer_struct(const void* src, void* dst);

void* allocate_tao_universe_pointer_struct_container();
void reallocate_tao_universe_pointer_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_universe_pointer_struct_container(void*) noexcept;
void access_tao_universe_pointer_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_super_universe_struct(int n, size_t* element_size);
void deallocate_fortran_tao_super_universe_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_super_universe_struct(const void* src, void* dst);

void* allocate_tao_super_universe_struct_container();
void reallocate_tao_super_universe_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_super_universe_struct_container(void*) noexcept;
void access_tao_super_universe_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_var_struct(int n, size_t* element_size);
void deallocate_fortran_tao_var_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_var_struct(const void* src, void* dst);

void* allocate_tao_var_struct_container();
void reallocate_tao_var_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_var_struct_container(void*) noexcept;
void access_tao_var_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_var_slave_struct(int n, size_t* element_size);
void deallocate_fortran_tao_var_slave_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_var_slave_struct(const void* src, void* dst);

void* allocate_tao_var_slave_struct_container();
void reallocate_tao_var_slave_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_var_slave_struct_container(void*) noexcept;
void access_tao_var_slave_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_lattice_struct(int n, size_t* element_size);
void deallocate_fortran_tao_lattice_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_lattice_struct(const void* src, void* dst);

void* allocate_tao_lattice_struct_container();
void reallocate_tao_lattice_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_lattice_struct_container(void*) noexcept;
void access_tao_lattice_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_beam_uni_struct(int n, size_t* element_size);
void deallocate_fortran_tao_beam_uni_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_beam_uni_struct(const void* src, void* dst);

void* allocate_tao_beam_uni_struct_container();
void reallocate_tao_beam_uni_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_beam_uni_struct_container(void*) noexcept;
void access_tao_beam_uni_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_dynamic_aperture_struct(int n, size_t* element_size);
void deallocate_fortran_tao_dynamic_aperture_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_dynamic_aperture_struct(const void* src, void* dst);

void* allocate_tao_dynamic_aperture_struct_container();
void reallocate_tao_dynamic_aperture_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_dynamic_aperture_struct_container(void*) noexcept;
void access_tao_dynamic_aperture_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_model_branch_struct(int n, size_t* element_size);
void deallocate_fortran_tao_model_branch_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_model_branch_struct(const void* src, void* dst);

void* allocate_tao_model_branch_struct_container();
void reallocate_tao_model_branch_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_model_branch_struct_container(void*) noexcept;
void access_tao_model_branch_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_spin_map_struct(int n, size_t* element_size);
void deallocate_fortran_tao_spin_map_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_spin_map_struct(const void* src, void* dst);

void* allocate_tao_spin_map_struct_container();
void reallocate_tao_spin_map_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_spin_map_struct_container(void*) noexcept;
void access_tao_spin_map_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_data_struct(int n, size_t* element_size);
void deallocate_fortran_tao_data_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_data_struct(const void* src, void* dst);

void* allocate_tao_data_struct_container();
void reallocate_tao_data_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_data_struct_container(void*) noexcept;
void access_tao_data_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_ping_scale_struct(int n, size_t* element_size);
void deallocate_fortran_tao_ping_scale_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_ping_scale_struct(const void* src, void* dst);

void* allocate_tao_ping_scale_struct_container();
void reallocate_tao_ping_scale_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_ping_scale_struct_container(void*) noexcept;
void access_tao_ping_scale_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_universe_calc_struct(int n, size_t* element_size);
void deallocate_fortran_tao_universe_calc_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_universe_calc_struct(const void* src, void* dst);

void* allocate_tao_universe_calc_struct_container();
void reallocate_tao_universe_calc_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_universe_calc_struct_container(void*) noexcept;
void access_tao_universe_calc_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_lat_ele_order_struct(int n, size_t* element_size);
void deallocate_fortran_lat_ele_order_struct(void* ptr, int n) noexcept;
void copy_fortran_lat_ele_order_struct(const void* src, void* dst);

void* allocate_lat_ele_order_struct_container();
void reallocate_lat_ele_order_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_lat_ele_order_struct_container(void*) noexcept;
void access_lat_ele_order_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_title_struct(int n, size_t* element_size);
void deallocate_fortran_tao_title_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_title_struct(const void* src, void* dst);

void* allocate_tao_title_struct_container();
void reallocate_tao_title_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_title_struct_container(void*) noexcept;
void access_tao_title_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_qp_rect_struct(int n, size_t* element_size);
void deallocate_fortran_qp_rect_struct(void* ptr, int n) noexcept;
void copy_fortran_qp_rect_struct(const void* src, void* dst);

void* allocate_qp_rect_struct_container();
void reallocate_qp_rect_struct_container_data(void*, int, size_t) noexcept;
void deallocate_qp_rect_struct_container(void*) noexcept;
void access_qp_rect_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_drawing_struct(int n, size_t* element_size);
void deallocate_fortran_tao_drawing_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_drawing_struct(const void* src, void* dst);

void* allocate_tao_drawing_struct_container();
void reallocate_tao_drawing_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_drawing_struct_container(void*) noexcept;
void access_tao_drawing_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_shape_pattern_struct(int n, size_t* element_size);
void deallocate_fortran_tao_shape_pattern_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_shape_pattern_struct(const void* src, void* dst);

void* allocate_tao_shape_pattern_struct_container();
void reallocate_tao_shape_pattern_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_shape_pattern_struct_container(void*) noexcept;
void access_tao_shape_pattern_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_shape_pattern_point_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_tao_shape_pattern_point_struct(
    void* ptr,
    int n) noexcept;
void copy_fortran_tao_shape_pattern_point_struct(const void* src, void* dst);

void* allocate_tao_shape_pattern_point_struct_container();
void reallocate_tao_shape_pattern_point_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_shape_pattern_point_struct_container(void*) noexcept;
void access_tao_shape_pattern_point_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_qp_axis_struct(int n, size_t* element_size);
void deallocate_fortran_qp_axis_struct(void* ptr, int n) noexcept;
void copy_fortran_qp_axis_struct(const void* src, void* dst);

void* allocate_qp_axis_struct_container();
void reallocate_qp_axis_struct_container_data(void*, int, size_t) noexcept;
void deallocate_qp_axis_struct_container(void*) noexcept;
void access_qp_axis_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_qp_legend_struct(int n, size_t* element_size);
void deallocate_fortran_qp_legend_struct(void* ptr, int n) noexcept;
void copy_fortran_qp_legend_struct(const void* src, void* dst);

void* allocate_qp_legend_struct_container();
void reallocate_qp_legend_struct_container_data(void*, int, size_t) noexcept;
void deallocate_qp_legend_struct_container(void*) noexcept;
void access_qp_legend_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_qp_point_struct(int n, size_t* element_size);
void deallocate_fortran_qp_point_struct(void* ptr, int n) noexcept;
void copy_fortran_qp_point_struct(const void* src, void* dst);

void* allocate_qp_point_struct_container();
void reallocate_qp_point_struct_container_data(void*, int, size_t) noexcept;
void deallocate_qp_point_struct_container(void*) noexcept;
void access_qp_point_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_qp_line_struct(int n, size_t* element_size);
void deallocate_fortran_qp_line_struct(void* ptr, int n) noexcept;
void copy_fortran_qp_line_struct(const void* src, void* dst);

void* allocate_qp_line_struct_container();
void reallocate_qp_line_struct_container_data(void*, int, size_t) noexcept;
void deallocate_qp_line_struct_container(void*) noexcept;
void access_qp_line_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_qp_symbol_struct(int n, size_t* element_size);
void deallocate_fortran_qp_symbol_struct(void* ptr, int n) noexcept;
void copy_fortran_qp_symbol_struct(const void* src, void* dst);

void* allocate_qp_symbol_struct_container();
void reallocate_qp_symbol_struct_container_data(void*, int, size_t) noexcept;
void deallocate_qp_symbol_struct_container(void*) noexcept;
void access_qp_symbol_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_floor_plan_struct(int n, size_t* element_size);
void deallocate_fortran_tao_floor_plan_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_floor_plan_struct(const void* src, void* dst);

void* allocate_tao_floor_plan_struct_container();
void reallocate_tao_floor_plan_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_floor_plan_struct_container(void*) noexcept;
void access_tao_floor_plan_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_v1_var_struct(int n, size_t* element_size);
void deallocate_fortran_tao_v1_var_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_v1_var_struct(const void* src, void* dst);

void* allocate_tao_v1_var_struct_container();
void reallocate_tao_v1_var_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_v1_var_struct_container(void*) noexcept;
void access_tao_v1_var_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_global_struct(int n, size_t* element_size);
void deallocate_fortran_tao_global_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_global_struct(const void* src, void* dst);

void* allocate_tao_global_struct_container();
void reallocate_tao_global_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_global_struct_container(void*) noexcept;
void access_tao_global_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_init_struct(int n, size_t* element_size);
void deallocate_fortran_tao_init_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_init_struct(const void* src, void* dst);

void* allocate_tao_init_struct_container();
void reallocate_tao_init_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_init_struct_container(void*) noexcept;
void access_tao_init_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_common_struct(int n, size_t* element_size);
void deallocate_fortran_tao_common_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_common_struct(const void* src, void* dst);

void* allocate_tao_common_struct_container();
void reallocate_tao_common_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_common_struct_container(void*) noexcept;
void access_tao_common_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_plot_page_struct(int n, size_t* element_size);
void deallocate_fortran_tao_plot_page_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_plot_page_struct(const void* src, void* dst);

void* allocate_tao_plot_page_struct_container();
void reallocate_tao_plot_page_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_plot_page_struct_container(void*) noexcept;
void access_tao_plot_page_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_building_wall_struct(int n, size_t* element_size);
void deallocate_fortran_tao_building_wall_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_building_wall_struct(const void* src, void* dst);

void* allocate_tao_building_wall_struct_container();
void reallocate_tao_building_wall_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_building_wall_struct_container(void*) noexcept;
void access_tao_building_wall_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_building_wall_orientation_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_tao_building_wall_orientation_struct(
    void* ptr,
    int n) noexcept;
void copy_fortran_tao_building_wall_orientation_struct(
    const void* src,
    void* dst);

void* allocate_tao_building_wall_orientation_struct_container();
void reallocate_tao_building_wall_orientation_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_building_wall_orientation_struct_container(void*) noexcept;
void access_tao_building_wall_orientation_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_building_wall_section_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_tao_building_wall_section_struct(
    void* ptr,
    int n) noexcept;
void copy_fortran_tao_building_wall_section_struct(const void* src, void* dst);

void* allocate_tao_building_wall_section_struct_container();
void reallocate_tao_building_wall_section_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_building_wall_section_struct_container(void*) noexcept;
void access_tao_building_wall_section_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_building_wall_point_struct(
    int n,
    size_t* element_size);
void deallocate_fortran_tao_building_wall_point_struct(
    void* ptr,
    int n) noexcept;
void copy_fortran_tao_building_wall_point_struct(const void* src, void* dst);

void* allocate_tao_building_wall_point_struct_container();
void reallocate_tao_building_wall_point_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_building_wall_point_struct_container(void*) noexcept;
void access_tao_building_wall_point_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_wave_struct(int n, size_t* element_size);
void deallocate_fortran_tao_wave_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_wave_struct(const void* src, void* dst);

void* allocate_tao_wave_struct_container();
void reallocate_tao_wave_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_wave_struct_container(void*) noexcept;
void access_tao_wave_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_wave_kick_pt_struct(int n, size_t* element_size);
void deallocate_fortran_tao_wave_kick_pt_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_wave_kick_pt_struct(const void* src, void* dst);

void* allocate_tao_wave_kick_pt_struct_container();
void reallocate_tao_wave_kick_pt_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_wave_kick_pt_struct_container(void*) noexcept;
void access_tao_wave_kick_pt_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_cmd_history_struct(int n, size_t* element_size);
void deallocate_fortran_tao_cmd_history_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_cmd_history_struct(const void* src, void* dst);

void* allocate_tao_cmd_history_struct_container();
void reallocate_tao_cmd_history_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_tao_cmd_history_struct_container(void*) noexcept;
void access_tao_cmd_history_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_tao_universe_struct(int n, size_t* element_size);
void deallocate_fortran_tao_universe_struct(void* ptr, int n) noexcept;
void copy_fortran_tao_universe_struct(const void* src, void* dst);

void* allocate_tao_universe_struct_container();
void reallocate_tao_universe_struct_container_data(void*, int, size_t) noexcept;
void deallocate_tao_universe_struct_container(void*) noexcept;
void access_tao_universe_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_mad_energy_struct(int n, size_t* element_size);
void deallocate_fortran_mad_energy_struct(void* ptr, int n) noexcept;
void copy_fortran_mad_energy_struct(const void* src, void* dst);

void* allocate_mad_energy_struct_container();
void reallocate_mad_energy_struct_container_data(void*, int, size_t) noexcept;
void deallocate_mad_energy_struct_container(void*) noexcept;
void access_mad_energy_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_mad_map_struct(int n, size_t* element_size);
void deallocate_fortran_mad_map_struct(void* ptr, int n) noexcept;
void copy_fortran_mad_map_struct(const void* src, void* dst);

void* allocate_mad_map_struct_container();
void reallocate_mad_map_struct_container_data(void*, int, size_t) noexcept;
void deallocate_mad_map_struct_container(void*) noexcept;
void access_mad_map_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_random_state_struct(int n, size_t* element_size);
void deallocate_fortran_random_state_struct(void* ptr, int n) noexcept;
void copy_fortran_random_state_struct(const void* src, void* dst);

void* allocate_random_state_struct_container();
void reallocate_random_state_struct_container_data(void*, int, size_t) noexcept;
void deallocate_random_state_struct_container(void*) noexcept;
void access_random_state_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bbu_stage_struct(int n, size_t* element_size);
void deallocate_fortran_bbu_stage_struct(void* ptr, int n) noexcept;
void copy_fortran_bbu_stage_struct(const void* src, void* dst);

void* allocate_bbu_stage_struct_container();
void reallocate_bbu_stage_struct_container_data(void*, int, size_t) noexcept;
void deallocate_bbu_stage_struct_container(void*) noexcept;
void access_bbu_stage_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bbu_beam_struct(int n, size_t* element_size);
void deallocate_fortran_bbu_beam_struct(void* ptr, int n) noexcept;
void copy_fortran_bbu_beam_struct(const void* src, void* dst);

void* allocate_bbu_beam_struct_container();
void reallocate_bbu_beam_struct_container_data(void*, int, size_t) noexcept;
void deallocate_bbu_beam_struct_container(void*) noexcept;
void access_bbu_beam_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_bbu_param_struct(int n, size_t* element_size);
void deallocate_fortran_bbu_param_struct(void* ptr, int n) noexcept;
void copy_fortran_bbu_param_struct(const void* src, void* dst);

void* allocate_bbu_param_struct_container();
void reallocate_bbu_param_struct_container_data(void*, int, size_t) noexcept;
void deallocate_bbu_param_struct_container(void*) noexcept;
void access_bbu_param_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_all_encompassing_struct(int n, size_t* element_size);
void deallocate_fortran_all_encompassing_struct(void* ptr, int n) noexcept;
void copy_fortran_all_encompassing_struct(const void* src, void* dst);

void* allocate_all_encompassing_struct_container();
void reallocate_all_encompassing_struct_container_data(
    void*,
    int,
    size_t) noexcept;
void deallocate_all_encompassing_struct_container(void*) noexcept;
void access_all_encompassing_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_test_sub_struct(int n, size_t* element_size);
void deallocate_fortran_test_sub_struct(void* ptr, int n) noexcept;
void copy_fortran_test_sub_struct(const void* src, void* dst);

void* allocate_test_sub_struct_container();
void reallocate_test_sub_struct_container_data(void*, int, size_t) noexcept;
void deallocate_test_sub_struct_container(void*) noexcept;
void access_test_sub_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);

void* allocate_fortran_test_sub_sub_struct(int n, size_t* element_size);
void deallocate_fortran_test_sub_sub_struct(void* ptr, int n) noexcept;
void copy_fortran_test_sub_sub_struct(const void* src, void* dst);

void* allocate_test_sub_sub_struct_container();
void reallocate_test_sub_sub_struct_container_data(void*, int, size_t) noexcept;
void deallocate_test_sub_sub_struct_container(void*) noexcept;
void access_test_sub_sub_struct_container(
    void* handle,
    void** data,
    int* lbound,
    int* size,
    size_t* elem_size,
    bool* alloc);
}

using RealAlloc1D = FAlloc1D<
    double,
    allocate_real_container,
    deallocate_real_container,
    reallocate_real_container_data,
    access_real_container>;

using Real16Alloc1D = FAlloc1D<
    long double,
    allocate_real16_container,
    deallocate_real16_container,
    reallocate_real16_container_data,
    access_real16_container>;

using IntAlloc1D = FAlloc1D<
    int,
    allocate_integer_container,
    deallocate_integer_container,
    reallocate_integer_container_data,
    access_integer_container>;

using Int8Alloc1D = FAlloc1D<
    int64_t,
    allocate_integer8_container,
    deallocate_integer8_container,
    reallocate_integer8_container_data,
    access_integer8_container>;

using BoolAlloc1D = FAlloc1D<
    bool,
    allocate_logical_container,
    deallocate_logical_container,
    reallocate_logical_container_data,
    access_logical_container>;

using ComplexAlloc1D = FAlloc1D<
    std::complex<double>,
    allocate_complex_container,
    deallocate_complex_container,
    reallocate_complex_container_data,
    access_complex_container>;

class SplineProxy;

using SplineProxyArray1D = FTypeArray1D<
    SplineProxy,
    allocate_fortran_spline_struct,
    deallocate_fortran_spline_struct>;
using SplineProxyArray2D = FTypeArray2D<SplineProxy>;
using SplineProxyArray3D = FTypeArray3D<SplineProxy>;

using SplineProxyAlloc1D = FTypeAlloc1D<
    SplineProxyArray1D,
    allocate_spline_struct_container,
    deallocate_spline_struct_container,
    reallocate_spline_struct_container_data,
    access_spline_struct_container>;

class SpinPolarProxy;

using SpinPolarProxyArray1D = FTypeArray1D<
    SpinPolarProxy,
    allocate_fortran_spin_polar_struct,
    deallocate_fortran_spin_polar_struct>;
using SpinPolarProxyArray2D = FTypeArray2D<SpinPolarProxy>;
using SpinPolarProxyArray3D = FTypeArray3D<SpinPolarProxy>;

using SpinPolarProxyAlloc1D = FTypeAlloc1D<
    SpinPolarProxyArray1D,
    allocate_spin_polar_struct_container,
    deallocate_spin_polar_struct_container,
    reallocate_spin_polar_struct_container_data,
    access_spin_polar_struct_container>;

class AcKickerTimeProxy;

using AcKickerTimeProxyArray1D = FTypeArray1D<
    AcKickerTimeProxy,
    allocate_fortran_ac_kicker_time_struct,
    deallocate_fortran_ac_kicker_time_struct>;
using AcKickerTimeProxyArray2D = FTypeArray2D<AcKickerTimeProxy>;
using AcKickerTimeProxyArray3D = FTypeArray3D<AcKickerTimeProxy>;

using AcKickerTimeProxyAlloc1D = FTypeAlloc1D<
    AcKickerTimeProxyArray1D,
    allocate_ac_kicker_time_struct_container,
    deallocate_ac_kicker_time_struct_container,
    reallocate_ac_kicker_time_struct_container_data,
    access_ac_kicker_time_struct_container>;

class AcKickerFreqProxy;

using AcKickerFreqProxyArray1D = FTypeArray1D<
    AcKickerFreqProxy,
    allocate_fortran_ac_kicker_freq_struct,
    deallocate_fortran_ac_kicker_freq_struct>;
using AcKickerFreqProxyArray2D = FTypeArray2D<AcKickerFreqProxy>;
using AcKickerFreqProxyArray3D = FTypeArray3D<AcKickerFreqProxy>;

using AcKickerFreqProxyAlloc1D = FTypeAlloc1D<
    AcKickerFreqProxyArray1D,
    allocate_ac_kicker_freq_struct_container,
    deallocate_ac_kicker_freq_struct_container,
    reallocate_ac_kicker_freq_struct_container_data,
    access_ac_kicker_freq_struct_container>;

class AcKickerProxy;

using AcKickerProxyArray1D = FTypeArray1D<
    AcKickerProxy,
    allocate_fortran_ac_kicker_struct,
    deallocate_fortran_ac_kicker_struct>;
using AcKickerProxyArray2D = FTypeArray2D<AcKickerProxy>;
using AcKickerProxyArray3D = FTypeArray3D<AcKickerProxy>;

using AcKickerProxyAlloc1D = FTypeAlloc1D<
    AcKickerProxyArray1D,
    allocate_ac_kicker_struct_container,
    deallocate_ac_kicker_struct_container,
    reallocate_ac_kicker_struct_container_data,
    access_ac_kicker_struct_container>;

class Interval1CoefProxy;

using Interval1CoefProxyArray1D = FTypeArray1D<
    Interval1CoefProxy,
    allocate_fortran_interval1_coef_struct,
    deallocate_fortran_interval1_coef_struct>;
using Interval1CoefProxyArray2D = FTypeArray2D<Interval1CoefProxy>;
using Interval1CoefProxyArray3D = FTypeArray3D<Interval1CoefProxy>;

using Interval1CoefProxyAlloc1D = FTypeAlloc1D<
    Interval1CoefProxyArray1D,
    allocate_interval1_coef_struct_container,
    deallocate_interval1_coef_struct_container,
    reallocate_interval1_coef_struct_container_data,
    access_interval1_coef_struct_container>;

class PhotonReflectTableProxy;

using PhotonReflectTableProxyArray1D = FTypeArray1D<
    PhotonReflectTableProxy,
    allocate_fortran_photon_reflect_table_struct,
    deallocate_fortran_photon_reflect_table_struct>;
using PhotonReflectTableProxyArray2D = FTypeArray2D<PhotonReflectTableProxy>;
using PhotonReflectTableProxyArray3D = FTypeArray3D<PhotonReflectTableProxy>;

using PhotonReflectTableProxyAlloc1D = FTypeAlloc1D<
    PhotonReflectTableProxyArray1D,
    allocate_photon_reflect_table_struct_container,
    deallocate_photon_reflect_table_struct_container,
    reallocate_photon_reflect_table_struct_container_data,
    access_photon_reflect_table_struct_container>;

class PhotonReflectSurfaceProxy;

using PhotonReflectSurfaceProxyArray1D = FTypeArray1D<
    PhotonReflectSurfaceProxy,
    allocate_fortran_photon_reflect_surface_struct,
    deallocate_fortran_photon_reflect_surface_struct>;
using PhotonReflectSurfaceProxyArray2D =
    FTypeArray2D<PhotonReflectSurfaceProxy>;
using PhotonReflectSurfaceProxyArray3D =
    FTypeArray3D<PhotonReflectSurfaceProxy>;

using PhotonReflectSurfaceProxyAlloc1D = FTypeAlloc1D<
    PhotonReflectSurfaceProxyArray1D,
    allocate_photon_reflect_surface_struct_container,
    deallocate_photon_reflect_surface_struct_container,
    reallocate_photon_reflect_surface_struct_container_data,
    access_photon_reflect_surface_struct_container>;

class CoordProxy;

using CoordProxyArray1D = FTypeArray1D<
    CoordProxy,
    allocate_fortran_coord_struct,
    deallocate_fortran_coord_struct>;
using CoordProxyArray2D = FTypeArray2D<CoordProxy>;
using CoordProxyArray3D = FTypeArray3D<CoordProxy>;

using CoordProxyAlloc1D = FTypeAlloc1D<
    CoordProxyArray1D,
    allocate_coord_struct_container,
    deallocate_coord_struct_container,
    reallocate_coord_struct_container_data,
    access_coord_struct_container>;

class CoordArrayProxy;

using CoordArrayProxyArray1D = FTypeArray1D<
    CoordArrayProxy,
    allocate_fortran_coord_array_struct,
    deallocate_fortran_coord_array_struct>;
using CoordArrayProxyArray2D = FTypeArray2D<CoordArrayProxy>;
using CoordArrayProxyArray3D = FTypeArray3D<CoordArrayProxy>;

using CoordArrayProxyAlloc1D = FTypeAlloc1D<
    CoordArrayProxyArray1D,
    allocate_coord_array_struct_container,
    deallocate_coord_array_struct_container,
    reallocate_coord_array_struct_container_data,
    access_coord_array_struct_container>;

class BpmPhaseCouplingProxy;

using BpmPhaseCouplingProxyArray1D = FTypeArray1D<
    BpmPhaseCouplingProxy,
    allocate_fortran_bpm_phase_coupling_struct,
    deallocate_fortran_bpm_phase_coupling_struct>;
using BpmPhaseCouplingProxyArray2D = FTypeArray2D<BpmPhaseCouplingProxy>;
using BpmPhaseCouplingProxyArray3D = FTypeArray3D<BpmPhaseCouplingProxy>;

using BpmPhaseCouplingProxyAlloc1D = FTypeAlloc1D<
    BpmPhaseCouplingProxyArray1D,
    allocate_bpm_phase_coupling_struct_container,
    deallocate_bpm_phase_coupling_struct_container,
    reallocate_bpm_phase_coupling_struct_container_data,
    access_bpm_phase_coupling_struct_container>;

class ExpressionAtomProxy;

using ExpressionAtomProxyArray1D = FTypeArray1D<
    ExpressionAtomProxy,
    allocate_fortran_expression_atom_struct,
    deallocate_fortran_expression_atom_struct>;
using ExpressionAtomProxyArray2D = FTypeArray2D<ExpressionAtomProxy>;
using ExpressionAtomProxyArray3D = FTypeArray3D<ExpressionAtomProxy>;

using ExpressionAtomProxyAlloc1D = FTypeAlloc1D<
    ExpressionAtomProxyArray1D,
    allocate_expression_atom_struct_container,
    deallocate_expression_atom_struct_container,
    reallocate_expression_atom_struct_container_data,
    access_expression_atom_struct_container>;

class WakeSrZLongProxy;

using WakeSrZLongProxyArray1D = FTypeArray1D<
    WakeSrZLongProxy,
    allocate_fortran_wake_sr_z_long_struct,
    deallocate_fortran_wake_sr_z_long_struct>;
using WakeSrZLongProxyArray2D = FTypeArray2D<WakeSrZLongProxy>;
using WakeSrZLongProxyArray3D = FTypeArray3D<WakeSrZLongProxy>;

using WakeSrZLongProxyAlloc1D = FTypeAlloc1D<
    WakeSrZLongProxyArray1D,
    allocate_wake_sr_z_long_struct_container,
    deallocate_wake_sr_z_long_struct_container,
    reallocate_wake_sr_z_long_struct_container_data,
    access_wake_sr_z_long_struct_container>;

class WakeSrModeProxy;

using WakeSrModeProxyArray1D = FTypeArray1D<
    WakeSrModeProxy,
    allocate_fortran_wake_sr_mode_struct,
    deallocate_fortran_wake_sr_mode_struct>;
using WakeSrModeProxyArray2D = FTypeArray2D<WakeSrModeProxy>;
using WakeSrModeProxyArray3D = FTypeArray3D<WakeSrModeProxy>;

using WakeSrModeProxyAlloc1D = FTypeAlloc1D<
    WakeSrModeProxyArray1D,
    allocate_wake_sr_mode_struct_container,
    deallocate_wake_sr_mode_struct_container,
    reallocate_wake_sr_mode_struct_container_data,
    access_wake_sr_mode_struct_container>;

class WakeSrProxy;

using WakeSrProxyArray1D = FTypeArray1D<
    WakeSrProxy,
    allocate_fortran_wake_sr_struct,
    deallocate_fortran_wake_sr_struct>;
using WakeSrProxyArray2D = FTypeArray2D<WakeSrProxy>;
using WakeSrProxyArray3D = FTypeArray3D<WakeSrProxy>;

using WakeSrProxyAlloc1D = FTypeAlloc1D<
    WakeSrProxyArray1D,
    allocate_wake_sr_struct_container,
    deallocate_wake_sr_struct_container,
    reallocate_wake_sr_struct_container_data,
    access_wake_sr_struct_container>;

class WakeLrModeProxy;

using WakeLrModeProxyArray1D = FTypeArray1D<
    WakeLrModeProxy,
    allocate_fortran_wake_lr_mode_struct,
    deallocate_fortran_wake_lr_mode_struct>;
using WakeLrModeProxyArray2D = FTypeArray2D<WakeLrModeProxy>;
using WakeLrModeProxyArray3D = FTypeArray3D<WakeLrModeProxy>;

using WakeLrModeProxyAlloc1D = FTypeAlloc1D<
    WakeLrModeProxyArray1D,
    allocate_wake_lr_mode_struct_container,
    deallocate_wake_lr_mode_struct_container,
    reallocate_wake_lr_mode_struct_container_data,
    access_wake_lr_mode_struct_container>;

class WakeLrProxy;

using WakeLrProxyArray1D = FTypeArray1D<
    WakeLrProxy,
    allocate_fortran_wake_lr_struct,
    deallocate_fortran_wake_lr_struct>;
using WakeLrProxyArray2D = FTypeArray2D<WakeLrProxy>;
using WakeLrProxyArray3D = FTypeArray3D<WakeLrProxy>;

using WakeLrProxyAlloc1D = FTypeAlloc1D<
    WakeLrProxyArray1D,
    allocate_wake_lr_struct_container,
    deallocate_wake_lr_struct_container,
    reallocate_wake_lr_struct_container_data,
    access_wake_lr_struct_container>;

class LatEleLocProxy;

using LatEleLocProxyArray1D = FTypeArray1D<
    LatEleLocProxy,
    allocate_fortran_lat_ele_loc_struct,
    deallocate_fortran_lat_ele_loc_struct>;
using LatEleLocProxyArray2D = FTypeArray2D<LatEleLocProxy>;
using LatEleLocProxyArray3D = FTypeArray3D<LatEleLocProxy>;

using LatEleLocProxyAlloc1D = FTypeAlloc1D<
    LatEleLocProxyArray1D,
    allocate_lat_ele_loc_struct_container,
    deallocate_lat_ele_loc_struct_container,
    reallocate_lat_ele_loc_struct_container_data,
    access_lat_ele_loc_struct_container>;

class WakeProxy;

using WakeProxyArray1D = FTypeArray1D<
    WakeProxy,
    allocate_fortran_wake_struct,
    deallocate_fortran_wake_struct>;
using WakeProxyArray2D = FTypeArray2D<WakeProxy>;
using WakeProxyArray3D = FTypeArray3D<WakeProxy>;

using WakeProxyAlloc1D = FTypeAlloc1D<
    WakeProxyArray1D,
    allocate_wake_struct_container,
    deallocate_wake_struct_container,
    reallocate_wake_struct_container_data,
    access_wake_struct_container>;

class TaylorTermProxy;

using TaylorTermProxyArray1D = FTypeArray1D<
    TaylorTermProxy,
    allocate_fortran_taylor_term_struct,
    deallocate_fortran_taylor_term_struct>;
using TaylorTermProxyArray2D = FTypeArray2D<TaylorTermProxy>;
using TaylorTermProxyArray3D = FTypeArray3D<TaylorTermProxy>;

using TaylorTermProxyAlloc1D = FTypeAlloc1D<
    TaylorTermProxyArray1D,
    allocate_taylor_term_struct_container,
    deallocate_taylor_term_struct_container,
    reallocate_taylor_term_struct_container_data,
    access_taylor_term_struct_container>;

class TaylorProxy;

using TaylorProxyArray1D = FTypeArray1D<
    TaylorProxy,
    allocate_fortran_taylor_struct,
    deallocate_fortran_taylor_struct>;
using TaylorProxyArray2D = FTypeArray2D<TaylorProxy>;
using TaylorProxyArray3D = FTypeArray3D<TaylorProxy>;

using TaylorProxyAlloc1D = FTypeAlloc1D<
    TaylorProxyArray1D,
    allocate_taylor_struct_container,
    deallocate_taylor_struct_container,
    reallocate_taylor_struct_container_data,
    access_taylor_struct_container>;

class EmTaylorTermProxy;

using EmTaylorTermProxyArray1D = FTypeArray1D<
    EmTaylorTermProxy,
    allocate_fortran_em_taylor_term_struct,
    deallocate_fortran_em_taylor_term_struct>;
using EmTaylorTermProxyArray2D = FTypeArray2D<EmTaylorTermProxy>;
using EmTaylorTermProxyArray3D = FTypeArray3D<EmTaylorTermProxy>;

using EmTaylorTermProxyAlloc1D = FTypeAlloc1D<
    EmTaylorTermProxyArray1D,
    allocate_em_taylor_term_struct_container,
    deallocate_em_taylor_term_struct_container,
    reallocate_em_taylor_term_struct_container_data,
    access_em_taylor_term_struct_container>;

class EmTaylorProxy;

using EmTaylorProxyArray1D = FTypeArray1D<
    EmTaylorProxy,
    allocate_fortran_em_taylor_struct,
    deallocate_fortran_em_taylor_struct>;
using EmTaylorProxyArray2D = FTypeArray2D<EmTaylorProxy>;
using EmTaylorProxyArray3D = FTypeArray3D<EmTaylorProxy>;

using EmTaylorProxyAlloc1D = FTypeAlloc1D<
    EmTaylorProxyArray1D,
    allocate_em_taylor_struct_container,
    deallocate_em_taylor_struct_container,
    reallocate_em_taylor_struct_container_data,
    access_em_taylor_struct_container>;

class CartesianMapTerm1Proxy;

using CartesianMapTerm1ProxyArray1D = FTypeArray1D<
    CartesianMapTerm1Proxy,
    allocate_fortran_cartesian_map_term1_struct,
    deallocate_fortran_cartesian_map_term1_struct>;
using CartesianMapTerm1ProxyArray2D = FTypeArray2D<CartesianMapTerm1Proxy>;
using CartesianMapTerm1ProxyArray3D = FTypeArray3D<CartesianMapTerm1Proxy>;

using CartesianMapTerm1ProxyAlloc1D = FTypeAlloc1D<
    CartesianMapTerm1ProxyArray1D,
    allocate_cartesian_map_term1_struct_container,
    deallocate_cartesian_map_term1_struct_container,
    reallocate_cartesian_map_term1_struct_container_data,
    access_cartesian_map_term1_struct_container>;

class CartesianMapTermProxy;

using CartesianMapTermProxyArray1D = FTypeArray1D<
    CartesianMapTermProxy,
    allocate_fortran_cartesian_map_term_struct,
    deallocate_fortran_cartesian_map_term_struct>;
using CartesianMapTermProxyArray2D = FTypeArray2D<CartesianMapTermProxy>;
using CartesianMapTermProxyArray3D = FTypeArray3D<CartesianMapTermProxy>;

using CartesianMapTermProxyAlloc1D = FTypeAlloc1D<
    CartesianMapTermProxyArray1D,
    allocate_cartesian_map_term_struct_container,
    deallocate_cartesian_map_term_struct_container,
    reallocate_cartesian_map_term_struct_container_data,
    access_cartesian_map_term_struct_container>;

class CartesianMapProxy;

using CartesianMapProxyArray1D = FTypeArray1D<
    CartesianMapProxy,
    allocate_fortran_cartesian_map_struct,
    deallocate_fortran_cartesian_map_struct>;
using CartesianMapProxyArray2D = FTypeArray2D<CartesianMapProxy>;
using CartesianMapProxyArray3D = FTypeArray3D<CartesianMapProxy>;

using CartesianMapProxyAlloc1D = FTypeAlloc1D<
    CartesianMapProxyArray1D,
    allocate_cartesian_map_struct_container,
    deallocate_cartesian_map_struct_container,
    reallocate_cartesian_map_struct_container_data,
    access_cartesian_map_struct_container>;

class CylindricalMapTerm1Proxy;

using CylindricalMapTerm1ProxyArray1D = FTypeArray1D<
    CylindricalMapTerm1Proxy,
    allocate_fortran_cylindrical_map_term1_struct,
    deallocate_fortran_cylindrical_map_term1_struct>;
using CylindricalMapTerm1ProxyArray2D = FTypeArray2D<CylindricalMapTerm1Proxy>;
using CylindricalMapTerm1ProxyArray3D = FTypeArray3D<CylindricalMapTerm1Proxy>;

using CylindricalMapTerm1ProxyAlloc1D = FTypeAlloc1D<
    CylindricalMapTerm1ProxyArray1D,
    allocate_cylindrical_map_term1_struct_container,
    deallocate_cylindrical_map_term1_struct_container,
    reallocate_cylindrical_map_term1_struct_container_data,
    access_cylindrical_map_term1_struct_container>;

class CylindricalMapTermProxy;

using CylindricalMapTermProxyArray1D = FTypeArray1D<
    CylindricalMapTermProxy,
    allocate_fortran_cylindrical_map_term_struct,
    deallocate_fortran_cylindrical_map_term_struct>;
using CylindricalMapTermProxyArray2D = FTypeArray2D<CylindricalMapTermProxy>;
using CylindricalMapTermProxyArray3D = FTypeArray3D<CylindricalMapTermProxy>;

using CylindricalMapTermProxyAlloc1D = FTypeAlloc1D<
    CylindricalMapTermProxyArray1D,
    allocate_cylindrical_map_term_struct_container,
    deallocate_cylindrical_map_term_struct_container,
    reallocate_cylindrical_map_term_struct_container_data,
    access_cylindrical_map_term_struct_container>;

class CylindricalMapProxy;

using CylindricalMapProxyArray1D = FTypeArray1D<
    CylindricalMapProxy,
    allocate_fortran_cylindrical_map_struct,
    deallocate_fortran_cylindrical_map_struct>;
using CylindricalMapProxyArray2D = FTypeArray2D<CylindricalMapProxy>;
using CylindricalMapProxyArray3D = FTypeArray3D<CylindricalMapProxy>;

using CylindricalMapProxyAlloc1D = FTypeAlloc1D<
    CylindricalMapProxyArray1D,
    allocate_cylindrical_map_struct_container,
    deallocate_cylindrical_map_struct_container,
    reallocate_cylindrical_map_struct_container_data,
    access_cylindrical_map_struct_container>;

class BicubicCmplxCoefProxy;

using BicubicCmplxCoefProxyArray1D = FTypeArray1D<
    BicubicCmplxCoefProxy,
    allocate_fortran_bicubic_cmplx_coef_struct,
    deallocate_fortran_bicubic_cmplx_coef_struct>;
using BicubicCmplxCoefProxyArray2D = FTypeArray2D<BicubicCmplxCoefProxy>;
using BicubicCmplxCoefProxyArray3D = FTypeArray3D<BicubicCmplxCoefProxy>;

using BicubicCmplxCoefProxyAlloc1D = FTypeAlloc1D<
    BicubicCmplxCoefProxyArray1D,
    allocate_bicubic_cmplx_coef_struct_container,
    deallocate_bicubic_cmplx_coef_struct_container,
    reallocate_bicubic_cmplx_coef_struct_container_data,
    access_bicubic_cmplx_coef_struct_container>;

class TricubicCmplxCoefProxy;

using TricubicCmplxCoefProxyArray1D = FTypeArray1D<
    TricubicCmplxCoefProxy,
    allocate_fortran_tricubic_cmplx_coef_struct,
    deallocate_fortran_tricubic_cmplx_coef_struct>;
using TricubicCmplxCoefProxyArray2D = FTypeArray2D<TricubicCmplxCoefProxy>;
using TricubicCmplxCoefProxyArray3D = FTypeArray3D<TricubicCmplxCoefProxy>;

using TricubicCmplxCoefProxyAlloc1D = FTypeAlloc1D<
    TricubicCmplxCoefProxyArray1D,
    allocate_tricubic_cmplx_coef_struct_container,
    deallocate_tricubic_cmplx_coef_struct_container,
    reallocate_tricubic_cmplx_coef_struct_container_data,
    access_tricubic_cmplx_coef_struct_container>;

class GridFieldPt1Proxy;

using GridFieldPt1ProxyArray1D = FTypeArray1D<
    GridFieldPt1Proxy,
    allocate_fortran_grid_field_pt1_struct,
    deallocate_fortran_grid_field_pt1_struct>;
using GridFieldPt1ProxyArray2D = FTypeArray2D<GridFieldPt1Proxy>;
using GridFieldPt1ProxyArray3D = FTypeArray3D<GridFieldPt1Proxy>;

using GridFieldPt1ProxyAlloc1D = FTypeAlloc1D<
    GridFieldPt1ProxyArray1D,
    allocate_grid_field_pt1_struct_container,
    deallocate_grid_field_pt1_struct_container,
    reallocate_grid_field_pt1_struct_container_data,
    access_grid_field_pt1_struct_container>;

class GridFieldPtProxy;

using GridFieldPtProxyArray1D = FTypeArray1D<
    GridFieldPtProxy,
    allocate_fortran_grid_field_pt_struct,
    deallocate_fortran_grid_field_pt_struct>;
using GridFieldPtProxyArray2D = FTypeArray2D<GridFieldPtProxy>;
using GridFieldPtProxyArray3D = FTypeArray3D<GridFieldPtProxy>;

using GridFieldPtProxyAlloc1D = FTypeAlloc1D<
    GridFieldPtProxyArray1D,
    allocate_grid_field_pt_struct_container,
    deallocate_grid_field_pt_struct_container,
    reallocate_grid_field_pt_struct_container_data,
    access_grid_field_pt_struct_container>;

class GridFieldProxy;

using GridFieldProxyArray1D = FTypeArray1D<
    GridFieldProxy,
    allocate_fortran_grid_field_struct,
    deallocate_fortran_grid_field_struct>;
using GridFieldProxyArray2D = FTypeArray2D<GridFieldProxy>;
using GridFieldProxyArray3D = FTypeArray3D<GridFieldProxy>;

using GridFieldProxyAlloc1D = FTypeAlloc1D<
    GridFieldProxyArray1D,
    allocate_grid_field_struct_container,
    deallocate_grid_field_struct_container,
    reallocate_grid_field_struct_container_data,
    access_grid_field_struct_container>;

class FloorPositionProxy;

using FloorPositionProxyArray1D = FTypeArray1D<
    FloorPositionProxy,
    allocate_fortran_floor_position_struct,
    deallocate_fortran_floor_position_struct>;
using FloorPositionProxyArray2D = FTypeArray2D<FloorPositionProxy>;
using FloorPositionProxyArray3D = FTypeArray3D<FloorPositionProxy>;

using FloorPositionProxyAlloc1D = FTypeAlloc1D<
    FloorPositionProxyArray1D,
    allocate_floor_position_struct_container,
    deallocate_floor_position_struct_container,
    reallocate_floor_position_struct_container_data,
    access_floor_position_struct_container>;

class HighEnergySpaceChargeProxy;

using HighEnergySpaceChargeProxyArray1D = FTypeArray1D<
    HighEnergySpaceChargeProxy,
    allocate_fortran_high_energy_space_charge_struct,
    deallocate_fortran_high_energy_space_charge_struct>;
using HighEnergySpaceChargeProxyArray2D =
    FTypeArray2D<HighEnergySpaceChargeProxy>;
using HighEnergySpaceChargeProxyArray3D =
    FTypeArray3D<HighEnergySpaceChargeProxy>;

using HighEnergySpaceChargeProxyAlloc1D = FTypeAlloc1D<
    HighEnergySpaceChargeProxyArray1D,
    allocate_high_energy_space_charge_struct_container,
    deallocate_high_energy_space_charge_struct_container,
    reallocate_high_energy_space_charge_struct_container_data,
    access_high_energy_space_charge_struct_container>;

class XyDispProxy;

using XyDispProxyArray1D = FTypeArray1D<
    XyDispProxy,
    allocate_fortran_xy_disp_struct,
    deallocate_fortran_xy_disp_struct>;
using XyDispProxyArray2D = FTypeArray2D<XyDispProxy>;
using XyDispProxyArray3D = FTypeArray3D<XyDispProxy>;

using XyDispProxyAlloc1D = FTypeAlloc1D<
    XyDispProxyArray1D,
    allocate_xy_disp_struct_container,
    deallocate_xy_disp_struct_container,
    reallocate_xy_disp_struct_container_data,
    access_xy_disp_struct_container>;

class TwissProxy;

using TwissProxyArray1D = FTypeArray1D<
    TwissProxy,
    allocate_fortran_twiss_struct,
    deallocate_fortran_twiss_struct>;
using TwissProxyArray2D = FTypeArray2D<TwissProxy>;
using TwissProxyArray3D = FTypeArray3D<TwissProxy>;

using TwissProxyAlloc1D = FTypeAlloc1D<
    TwissProxyArray1D,
    allocate_twiss_struct_container,
    deallocate_twiss_struct_container,
    reallocate_twiss_struct_container_data,
    access_twiss_struct_container>;

class Mode3Proxy;

using Mode3ProxyArray1D = FTypeArray1D<
    Mode3Proxy,
    allocate_fortran_mode3_struct,
    deallocate_fortran_mode3_struct>;
using Mode3ProxyArray2D = FTypeArray2D<Mode3Proxy>;
using Mode3ProxyArray3D = FTypeArray3D<Mode3Proxy>;

using Mode3ProxyAlloc1D = FTypeAlloc1D<
    Mode3ProxyArray1D,
    allocate_mode3_struct_container,
    deallocate_mode3_struct_container,
    reallocate_mode3_struct_container_data,
    access_mode3_struct_container>;

class BookkeepingStateProxy;

using BookkeepingStateProxyArray1D = FTypeArray1D<
    BookkeepingStateProxy,
    allocate_fortran_bookkeeping_state_struct,
    deallocate_fortran_bookkeeping_state_struct>;
using BookkeepingStateProxyArray2D = FTypeArray2D<BookkeepingStateProxy>;
using BookkeepingStateProxyArray3D = FTypeArray3D<BookkeepingStateProxy>;

using BookkeepingStateProxyAlloc1D = FTypeAlloc1D<
    BookkeepingStateProxyArray1D,
    allocate_bookkeeping_state_struct_container,
    deallocate_bookkeeping_state_struct_container,
    reallocate_bookkeeping_state_struct_container_data,
    access_bookkeeping_state_struct_container>;

class RadMapProxy;

using RadMapProxyArray1D = FTypeArray1D<
    RadMapProxy,
    allocate_fortran_rad_map_struct,
    deallocate_fortran_rad_map_struct>;
using RadMapProxyArray2D = FTypeArray2D<RadMapProxy>;
using RadMapProxyArray3D = FTypeArray3D<RadMapProxy>;

using RadMapProxyAlloc1D = FTypeAlloc1D<
    RadMapProxyArray1D,
    allocate_rad_map_struct_container,
    deallocate_rad_map_struct_container,
    reallocate_rad_map_struct_container_data,
    access_rad_map_struct_container>;

class RadMapEleProxy;

using RadMapEleProxyArray1D = FTypeArray1D<
    RadMapEleProxy,
    allocate_fortran_rad_map_ele_struct,
    deallocate_fortran_rad_map_ele_struct>;
using RadMapEleProxyArray2D = FTypeArray2D<RadMapEleProxy>;
using RadMapEleProxyArray3D = FTypeArray3D<RadMapEleProxy>;

using RadMapEleProxyAlloc1D = FTypeAlloc1D<
    RadMapEleProxyArray1D,
    allocate_rad_map_ele_struct_container,
    deallocate_rad_map_ele_struct_container,
    reallocate_rad_map_ele_struct_container_data,
    access_rad_map_ele_struct_container>;

class GenGrad1Proxy;

using GenGrad1ProxyArray1D = FTypeArray1D<
    GenGrad1Proxy,
    allocate_fortran_gen_grad1_struct,
    deallocate_fortran_gen_grad1_struct>;
using GenGrad1ProxyArray2D = FTypeArray2D<GenGrad1Proxy>;
using GenGrad1ProxyArray3D = FTypeArray3D<GenGrad1Proxy>;

using GenGrad1ProxyAlloc1D = FTypeAlloc1D<
    GenGrad1ProxyArray1D,
    allocate_gen_grad1_struct_container,
    deallocate_gen_grad1_struct_container,
    reallocate_gen_grad1_struct_container_data,
    access_gen_grad1_struct_container>;

class GenGradMapProxy;

using GenGradMapProxyArray1D = FTypeArray1D<
    GenGradMapProxy,
    allocate_fortran_gen_grad_map_struct,
    deallocate_fortran_gen_grad_map_struct>;
using GenGradMapProxyArray2D = FTypeArray2D<GenGradMapProxy>;
using GenGradMapProxyArray3D = FTypeArray3D<GenGradMapProxy>;

using GenGradMapProxyAlloc1D = FTypeAlloc1D<
    GenGradMapProxyArray1D,
    allocate_gen_grad_map_struct_container,
    deallocate_gen_grad_map_struct_container,
    reallocate_gen_grad_map_struct_container_data,
    access_gen_grad_map_struct_container>;

class SurfaceSegmentedPtProxy;

using SurfaceSegmentedPtProxyArray1D = FTypeArray1D<
    SurfaceSegmentedPtProxy,
    allocate_fortran_surface_segmented_pt_struct,
    deallocate_fortran_surface_segmented_pt_struct>;
using SurfaceSegmentedPtProxyArray2D = FTypeArray2D<SurfaceSegmentedPtProxy>;
using SurfaceSegmentedPtProxyArray3D = FTypeArray3D<SurfaceSegmentedPtProxy>;

using SurfaceSegmentedPtProxyAlloc1D = FTypeAlloc1D<
    SurfaceSegmentedPtProxyArray1D,
    allocate_surface_segmented_pt_struct_container,
    deallocate_surface_segmented_pt_struct_container,
    reallocate_surface_segmented_pt_struct_container_data,
    access_surface_segmented_pt_struct_container>;

class SurfaceSegmentedProxy;

using SurfaceSegmentedProxyArray1D = FTypeArray1D<
    SurfaceSegmentedProxy,
    allocate_fortran_surface_segmented_struct,
    deallocate_fortran_surface_segmented_struct>;
using SurfaceSegmentedProxyArray2D = FTypeArray2D<SurfaceSegmentedProxy>;
using SurfaceSegmentedProxyArray3D = FTypeArray3D<SurfaceSegmentedProxy>;

using SurfaceSegmentedProxyAlloc1D = FTypeAlloc1D<
    SurfaceSegmentedProxyArray1D,
    allocate_surface_segmented_struct_container,
    deallocate_surface_segmented_struct_container,
    reallocate_surface_segmented_struct_container_data,
    access_surface_segmented_struct_container>;

class SurfaceHMisalignPtProxy;

using SurfaceHMisalignPtProxyArray1D = FTypeArray1D<
    SurfaceHMisalignPtProxy,
    allocate_fortran_surface_h_misalign_pt_struct,
    deallocate_fortran_surface_h_misalign_pt_struct>;
using SurfaceHMisalignPtProxyArray2D = FTypeArray2D<SurfaceHMisalignPtProxy>;
using SurfaceHMisalignPtProxyArray3D = FTypeArray3D<SurfaceHMisalignPtProxy>;

using SurfaceHMisalignPtProxyAlloc1D = FTypeAlloc1D<
    SurfaceHMisalignPtProxyArray1D,
    allocate_surface_h_misalign_pt_struct_container,
    deallocate_surface_h_misalign_pt_struct_container,
    reallocate_surface_h_misalign_pt_struct_container_data,
    access_surface_h_misalign_pt_struct_container>;

class SurfaceHMisalignProxy;

using SurfaceHMisalignProxyArray1D = FTypeArray1D<
    SurfaceHMisalignProxy,
    allocate_fortran_surface_h_misalign_struct,
    deallocate_fortran_surface_h_misalign_struct>;
using SurfaceHMisalignProxyArray2D = FTypeArray2D<SurfaceHMisalignProxy>;
using SurfaceHMisalignProxyArray3D = FTypeArray3D<SurfaceHMisalignProxy>;

using SurfaceHMisalignProxyAlloc1D = FTypeAlloc1D<
    SurfaceHMisalignProxyArray1D,
    allocate_surface_h_misalign_struct_container,
    deallocate_surface_h_misalign_struct_container,
    reallocate_surface_h_misalign_struct_container_data,
    access_surface_h_misalign_struct_container>;

class SurfaceDisplacementPtProxy;

using SurfaceDisplacementPtProxyArray1D = FTypeArray1D<
    SurfaceDisplacementPtProxy,
    allocate_fortran_surface_displacement_pt_struct,
    deallocate_fortran_surface_displacement_pt_struct>;
using SurfaceDisplacementPtProxyArray2D =
    FTypeArray2D<SurfaceDisplacementPtProxy>;
using SurfaceDisplacementPtProxyArray3D =
    FTypeArray3D<SurfaceDisplacementPtProxy>;

using SurfaceDisplacementPtProxyAlloc1D = FTypeAlloc1D<
    SurfaceDisplacementPtProxyArray1D,
    allocate_surface_displacement_pt_struct_container,
    deallocate_surface_displacement_pt_struct_container,
    reallocate_surface_displacement_pt_struct_container_data,
    access_surface_displacement_pt_struct_container>;

class SurfaceDisplacementProxy;

using SurfaceDisplacementProxyArray1D = FTypeArray1D<
    SurfaceDisplacementProxy,
    allocate_fortran_surface_displacement_struct,
    deallocate_fortran_surface_displacement_struct>;
using SurfaceDisplacementProxyArray2D = FTypeArray2D<SurfaceDisplacementProxy>;
using SurfaceDisplacementProxyArray3D = FTypeArray3D<SurfaceDisplacementProxy>;

using SurfaceDisplacementProxyAlloc1D = FTypeAlloc1D<
    SurfaceDisplacementProxyArray1D,
    allocate_surface_displacement_struct_container,
    deallocate_surface_displacement_struct_container,
    reallocate_surface_displacement_struct_container_data,
    access_surface_displacement_struct_container>;

class TargetPointProxy;

using TargetPointProxyArray1D = FTypeArray1D<
    TargetPointProxy,
    allocate_fortran_target_point_struct,
    deallocate_fortran_target_point_struct>;
using TargetPointProxyArray2D = FTypeArray2D<TargetPointProxy>;
using TargetPointProxyArray3D = FTypeArray3D<TargetPointProxy>;

using TargetPointProxyAlloc1D = FTypeAlloc1D<
    TargetPointProxyArray1D,
    allocate_target_point_struct_container,
    deallocate_target_point_struct_container,
    reallocate_target_point_struct_container_data,
    access_target_point_struct_container>;

class SurfaceCurvatureProxy;

using SurfaceCurvatureProxyArray1D = FTypeArray1D<
    SurfaceCurvatureProxy,
    allocate_fortran_surface_curvature_struct,
    deallocate_fortran_surface_curvature_struct>;
using SurfaceCurvatureProxyArray2D = FTypeArray2D<SurfaceCurvatureProxy>;
using SurfaceCurvatureProxyArray3D = FTypeArray3D<SurfaceCurvatureProxy>;

using SurfaceCurvatureProxyAlloc1D = FTypeAlloc1D<
    SurfaceCurvatureProxyArray1D,
    allocate_surface_curvature_struct_container,
    deallocate_surface_curvature_struct_container,
    reallocate_surface_curvature_struct_container_data,
    access_surface_curvature_struct_container>;

class PhotonTargetProxy;

using PhotonTargetProxyArray1D = FTypeArray1D<
    PhotonTargetProxy,
    allocate_fortran_photon_target_struct,
    deallocate_fortran_photon_target_struct>;
using PhotonTargetProxyArray2D = FTypeArray2D<PhotonTargetProxy>;
using PhotonTargetProxyArray3D = FTypeArray3D<PhotonTargetProxy>;

using PhotonTargetProxyAlloc1D = FTypeAlloc1D<
    PhotonTargetProxyArray1D,
    allocate_photon_target_struct_container,
    deallocate_photon_target_struct_container,
    reallocate_photon_target_struct_container_data,
    access_photon_target_struct_container>;

class PhotonMaterialProxy;

using PhotonMaterialProxyArray1D = FTypeArray1D<
    PhotonMaterialProxy,
    allocate_fortran_photon_material_struct,
    deallocate_fortran_photon_material_struct>;
using PhotonMaterialProxyArray2D = FTypeArray2D<PhotonMaterialProxy>;
using PhotonMaterialProxyArray3D = FTypeArray3D<PhotonMaterialProxy>;

using PhotonMaterialProxyAlloc1D = FTypeAlloc1D<
    PhotonMaterialProxyArray1D,
    allocate_photon_material_struct_container,
    deallocate_photon_material_struct_container,
    reallocate_photon_material_struct_container_data,
    access_photon_material_struct_container>;

class PixelPtProxy;

using PixelPtProxyArray1D = FTypeArray1D<
    PixelPtProxy,
    allocate_fortran_pixel_pt_struct,
    deallocate_fortran_pixel_pt_struct>;
using PixelPtProxyArray2D = FTypeArray2D<PixelPtProxy>;
using PixelPtProxyArray3D = FTypeArray3D<PixelPtProxy>;

using PixelPtProxyAlloc1D = FTypeAlloc1D<
    PixelPtProxyArray1D,
    allocate_pixel_pt_struct_container,
    deallocate_pixel_pt_struct_container,
    reallocate_pixel_pt_struct_container_data,
    access_pixel_pt_struct_container>;

class PixelDetecProxy;

using PixelDetecProxyArray1D = FTypeArray1D<
    PixelDetecProxy,
    allocate_fortran_pixel_detec_struct,
    deallocate_fortran_pixel_detec_struct>;
using PixelDetecProxyArray2D = FTypeArray2D<PixelDetecProxy>;
using PixelDetecProxyArray3D = FTypeArray3D<PixelDetecProxy>;

using PixelDetecProxyAlloc1D = FTypeAlloc1D<
    PixelDetecProxyArray1D,
    allocate_pixel_detec_struct_container,
    deallocate_pixel_detec_struct_container,
    reallocate_pixel_detec_struct_container_data,
    access_pixel_detec_struct_container>;

class PhotonElementProxy;

using PhotonElementProxyArray1D = FTypeArray1D<
    PhotonElementProxy,
    allocate_fortran_photon_element_struct,
    deallocate_fortran_photon_element_struct>;
using PhotonElementProxyArray2D = FTypeArray2D<PhotonElementProxy>;
using PhotonElementProxyArray3D = FTypeArray3D<PhotonElementProxy>;

using PhotonElementProxyAlloc1D = FTypeAlloc1D<
    PhotonElementProxyArray1D,
    allocate_photon_element_struct_container,
    deallocate_photon_element_struct_container,
    reallocate_photon_element_struct_container_data,
    access_photon_element_struct_container>;

class Wall3dVertexProxy;

using Wall3dVertexProxyArray1D = FTypeArray1D<
    Wall3dVertexProxy,
    allocate_fortran_wall3d_vertex_struct,
    deallocate_fortran_wall3d_vertex_struct>;
using Wall3dVertexProxyArray2D = FTypeArray2D<Wall3dVertexProxy>;
using Wall3dVertexProxyArray3D = FTypeArray3D<Wall3dVertexProxy>;

using Wall3dVertexProxyAlloc1D = FTypeAlloc1D<
    Wall3dVertexProxyArray1D,
    allocate_wall3d_vertex_struct_container,
    deallocate_wall3d_vertex_struct_container,
    reallocate_wall3d_vertex_struct_container_data,
    access_wall3d_vertex_struct_container>;

class Wall3dSectionProxy;

using Wall3dSectionProxyArray1D = FTypeArray1D<
    Wall3dSectionProxy,
    allocate_fortran_wall3d_section_struct,
    deallocate_fortran_wall3d_section_struct>;
using Wall3dSectionProxyArray2D = FTypeArray2D<Wall3dSectionProxy>;
using Wall3dSectionProxyArray3D = FTypeArray3D<Wall3dSectionProxy>;

using Wall3dSectionProxyAlloc1D = FTypeAlloc1D<
    Wall3dSectionProxyArray1D,
    allocate_wall3d_section_struct_container,
    deallocate_wall3d_section_struct_container,
    reallocate_wall3d_section_struct_container_data,
    access_wall3d_section_struct_container>;

class Wall3dProxy;

using Wall3dProxyArray1D = FTypeArray1D<
    Wall3dProxy,
    allocate_fortran_wall3d_struct,
    deallocate_fortran_wall3d_struct>;
using Wall3dProxyArray2D = FTypeArray2D<Wall3dProxy>;
using Wall3dProxyArray3D = FTypeArray3D<Wall3dProxy>;

using Wall3dProxyAlloc1D = FTypeAlloc1D<
    Wall3dProxyArray1D,
    allocate_wall3d_struct_container,
    deallocate_wall3d_struct_container,
    reallocate_wall3d_struct_container_data,
    access_wall3d_struct_container>;

class RamperLordProxy;

using RamperLordProxyArray1D = FTypeArray1D<
    RamperLordProxy,
    allocate_fortran_ramper_lord_struct,
    deallocate_fortran_ramper_lord_struct>;
using RamperLordProxyArray2D = FTypeArray2D<RamperLordProxy>;
using RamperLordProxyArray3D = FTypeArray3D<RamperLordProxy>;

using RamperLordProxyAlloc1D = FTypeAlloc1D<
    RamperLordProxyArray1D,
    allocate_ramper_lord_struct_container,
    deallocate_ramper_lord_struct_container,
    reallocate_ramper_lord_struct_container_data,
    access_ramper_lord_struct_container>;

class ControlProxy;

using ControlProxyArray1D = FTypeArray1D<
    ControlProxy,
    allocate_fortran_control_struct,
    deallocate_fortran_control_struct>;
using ControlProxyArray2D = FTypeArray2D<ControlProxy>;
using ControlProxyArray3D = FTypeArray3D<ControlProxy>;

using ControlProxyAlloc1D = FTypeAlloc1D<
    ControlProxyArray1D,
    allocate_control_struct_container,
    deallocate_control_struct_container,
    reallocate_control_struct_container_data,
    access_control_struct_container>;

class ControlVar1Proxy;

using ControlVar1ProxyArray1D = FTypeArray1D<
    ControlVar1Proxy,
    allocate_fortran_control_var1_struct,
    deallocate_fortran_control_var1_struct>;
using ControlVar1ProxyArray2D = FTypeArray2D<ControlVar1Proxy>;
using ControlVar1ProxyArray3D = FTypeArray3D<ControlVar1Proxy>;

using ControlVar1ProxyAlloc1D = FTypeAlloc1D<
    ControlVar1ProxyArray1D,
    allocate_control_var1_struct_container,
    deallocate_control_var1_struct_container,
    reallocate_control_var1_struct_container_data,
    access_control_var1_struct_container>;

class ControlRamp1Proxy;

using ControlRamp1ProxyArray1D = FTypeArray1D<
    ControlRamp1Proxy,
    allocate_fortran_control_ramp1_struct,
    deallocate_fortran_control_ramp1_struct>;
using ControlRamp1ProxyArray2D = FTypeArray2D<ControlRamp1Proxy>;
using ControlRamp1ProxyArray3D = FTypeArray3D<ControlRamp1Proxy>;

using ControlRamp1ProxyAlloc1D = FTypeAlloc1D<
    ControlRamp1ProxyArray1D,
    allocate_control_ramp1_struct_container,
    deallocate_control_ramp1_struct_container,
    reallocate_control_ramp1_struct_container_data,
    access_control_ramp1_struct_container>;

class ControllerProxy;

using ControllerProxyArray1D = FTypeArray1D<
    ControllerProxy,
    allocate_fortran_controller_struct,
    deallocate_fortran_controller_struct>;
using ControllerProxyArray2D = FTypeArray2D<ControllerProxy>;
using ControllerProxyArray3D = FTypeArray3D<ControllerProxy>;

using ControllerProxyAlloc1D = FTypeAlloc1D<
    ControllerProxyArray1D,
    allocate_controller_struct_container,
    deallocate_controller_struct_container,
    reallocate_controller_struct_container_data,
    access_controller_struct_container>;

class EllipseBeamInitProxy;

using EllipseBeamInitProxyArray1D = FTypeArray1D<
    EllipseBeamInitProxy,
    allocate_fortran_ellipse_beam_init_struct,
    deallocate_fortran_ellipse_beam_init_struct>;
using EllipseBeamInitProxyArray2D = FTypeArray2D<EllipseBeamInitProxy>;
using EllipseBeamInitProxyArray3D = FTypeArray3D<EllipseBeamInitProxy>;

using EllipseBeamInitProxyAlloc1D = FTypeAlloc1D<
    EllipseBeamInitProxyArray1D,
    allocate_ellipse_beam_init_struct_container,
    deallocate_ellipse_beam_init_struct_container,
    reallocate_ellipse_beam_init_struct_container_data,
    access_ellipse_beam_init_struct_container>;

class KvBeamInitProxy;

using KvBeamInitProxyArray1D = FTypeArray1D<
    KvBeamInitProxy,
    allocate_fortran_kv_beam_init_struct,
    deallocate_fortran_kv_beam_init_struct>;
using KvBeamInitProxyArray2D = FTypeArray2D<KvBeamInitProxy>;
using KvBeamInitProxyArray3D = FTypeArray3D<KvBeamInitProxy>;

using KvBeamInitProxyAlloc1D = FTypeAlloc1D<
    KvBeamInitProxyArray1D,
    allocate_kv_beam_init_struct_container,
    deallocate_kv_beam_init_struct_container,
    reallocate_kv_beam_init_struct_container_data,
    access_kv_beam_init_struct_container>;

class GridBeamInitProxy;

using GridBeamInitProxyArray1D = FTypeArray1D<
    GridBeamInitProxy,
    allocate_fortran_grid_beam_init_struct,
    deallocate_fortran_grid_beam_init_struct>;
using GridBeamInitProxyArray2D = FTypeArray2D<GridBeamInitProxy>;
using GridBeamInitProxyArray3D = FTypeArray3D<GridBeamInitProxy>;

using GridBeamInitProxyAlloc1D = FTypeAlloc1D<
    GridBeamInitProxyArray1D,
    allocate_grid_beam_init_struct_container,
    deallocate_grid_beam_init_struct_container,
    reallocate_grid_beam_init_struct_container_data,
    access_grid_beam_init_struct_container>;

class BeamInitProxy;

using BeamInitProxyArray1D = FTypeArray1D<
    BeamInitProxy,
    allocate_fortran_beam_init_struct,
    deallocate_fortran_beam_init_struct>;
using BeamInitProxyArray2D = FTypeArray2D<BeamInitProxy>;
using BeamInitProxyArray3D = FTypeArray3D<BeamInitProxy>;

using BeamInitProxyAlloc1D = FTypeAlloc1D<
    BeamInitProxyArray1D,
    allocate_beam_init_struct_container,
    deallocate_beam_init_struct_container,
    reallocate_beam_init_struct_container_data,
    access_beam_init_struct_container>;

class LatParamProxy;

using LatParamProxyArray1D = FTypeArray1D<
    LatParamProxy,
    allocate_fortran_lat_param_struct,
    deallocate_fortran_lat_param_struct>;
using LatParamProxyArray2D = FTypeArray2D<LatParamProxy>;
using LatParamProxyArray3D = FTypeArray3D<LatParamProxy>;

using LatParamProxyAlloc1D = FTypeAlloc1D<
    LatParamProxyArray1D,
    allocate_lat_param_struct_container,
    deallocate_lat_param_struct_container,
    reallocate_lat_param_struct_container_data,
    access_lat_param_struct_container>;

class ModeInfoProxy;

using ModeInfoProxyArray1D = FTypeArray1D<
    ModeInfoProxy,
    allocate_fortran_mode_info_struct,
    deallocate_fortran_mode_info_struct>;
using ModeInfoProxyArray2D = FTypeArray2D<ModeInfoProxy>;
using ModeInfoProxyArray3D = FTypeArray3D<ModeInfoProxy>;

using ModeInfoProxyAlloc1D = FTypeAlloc1D<
    ModeInfoProxyArray1D,
    allocate_mode_info_struct_container,
    deallocate_mode_info_struct_container,
    reallocate_mode_info_struct_container_data,
    access_mode_info_struct_container>;

class PreTrackerProxy;

using PreTrackerProxyArray1D = FTypeArray1D<
    PreTrackerProxy,
    allocate_fortran_pre_tracker_struct,
    deallocate_fortran_pre_tracker_struct>;
using PreTrackerProxyArray2D = FTypeArray2D<PreTrackerProxy>;
using PreTrackerProxyArray3D = FTypeArray3D<PreTrackerProxy>;

using PreTrackerProxyAlloc1D = FTypeAlloc1D<
    PreTrackerProxyArray1D,
    allocate_pre_tracker_struct_container,
    deallocate_pre_tracker_struct_container,
    reallocate_pre_tracker_struct_container_data,
    access_pre_tracker_struct_container>;

class AnormalModeProxy;

using AnormalModeProxyArray1D = FTypeArray1D<
    AnormalModeProxy,
    allocate_fortran_anormal_mode_struct,
    deallocate_fortran_anormal_mode_struct>;
using AnormalModeProxyArray2D = FTypeArray2D<AnormalModeProxy>;
using AnormalModeProxyArray3D = FTypeArray3D<AnormalModeProxy>;

using AnormalModeProxyAlloc1D = FTypeAlloc1D<
    AnormalModeProxyArray1D,
    allocate_anormal_mode_struct_container,
    deallocate_anormal_mode_struct_container,
    reallocate_anormal_mode_struct_container_data,
    access_anormal_mode_struct_container>;

class LinacNormalModeProxy;

using LinacNormalModeProxyArray1D = FTypeArray1D<
    LinacNormalModeProxy,
    allocate_fortran_linac_normal_mode_struct,
    deallocate_fortran_linac_normal_mode_struct>;
using LinacNormalModeProxyArray2D = FTypeArray2D<LinacNormalModeProxy>;
using LinacNormalModeProxyArray3D = FTypeArray3D<LinacNormalModeProxy>;

using LinacNormalModeProxyAlloc1D = FTypeAlloc1D<
    LinacNormalModeProxyArray1D,
    allocate_linac_normal_mode_struct_container,
    deallocate_linac_normal_mode_struct_container,
    reallocate_linac_normal_mode_struct_container_data,
    access_linac_normal_mode_struct_container>;

class NormalModesProxy;

using NormalModesProxyArray1D = FTypeArray1D<
    NormalModesProxy,
    allocate_fortran_normal_modes_struct,
    deallocate_fortran_normal_modes_struct>;
using NormalModesProxyArray2D = FTypeArray2D<NormalModesProxy>;
using NormalModesProxyArray3D = FTypeArray3D<NormalModesProxy>;

using NormalModesProxyAlloc1D = FTypeAlloc1D<
    NormalModesProxyArray1D,
    allocate_normal_modes_struct_container,
    deallocate_normal_modes_struct_container,
    reallocate_normal_modes_struct_container_data,
    access_normal_modes_struct_container>;

class EmFieldProxy;

using EmFieldProxyArray1D = FTypeArray1D<
    EmFieldProxy,
    allocate_fortran_em_field_struct,
    deallocate_fortran_em_field_struct>;
using EmFieldProxyArray2D = FTypeArray2D<EmFieldProxy>;
using EmFieldProxyArray3D = FTypeArray3D<EmFieldProxy>;

using EmFieldProxyAlloc1D = FTypeAlloc1D<
    EmFieldProxyArray1D,
    allocate_em_field_struct_container,
    deallocate_em_field_struct_container,
    reallocate_em_field_struct_container_data,
    access_em_field_struct_container>;

class StrongBeamProxy;

using StrongBeamProxyArray1D = FTypeArray1D<
    StrongBeamProxy,
    allocate_fortran_strong_beam_struct,
    deallocate_fortran_strong_beam_struct>;
using StrongBeamProxyArray2D = FTypeArray2D<StrongBeamProxy>;
using StrongBeamProxyArray3D = FTypeArray3D<StrongBeamProxy>;

using StrongBeamProxyAlloc1D = FTypeAlloc1D<
    StrongBeamProxyArray1D,
    allocate_strong_beam_struct_container,
    deallocate_strong_beam_struct_container,
    reallocate_strong_beam_struct_container_data,
    access_strong_beam_struct_container>;

class TrackPointProxy;

using TrackPointProxyArray1D = FTypeArray1D<
    TrackPointProxy,
    allocate_fortran_track_point_struct,
    deallocate_fortran_track_point_struct>;
using TrackPointProxyArray2D = FTypeArray2D<TrackPointProxy>;
using TrackPointProxyArray3D = FTypeArray3D<TrackPointProxy>;

using TrackPointProxyAlloc1D = FTypeAlloc1D<
    TrackPointProxyArray1D,
    allocate_track_point_struct_container,
    deallocate_track_point_struct_container,
    reallocate_track_point_struct_container_data,
    access_track_point_struct_container>;

class TrackProxy;

using TrackProxyArray1D = FTypeArray1D<
    TrackProxy,
    allocate_fortran_track_struct,
    deallocate_fortran_track_struct>;
using TrackProxyArray2D = FTypeArray2D<TrackProxy>;
using TrackProxyArray3D = FTypeArray3D<TrackProxy>;

using TrackProxyAlloc1D = FTypeAlloc1D<
    TrackProxyArray1D,
    allocate_track_struct_container,
    deallocate_track_struct_container,
    reallocate_track_struct_container_data,
    access_track_struct_container>;

class SpaceChargeCommonProxy;

using SpaceChargeCommonProxyArray1D = FTypeArray1D<
    SpaceChargeCommonProxy,
    allocate_fortran_space_charge_common_struct,
    deallocate_fortran_space_charge_common_struct>;
using SpaceChargeCommonProxyArray2D = FTypeArray2D<SpaceChargeCommonProxy>;
using SpaceChargeCommonProxyArray3D = FTypeArray3D<SpaceChargeCommonProxy>;

using SpaceChargeCommonProxyAlloc1D = FTypeAlloc1D<
    SpaceChargeCommonProxyArray1D,
    allocate_space_charge_common_struct_container,
    deallocate_space_charge_common_struct_container,
    reallocate_space_charge_common_struct_container_data,
    access_space_charge_common_struct_container>;

class BmadCommonProxy;

using BmadCommonProxyArray1D = FTypeArray1D<
    BmadCommonProxy,
    allocate_fortran_bmad_common_struct,
    deallocate_fortran_bmad_common_struct>;
using BmadCommonProxyArray2D = FTypeArray2D<BmadCommonProxy>;
using BmadCommonProxyArray3D = FTypeArray3D<BmadCommonProxy>;

using BmadCommonProxyAlloc1D = FTypeAlloc1D<
    BmadCommonProxyArray1D,
    allocate_bmad_common_struct_container,
    deallocate_bmad_common_struct_container,
    reallocate_bmad_common_struct_container_data,
    access_bmad_common_struct_container>;

class RadInt1Proxy;

using RadInt1ProxyArray1D = FTypeArray1D<
    RadInt1Proxy,
    allocate_fortran_rad_int1_struct,
    deallocate_fortran_rad_int1_struct>;
using RadInt1ProxyArray2D = FTypeArray2D<RadInt1Proxy>;
using RadInt1ProxyArray3D = FTypeArray3D<RadInt1Proxy>;

using RadInt1ProxyAlloc1D = FTypeAlloc1D<
    RadInt1ProxyArray1D,
    allocate_rad_int1_struct_container,
    deallocate_rad_int1_struct_container,
    reallocate_rad_int1_struct_container_data,
    access_rad_int1_struct_container>;

class RadIntBranchProxy;

using RadIntBranchProxyArray1D = FTypeArray1D<
    RadIntBranchProxy,
    allocate_fortran_rad_int_branch_struct,
    deallocate_fortran_rad_int_branch_struct>;
using RadIntBranchProxyArray2D = FTypeArray2D<RadIntBranchProxy>;
using RadIntBranchProxyArray3D = FTypeArray3D<RadIntBranchProxy>;

using RadIntBranchProxyAlloc1D = FTypeAlloc1D<
    RadIntBranchProxyArray1D,
    allocate_rad_int_branch_struct_container,
    deallocate_rad_int_branch_struct_container,
    reallocate_rad_int_branch_struct_container_data,
    access_rad_int_branch_struct_container>;

class RadIntAllEleProxy;

using RadIntAllEleProxyArray1D = FTypeArray1D<
    RadIntAllEleProxy,
    allocate_fortran_rad_int_all_ele_struct,
    deallocate_fortran_rad_int_all_ele_struct>;
using RadIntAllEleProxyArray2D = FTypeArray2D<RadIntAllEleProxy>;
using RadIntAllEleProxyArray3D = FTypeArray3D<RadIntAllEleProxy>;

using RadIntAllEleProxyAlloc1D = FTypeAlloc1D<
    RadIntAllEleProxyArray1D,
    allocate_rad_int_all_ele_struct_container,
    deallocate_rad_int_all_ele_struct_container,
    reallocate_rad_int_all_ele_struct_container_data,
    access_rad_int_all_ele_struct_container>;

class RfStairStepProxy;

using RfStairStepProxyArray1D = FTypeArray1D<
    RfStairStepProxy,
    allocate_fortran_rf_stair_step_struct,
    deallocate_fortran_rf_stair_step_struct>;
using RfStairStepProxyArray2D = FTypeArray2D<RfStairStepProxy>;
using RfStairStepProxyArray3D = FTypeArray3D<RfStairStepProxy>;

using RfStairStepProxyAlloc1D = FTypeAlloc1D<
    RfStairStepProxyArray1D,
    allocate_rf_stair_step_struct_container,
    deallocate_rf_stair_step_struct_container,
    reallocate_rf_stair_step_struct_container_data,
    access_rf_stair_step_struct_container>;

class RfEleProxy;

using RfEleProxyArray1D = FTypeArray1D<
    RfEleProxy,
    allocate_fortran_rf_ele_struct,
    deallocate_fortran_rf_ele_struct>;
using RfEleProxyArray2D = FTypeArray2D<RfEleProxy>;
using RfEleProxyArray3D = FTypeArray3D<RfEleProxy>;

using RfEleProxyAlloc1D = FTypeAlloc1D<
    RfEleProxyArray1D,
    allocate_rf_ele_struct_container,
    deallocate_rf_ele_struct_container,
    reallocate_rf_ele_struct_container_data,
    access_rf_ele_struct_container>;

class EleProxy;

using EleProxyArray1D = FTypeArray1D<
    EleProxy,
    allocate_fortran_ele_struct,
    deallocate_fortran_ele_struct>;
using EleProxyArray2D = FTypeArray2D<EleProxy>;
using EleProxyArray3D = FTypeArray3D<EleProxy>;

using EleProxyAlloc1D = FTypeAlloc1D<
    EleProxyArray1D,
    allocate_ele_struct_container,
    deallocate_ele_struct_container,
    reallocate_ele_struct_container_data,
    access_ele_struct_container>;

class ComplexTaylorTermProxy;

using ComplexTaylorTermProxyArray1D = FTypeArray1D<
    ComplexTaylorTermProxy,
    allocate_fortran_complex_taylor_term_struct,
    deallocate_fortran_complex_taylor_term_struct>;
using ComplexTaylorTermProxyArray2D = FTypeArray2D<ComplexTaylorTermProxy>;
using ComplexTaylorTermProxyArray3D = FTypeArray3D<ComplexTaylorTermProxy>;

using ComplexTaylorTermProxyAlloc1D = FTypeAlloc1D<
    ComplexTaylorTermProxyArray1D,
    allocate_complex_taylor_term_struct_container,
    deallocate_complex_taylor_term_struct_container,
    reallocate_complex_taylor_term_struct_container_data,
    access_complex_taylor_term_struct_container>;

class ComplexTaylorProxy;

using ComplexTaylorProxyArray1D = FTypeArray1D<
    ComplexTaylorProxy,
    allocate_fortran_complex_taylor_struct,
    deallocate_fortran_complex_taylor_struct>;
using ComplexTaylorProxyArray2D = FTypeArray2D<ComplexTaylorProxy>;
using ComplexTaylorProxyArray3D = FTypeArray3D<ComplexTaylorProxy>;

using ComplexTaylorProxyAlloc1D = FTypeAlloc1D<
    ComplexTaylorProxyArray1D,
    allocate_complex_taylor_struct_container,
    deallocate_complex_taylor_struct_container,
    reallocate_complex_taylor_struct_container_data,
    access_complex_taylor_struct_container>;

class BranchProxy;

using BranchProxyArray1D = FTypeArray1D<
    BranchProxy,
    allocate_fortran_branch_struct,
    deallocate_fortran_branch_struct>;
using BranchProxyArray2D = FTypeArray2D<BranchProxy>;
using BranchProxyArray3D = FTypeArray3D<BranchProxy>;

using BranchProxyAlloc1D = FTypeAlloc1D<
    BranchProxyArray1D,
    allocate_branch_struct_container,
    deallocate_branch_struct_container,
    reallocate_branch_struct_container_data,
    access_branch_struct_container>;

class LatProxy;

using LatProxyArray1D = FTypeArray1D<
    LatProxy,
    allocate_fortran_lat_struct,
    deallocate_fortran_lat_struct>;
using LatProxyArray2D = FTypeArray2D<LatProxy>;
using LatProxyArray3D = FTypeArray3D<LatProxy>;

using LatProxyAlloc1D = FTypeAlloc1D<
    LatProxyArray1D,
    allocate_lat_struct_container,
    deallocate_lat_struct_container,
    reallocate_lat_struct_container_data,
    access_lat_struct_container>;

class BunchProxy;

using BunchProxyArray1D = FTypeArray1D<
    BunchProxy,
    allocate_fortran_bunch_struct,
    deallocate_fortran_bunch_struct>;
using BunchProxyArray2D = FTypeArray2D<BunchProxy>;
using BunchProxyArray3D = FTypeArray3D<BunchProxy>;

using BunchProxyAlloc1D = FTypeAlloc1D<
    BunchProxyArray1D,
    allocate_bunch_struct_container,
    deallocate_bunch_struct_container,
    reallocate_bunch_struct_container_data,
    access_bunch_struct_container>;

class BunchParamsProxy;

using BunchParamsProxyArray1D = FTypeArray1D<
    BunchParamsProxy,
    allocate_fortran_bunch_params_struct,
    deallocate_fortran_bunch_params_struct>;
using BunchParamsProxyArray2D = FTypeArray2D<BunchParamsProxy>;
using BunchParamsProxyArray3D = FTypeArray3D<BunchParamsProxy>;

using BunchParamsProxyAlloc1D = FTypeAlloc1D<
    BunchParamsProxyArray1D,
    allocate_bunch_params_struct_container,
    deallocate_bunch_params_struct_container,
    reallocate_bunch_params_struct_container_data,
    access_bunch_params_struct_container>;

class BeamProxy;

using BeamProxyArray1D = FTypeArray1D<
    BeamProxy,
    allocate_fortran_beam_struct,
    deallocate_fortran_beam_struct>;
using BeamProxyArray2D = FTypeArray2D<BeamProxy>;
using BeamProxyArray3D = FTypeArray3D<BeamProxy>;

using BeamProxyAlloc1D = FTypeAlloc1D<
    BeamProxyArray1D,
    allocate_beam_struct_container,
    deallocate_beam_struct_container,
    reallocate_beam_struct_container_data,
    access_beam_struct_container>;

class AperturePointProxy;

using AperturePointProxyArray1D = FTypeArray1D<
    AperturePointProxy,
    allocate_fortran_aperture_point_struct,
    deallocate_fortran_aperture_point_struct>;
using AperturePointProxyArray2D = FTypeArray2D<AperturePointProxy>;
using AperturePointProxyArray3D = FTypeArray3D<AperturePointProxy>;

using AperturePointProxyAlloc1D = FTypeAlloc1D<
    AperturePointProxyArray1D,
    allocate_aperture_point_struct_container,
    deallocate_aperture_point_struct_container,
    reallocate_aperture_point_struct_container_data,
    access_aperture_point_struct_container>;

class ApertureParamProxy;

using ApertureParamProxyArray1D = FTypeArray1D<
    ApertureParamProxy,
    allocate_fortran_aperture_param_struct,
    deallocate_fortran_aperture_param_struct>;
using ApertureParamProxyArray2D = FTypeArray2D<ApertureParamProxy>;
using ApertureParamProxyArray3D = FTypeArray3D<ApertureParamProxy>;

using ApertureParamProxyAlloc1D = FTypeAlloc1D<
    ApertureParamProxyArray1D,
    allocate_aperture_param_struct_container,
    deallocate_aperture_param_struct_container,
    reallocate_aperture_param_struct_container_data,
    access_aperture_param_struct_container>;

class ApertureScanProxy;

using ApertureScanProxyArray1D = FTypeArray1D<
    ApertureScanProxy,
    allocate_fortran_aperture_scan_struct,
    deallocate_fortran_aperture_scan_struct>;
using ApertureScanProxyArray2D = FTypeArray2D<ApertureScanProxy>;
using ApertureScanProxyArray3D = FTypeArray3D<ApertureScanProxy>;

using ApertureScanProxyAlloc1D = FTypeAlloc1D<
    ApertureScanProxyArray1D,
    allocate_aperture_scan_struct_container,
    deallocate_aperture_scan_struct_container,
    reallocate_aperture_scan_struct_container_data,
    access_aperture_scan_struct_container>;

class ElePointerProxy;

using ElePointerProxyArray1D = FTypeArray1D<
    ElePointerProxy,
    allocate_fortran_ele_pointer_struct,
    deallocate_fortran_ele_pointer_struct>;
using ElePointerProxyArray2D = FTypeArray2D<ElePointerProxy>;
using ElePointerProxyArray3D = FTypeArray3D<ElePointerProxy>;

using ElePointerProxyAlloc1D = FTypeAlloc1D<
    ElePointerProxyArray1D,
    allocate_ele_pointer_struct_container,
    deallocate_ele_pointer_struct_container,
    reallocate_ele_pointer_struct_container_data,
    access_ele_pointer_struct_container>;

class ExpressionTreeProxy;

using ExpressionTreeProxyArray1D = FTypeArray1D<
    ExpressionTreeProxy,
    allocate_fortran_expression_tree_struct,
    deallocate_fortran_expression_tree_struct>;
using ExpressionTreeProxyArray2D = FTypeArray2D<ExpressionTreeProxy>;
using ExpressionTreeProxyArray3D = FTypeArray3D<ExpressionTreeProxy>;

using ExpressionTreeProxyAlloc1D = FTypeAlloc1D<
    ExpressionTreeProxyArray1D,
    allocate_expression_tree_struct_container,
    deallocate_expression_tree_struct_container,
    reallocate_expression_tree_struct_container_data,
    access_expression_tree_struct_container>;

class NametableProxy;

using NametableProxyArray1D = FTypeArray1D<
    NametableProxy,
    allocate_fortran_nametable_struct,
    deallocate_fortran_nametable_struct>;
using NametableProxyArray2D = FTypeArray2D<NametableProxy>;
using NametableProxyArray3D = FTypeArray3D<NametableProxy>;

using NametableProxyAlloc1D = FTypeAlloc1D<
    NametableProxyArray1D,
    allocate_nametable_struct_container,
    deallocate_nametable_struct_container,
    reallocate_nametable_struct_container_data,
    access_nametable_struct_container>;

class TaoSpinDnDpzProxy;

using TaoSpinDnDpzProxyArray1D = FTypeArray1D<
    TaoSpinDnDpzProxy,
    allocate_fortran_tao_spin_dn_dpz_struct,
    deallocate_fortran_tao_spin_dn_dpz_struct>;
using TaoSpinDnDpzProxyArray2D = FTypeArray2D<TaoSpinDnDpzProxy>;
using TaoSpinDnDpzProxyArray3D = FTypeArray3D<TaoSpinDnDpzProxy>;

using TaoSpinDnDpzProxyAlloc1D = FTypeAlloc1D<
    TaoSpinDnDpzProxyArray1D,
    allocate_tao_spin_dn_dpz_struct_container,
    deallocate_tao_spin_dn_dpz_struct_container,
    reallocate_tao_spin_dn_dpz_struct_container_data,
    access_tao_spin_dn_dpz_struct_container>;

class ResonanceHProxy;

using ResonanceHProxyArray1D = FTypeArray1D<
    ResonanceHProxy,
    allocate_fortran_resonance_h_struct,
    deallocate_fortran_resonance_h_struct>;
using ResonanceHProxyArray2D = FTypeArray2D<ResonanceHProxy>;
using ResonanceHProxyArray3D = FTypeArray3D<ResonanceHProxy>;

using ResonanceHProxyAlloc1D = FTypeAlloc1D<
    ResonanceHProxyArray1D,
    allocate_resonance_h_struct_container,
    deallocate_resonance_h_struct_container,
    reallocate_resonance_h_struct_container_data,
    access_resonance_h_struct_container>;

class SpinOrbitMap1Proxy;

using SpinOrbitMap1ProxyArray1D = FTypeArray1D<
    SpinOrbitMap1Proxy,
    allocate_fortran_spin_orbit_map1_struct,
    deallocate_fortran_spin_orbit_map1_struct>;
using SpinOrbitMap1ProxyArray2D = FTypeArray2D<SpinOrbitMap1Proxy>;
using SpinOrbitMap1ProxyArray3D = FTypeArray3D<SpinOrbitMap1Proxy>;

using SpinOrbitMap1ProxyAlloc1D = FTypeAlloc1D<
    SpinOrbitMap1ProxyArray1D,
    allocate_spin_orbit_map1_struct_container,
    deallocate_spin_orbit_map1_struct_container,
    reallocate_spin_orbit_map1_struct_container_data,
    access_spin_orbit_map1_struct_container>;

class SpinAxisProxy;

using SpinAxisProxyArray1D = FTypeArray1D<
    SpinAxisProxy,
    allocate_fortran_spin_axis_struct,
    deallocate_fortran_spin_axis_struct>;
using SpinAxisProxyArray2D = FTypeArray2D<SpinAxisProxy>;
using SpinAxisProxyArray3D = FTypeArray3D<SpinAxisProxy>;

using SpinAxisProxyAlloc1D = FTypeAlloc1D<
    SpinAxisProxyArray1D,
    allocate_spin_axis_struct_container,
    deallocate_spin_axis_struct_container,
    reallocate_spin_axis_struct_container_data,
    access_spin_axis_struct_container>;

class PtcNormalFormProxy;

using PtcNormalFormProxyArray1D = FTypeArray1D<
    PtcNormalFormProxy,
    allocate_fortran_ptc_normal_form_struct,
    deallocate_fortran_ptc_normal_form_struct>;
using PtcNormalFormProxyArray2D = FTypeArray2D<PtcNormalFormProxy>;
using PtcNormalFormProxyArray3D = FTypeArray3D<PtcNormalFormProxy>;

using PtcNormalFormProxyAlloc1D = FTypeAlloc1D<
    PtcNormalFormProxyArray1D,
    allocate_ptc_normal_form_struct_container,
    deallocate_ptc_normal_form_struct_container,
    reallocate_ptc_normal_form_struct_container_data,
    access_ptc_normal_form_struct_container>;

class BmadNormalFormProxy;

using BmadNormalFormProxyArray1D = FTypeArray1D<
    BmadNormalFormProxy,
    allocate_fortran_bmad_normal_form_struct,
    deallocate_fortran_bmad_normal_form_struct>;
using BmadNormalFormProxyArray2D = FTypeArray2D<BmadNormalFormProxy>;
using BmadNormalFormProxyArray3D = FTypeArray3D<BmadNormalFormProxy>;

using BmadNormalFormProxyAlloc1D = FTypeAlloc1D<
    BmadNormalFormProxyArray1D,
    allocate_bmad_normal_form_struct_container,
    deallocate_bmad_normal_form_struct_container,
    reallocate_bmad_normal_form_struct_container_data,
    access_bmad_normal_form_struct_container>;

class BunchTrackProxy;

using BunchTrackProxyArray1D = FTypeArray1D<
    BunchTrackProxy,
    allocate_fortran_bunch_track_struct,
    deallocate_fortran_bunch_track_struct>;
using BunchTrackProxyArray2D = FTypeArray2D<BunchTrackProxy>;
using BunchTrackProxyArray3D = FTypeArray3D<BunchTrackProxy>;

using BunchTrackProxyAlloc1D = FTypeAlloc1D<
    BunchTrackProxyArray1D,
    allocate_bunch_track_struct_container,
    deallocate_bunch_track_struct_container,
    reallocate_bunch_track_struct_container_data,
    access_bunch_track_struct_container>;

class SummationRdtProxy;

using SummationRdtProxyArray1D = FTypeArray1D<
    SummationRdtProxy,
    allocate_fortran_summation_rdt_struct,
    deallocate_fortran_summation_rdt_struct>;
using SummationRdtProxyArray2D = FTypeArray2D<SummationRdtProxy>;
using SummationRdtProxyArray3D = FTypeArray3D<SummationRdtProxy>;

using SummationRdtProxyAlloc1D = FTypeAlloc1D<
    SummationRdtProxyArray1D,
    allocate_summation_rdt_struct_container,
    deallocate_summation_rdt_struct_container,
    reallocate_summation_rdt_struct_container_data,
    access_summation_rdt_struct_container>;

class TaoEleShapeProxy;

using TaoEleShapeProxyArray1D = FTypeArray1D<
    TaoEleShapeProxy,
    allocate_fortran_tao_ele_shape_struct,
    deallocate_fortran_tao_ele_shape_struct>;
using TaoEleShapeProxyArray2D = FTypeArray2D<TaoEleShapeProxy>;
using TaoEleShapeProxyArray3D = FTypeArray3D<TaoEleShapeProxy>;

using TaoEleShapeProxyAlloc1D = FTypeAlloc1D<
    TaoEleShapeProxyArray1D,
    allocate_tao_ele_shape_struct_container,
    deallocate_tao_ele_shape_struct_container,
    reallocate_tao_ele_shape_struct_container_data,
    access_tao_ele_shape_struct_container>;

class TaoElePointerProxy;

using TaoElePointerProxyArray1D = FTypeArray1D<
    TaoElePointerProxy,
    allocate_fortran_tao_ele_pointer_struct,
    deallocate_fortran_tao_ele_pointer_struct>;
using TaoElePointerProxyArray2D = FTypeArray2D<TaoElePointerProxy>;
using TaoElePointerProxyArray3D = FTypeArray3D<TaoElePointerProxy>;

using TaoElePointerProxyAlloc1D = FTypeAlloc1D<
    TaoElePointerProxyArray1D,
    allocate_tao_ele_pointer_struct_container,
    deallocate_tao_ele_pointer_struct_container,
    reallocate_tao_ele_pointer_struct_container_data,
    access_tao_ele_pointer_struct_container>;

class TaoCurveProxy;

using TaoCurveProxyArray1D = FTypeArray1D<
    TaoCurveProxy,
    allocate_fortran_tao_curve_struct,
    deallocate_fortran_tao_curve_struct>;
using TaoCurveProxyArray2D = FTypeArray2D<TaoCurveProxy>;
using TaoCurveProxyArray3D = FTypeArray3D<TaoCurveProxy>;

using TaoCurveProxyAlloc1D = FTypeAlloc1D<
    TaoCurveProxyArray1D,
    allocate_tao_curve_struct_container,
    deallocate_tao_curve_struct_container,
    reallocate_tao_curve_struct_container_data,
    access_tao_curve_struct_container>;

class TaoCurveColorProxy;

using TaoCurveColorProxyArray1D = FTypeArray1D<
    TaoCurveColorProxy,
    allocate_fortran_tao_curve_color_struct,
    deallocate_fortran_tao_curve_color_struct>;
using TaoCurveColorProxyArray2D = FTypeArray2D<TaoCurveColorProxy>;
using TaoCurveColorProxyArray3D = FTypeArray3D<TaoCurveColorProxy>;

using TaoCurveColorProxyAlloc1D = FTypeAlloc1D<
    TaoCurveColorProxyArray1D,
    allocate_tao_curve_color_struct_container,
    deallocate_tao_curve_color_struct_container,
    reallocate_tao_curve_color_struct_container_data,
    access_tao_curve_color_struct_container>;

class TaoCurveOrbitProxy;

using TaoCurveOrbitProxyArray1D = FTypeArray1D<
    TaoCurveOrbitProxy,
    allocate_fortran_tao_curve_orbit_struct,
    deallocate_fortran_tao_curve_orbit_struct>;
using TaoCurveOrbitProxyArray2D = FTypeArray2D<TaoCurveOrbitProxy>;
using TaoCurveOrbitProxyArray3D = FTypeArray3D<TaoCurveOrbitProxy>;

using TaoCurveOrbitProxyAlloc1D = FTypeAlloc1D<
    TaoCurveOrbitProxyArray1D,
    allocate_tao_curve_orbit_struct_container,
    deallocate_tao_curve_orbit_struct_container,
    reallocate_tao_curve_orbit_struct_container_data,
    access_tao_curve_orbit_struct_container>;

class TaoHistogramProxy;

using TaoHistogramProxyArray1D = FTypeArray1D<
    TaoHistogramProxy,
    allocate_fortran_tao_histogram_struct,
    deallocate_fortran_tao_histogram_struct>;
using TaoHistogramProxyArray2D = FTypeArray2D<TaoHistogramProxy>;
using TaoHistogramProxyArray3D = FTypeArray3D<TaoHistogramProxy>;

using TaoHistogramProxyAlloc1D = FTypeAlloc1D<
    TaoHistogramProxyArray1D,
    allocate_tao_histogram_struct_container,
    deallocate_tao_histogram_struct_container,
    reallocate_tao_histogram_struct_container_data,
    access_tao_histogram_struct_container>;

class LatEleOrder1Proxy;

using LatEleOrder1ProxyArray1D = FTypeArray1D<
    LatEleOrder1Proxy,
    allocate_fortran_lat_ele_order1_struct,
    deallocate_fortran_lat_ele_order1_struct>;
using LatEleOrder1ProxyArray2D = FTypeArray2D<LatEleOrder1Proxy>;
using LatEleOrder1ProxyArray3D = FTypeArray3D<LatEleOrder1Proxy>;

using LatEleOrder1ProxyAlloc1D = FTypeAlloc1D<
    LatEleOrder1ProxyArray1D,
    allocate_lat_ele_order1_struct_container,
    deallocate_lat_ele_order1_struct_container,
    reallocate_lat_ele_order1_struct_container_data,
    access_lat_ele_order1_struct_container>;

class LatEleOrderArrayProxy;

using LatEleOrderArrayProxyArray1D = FTypeArray1D<
    LatEleOrderArrayProxy,
    allocate_fortran_lat_ele_order_array_struct,
    deallocate_fortran_lat_ele_order_array_struct>;
using LatEleOrderArrayProxyArray2D = FTypeArray2D<LatEleOrderArrayProxy>;
using LatEleOrderArrayProxyArray3D = FTypeArray3D<LatEleOrderArrayProxy>;

using LatEleOrderArrayProxyAlloc1D = FTypeAlloc1D<
    LatEleOrderArrayProxyArray1D,
    allocate_lat_ele_order_array_struct_container,
    deallocate_lat_ele_order_array_struct_container,
    reallocate_lat_ele_order_array_struct_container_data,
    access_lat_ele_order_array_struct_container>;

class TaoLatSigmaProxy;

using TaoLatSigmaProxyArray1D = FTypeArray1D<
    TaoLatSigmaProxy,
    allocate_fortran_tao_lat_sigma_struct,
    deallocate_fortran_tao_lat_sigma_struct>;
using TaoLatSigmaProxyArray2D = FTypeArray2D<TaoLatSigmaProxy>;
using TaoLatSigmaProxyArray3D = FTypeArray3D<TaoLatSigmaProxy>;

using TaoLatSigmaProxyAlloc1D = FTypeAlloc1D<
    TaoLatSigmaProxyArray1D,
    allocate_tao_lat_sigma_struct_container,
    deallocate_tao_lat_sigma_struct_container,
    reallocate_tao_lat_sigma_struct_container_data,
    access_tao_lat_sigma_struct_container>;

class TaoSpinEleProxy;

using TaoSpinEleProxyArray1D = FTypeArray1D<
    TaoSpinEleProxy,
    allocate_fortran_tao_spin_ele_struct,
    deallocate_fortran_tao_spin_ele_struct>;
using TaoSpinEleProxyArray2D = FTypeArray2D<TaoSpinEleProxy>;
using TaoSpinEleProxyArray3D = FTypeArray3D<TaoSpinEleProxy>;

using TaoSpinEleProxyAlloc1D = FTypeAlloc1D<
    TaoSpinEleProxyArray1D,
    allocate_tao_spin_ele_struct_container,
    deallocate_tao_spin_ele_struct_container,
    reallocate_tao_spin_ele_struct_container_data,
    access_tao_spin_ele_struct_container>;

class TaoPlotCacheProxy;

using TaoPlotCacheProxyArray1D = FTypeArray1D<
    TaoPlotCacheProxy,
    allocate_fortran_tao_plot_cache_struct,
    deallocate_fortran_tao_plot_cache_struct>;
using TaoPlotCacheProxyArray2D = FTypeArray2D<TaoPlotCacheProxy>;
using TaoPlotCacheProxyArray3D = FTypeArray3D<TaoPlotCacheProxy>;

using TaoPlotCacheProxyAlloc1D = FTypeAlloc1D<
    TaoPlotCacheProxyArray1D,
    allocate_tao_plot_cache_struct_container,
    deallocate_tao_plot_cache_struct_container,
    reallocate_tao_plot_cache_struct_container_data,
    access_tao_plot_cache_struct_container>;

class TaoSpinPolarizationProxy;

using TaoSpinPolarizationProxyArray1D = FTypeArray1D<
    TaoSpinPolarizationProxy,
    allocate_fortran_tao_spin_polarization_struct,
    deallocate_fortran_tao_spin_polarization_struct>;
using TaoSpinPolarizationProxyArray2D = FTypeArray2D<TaoSpinPolarizationProxy>;
using TaoSpinPolarizationProxyArray3D = FTypeArray3D<TaoSpinPolarizationProxy>;

using TaoSpinPolarizationProxyAlloc1D = FTypeAlloc1D<
    TaoSpinPolarizationProxyArray1D,
    allocate_tao_spin_polarization_struct_container,
    deallocate_tao_spin_polarization_struct_container,
    reallocate_tao_spin_polarization_struct_container_data,
    access_tao_spin_polarization_struct_container>;

class TaoLatticeBranchProxy;

using TaoLatticeBranchProxyArray1D = FTypeArray1D<
    TaoLatticeBranchProxy,
    allocate_fortran_tao_lattice_branch_struct,
    deallocate_fortran_tao_lattice_branch_struct>;
using TaoLatticeBranchProxyArray2D = FTypeArray2D<TaoLatticeBranchProxy>;
using TaoLatticeBranchProxyArray3D = FTypeArray3D<TaoLatticeBranchProxy>;

using TaoLatticeBranchProxyAlloc1D = FTypeAlloc1D<
    TaoLatticeBranchProxyArray1D,
    allocate_tao_lattice_branch_struct_container,
    deallocate_tao_lattice_branch_struct_container,
    reallocate_tao_lattice_branch_struct_container_data,
    access_tao_lattice_branch_struct_container>;

class TaoModelElementProxy;

using TaoModelElementProxyArray1D = FTypeArray1D<
    TaoModelElementProxy,
    allocate_fortran_tao_model_element_struct,
    deallocate_fortran_tao_model_element_struct>;
using TaoModelElementProxyArray2D = FTypeArray2D<TaoModelElementProxy>;
using TaoModelElementProxyArray3D = FTypeArray3D<TaoModelElementProxy>;

using TaoModelElementProxyAlloc1D = FTypeAlloc1D<
    TaoModelElementProxyArray1D,
    allocate_tao_model_element_struct_container,
    deallocate_tao_model_element_struct_container,
    reallocate_tao_model_element_struct_container_data,
    access_tao_model_element_struct_container>;

class TaoBeamBranchProxy;

using TaoBeamBranchProxyArray1D = FTypeArray1D<
    TaoBeamBranchProxy,
    allocate_fortran_tao_beam_branch_struct,
    deallocate_fortran_tao_beam_branch_struct>;
using TaoBeamBranchProxyArray2D = FTypeArray2D<TaoBeamBranchProxy>;
using TaoBeamBranchProxyArray3D = FTypeArray3D<TaoBeamBranchProxy>;

using TaoBeamBranchProxyAlloc1D = FTypeAlloc1D<
    TaoBeamBranchProxyArray1D,
    allocate_tao_beam_branch_struct_container,
    deallocate_tao_beam_branch_struct_container,
    reallocate_tao_beam_branch_struct_container_data,
    access_tao_beam_branch_struct_container>;

class TaoD1DataProxy;

using TaoD1DataProxyArray1D = FTypeArray1D<
    TaoD1DataProxy,
    allocate_fortran_tao_d1_data_struct,
    deallocate_fortran_tao_d1_data_struct>;
using TaoD1DataProxyArray2D = FTypeArray2D<TaoD1DataProxy>;
using TaoD1DataProxyArray3D = FTypeArray3D<TaoD1DataProxy>;

using TaoD1DataProxyAlloc1D = FTypeAlloc1D<
    TaoD1DataProxyArray1D,
    allocate_tao_d1_data_struct_container,
    deallocate_tao_d1_data_struct_container,
    reallocate_tao_d1_data_struct_container_data,
    access_tao_d1_data_struct_container>;

class TaoD2DataProxy;

using TaoD2DataProxyArray1D = FTypeArray1D<
    TaoD2DataProxy,
    allocate_fortran_tao_d2_data_struct,
    deallocate_fortran_tao_d2_data_struct>;
using TaoD2DataProxyArray2D = FTypeArray2D<TaoD2DataProxy>;
using TaoD2DataProxyArray3D = FTypeArray3D<TaoD2DataProxy>;

using TaoD2DataProxyAlloc1D = FTypeAlloc1D<
    TaoD2DataProxyArray1D,
    allocate_tao_d2_data_struct_container,
    deallocate_tao_d2_data_struct_container,
    reallocate_tao_d2_data_struct_container_data,
    access_tao_d2_data_struct_container>;

class TaoDataVarComponentProxy;

using TaoDataVarComponentProxyArray1D = FTypeArray1D<
    TaoDataVarComponentProxy,
    allocate_fortran_tao_data_var_component_struct,
    deallocate_fortran_tao_data_var_component_struct>;
using TaoDataVarComponentProxyArray2D = FTypeArray2D<TaoDataVarComponentProxy>;
using TaoDataVarComponentProxyArray3D = FTypeArray3D<TaoDataVarComponentProxy>;

using TaoDataVarComponentProxyAlloc1D = FTypeAlloc1D<
    TaoDataVarComponentProxyArray1D,
    allocate_tao_data_var_component_struct_container,
    deallocate_tao_data_var_component_struct_container,
    reallocate_tao_data_var_component_struct_container_data,
    access_tao_data_var_component_struct_container>;

class TaoGraphProxy;

using TaoGraphProxyArray1D = FTypeArray1D<
    TaoGraphProxy,
    allocate_fortran_tao_graph_struct,
    deallocate_fortran_tao_graph_struct>;
using TaoGraphProxyArray2D = FTypeArray2D<TaoGraphProxy>;
using TaoGraphProxyArray3D = FTypeArray3D<TaoGraphProxy>;

using TaoGraphProxyAlloc1D = FTypeAlloc1D<
    TaoGraphProxyArray1D,
    allocate_tao_graph_struct_container,
    deallocate_tao_graph_struct_container,
    reallocate_tao_graph_struct_container_data,
    access_tao_graph_struct_container>;

class TaoPlotProxy;

using TaoPlotProxyArray1D = FTypeArray1D<
    TaoPlotProxy,
    allocate_fortran_tao_plot_struct,
    deallocate_fortran_tao_plot_struct>;
using TaoPlotProxyArray2D = FTypeArray2D<TaoPlotProxy>;
using TaoPlotProxyArray3D = FTypeArray3D<TaoPlotProxy>;

using TaoPlotProxyAlloc1D = FTypeAlloc1D<
    TaoPlotProxyArray1D,
    allocate_tao_plot_struct_container,
    deallocate_tao_plot_struct_container,
    reallocate_tao_plot_struct_container_data,
    access_tao_plot_struct_container>;

class TaoPlotRegionProxy;

using TaoPlotRegionProxyArray1D = FTypeArray1D<
    TaoPlotRegionProxy,
    allocate_fortran_tao_plot_region_struct,
    deallocate_fortran_tao_plot_region_struct>;
using TaoPlotRegionProxyArray2D = FTypeArray2D<TaoPlotRegionProxy>;
using TaoPlotRegionProxyArray3D = FTypeArray3D<TaoPlotRegionProxy>;

using TaoPlotRegionProxyAlloc1D = FTypeAlloc1D<
    TaoPlotRegionProxyArray1D,
    allocate_tao_plot_region_struct_container,
    deallocate_tao_plot_region_struct_container,
    reallocate_tao_plot_region_struct_container_data,
    access_tao_plot_region_struct_container>;

class TaoUniversePointerProxy;

using TaoUniversePointerProxyArray1D = FTypeArray1D<
    TaoUniversePointerProxy,
    allocate_fortran_tao_universe_pointer_struct,
    deallocate_fortran_tao_universe_pointer_struct>;
using TaoUniversePointerProxyArray2D = FTypeArray2D<TaoUniversePointerProxy>;
using TaoUniversePointerProxyArray3D = FTypeArray3D<TaoUniversePointerProxy>;

using TaoUniversePointerProxyAlloc1D = FTypeAlloc1D<
    TaoUniversePointerProxyArray1D,
    allocate_tao_universe_pointer_struct_container,
    deallocate_tao_universe_pointer_struct_container,
    reallocate_tao_universe_pointer_struct_container_data,
    access_tao_universe_pointer_struct_container>;

class TaoSuperUniverseProxy;

using TaoSuperUniverseProxyArray1D = FTypeArray1D<
    TaoSuperUniverseProxy,
    allocate_fortran_tao_super_universe_struct,
    deallocate_fortran_tao_super_universe_struct>;
using TaoSuperUniverseProxyArray2D = FTypeArray2D<TaoSuperUniverseProxy>;
using TaoSuperUniverseProxyArray3D = FTypeArray3D<TaoSuperUniverseProxy>;

using TaoSuperUniverseProxyAlloc1D = FTypeAlloc1D<
    TaoSuperUniverseProxyArray1D,
    allocate_tao_super_universe_struct_container,
    deallocate_tao_super_universe_struct_container,
    reallocate_tao_super_universe_struct_container_data,
    access_tao_super_universe_struct_container>;

class TaoVarProxy;

using TaoVarProxyArray1D = FTypeArray1D<
    TaoVarProxy,
    allocate_fortran_tao_var_struct,
    deallocate_fortran_tao_var_struct>;
using TaoVarProxyArray2D = FTypeArray2D<TaoVarProxy>;
using TaoVarProxyArray3D = FTypeArray3D<TaoVarProxy>;

using TaoVarProxyAlloc1D = FTypeAlloc1D<
    TaoVarProxyArray1D,
    allocate_tao_var_struct_container,
    deallocate_tao_var_struct_container,
    reallocate_tao_var_struct_container_data,
    access_tao_var_struct_container>;

class TaoVarSlaveProxy;

using TaoVarSlaveProxyArray1D = FTypeArray1D<
    TaoVarSlaveProxy,
    allocate_fortran_tao_var_slave_struct,
    deallocate_fortran_tao_var_slave_struct>;
using TaoVarSlaveProxyArray2D = FTypeArray2D<TaoVarSlaveProxy>;
using TaoVarSlaveProxyArray3D = FTypeArray3D<TaoVarSlaveProxy>;

using TaoVarSlaveProxyAlloc1D = FTypeAlloc1D<
    TaoVarSlaveProxyArray1D,
    allocate_tao_var_slave_struct_container,
    deallocate_tao_var_slave_struct_container,
    reallocate_tao_var_slave_struct_container_data,
    access_tao_var_slave_struct_container>;

class TaoLatticeProxy;

using TaoLatticeProxyArray1D = FTypeArray1D<
    TaoLatticeProxy,
    allocate_fortran_tao_lattice_struct,
    deallocate_fortran_tao_lattice_struct>;
using TaoLatticeProxyArray2D = FTypeArray2D<TaoLatticeProxy>;
using TaoLatticeProxyArray3D = FTypeArray3D<TaoLatticeProxy>;

using TaoLatticeProxyAlloc1D = FTypeAlloc1D<
    TaoLatticeProxyArray1D,
    allocate_tao_lattice_struct_container,
    deallocate_tao_lattice_struct_container,
    reallocate_tao_lattice_struct_container_data,
    access_tao_lattice_struct_container>;

class TaoBeamUniProxy;

using TaoBeamUniProxyArray1D = FTypeArray1D<
    TaoBeamUniProxy,
    allocate_fortran_tao_beam_uni_struct,
    deallocate_fortran_tao_beam_uni_struct>;
using TaoBeamUniProxyArray2D = FTypeArray2D<TaoBeamUniProxy>;
using TaoBeamUniProxyArray3D = FTypeArray3D<TaoBeamUniProxy>;

using TaoBeamUniProxyAlloc1D = FTypeAlloc1D<
    TaoBeamUniProxyArray1D,
    allocate_tao_beam_uni_struct_container,
    deallocate_tao_beam_uni_struct_container,
    reallocate_tao_beam_uni_struct_container_data,
    access_tao_beam_uni_struct_container>;

class TaoDynamicApertureProxy;

using TaoDynamicApertureProxyArray1D = FTypeArray1D<
    TaoDynamicApertureProxy,
    allocate_fortran_tao_dynamic_aperture_struct,
    deallocate_fortran_tao_dynamic_aperture_struct>;
using TaoDynamicApertureProxyArray2D = FTypeArray2D<TaoDynamicApertureProxy>;
using TaoDynamicApertureProxyArray3D = FTypeArray3D<TaoDynamicApertureProxy>;

using TaoDynamicApertureProxyAlloc1D = FTypeAlloc1D<
    TaoDynamicApertureProxyArray1D,
    allocate_tao_dynamic_aperture_struct_container,
    deallocate_tao_dynamic_aperture_struct_container,
    reallocate_tao_dynamic_aperture_struct_container_data,
    access_tao_dynamic_aperture_struct_container>;

class TaoModelBranchProxy;

using TaoModelBranchProxyArray1D = FTypeArray1D<
    TaoModelBranchProxy,
    allocate_fortran_tao_model_branch_struct,
    deallocate_fortran_tao_model_branch_struct>;
using TaoModelBranchProxyArray2D = FTypeArray2D<TaoModelBranchProxy>;
using TaoModelBranchProxyArray3D = FTypeArray3D<TaoModelBranchProxy>;

using TaoModelBranchProxyAlloc1D = FTypeAlloc1D<
    TaoModelBranchProxyArray1D,
    allocate_tao_model_branch_struct_container,
    deallocate_tao_model_branch_struct_container,
    reallocate_tao_model_branch_struct_container_data,
    access_tao_model_branch_struct_container>;

class TaoSpinMapProxy;

using TaoSpinMapProxyArray1D = FTypeArray1D<
    TaoSpinMapProxy,
    allocate_fortran_tao_spin_map_struct,
    deallocate_fortran_tao_spin_map_struct>;
using TaoSpinMapProxyArray2D = FTypeArray2D<TaoSpinMapProxy>;
using TaoSpinMapProxyArray3D = FTypeArray3D<TaoSpinMapProxy>;

using TaoSpinMapProxyAlloc1D = FTypeAlloc1D<
    TaoSpinMapProxyArray1D,
    allocate_tao_spin_map_struct_container,
    deallocate_tao_spin_map_struct_container,
    reallocate_tao_spin_map_struct_container_data,
    access_tao_spin_map_struct_container>;

class TaoDataProxy;

using TaoDataProxyArray1D = FTypeArray1D<
    TaoDataProxy,
    allocate_fortran_tao_data_struct,
    deallocate_fortran_tao_data_struct>;
using TaoDataProxyArray2D = FTypeArray2D<TaoDataProxy>;
using TaoDataProxyArray3D = FTypeArray3D<TaoDataProxy>;

using TaoDataProxyAlloc1D = FTypeAlloc1D<
    TaoDataProxyArray1D,
    allocate_tao_data_struct_container,
    deallocate_tao_data_struct_container,
    reallocate_tao_data_struct_container_data,
    access_tao_data_struct_container>;

class TaoPingScaleProxy;

using TaoPingScaleProxyArray1D = FTypeArray1D<
    TaoPingScaleProxy,
    allocate_fortran_tao_ping_scale_struct,
    deallocate_fortran_tao_ping_scale_struct>;
using TaoPingScaleProxyArray2D = FTypeArray2D<TaoPingScaleProxy>;
using TaoPingScaleProxyArray3D = FTypeArray3D<TaoPingScaleProxy>;

using TaoPingScaleProxyAlloc1D = FTypeAlloc1D<
    TaoPingScaleProxyArray1D,
    allocate_tao_ping_scale_struct_container,
    deallocate_tao_ping_scale_struct_container,
    reallocate_tao_ping_scale_struct_container_data,
    access_tao_ping_scale_struct_container>;

class TaoUniverseCalcProxy;

using TaoUniverseCalcProxyArray1D = FTypeArray1D<
    TaoUniverseCalcProxy,
    allocate_fortran_tao_universe_calc_struct,
    deallocate_fortran_tao_universe_calc_struct>;
using TaoUniverseCalcProxyArray2D = FTypeArray2D<TaoUniverseCalcProxy>;
using TaoUniverseCalcProxyArray3D = FTypeArray3D<TaoUniverseCalcProxy>;

using TaoUniverseCalcProxyAlloc1D = FTypeAlloc1D<
    TaoUniverseCalcProxyArray1D,
    allocate_tao_universe_calc_struct_container,
    deallocate_tao_universe_calc_struct_container,
    reallocate_tao_universe_calc_struct_container_data,
    access_tao_universe_calc_struct_container>;

class LatEleOrderProxy;

using LatEleOrderProxyArray1D = FTypeArray1D<
    LatEleOrderProxy,
    allocate_fortran_lat_ele_order_struct,
    deallocate_fortran_lat_ele_order_struct>;
using LatEleOrderProxyArray2D = FTypeArray2D<LatEleOrderProxy>;
using LatEleOrderProxyArray3D = FTypeArray3D<LatEleOrderProxy>;

using LatEleOrderProxyAlloc1D = FTypeAlloc1D<
    LatEleOrderProxyArray1D,
    allocate_lat_ele_order_struct_container,
    deallocate_lat_ele_order_struct_container,
    reallocate_lat_ele_order_struct_container_data,
    access_lat_ele_order_struct_container>;

class TaoTitleProxy;

using TaoTitleProxyArray1D = FTypeArray1D<
    TaoTitleProxy,
    allocate_fortran_tao_title_struct,
    deallocate_fortran_tao_title_struct>;
using TaoTitleProxyArray2D = FTypeArray2D<TaoTitleProxy>;
using TaoTitleProxyArray3D = FTypeArray3D<TaoTitleProxy>;

using TaoTitleProxyAlloc1D = FTypeAlloc1D<
    TaoTitleProxyArray1D,
    allocate_tao_title_struct_container,
    deallocate_tao_title_struct_container,
    reallocate_tao_title_struct_container_data,
    access_tao_title_struct_container>;

class QpRectProxy;

using QpRectProxyArray1D = FTypeArray1D<
    QpRectProxy,
    allocate_fortran_qp_rect_struct,
    deallocate_fortran_qp_rect_struct>;
using QpRectProxyArray2D = FTypeArray2D<QpRectProxy>;
using QpRectProxyArray3D = FTypeArray3D<QpRectProxy>;

using QpRectProxyAlloc1D = FTypeAlloc1D<
    QpRectProxyArray1D,
    allocate_qp_rect_struct_container,
    deallocate_qp_rect_struct_container,
    reallocate_qp_rect_struct_container_data,
    access_qp_rect_struct_container>;

class TaoDrawingProxy;

using TaoDrawingProxyArray1D = FTypeArray1D<
    TaoDrawingProxy,
    allocate_fortran_tao_drawing_struct,
    deallocate_fortran_tao_drawing_struct>;
using TaoDrawingProxyArray2D = FTypeArray2D<TaoDrawingProxy>;
using TaoDrawingProxyArray3D = FTypeArray3D<TaoDrawingProxy>;

using TaoDrawingProxyAlloc1D = FTypeAlloc1D<
    TaoDrawingProxyArray1D,
    allocate_tao_drawing_struct_container,
    deallocate_tao_drawing_struct_container,
    reallocate_tao_drawing_struct_container_data,
    access_tao_drawing_struct_container>;

class TaoShapePatternProxy;

using TaoShapePatternProxyArray1D = FTypeArray1D<
    TaoShapePatternProxy,
    allocate_fortran_tao_shape_pattern_struct,
    deallocate_fortran_tao_shape_pattern_struct>;
using TaoShapePatternProxyArray2D = FTypeArray2D<TaoShapePatternProxy>;
using TaoShapePatternProxyArray3D = FTypeArray3D<TaoShapePatternProxy>;

using TaoShapePatternProxyAlloc1D = FTypeAlloc1D<
    TaoShapePatternProxyArray1D,
    allocate_tao_shape_pattern_struct_container,
    deallocate_tao_shape_pattern_struct_container,
    reallocate_tao_shape_pattern_struct_container_data,
    access_tao_shape_pattern_struct_container>;

class TaoShapePatternPointProxy;

using TaoShapePatternPointProxyArray1D = FTypeArray1D<
    TaoShapePatternPointProxy,
    allocate_fortran_tao_shape_pattern_point_struct,
    deallocate_fortran_tao_shape_pattern_point_struct>;
using TaoShapePatternPointProxyArray2D =
    FTypeArray2D<TaoShapePatternPointProxy>;
using TaoShapePatternPointProxyArray3D =
    FTypeArray3D<TaoShapePatternPointProxy>;

using TaoShapePatternPointProxyAlloc1D = FTypeAlloc1D<
    TaoShapePatternPointProxyArray1D,
    allocate_tao_shape_pattern_point_struct_container,
    deallocate_tao_shape_pattern_point_struct_container,
    reallocate_tao_shape_pattern_point_struct_container_data,
    access_tao_shape_pattern_point_struct_container>;

class QpAxisProxy;

using QpAxisProxyArray1D = FTypeArray1D<
    QpAxisProxy,
    allocate_fortran_qp_axis_struct,
    deallocate_fortran_qp_axis_struct>;
using QpAxisProxyArray2D = FTypeArray2D<QpAxisProxy>;
using QpAxisProxyArray3D = FTypeArray3D<QpAxisProxy>;

using QpAxisProxyAlloc1D = FTypeAlloc1D<
    QpAxisProxyArray1D,
    allocate_qp_axis_struct_container,
    deallocate_qp_axis_struct_container,
    reallocate_qp_axis_struct_container_data,
    access_qp_axis_struct_container>;

class QpLegendProxy;

using QpLegendProxyArray1D = FTypeArray1D<
    QpLegendProxy,
    allocate_fortran_qp_legend_struct,
    deallocate_fortran_qp_legend_struct>;
using QpLegendProxyArray2D = FTypeArray2D<QpLegendProxy>;
using QpLegendProxyArray3D = FTypeArray3D<QpLegendProxy>;

using QpLegendProxyAlloc1D = FTypeAlloc1D<
    QpLegendProxyArray1D,
    allocate_qp_legend_struct_container,
    deallocate_qp_legend_struct_container,
    reallocate_qp_legend_struct_container_data,
    access_qp_legend_struct_container>;

class QpPointProxy;

using QpPointProxyArray1D = FTypeArray1D<
    QpPointProxy,
    allocate_fortran_qp_point_struct,
    deallocate_fortran_qp_point_struct>;
using QpPointProxyArray2D = FTypeArray2D<QpPointProxy>;
using QpPointProxyArray3D = FTypeArray3D<QpPointProxy>;

using QpPointProxyAlloc1D = FTypeAlloc1D<
    QpPointProxyArray1D,
    allocate_qp_point_struct_container,
    deallocate_qp_point_struct_container,
    reallocate_qp_point_struct_container_data,
    access_qp_point_struct_container>;

class QpLineProxy;

using QpLineProxyArray1D = FTypeArray1D<
    QpLineProxy,
    allocate_fortran_qp_line_struct,
    deallocate_fortran_qp_line_struct>;
using QpLineProxyArray2D = FTypeArray2D<QpLineProxy>;
using QpLineProxyArray3D = FTypeArray3D<QpLineProxy>;

using QpLineProxyAlloc1D = FTypeAlloc1D<
    QpLineProxyArray1D,
    allocate_qp_line_struct_container,
    deallocate_qp_line_struct_container,
    reallocate_qp_line_struct_container_data,
    access_qp_line_struct_container>;

class QpSymbolProxy;

using QpSymbolProxyArray1D = FTypeArray1D<
    QpSymbolProxy,
    allocate_fortran_qp_symbol_struct,
    deallocate_fortran_qp_symbol_struct>;
using QpSymbolProxyArray2D = FTypeArray2D<QpSymbolProxy>;
using QpSymbolProxyArray3D = FTypeArray3D<QpSymbolProxy>;

using QpSymbolProxyAlloc1D = FTypeAlloc1D<
    QpSymbolProxyArray1D,
    allocate_qp_symbol_struct_container,
    deallocate_qp_symbol_struct_container,
    reallocate_qp_symbol_struct_container_data,
    access_qp_symbol_struct_container>;

class TaoFloorPlanProxy;

using TaoFloorPlanProxyArray1D = FTypeArray1D<
    TaoFloorPlanProxy,
    allocate_fortran_tao_floor_plan_struct,
    deallocate_fortran_tao_floor_plan_struct>;
using TaoFloorPlanProxyArray2D = FTypeArray2D<TaoFloorPlanProxy>;
using TaoFloorPlanProxyArray3D = FTypeArray3D<TaoFloorPlanProxy>;

using TaoFloorPlanProxyAlloc1D = FTypeAlloc1D<
    TaoFloorPlanProxyArray1D,
    allocate_tao_floor_plan_struct_container,
    deallocate_tao_floor_plan_struct_container,
    reallocate_tao_floor_plan_struct_container_data,
    access_tao_floor_plan_struct_container>;

class TaoV1VarProxy;

using TaoV1VarProxyArray1D = FTypeArray1D<
    TaoV1VarProxy,
    allocate_fortran_tao_v1_var_struct,
    deallocate_fortran_tao_v1_var_struct>;
using TaoV1VarProxyArray2D = FTypeArray2D<TaoV1VarProxy>;
using TaoV1VarProxyArray3D = FTypeArray3D<TaoV1VarProxy>;

using TaoV1VarProxyAlloc1D = FTypeAlloc1D<
    TaoV1VarProxyArray1D,
    allocate_tao_v1_var_struct_container,
    deallocate_tao_v1_var_struct_container,
    reallocate_tao_v1_var_struct_container_data,
    access_tao_v1_var_struct_container>;

class TaoGlobalProxy;

using TaoGlobalProxyArray1D = FTypeArray1D<
    TaoGlobalProxy,
    allocate_fortran_tao_global_struct,
    deallocate_fortran_tao_global_struct>;
using TaoGlobalProxyArray2D = FTypeArray2D<TaoGlobalProxy>;
using TaoGlobalProxyArray3D = FTypeArray3D<TaoGlobalProxy>;

using TaoGlobalProxyAlloc1D = FTypeAlloc1D<
    TaoGlobalProxyArray1D,
    allocate_tao_global_struct_container,
    deallocate_tao_global_struct_container,
    reallocate_tao_global_struct_container_data,
    access_tao_global_struct_container>;

class TaoInitProxy;

using TaoInitProxyArray1D = FTypeArray1D<
    TaoInitProxy,
    allocate_fortran_tao_init_struct,
    deallocate_fortran_tao_init_struct>;
using TaoInitProxyArray2D = FTypeArray2D<TaoInitProxy>;
using TaoInitProxyArray3D = FTypeArray3D<TaoInitProxy>;

using TaoInitProxyAlloc1D = FTypeAlloc1D<
    TaoInitProxyArray1D,
    allocate_tao_init_struct_container,
    deallocate_tao_init_struct_container,
    reallocate_tao_init_struct_container_data,
    access_tao_init_struct_container>;

class TaoCommonProxy;

using TaoCommonProxyArray1D = FTypeArray1D<
    TaoCommonProxy,
    allocate_fortran_tao_common_struct,
    deallocate_fortran_tao_common_struct>;
using TaoCommonProxyArray2D = FTypeArray2D<TaoCommonProxy>;
using TaoCommonProxyArray3D = FTypeArray3D<TaoCommonProxy>;

using TaoCommonProxyAlloc1D = FTypeAlloc1D<
    TaoCommonProxyArray1D,
    allocate_tao_common_struct_container,
    deallocate_tao_common_struct_container,
    reallocate_tao_common_struct_container_data,
    access_tao_common_struct_container>;

class TaoPlotPageProxy;

using TaoPlotPageProxyArray1D = FTypeArray1D<
    TaoPlotPageProxy,
    allocate_fortran_tao_plot_page_struct,
    deallocate_fortran_tao_plot_page_struct>;
using TaoPlotPageProxyArray2D = FTypeArray2D<TaoPlotPageProxy>;
using TaoPlotPageProxyArray3D = FTypeArray3D<TaoPlotPageProxy>;

using TaoPlotPageProxyAlloc1D = FTypeAlloc1D<
    TaoPlotPageProxyArray1D,
    allocate_tao_plot_page_struct_container,
    deallocate_tao_plot_page_struct_container,
    reallocate_tao_plot_page_struct_container_data,
    access_tao_plot_page_struct_container>;

class TaoBuildingWallProxy;

using TaoBuildingWallProxyArray1D = FTypeArray1D<
    TaoBuildingWallProxy,
    allocate_fortran_tao_building_wall_struct,
    deallocate_fortran_tao_building_wall_struct>;
using TaoBuildingWallProxyArray2D = FTypeArray2D<TaoBuildingWallProxy>;
using TaoBuildingWallProxyArray3D = FTypeArray3D<TaoBuildingWallProxy>;

using TaoBuildingWallProxyAlloc1D = FTypeAlloc1D<
    TaoBuildingWallProxyArray1D,
    allocate_tao_building_wall_struct_container,
    deallocate_tao_building_wall_struct_container,
    reallocate_tao_building_wall_struct_container_data,
    access_tao_building_wall_struct_container>;

class TaoBuildingWallOrientationProxy;

using TaoBuildingWallOrientationProxyArray1D = FTypeArray1D<
    TaoBuildingWallOrientationProxy,
    allocate_fortran_tao_building_wall_orientation_struct,
    deallocate_fortran_tao_building_wall_orientation_struct>;
using TaoBuildingWallOrientationProxyArray2D =
    FTypeArray2D<TaoBuildingWallOrientationProxy>;
using TaoBuildingWallOrientationProxyArray3D =
    FTypeArray3D<TaoBuildingWallOrientationProxy>;

using TaoBuildingWallOrientationProxyAlloc1D = FTypeAlloc1D<
    TaoBuildingWallOrientationProxyArray1D,
    allocate_tao_building_wall_orientation_struct_container,
    deallocate_tao_building_wall_orientation_struct_container,
    reallocate_tao_building_wall_orientation_struct_container_data,
    access_tao_building_wall_orientation_struct_container>;

class TaoBuildingWallSectionProxy;

using TaoBuildingWallSectionProxyArray1D = FTypeArray1D<
    TaoBuildingWallSectionProxy,
    allocate_fortran_tao_building_wall_section_struct,
    deallocate_fortran_tao_building_wall_section_struct>;
using TaoBuildingWallSectionProxyArray2D =
    FTypeArray2D<TaoBuildingWallSectionProxy>;
using TaoBuildingWallSectionProxyArray3D =
    FTypeArray3D<TaoBuildingWallSectionProxy>;

using TaoBuildingWallSectionProxyAlloc1D = FTypeAlloc1D<
    TaoBuildingWallSectionProxyArray1D,
    allocate_tao_building_wall_section_struct_container,
    deallocate_tao_building_wall_section_struct_container,
    reallocate_tao_building_wall_section_struct_container_data,
    access_tao_building_wall_section_struct_container>;

class TaoBuildingWallPointProxy;

using TaoBuildingWallPointProxyArray1D = FTypeArray1D<
    TaoBuildingWallPointProxy,
    allocate_fortran_tao_building_wall_point_struct,
    deallocate_fortran_tao_building_wall_point_struct>;
using TaoBuildingWallPointProxyArray2D =
    FTypeArray2D<TaoBuildingWallPointProxy>;
using TaoBuildingWallPointProxyArray3D =
    FTypeArray3D<TaoBuildingWallPointProxy>;

using TaoBuildingWallPointProxyAlloc1D = FTypeAlloc1D<
    TaoBuildingWallPointProxyArray1D,
    allocate_tao_building_wall_point_struct_container,
    deallocate_tao_building_wall_point_struct_container,
    reallocate_tao_building_wall_point_struct_container_data,
    access_tao_building_wall_point_struct_container>;

class TaoWaveProxy;

using TaoWaveProxyArray1D = FTypeArray1D<
    TaoWaveProxy,
    allocate_fortran_tao_wave_struct,
    deallocate_fortran_tao_wave_struct>;
using TaoWaveProxyArray2D = FTypeArray2D<TaoWaveProxy>;
using TaoWaveProxyArray3D = FTypeArray3D<TaoWaveProxy>;

using TaoWaveProxyAlloc1D = FTypeAlloc1D<
    TaoWaveProxyArray1D,
    allocate_tao_wave_struct_container,
    deallocate_tao_wave_struct_container,
    reallocate_tao_wave_struct_container_data,
    access_tao_wave_struct_container>;

class TaoWaveKickPtProxy;

using TaoWaveKickPtProxyArray1D = FTypeArray1D<
    TaoWaveKickPtProxy,
    allocate_fortran_tao_wave_kick_pt_struct,
    deallocate_fortran_tao_wave_kick_pt_struct>;
using TaoWaveKickPtProxyArray2D = FTypeArray2D<TaoWaveKickPtProxy>;
using TaoWaveKickPtProxyArray3D = FTypeArray3D<TaoWaveKickPtProxy>;

using TaoWaveKickPtProxyAlloc1D = FTypeAlloc1D<
    TaoWaveKickPtProxyArray1D,
    allocate_tao_wave_kick_pt_struct_container,
    deallocate_tao_wave_kick_pt_struct_container,
    reallocate_tao_wave_kick_pt_struct_container_data,
    access_tao_wave_kick_pt_struct_container>;

class TaoCmdHistoryProxy;

using TaoCmdHistoryProxyArray1D = FTypeArray1D<
    TaoCmdHistoryProxy,
    allocate_fortran_tao_cmd_history_struct,
    deallocate_fortran_tao_cmd_history_struct>;
using TaoCmdHistoryProxyArray2D = FTypeArray2D<TaoCmdHistoryProxy>;
using TaoCmdHistoryProxyArray3D = FTypeArray3D<TaoCmdHistoryProxy>;

using TaoCmdHistoryProxyAlloc1D = FTypeAlloc1D<
    TaoCmdHistoryProxyArray1D,
    allocate_tao_cmd_history_struct_container,
    deallocate_tao_cmd_history_struct_container,
    reallocate_tao_cmd_history_struct_container_data,
    access_tao_cmd_history_struct_container>;

class TaoUniverseProxy;

using TaoUniverseProxyArray1D = FTypeArray1D<
    TaoUniverseProxy,
    allocate_fortran_tao_universe_struct,
    deallocate_fortran_tao_universe_struct>;
using TaoUniverseProxyArray2D = FTypeArray2D<TaoUniverseProxy>;
using TaoUniverseProxyArray3D = FTypeArray3D<TaoUniverseProxy>;

using TaoUniverseProxyAlloc1D = FTypeAlloc1D<
    TaoUniverseProxyArray1D,
    allocate_tao_universe_struct_container,
    deallocate_tao_universe_struct_container,
    reallocate_tao_universe_struct_container_data,
    access_tao_universe_struct_container>;

class MadEnergyProxy;

using MadEnergyProxyArray1D = FTypeArray1D<
    MadEnergyProxy,
    allocate_fortran_mad_energy_struct,
    deallocate_fortran_mad_energy_struct>;
using MadEnergyProxyArray2D = FTypeArray2D<MadEnergyProxy>;
using MadEnergyProxyArray3D = FTypeArray3D<MadEnergyProxy>;

using MadEnergyProxyAlloc1D = FTypeAlloc1D<
    MadEnergyProxyArray1D,
    allocate_mad_energy_struct_container,
    deallocate_mad_energy_struct_container,
    reallocate_mad_energy_struct_container_data,
    access_mad_energy_struct_container>;

class MadMapProxy;

using MadMapProxyArray1D = FTypeArray1D<
    MadMapProxy,
    allocate_fortran_mad_map_struct,
    deallocate_fortran_mad_map_struct>;
using MadMapProxyArray2D = FTypeArray2D<MadMapProxy>;
using MadMapProxyArray3D = FTypeArray3D<MadMapProxy>;

using MadMapProxyAlloc1D = FTypeAlloc1D<
    MadMapProxyArray1D,
    allocate_mad_map_struct_container,
    deallocate_mad_map_struct_container,
    reallocate_mad_map_struct_container_data,
    access_mad_map_struct_container>;

class RandomStateProxy;

using RandomStateProxyArray1D = FTypeArray1D<
    RandomStateProxy,
    allocate_fortran_random_state_struct,
    deallocate_fortran_random_state_struct>;
using RandomStateProxyArray2D = FTypeArray2D<RandomStateProxy>;
using RandomStateProxyArray3D = FTypeArray3D<RandomStateProxy>;

using RandomStateProxyAlloc1D = FTypeAlloc1D<
    RandomStateProxyArray1D,
    allocate_random_state_struct_container,
    deallocate_random_state_struct_container,
    reallocate_random_state_struct_container_data,
    access_random_state_struct_container>;

class BbuStageProxy;

using BbuStageProxyArray1D = FTypeArray1D<
    BbuStageProxy,
    allocate_fortran_bbu_stage_struct,
    deallocate_fortran_bbu_stage_struct>;
using BbuStageProxyArray2D = FTypeArray2D<BbuStageProxy>;
using BbuStageProxyArray3D = FTypeArray3D<BbuStageProxy>;

using BbuStageProxyAlloc1D = FTypeAlloc1D<
    BbuStageProxyArray1D,
    allocate_bbu_stage_struct_container,
    deallocate_bbu_stage_struct_container,
    reallocate_bbu_stage_struct_container_data,
    access_bbu_stage_struct_container>;

class BbuBeamProxy;

using BbuBeamProxyArray1D = FTypeArray1D<
    BbuBeamProxy,
    allocate_fortran_bbu_beam_struct,
    deallocate_fortran_bbu_beam_struct>;
using BbuBeamProxyArray2D = FTypeArray2D<BbuBeamProxy>;
using BbuBeamProxyArray3D = FTypeArray3D<BbuBeamProxy>;

using BbuBeamProxyAlloc1D = FTypeAlloc1D<
    BbuBeamProxyArray1D,
    allocate_bbu_beam_struct_container,
    deallocate_bbu_beam_struct_container,
    reallocate_bbu_beam_struct_container_data,
    access_bbu_beam_struct_container>;

class BbuParamProxy;

using BbuParamProxyArray1D = FTypeArray1D<
    BbuParamProxy,
    allocate_fortran_bbu_param_struct,
    deallocate_fortran_bbu_param_struct>;
using BbuParamProxyArray2D = FTypeArray2D<BbuParamProxy>;
using BbuParamProxyArray3D = FTypeArray3D<BbuParamProxy>;

using BbuParamProxyAlloc1D = FTypeAlloc1D<
    BbuParamProxyArray1D,
    allocate_bbu_param_struct_container,
    deallocate_bbu_param_struct_container,
    reallocate_bbu_param_struct_container_data,
    access_bbu_param_struct_container>;

class AllEncompassingProxy;

using AllEncompassingProxyArray1D = FTypeArray1D<
    AllEncompassingProxy,
    allocate_fortran_all_encompassing_struct,
    deallocate_fortran_all_encompassing_struct>;
using AllEncompassingProxyArray2D = FTypeArray2D<AllEncompassingProxy>;
using AllEncompassingProxyArray3D = FTypeArray3D<AllEncompassingProxy>;

using AllEncompassingProxyAlloc1D = FTypeAlloc1D<
    AllEncompassingProxyArray1D,
    allocate_all_encompassing_struct_container,
    deallocate_all_encompassing_struct_container,
    reallocate_all_encompassing_struct_container_data,
    access_all_encompassing_struct_container>;

class TestSubProxy;

using TestSubProxyArray1D = FTypeArray1D<
    TestSubProxy,
    allocate_fortran_test_sub_struct,
    deallocate_fortran_test_sub_struct>;
using TestSubProxyArray2D = FTypeArray2D<TestSubProxy>;
using TestSubProxyArray3D = FTypeArray3D<TestSubProxy>;

using TestSubProxyAlloc1D = FTypeAlloc1D<
    TestSubProxyArray1D,
    allocate_test_sub_struct_container,
    deallocate_test_sub_struct_container,
    reallocate_test_sub_struct_container_data,
    access_test_sub_struct_container>;

class TestSubSubProxy;

using TestSubSubProxyArray1D = FTypeArray1D<
    TestSubSubProxy,
    allocate_fortran_test_sub_sub_struct,
    deallocate_fortran_test_sub_sub_struct>;
using TestSubSubProxyArray2D = FTypeArray2D<TestSubSubProxy>;
using TestSubSubProxyArray3D = FTypeArray3D<TestSubSubProxy>;

using TestSubSubProxyAlloc1D = FTypeAlloc1D<
    TestSubSubProxyArray1D,
    allocate_test_sub_sub_struct_container,
    deallocate_test_sub_sub_struct_container,
    reallocate_test_sub_sub_struct_container_data,
    access_test_sub_sub_struct_container>;

template <>
struct FortranTraits<SplineProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_spline_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_spline_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_spline_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "spline_struct";
  }
};

class SplineProxy : public FortranProxy<SplineProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double x0() const; // 0D_NOT_real
  void set_x0(double value);
  double y0() const; // 0D_NOT_real
  void set_y0(double value);
  double x1() const; // 0D_NOT_real
  void set_x1(double value);
  FArray1D<double> coef() const; // 1D_NOT_real
};

template <>
struct FortranTraits<SpinPolarProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_spin_polar_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_spin_polar_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_spin_polar_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "spin_polar_struct";
  }
};

class SpinPolarProxy : public FortranProxy<SpinPolarProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double polarization() const; // 0D_NOT_real
  void set_polarization(double value);
  double theta() const; // 0D_NOT_real
  void set_theta(double value);
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
  double xi() const; // 0D_NOT_real
  void set_xi(double value);
};

template <>
struct FortranTraits<AcKickerTimeProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_ac_kicker_time_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_ac_kicker_time_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_ac_kicker_time_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "ac_kicker_time_struct";
  }
};

class AcKickerTimeProxy : public FortranProxy<AcKickerTimeProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double amp() const; // 0D_NOT_real
  void set_amp(double value);
  double time() const; // 0D_NOT_real
  void set_time(double value);
  SplineProxy spline() const; // 0D_NOT_type
  void set_spline(const SplineProxy& src);
};

template <>
struct FortranTraits<AcKickerFreqProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_ac_kicker_freq_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_ac_kicker_freq_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_ac_kicker_freq_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "ac_kicker_freq_struct";
  }
};

class AcKickerFreqProxy : public FortranProxy<AcKickerFreqProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double f() const; // 0D_NOT_real
  void set_f(double value);
  double amp() const; // 0D_NOT_real
  void set_amp(double value);
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
};

template <>
struct FortranTraits<AcKickerProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_ac_kicker_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_ac_kicker_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_ac_kicker_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "ac_kicker_struct";
  }
};

class AcKickerProxy : public FortranProxy<AcKickerProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  AcKickerTimeProxyArray1D amp_vs_time() const; // 1D_ALLOC_type
  AcKickerFreqProxyArray1D frequency() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<Interval1CoefProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_interval1_coef_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_interval1_coef_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_interval1_coef_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "interval1_coef_struct";
  }
};

class Interval1CoefProxy : public FortranProxy<Interval1CoefProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double c0() const; // 0D_NOT_real
  void set_c0(double value);
  double c1() const; // 0D_NOT_real
  void set_c1(double value);
  double n_exp() const; // 0D_NOT_real
  void set_n_exp(double value);
};

template <>
struct FortranTraits<PhotonReflectTableProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_photon_reflect_table_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_photon_reflect_table_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_photon_reflect_table_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "photon_reflect_table_struct";
  }
};

class PhotonReflectTableProxy : public FortranProxy<PhotonReflectTableProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> angle() const; // 1D_ALLOC_real
  FArray1D<double> energy() const; // 1D_ALLOC_real
  Interval1CoefProxyArray1D int1() const; // 1D_ALLOC_type
  FArray2D<double> p_reflect() const; // 2D_ALLOC_real
  double max_energy() const; // 0D_NOT_real
  void set_max_energy(double value);
  FArray1D<double> p_reflect_scratch() const; // 1D_ALLOC_real
  FArray1D<double> bragg_angle() const; // 1D_ALLOC_real
};

template <>
struct FortranTraits<PhotonReflectSurfaceProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_photon_reflect_surface_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_photon_reflect_surface_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_photon_reflect_surface_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "photon_reflect_surface_struct";
  }
};

class PhotonReflectSurfaceProxy
    : public FortranProxy<PhotonReflectSurfaceProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  std::string description() const; // 0D_NOT_character
  void set_description(const std::string& value);
  std::string reflectivity_file() const; // 0D_NOT_character
  void set_reflectivity_file(const std::string& value);
  PhotonReflectTableProxyArray1D table() const; // 1D_ALLOC_type
  double surface_roughness_rms() const; // 0D_NOT_real
  void set_surface_roughness_rms(double value);
  double roughness_correlation_len() const; // 0D_NOT_real
  void set_roughness_correlation_len(double value);
  int ix_surface() const; // 0D_NOT_integer
  void set_ix_surface(int value);
};

template <>
struct FortranTraits<CoordProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_coord_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_coord_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_coord_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "coord_struct";
  }
};

class CoordProxy : public FortranProxy<CoordProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> vec() const; // 1D_NOT_real
  double s() const; // 0D_NOT_real
  void set_s(double value);
  long double t() const; // 0D_NOT_real16
  void set_t(long double value);
  FArray1D<double> spin() const; // 1D_NOT_real
  FArray1D<double> field() const; // 1D_NOT_real
  FArray1D<double> phase() const; // 1D_NOT_real
  double charge() const; // 0D_NOT_real
  void set_charge(double value);
  double dt_ref() const; // 0D_NOT_real
  void set_dt_ref(double value);
  double r() const; // 0D_NOT_real
  void set_r(double value);
  double p0c() const; // 0D_NOT_real
  void set_p0c(double value);
  double E_potential() const; // 0D_NOT_real
  void set_E_potential(double value);
  double beta() const; // 0D_NOT_real
  void set_beta(double value);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_turn() const; // 0D_NOT_integer
  void set_ix_turn(int value);
  int ix_user() const; // 0D_NOT_integer
  void set_ix_user(int value);
  int state() const; // 0D_NOT_integer
  void set_state(int value);
  int direction() const; // 0D_NOT_integer
  void set_direction(int value);
  int time_dir() const; // 0D_NOT_integer
  void set_time_dir(int value);
  int species() const; // 0D_NOT_integer
  void set_species(int value);
  int location() const; // 0D_NOT_integer
  void set_location(int value);
};

template <>
struct FortranTraits<CoordArrayProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_coord_array_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_coord_array_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_coord_array_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "coord_array_struct";
  }
};

class CoordArrayProxy : public FortranProxy<CoordArrayProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  CoordProxyArray1D orbit() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<BpmPhaseCouplingProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bpm_phase_coupling_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bpm_phase_coupling_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bpm_phase_coupling_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bpm_phase_coupling_struct";
  }
};

class BpmPhaseCouplingProxy : public FortranProxy<BpmPhaseCouplingProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double K_22a() const; // 0D_NOT_real
  void set_K_22a(double value);
  double K_12a() const; // 0D_NOT_real
  void set_K_12a(double value);
  double K_11b() const; // 0D_NOT_real
  void set_K_11b(double value);
  double K_12b() const; // 0D_NOT_real
  void set_K_12b(double value);
  double Cbar22_a() const; // 0D_NOT_real
  void set_Cbar22_a(double value);
  double Cbar12_a() const; // 0D_NOT_real
  void set_Cbar12_a(double value);
  double Cbar11_b() const; // 0D_NOT_real
  void set_Cbar11_b(double value);
  double Cbar12_b() const; // 0D_NOT_real
  void set_Cbar12_b(double value);
  double phi_a() const; // 0D_NOT_real
  void set_phi_a(double value);
  double phi_b() const; // 0D_NOT_real
  void set_phi_b(double value);
};

template <>
struct FortranTraits<ExpressionAtomProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_expression_atom_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_expression_atom_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_expression_atom_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "expression_atom_struct";
  }
};

class ExpressionAtomProxy : public FortranProxy<ExpressionAtomProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  int type() const; // 0D_NOT_integer
  void set_type(int value);
  double value() const; // 0D_NOT_real
  void set_value(double value);
};

template <>
struct FortranTraits<WakeSrZLongProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_wake_sr_z_long_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_wake_sr_z_long_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_wake_sr_z_long_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "wake_sr_z_long_struct";
  }
};

class WakeSrZLongProxy : public FortranProxy<WakeSrZLongProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> w() const; // 1D_ALLOC_real
  FArray1D<std::complex<double>> fw() const; // 1D_ALLOC_complex
  FArray1D<std::complex<double>> fbunch() const; // 1D_ALLOC_complex
  FArray1D<std::complex<double>> w_out() const; // 1D_ALLOC_complex
  double dz() const; // 0D_NOT_real
  void set_dz(double value);
  double z0() const; // 0D_NOT_real
  void set_z0(double value);
  double smoothing_sigma() const; // 0D_NOT_real
  void set_smoothing_sigma(double value);
  int position_dependence() const; // 0D_NOT_integer
  void set_position_dependence(int value);
  bool time_based() const; // 0D_NOT_logical
  void set_time_based(bool value);
};

template <>
struct FortranTraits<WakeSrModeProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_wake_sr_mode_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_wake_sr_mode_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_wake_sr_mode_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "wake_sr_mode_struct";
  }
};

class WakeSrModeProxy : public FortranProxy<WakeSrModeProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double amp() const; // 0D_NOT_real
  void set_amp(double value);
  double damp() const; // 0D_NOT_real
  void set_damp(double value);
  double k() const; // 0D_NOT_real
  void set_k(double value);
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
  double b_sin() const; // 0D_NOT_real
  void set_b_sin(double value);
  double b_cos() const; // 0D_NOT_real
  void set_b_cos(double value);
  double a_sin() const; // 0D_NOT_real
  void set_a_sin(double value);
  double a_cos() const; // 0D_NOT_real
  void set_a_cos(double value);
  int polarization() const; // 0D_NOT_integer
  void set_polarization(int value);
  int position_dependence() const; // 0D_NOT_integer
  void set_position_dependence(int value);
};

template <>
struct FortranTraits<WakeSrProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_wake_sr_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_wake_sr_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_wake_sr_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "wake_sr_struct";
  }
};

class WakeSrProxy : public FortranProxy<WakeSrProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string& value);
  WakeSrZLongProxy z_long() const; // 0D_NOT_type
  void set_z_long(const WakeSrZLongProxy& src);
  WakeSrModeProxyArray1D long_wake() const; // 1D_ALLOC_type
  WakeSrModeProxyArray1D trans_wake() const; // 1D_ALLOC_type
  double z_ref_long() const; // 0D_NOT_real
  void set_z_ref_long(double value);
  double z_ref_trans() const; // 0D_NOT_real
  void set_z_ref_trans(double value);
  double z_max() const; // 0D_NOT_real
  void set_z_max(double value);
  double amp_scale() const; // 0D_NOT_real
  void set_amp_scale(double value);
  double z_scale() const; // 0D_NOT_real
  void set_z_scale(double value);
  bool scale_with_length() const; // 0D_NOT_logical
  void set_scale_with_length(bool value);
};

template <>
struct FortranTraits<WakeLrModeProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_wake_lr_mode_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_wake_lr_mode_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_wake_lr_mode_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "wake_lr_mode_struct";
  }
};

class WakeLrModeProxy : public FortranProxy<WakeLrModeProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double freq() const; // 0D_NOT_real
  void set_freq(double value);
  double freq_in() const; // 0D_NOT_real
  void set_freq_in(double value);
  double R_over_Q() const; // 0D_NOT_real
  void set_R_over_Q(double value);
  double Q() const; // 0D_NOT_real
  void set_Q(double value);
  double damp() const; // 0D_NOT_real
  void set_damp(double value);
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
  double angle() const; // 0D_NOT_real
  void set_angle(double value);
  double b_sin() const; // 0D_NOT_real
  void set_b_sin(double value);
  double b_cos() const; // 0D_NOT_real
  void set_b_cos(double value);
  double a_sin() const; // 0D_NOT_real
  void set_a_sin(double value);
  double a_cos() const; // 0D_NOT_real
  void set_a_cos(double value);
  int m() const; // 0D_NOT_integer
  void set_m(int value);
  bool polarized() const; // 0D_NOT_logical
  void set_polarized(bool value);
};

template <>
struct FortranTraits<WakeLrProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_wake_lr_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_wake_lr_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_wake_lr_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "wake_lr_struct";
  }
};

class WakeLrProxy : public FortranProxy<WakeLrProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string& value);
  WakeLrModeProxyArray1D mode() const; // 1D_ALLOC_type
  double t_ref() const; // 0D_NOT_real
  void set_t_ref(double value);
  double freq_spread() const; // 0D_NOT_real
  void set_freq_spread(double value);
  double amp_scale() const; // 0D_NOT_real
  void set_amp_scale(double value);
  double time_scale() const; // 0D_NOT_real
  void set_time_scale(double value);
  bool self_wake_on() const; // 0D_NOT_logical
  void set_self_wake_on(bool value);
};

template <>
struct FortranTraits<LatEleLocProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_lat_ele_loc_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_lat_ele_loc_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_lat_ele_loc_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "lat_ele_loc_struct";
  }
};

class LatEleLocProxy : public FortranProxy<LatEleLocProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
};

template <>
struct FortranTraits<WakeProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_wake_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_wake_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_wake_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "wake_struct";
  }
};

class WakeProxy : public FortranProxy<WakeProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  WakeSrProxy sr() const; // 0D_NOT_type
  void set_sr(const WakeSrProxy& src);
  WakeLrProxy lr() const; // 0D_NOT_type
  void set_lr(const WakeLrProxy& src);
};

template <>
struct FortranTraits<TaylorTermProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_taylor_term_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_taylor_term_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_taylor_term_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "taylor_term_struct";
  }
};

class TaylorTermProxy : public FortranProxy<TaylorTermProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double coef() const; // 0D_NOT_real
  void set_coef(double value);
  FArray1D<int> expn() const; // 1D_NOT_integer
};

template <>
struct FortranTraits<TaylorProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_taylor_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_taylor_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_taylor_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "taylor_struct";
  }
};

class TaylorProxy : public FortranProxy<TaylorProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double ref() const; // 0D_NOT_real
  void set_ref(double value);
  TaylorTermProxyArray1D term() const; // 1D_PTR_type
};

template <>
struct FortranTraits<EmTaylorTermProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_em_taylor_term_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_em_taylor_term_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_em_taylor_term_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "em_taylor_term_struct";
  }
};

class EmTaylorTermProxy : public FortranProxy<EmTaylorTermProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double coef() const; // 0D_NOT_real
  void set_coef(double value);
  FArray1D<int> expn() const; // 1D_NOT_integer
};

template <>
struct FortranTraits<EmTaylorProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_em_taylor_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_em_taylor_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_em_taylor_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "em_taylor_struct";
  }
};

class EmTaylorProxy : public FortranProxy<EmTaylorProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double ref() const; // 0D_NOT_real
  void set_ref(double value);
  EmTaylorTermProxyArray1D term() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<CartesianMapTerm1Proxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_cartesian_map_term1_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_cartesian_map_term1_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_cartesian_map_term1_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "cartesian_map_term1_struct";
  }
};

class CartesianMapTerm1Proxy : public FortranProxy<CartesianMapTerm1Proxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double coef() const; // 0D_NOT_real
  void set_coef(double value);
  double kx() const; // 0D_NOT_real
  void set_kx(double value);
  double ky() const; // 0D_NOT_real
  void set_ky(double value);
  double kz() const; // 0D_NOT_real
  void set_kz(double value);
  double x0() const; // 0D_NOT_real
  void set_x0(double value);
  double y0() const; // 0D_NOT_real
  void set_y0(double value);
  double phi_z() const; // 0D_NOT_real
  void set_phi_z(double value);
  int family() const; // 0D_NOT_integer
  void set_family(int value);
  int form() const; // 0D_NOT_integer
  void set_form(int value);
};

template <>
struct FortranTraits<CartesianMapTermProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_cartesian_map_term_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_cartesian_map_term_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_cartesian_map_term_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "cartesian_map_term_struct";
  }
};

class CartesianMapTermProxy : public FortranProxy<CartesianMapTermProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string& value);
  int n_link() const; // 0D_NOT_integer
  void set_n_link(int value);
  CartesianMapTerm1ProxyArray1D term() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<CartesianMapProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_cartesian_map_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_cartesian_map_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_cartesian_map_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "cartesian_map_struct";
  }
};

class CartesianMapProxy : public FortranProxy<CartesianMapProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double field_scale() const; // 0D_NOT_real
  void set_field_scale(double value);
  FArray1D<double> r0() const; // 1D_NOT_real
  int master_parameter() const; // 0D_NOT_integer
  void set_master_parameter(int value);
  int ele_anchor_pt() const; // 0D_NOT_integer
  void set_ele_anchor_pt(int value);
  int field_type() const; // 0D_NOT_integer
  void set_field_type(int value);
  std::optional<CartesianMapTermProxy> ptr() const; // 0D_PTR_type
  void set_ptr(const CartesianMapTermProxy& src);
};

template <>
struct FortranTraits<CylindricalMapTerm1Proxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_cylindrical_map_term1_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_cylindrical_map_term1_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_cylindrical_map_term1_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "cylindrical_map_term1_struct";
  }
};

class CylindricalMapTerm1Proxy : public FortranProxy<CylindricalMapTerm1Proxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::complex<double> e_coef() const; // 0D_NOT_complex
  void set_e_coef(std::complex<double> value);
  std::complex<double> b_coef() const; // 0D_NOT_complex
  void set_b_coef(std::complex<double> value);
};

template <>
struct FortranTraits<CylindricalMapTermProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_cylindrical_map_term_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_cylindrical_map_term_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_cylindrical_map_term_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "cylindrical_map_term_struct";
  }
};

class CylindricalMapTermProxy : public FortranProxy<CylindricalMapTermProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string& value);
  int n_link() const; // 0D_NOT_integer
  void set_n_link(int value);
  CylindricalMapTerm1ProxyArray1D term() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<CylindricalMapProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_cylindrical_map_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_cylindrical_map_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_cylindrical_map_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "cylindrical_map_struct";
  }
};

class CylindricalMapProxy : public FortranProxy<CylindricalMapProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int m() const; // 0D_NOT_integer
  void set_m(int value);
  int harmonic() const; // 0D_NOT_integer
  void set_harmonic(int value);
  double phi0_fieldmap() const; // 0D_NOT_real
  void set_phi0_fieldmap(double value);
  double theta0_azimuth() const; // 0D_NOT_real
  void set_theta0_azimuth(double value);
  double field_scale() const; // 0D_NOT_real
  void set_field_scale(double value);
  int master_parameter() const; // 0D_NOT_integer
  void set_master_parameter(int value);
  int ele_anchor_pt() const; // 0D_NOT_integer
  void set_ele_anchor_pt(int value);
  double dz() const; // 0D_NOT_real
  void set_dz(double value);
  FArray1D<double> r0() const; // 1D_NOT_real
  std::optional<CylindricalMapTermProxy> ptr() const; // 0D_PTR_type
  void set_ptr(const CylindricalMapTermProxy& src);
};

template <>
struct FortranTraits<BicubicCmplxCoefProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bicubic_cmplx_coef_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bicubic_cmplx_coef_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bicubic_cmplx_coef_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bicubic_cmplx_coef_struct";
  }
};

class BicubicCmplxCoefProxy : public FortranProxy<BicubicCmplxCoefProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray2D<std::complex<double>> coef() const; // 2D_NOT_complex
  FArray1D<int> i_box() const; // 1D_NOT_integer
};

template <>
struct FortranTraits<TricubicCmplxCoefProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tricubic_cmplx_coef_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tricubic_cmplx_coef_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tricubic_cmplx_coef_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tricubic_cmplx_coef_struct";
  }
};

class TricubicCmplxCoefProxy : public FortranProxy<TricubicCmplxCoefProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray3D<std::complex<double>> coef() const; // 3D_NOT_complex
  FArray1D<int> i_box() const; // 1D_NOT_integer
};

template <>
struct FortranTraits<GridFieldPt1Proxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_grid_field_pt1_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_grid_field_pt1_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_grid_field_pt1_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "grid_field_pt1_struct";
  }
};

class GridFieldPt1Proxy : public FortranProxy<GridFieldPt1Proxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<std::complex<double>> E() const; // 1D_NOT_complex
  FArray1D<std::complex<double>> B() const; // 1D_NOT_complex
};

template <>
struct FortranTraits<GridFieldPtProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_grid_field_pt_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_grid_field_pt_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_grid_field_pt_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "grid_field_pt_struct";
  }
};

class GridFieldPtProxy : public FortranProxy<GridFieldPtProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string& value);
  int n_link() const; // 0D_NOT_integer
  void set_n_link(int value);
  GridFieldPt1ProxyArray3D pt() const; // 3D_ALLOC_type
};

template <>
struct FortranTraits<GridFieldProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_grid_field_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_grid_field_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_grid_field_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "grid_field_struct";
  }
};

class GridFieldProxy : public FortranProxy<GridFieldProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int geometry() const; // 0D_NOT_integer
  void set_geometry(int value);
  int harmonic() const; // 0D_NOT_integer
  void set_harmonic(int value);
  double phi0_fieldmap() const; // 0D_NOT_real
  void set_phi0_fieldmap(double value);
  double field_scale() const; // 0D_NOT_real
  void set_field_scale(double value);
  int field_type() const; // 0D_NOT_integer
  void set_field_type(int value);
  int master_parameter() const; // 0D_NOT_integer
  void set_master_parameter(int value);
  int ele_anchor_pt() const; // 0D_NOT_integer
  void set_ele_anchor_pt(int value);
  int interpolation_order() const; // 0D_NOT_integer
  void set_interpolation_order(int value);
  FArray1D<double> dr() const; // 1D_NOT_real
  FArray1D<double> r0() const; // 1D_NOT_real
  bool curved_ref_frame() const; // 0D_NOT_logical
  void set_curved_ref_frame(bool value);
  std::optional<GridFieldPtProxy> ptr() const; // 0D_PTR_type
  void set_ptr(const GridFieldPtProxy& src);
  BicubicCmplxCoefProxyArray3D bi_coef() const; // 3D_NOT_type
  TricubicCmplxCoefProxyArray3D tri_coef() const; // 3D_NOT_type
};

template <>
struct FortranTraits<FloorPositionProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_floor_position_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_floor_position_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_floor_position_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "floor_position_struct";
  }
};

class FloorPositionProxy : public FortranProxy<FloorPositionProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> r() const; // 1D_NOT_real
  FArray2D<double> w() const; // 2D_NOT_real
  double theta() const; // 0D_NOT_real
  void set_theta(double value);
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
  double psi() const; // 0D_NOT_real
  void set_psi(double value);
};

template <>
struct FortranTraits<HighEnergySpaceChargeProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_high_energy_space_charge_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_high_energy_space_charge_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_high_energy_space_charge_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "high_energy_space_charge_struct";
  }
};

class HighEnergySpaceChargeProxy
    : public FortranProxy<HighEnergySpaceChargeProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  CoordProxy closed_orb() const; // 0D_NOT_type
  void set_closed_orb(const CoordProxy& src);
  double kick_const() const; // 0D_NOT_real
  void set_kick_const(double value);
  double sig_x() const; // 0D_NOT_real
  void set_sig_x(double value);
  double sig_y() const; // 0D_NOT_real
  void set_sig_y(double value);
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
  double sin_phi() const; // 0D_NOT_real
  void set_sin_phi(double value);
  double cos_phi() const; // 0D_NOT_real
  void set_cos_phi(double value);
  double sig_z() const; // 0D_NOT_real
  void set_sig_z(double value);
};

template <>
struct FortranTraits<XyDispProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_xy_disp_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_xy_disp_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_xy_disp_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "xy_disp_struct";
  }
};

class XyDispProxy : public FortranProxy<XyDispProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double eta() const; // 0D_NOT_real
  void set_eta(double value);
  double etap() const; // 0D_NOT_real
  void set_etap(double value);
  double deta_ds() const; // 0D_NOT_real
  void set_deta_ds(double value);
  double sigma() const; // 0D_NOT_real
  void set_sigma(double value);
  double deta_dpz() const; // 0D_NOT_real
  void set_deta_dpz(double value);
  double detap_dpz() const; // 0D_NOT_real
  void set_detap_dpz(double value);
};

template <>
struct FortranTraits<TwissProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_twiss_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_twiss_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_twiss_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "twiss_struct";
  }
};

class TwissProxy : public FortranProxy<TwissProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double beta() const; // 0D_NOT_real
  void set_beta(double value);
  double alpha() const; // 0D_NOT_real
  void set_alpha(double value);
  double gamma() const; // 0D_NOT_real
  void set_gamma(double value);
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
  double eta() const; // 0D_NOT_real
  void set_eta(double value);
  double etap() const; // 0D_NOT_real
  void set_etap(double value);
  double deta_ds() const; // 0D_NOT_real
  void set_deta_ds(double value);
  double sigma() const; // 0D_NOT_real
  void set_sigma(double value);
  double sigma_p() const; // 0D_NOT_real
  void set_sigma_p(double value);
  double emit() const; // 0D_NOT_real
  void set_emit(double value);
  double norm_emit() const; // 0D_NOT_real
  void set_norm_emit(double value);
  double chrom() const; // 0D_NOT_real
  void set_chrom(double value);
  double dbeta_dpz() const; // 0D_NOT_real
  void set_dbeta_dpz(double value);
  double dalpha_dpz() const; // 0D_NOT_real
  void set_dalpha_dpz(double value);
  double deta_dpz() const; // 0D_NOT_real
  void set_deta_dpz(double value);
  double detap_dpz() const; // 0D_NOT_real
  void set_detap_dpz(double value);
};

template <>
struct FortranTraits<Mode3Proxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_mode3_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_mode3_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_mode3_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "mode3_struct";
  }
};

class Mode3Proxy : public FortranProxy<Mode3Proxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray2D<double> v() const; // 2D_NOT_real
  TwissProxy a() const; // 0D_NOT_type
  void set_a(const TwissProxy& src);
  TwissProxy b() const; // 0D_NOT_type
  void set_b(const TwissProxy& src);
  TwissProxy c() const; // 0D_NOT_type
  void set_c(const TwissProxy& src);
  TwissProxy x() const; // 0D_NOT_type
  void set_x(const TwissProxy& src);
  TwissProxy y() const; // 0D_NOT_type
  void set_y(const TwissProxy& src);
};

template <>
struct FortranTraits<BookkeepingStateProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bookkeeping_state_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bookkeeping_state_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bookkeeping_state_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bookkeeping_state_struct";
  }
};

class BookkeepingStateProxy : public FortranProxy<BookkeepingStateProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int attributes() const; // 0D_NOT_integer
  void set_attributes(int value);
  int control() const; // 0D_NOT_integer
  void set_control(int value);
  int floor_position() const; // 0D_NOT_integer
  void set_floor_position(int value);
  int s_position() const; // 0D_NOT_integer
  void set_s_position(int value);
  int ref_energy() const; // 0D_NOT_integer
  void set_ref_energy(int value);
  int mat6() const; // 0D_NOT_integer
  void set_mat6(int value);
  int rad_int() const; // 0D_NOT_integer
  void set_rad_int(int value);
  int ptc() const; // 0D_NOT_integer
  void set_ptc(int value);
  bool has_misalign() const; // 0D_NOT_logical
  void set_has_misalign(bool value);
};

template <>
struct FortranTraits<RadMapProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_rad_map_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_rad_map_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_rad_map_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "rad_map_struct";
  }
};

class RadMapProxy : public FortranProxy<RadMapProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> ref_orb() const; // 1D_NOT_real
  FArray2D<double> damp_dmat() const; // 2D_NOT_real
  FArray1D<double> xfer_damp_vec() const; // 1D_NOT_real
  FArray2D<double> xfer_damp_mat() const; // 2D_NOT_real
  FArray2D<double> stoc_mat() const; // 2D_NOT_real
};

template <>
struct FortranTraits<RadMapEleProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_rad_map_ele_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_rad_map_ele_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_rad_map_ele_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "rad_map_ele_struct";
  }
};

class RadMapEleProxy : public FortranProxy<RadMapEleProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  RadMapProxy rm0() const; // 0D_NOT_type
  void set_rm0(const RadMapProxy& src);
  RadMapProxy rm1() const; // 0D_NOT_type
  void set_rm1(const RadMapProxy& src);
  bool stale() const; // 0D_NOT_logical
  void set_stale(bool value);
};

template <>
struct FortranTraits<GenGrad1Proxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_gen_grad1_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_gen_grad1_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_gen_grad1_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "gen_grad1_struct";
  }
};

class GenGrad1Proxy : public FortranProxy<GenGrad1Proxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int m() const; // 0D_NOT_integer
  void set_m(int value);
  int sincos() const; // 0D_NOT_integer
  void set_sincos(int value);
  int n_deriv_max() const; // 0D_NOT_integer
  void set_n_deriv_max(int value);
  FArray2D<double> deriv() const; // 2D_ALLOC_real
};

template <>
struct FortranTraits<GenGradMapProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_gen_grad_map_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_gen_grad_map_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_gen_grad_map_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "gen_grad_map_struct";
  }
};

class GenGradMapProxy : public FortranProxy<GenGradMapProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string& value);
  GenGrad1ProxyArray1D gg() const; // 1D_ALLOC_type
  int ele_anchor_pt() const; // 0D_NOT_integer
  void set_ele_anchor_pt(int value);
  int field_type() const; // 0D_NOT_integer
  void set_field_type(int value);
  int iz0() const; // 0D_NOT_integer
  void set_iz0(int value);
  int iz1() const; // 0D_NOT_integer
  void set_iz1(int value);
  double dz() const; // 0D_NOT_real
  void set_dz(double value);
  FArray1D<double> r0() const; // 1D_NOT_real
  double field_scale() const; // 0D_NOT_real
  void set_field_scale(double value);
  int master_parameter() const; // 0D_NOT_integer
  void set_master_parameter(int value);
  bool curved_ref_frame() const; // 0D_NOT_logical
  void set_curved_ref_frame(bool value);
};

template <>
struct FortranTraits<SurfaceSegmentedPtProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_surface_segmented_pt_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_surface_segmented_pt_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_surface_segmented_pt_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "surface_segmented_pt_struct";
  }
};

class SurfaceSegmentedPtProxy : public FortranProxy<SurfaceSegmentedPtProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double x0() const; // 0D_NOT_real
  void set_x0(double value);
  double y0() const; // 0D_NOT_real
  void set_y0(double value);
  double z0() const; // 0D_NOT_real
  void set_z0(double value);
  double dz_dx() const; // 0D_NOT_real
  void set_dz_dx(double value);
  double dz_dy() const; // 0D_NOT_real
  void set_dz_dy(double value);
};

template <>
struct FortranTraits<SurfaceSegmentedProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_surface_segmented_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_surface_segmented_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_surface_segmented_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "surface_segmented_struct";
  }
};

class SurfaceSegmentedProxy : public FortranProxy<SurfaceSegmentedProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool active() const; // 0D_NOT_logical
  void set_active(bool value);
  FArray1D<double> dr() const; // 1D_NOT_real
  FArray1D<double> r0() const; // 1D_NOT_real
  SurfaceSegmentedPtProxyArray2D pt() const; // 2D_ALLOC_type
};

template <>
struct FortranTraits<SurfaceHMisalignPtProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_surface_h_misalign_pt_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_surface_h_misalign_pt_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_surface_h_misalign_pt_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "surface_h_misalign_pt_struct";
  }
};

class SurfaceHMisalignPtProxy : public FortranProxy<SurfaceHMisalignPtProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double x0() const; // 0D_NOT_real
  void set_x0(double value);
  double y0() const; // 0D_NOT_real
  void set_y0(double value);
  double rot_y() const; // 0D_NOT_real
  void set_rot_y(double value);
  double rot_t() const; // 0D_NOT_real
  void set_rot_t(double value);
  double rot_y_rms() const; // 0D_NOT_real
  void set_rot_y_rms(double value);
  double rot_t_rms() const; // 0D_NOT_real
  void set_rot_t_rms(double value);
};

template <>
struct FortranTraits<SurfaceHMisalignProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_surface_h_misalign_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_surface_h_misalign_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_surface_h_misalign_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "surface_h_misalign_struct";
  }
};

class SurfaceHMisalignProxy : public FortranProxy<SurfaceHMisalignProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool active() const; // 0D_NOT_logical
  void set_active(bool value);
  FArray1D<double> dr() const; // 1D_NOT_real
  FArray1D<double> r0() const; // 1D_NOT_real
  SurfaceHMisalignPtProxyArray2D pt() const; // 2D_ALLOC_type
};

template <>
struct FortranTraits<SurfaceDisplacementPtProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_surface_displacement_pt_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_surface_displacement_pt_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_surface_displacement_pt_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "surface_displacement_pt_struct";
  }
};

class SurfaceDisplacementPtProxy
    : public FortranProxy<SurfaceDisplacementPtProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double x0() const; // 0D_NOT_real
  void set_x0(double value);
  double y0() const; // 0D_NOT_real
  void set_y0(double value);
  double z0() const; // 0D_NOT_real
  void set_z0(double value);
  double dz_dx() const; // 0D_NOT_real
  void set_dz_dx(double value);
  double dz_dy() const; // 0D_NOT_real
  void set_dz_dy(double value);
  double d2z_dxdy() const; // 0D_NOT_real
  void set_d2z_dxdy(double value);
};

template <>
struct FortranTraits<SurfaceDisplacementProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_surface_displacement_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_surface_displacement_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_surface_displacement_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "surface_displacement_struct";
  }
};

class SurfaceDisplacementProxy : public FortranProxy<SurfaceDisplacementProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool active() const; // 0D_NOT_logical
  void set_active(bool value);
  FArray1D<double> dr() const; // 1D_NOT_real
  FArray1D<double> r0() const; // 1D_NOT_real
  SurfaceDisplacementPtProxyArray2D pt() const; // 2D_ALLOC_type
};

template <>
struct FortranTraits<TargetPointProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_target_point_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_target_point_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_target_point_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "target_point_struct";
  }
};

class TargetPointProxy : public FortranProxy<TargetPointProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> r() const; // 1D_NOT_real
};

template <>
struct FortranTraits<SurfaceCurvatureProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_surface_curvature_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_surface_curvature_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_surface_curvature_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "surface_curvature_struct";
  }
};

class SurfaceCurvatureProxy : public FortranProxy<SurfaceCurvatureProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray2D<double> xy() const; // 2D_NOT_real
  double spherical() const; // 0D_NOT_real
  void set_spherical(double value);
  FArray1D<double> elliptical() const; // 1D_NOT_real
  bool has_curvature() const; // 0D_NOT_logical
  void set_has_curvature(bool value);
};

template <>
struct FortranTraits<PhotonTargetProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_photon_target_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_photon_target_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_photon_target_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "photon_target_struct";
  }
};

class PhotonTargetProxy : public FortranProxy<PhotonTargetProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int type() const; // 0D_NOT_integer
  void set_type(int value);
  int n_corner() const; // 0D_NOT_integer
  void set_n_corner(int value);
  LatEleLocProxy ele_loc() const; // 0D_NOT_type
  void set_ele_loc(const LatEleLocProxy& src);
  TargetPointProxyArray1D corner() const; // 1D_NOT_type
  TargetPointProxy center() const; // 0D_NOT_type
  void set_center(const TargetPointProxy& src);
};

template <>
struct FortranTraits<PhotonMaterialProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_photon_material_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_photon_material_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_photon_material_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "photon_material_struct";
  }
};

class PhotonMaterialProxy : public FortranProxy<PhotonMaterialProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::complex<double> f0_m1() const; // 0D_NOT_complex
  void set_f0_m1(std::complex<double> value);
  std::complex<double> f0_m2() const; // 0D_NOT_complex
  void set_f0_m2(std::complex<double> value);
  std::complex<double> f_0() const; // 0D_NOT_complex
  void set_f_0(std::complex<double> value);
  std::complex<double> f_h() const; // 0D_NOT_complex
  void set_f_h(std::complex<double> value);
  std::complex<double> f_hbar() const; // 0D_NOT_complex
  void set_f_hbar(std::complex<double> value);
  std::complex<double> f_hkl() const; // 0D_NOT_complex
  void set_f_hkl(std::complex<double> value);
  FArray1D<double> h_norm() const; // 1D_NOT_real
  FArray1D<double> l_ref() const; // 1D_NOT_real
};

template <>
struct FortranTraits<PixelPtProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_pixel_pt_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_pixel_pt_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_pixel_pt_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "pixel_pt_struct";
  }
};

class PixelPtProxy : public FortranProxy<PixelPtProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int64_t n_photon() const; // 0D_NOT_integer8
  void set_n_photon(int64_t value);
  std::complex<double> E_x() const; // 0D_NOT_complex
  void set_E_x(std::complex<double> value);
  std::complex<double> E_y() const; // 0D_NOT_complex
  void set_E_y(std::complex<double> value);
  double intensity_x() const; // 0D_NOT_real
  void set_intensity_x(double value);
  double intensity_y() const; // 0D_NOT_real
  void set_intensity_y(double value);
  double intensity() const; // 0D_NOT_real
  void set_intensity(double value);
  FArray1D<double> orbit() const; // 1D_NOT_real
  FArray1D<double> orbit_rms() const; // 1D_NOT_real
  FArray1D<double> init_orbit() const; // 1D_NOT_real
  FArray1D<double> init_orbit_rms() const; // 1D_NOT_real
};

template <>
struct FortranTraits<PixelDetecProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_pixel_detec_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_pixel_detec_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_pixel_detec_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "pixel_detec_struct";
  }
};

class PixelDetecProxy : public FortranProxy<PixelDetecProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> dr() const; // 1D_NOT_real
  FArray1D<double> r0() const; // 1D_NOT_real
  int64_t n_track_tot() const; // 0D_NOT_integer8
  void set_n_track_tot(int64_t value);
  int64_t n_hit_detec() const; // 0D_NOT_integer8
  void set_n_hit_detec(int64_t value);
  int64_t n_hit_pixel() const; // 0D_NOT_integer8
  void set_n_hit_pixel(int64_t value);
  PixelPtProxyArray2D pt() const; // 2D_ALLOC_type
};

template <>
struct FortranTraits<PhotonElementProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_photon_element_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_photon_element_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_photon_element_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "photon_element_struct";
  }
};

class PhotonElementProxy : public FortranProxy<PhotonElementProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  SurfaceCurvatureProxy curvature() const; // 0D_NOT_type
  void set_curvature(const SurfaceCurvatureProxy& src);
  PhotonTargetProxy target() const; // 0D_NOT_type
  void set_target(const PhotonTargetProxy& src);
  PhotonMaterialProxy material() const; // 0D_NOT_type
  void set_material(const PhotonMaterialProxy& src);
  SurfaceSegmentedProxy segmented() const; // 0D_NOT_type
  void set_segmented(const SurfaceSegmentedProxy& src);
  SurfaceHMisalignProxy h_misalign() const; // 0D_NOT_type
  void set_h_misalign(const SurfaceHMisalignProxy& src);
  SurfaceDisplacementProxy displacement() const; // 0D_NOT_type
  void set_displacement(const SurfaceDisplacementProxy& src);
  PixelDetecProxy pixel() const; // 0D_NOT_type
  void set_pixel(const PixelDetecProxy& src);
  int reflectivity_table_type() const; // 0D_NOT_integer
  void set_reflectivity_table_type(int value);
  PhotonReflectTableProxy reflectivity_table_sigma() const; // 0D_NOT_type
  void set_reflectivity_table_sigma(const PhotonReflectTableProxy& src);
  PhotonReflectTableProxy reflectivity_table_pi() const; // 0D_NOT_type
  void set_reflectivity_table_pi(const PhotonReflectTableProxy& src);
  SplineProxyArray1D init_energy_prob() const; // 1D_ALLOC_type
  FArray1D<double> integrated_init_energy_prob() const; // 1D_ALLOC_real
};

template <>
struct FortranTraits<Wall3dVertexProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_wall3d_vertex_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_wall3d_vertex_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_wall3d_vertex_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "wall3d_vertex_struct";
  }
};

class Wall3dVertexProxy : public FortranProxy<Wall3dVertexProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double x() const; // 0D_NOT_real
  void set_x(double value);
  double y() const; // 0D_NOT_real
  void set_y(double value);
  double radius_x() const; // 0D_NOT_real
  void set_radius_x(double value);
  double radius_y() const; // 0D_NOT_real
  void set_radius_y(double value);
  double tilt() const; // 0D_NOT_real
  void set_tilt(double value);
  double angle() const; // 0D_NOT_real
  void set_angle(double value);
  double x0() const; // 0D_NOT_real
  void set_x0(double value);
  double y0() const; // 0D_NOT_real
  void set_y0(double value);
  int type() const; // 0D_NOT_integer
  void set_type(int value);
};

template <>
struct FortranTraits<Wall3dSectionProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_wall3d_section_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_wall3d_section_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_wall3d_section_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "wall3d_section_struct";
  }
};

class Wall3dSectionProxy : public FortranProxy<Wall3dSectionProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  std::string material() const; // 0D_NOT_character
  void set_material(const std::string& value);
  Wall3dVertexProxyArray1D v() const; // 1D_ALLOC_type
  std::optional<PhotonReflectSurfaceProxy> surface() const; // 0D_PTR_type
  void set_surface(const PhotonReflectSurfaceProxy& src);
  int type() const; // 0D_NOT_integer
  void set_type(int value);
  int n_vertex_input() const; // 0D_NOT_integer
  void set_n_vertex_input(int value);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int vertices_state() const; // 0D_NOT_integer
  void set_vertices_state(int value);
  bool patch_in_region() const; // 0D_NOT_logical
  void set_patch_in_region(bool value);
  double thickness() const; // 0D_NOT_real
  void set_thickness(double value);
  double s() const; // 0D_NOT_real
  void set_s(double value);
  FArray1D<double> r0() const; // 1D_NOT_real
  double dx0_ds() const; // 0D_NOT_real
  void set_dx0_ds(double value);
  double dy0_ds() const; // 0D_NOT_real
  void set_dy0_ds(double value);
  FArray1D<double> x0_coef() const; // 1D_NOT_real
  FArray1D<double> y0_coef() const; // 1D_NOT_real
  double dr_ds() const; // 0D_NOT_real
  void set_dr_ds(double value);
  FArray1D<double> p1_coef() const; // 1D_NOT_real
  FArray1D<double> p2_coef() const; // 1D_NOT_real
};

template <>
struct FortranTraits<Wall3dProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_wall3d_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_wall3d_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_wall3d_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "wall3d_struct";
  }
};

class Wall3dProxy : public FortranProxy<Wall3dProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  int type() const; // 0D_NOT_integer
  void set_type(int value);
  int ix_wall3d() const; // 0D_NOT_integer
  void set_ix_wall3d(int value);
  int n_link() const; // 0D_NOT_integer
  void set_n_link(int value);
  double thickness() const; // 0D_NOT_real
  void set_thickness(double value);
  std::string clear_material() const; // 0D_NOT_character
  void set_clear_material(const std::string& value);
  std::string opaque_material() const; // 0D_NOT_character
  void set_opaque_material(const std::string& value);
  bool superimpose() const; // 0D_NOT_logical
  void set_superimpose(bool value);
  int ele_anchor_pt() const; // 0D_NOT_integer
  void set_ele_anchor_pt(int value);
  Wall3dSectionProxyArray1D section() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<RamperLordProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_ramper_lord_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_ramper_lord_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_ramper_lord_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "ramper_lord_struct";
  }
};

class RamperLordProxy : public FortranProxy<RamperLordProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_con() const; // 0D_NOT_integer
  void set_ix_con(int value);
  double* attrib_ptr() const; // 0D_PTR_real
  void set_attrib_ptr(double value);
};

template <>
struct FortranTraits<ControlProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_control_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_control_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_control_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "control_struct";
  }
};

class ControlProxy : public FortranProxy<ControlProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double value() const; // 0D_NOT_real
  void set_value(double value);
  FArray1D<double> y_knot() const; // 1D_ALLOC_real
  ExpressionAtomProxyArray1D stack() const; // 1D_ALLOC_type
  LatEleLocProxy slave() const; // 0D_NOT_type
  void set_slave(const LatEleLocProxy& src);
  LatEleLocProxy lord() const; // 0D_NOT_type
  void set_lord(const LatEleLocProxy& src);
  std::string slave_name() const; // 0D_NOT_character
  void set_slave_name(const std::string& value);
  std::string attribute() const; // 0D_NOT_character
  void set_attribute(const std::string& value);
  int ix_attrib() const; // 0D_NOT_integer
  void set_ix_attrib(int value);
};

template <>
struct FortranTraits<ControlVar1Proxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_control_var1_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_control_var1_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_control_var1_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "control_var1_struct";
  }
};

class ControlVar1Proxy : public FortranProxy<ControlVar1Proxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  double value() const; // 0D_NOT_real
  void set_value(double value);
  double old_value() const; // 0D_NOT_real
  void set_old_value(double value);
};

template <>
struct FortranTraits<ControlRamp1Proxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_control_ramp1_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_control_ramp1_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_control_ramp1_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "control_ramp1_struct";
  }
};

class ControlRamp1Proxy : public FortranProxy<ControlRamp1Proxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> y_knot() const; // 1D_ALLOC_real
  ExpressionAtomProxyArray1D stack() const; // 1D_ALLOC_type
  std::string attribute() const; // 0D_NOT_character
  void set_attribute(const std::string& value);
  std::string slave_name() const; // 0D_NOT_character
  void set_slave_name(const std::string& value);
  bool is_controller() const; // 0D_NOT_logical
  void set_is_controller(bool value);
};

template <>
struct FortranTraits<ControllerProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_controller_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_controller_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_controller_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "controller_struct";
  }
};

class ControllerProxy : public FortranProxy<ControllerProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  ControlVar1ProxyArray1D var() const; // 1D_ALLOC_type
  ControlRamp1ProxyArray1D ramp() const; // 1D_ALLOC_type
  RamperLordProxyArray1D ramper_lord() const; // 1D_ALLOC_type
  FArray1D<double> x_knot() const; // 1D_ALLOC_real
};

template <>
struct FortranTraits<EllipseBeamInitProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_ellipse_beam_init_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_ellipse_beam_init_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_ellipse_beam_init_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "ellipse_beam_init_struct";
  }
};

class EllipseBeamInitProxy : public FortranProxy<EllipseBeamInitProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int part_per_ellipse() const; // 0D_NOT_integer
  void set_part_per_ellipse(int value);
  int n_ellipse() const; // 0D_NOT_integer
  void set_n_ellipse(int value);
  double sigma_cutoff() const; // 0D_NOT_real
  void set_sigma_cutoff(double value);
};

template <>
struct FortranTraits<KvBeamInitProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_kv_beam_init_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_kv_beam_init_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_kv_beam_init_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "kv_beam_init_struct";
  }
};

class KvBeamInitProxy : public FortranProxy<KvBeamInitProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<int> part_per_phi() const; // 1D_NOT_integer
  int n_I2() const; // 0D_NOT_integer
  void set_n_I2(int value);
  double A() const; // 0D_NOT_real
  void set_A(double value);
};

template <>
struct FortranTraits<GridBeamInitProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_grid_beam_init_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_grid_beam_init_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_grid_beam_init_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "grid_beam_init_struct";
  }
};

class GridBeamInitProxy : public FortranProxy<GridBeamInitProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int n_x() const; // 0D_NOT_integer
  void set_n_x(int value);
  int n_px() const; // 0D_NOT_integer
  void set_n_px(int value);
  double x_min() const; // 0D_NOT_real
  void set_x_min(double value);
  double x_max() const; // 0D_NOT_real
  void set_x_max(double value);
  double px_min() const; // 0D_NOT_real
  void set_px_min(double value);
  double px_max() const; // 0D_NOT_real
  void set_px_max(double value);
};

template <>
struct FortranTraits<BeamInitProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_beam_init_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_beam_init_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_beam_init_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "beam_init_struct";
  }
};

class BeamInitProxy : public FortranProxy<BeamInitProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string position_file() const; // 0D_NOT_character
  void set_position_file(const std::string& value);
  FCharArray1D distribution_type() const; // 1D_NOT_character
  FArray1D<double> spin() const; // 1D_NOT_real
  EllipseBeamInitProxyArray1D ellipse() const; // 1D_NOT_type
  KvBeamInitProxy KV() const; // 0D_NOT_type
  void set_KV(const KvBeamInitProxy& src);
  GridBeamInitProxyArray1D grid() const; // 1D_NOT_type
  FArray1D<double> center_jitter() const; // 1D_NOT_real
  FArray1D<double> emit_jitter() const; // 1D_NOT_real
  double sig_z_jitter() const; // 0D_NOT_real
  void set_sig_z_jitter(double value);
  double sig_pz_jitter() const; // 0D_NOT_real
  void set_sig_pz_jitter(double value);
  int n_particle() const; // 0D_NOT_integer
  void set_n_particle(int value);
  bool renorm_center() const; // 0D_NOT_logical
  void set_renorm_center(bool value);
  bool renorm_sigma() const; // 0D_NOT_logical
  void set_renorm_sigma(bool value);
  std::string random_engine() const; // 0D_NOT_character
  void set_random_engine(const std::string& value);
  std::string random_gauss_converter() const; // 0D_NOT_character
  void set_random_gauss_converter(const std::string& value);
  double random_sigma_cutoff() const; // 0D_NOT_real
  void set_random_sigma_cutoff(double value);
  double a_norm_emit() const; // 0D_NOT_real
  void set_a_norm_emit(double value);
  double b_norm_emit() const; // 0D_NOT_real
  void set_b_norm_emit(double value);
  double a_emit() const; // 0D_NOT_real
  void set_a_emit(double value);
  double b_emit() const; // 0D_NOT_real
  void set_b_emit(double value);
  double dPz_dz() const; // 0D_NOT_real
  void set_dPz_dz(double value);
  FArray1D<double> center() const; // 1D_NOT_real
  double t_offset() const; // 0D_NOT_real
  void set_t_offset(double value);
  double dt_bunch() const; // 0D_NOT_real
  void set_dt_bunch(double value);
  double sig_z() const; // 0D_NOT_real
  void set_sig_z(double value);
  double sig_pz() const; // 0D_NOT_real
  void set_sig_pz(double value);
  double bunch_charge() const; // 0D_NOT_real
  void set_bunch_charge(double value);
  int n_bunch() const; // 0D_NOT_integer
  void set_n_bunch(int value);
  int ix_turn() const; // 0D_NOT_integer
  void set_ix_turn(int value);
  std::string species() const; // 0D_NOT_character
  void set_species(const std::string& value);
  bool full_6D_coupling_calc() const; // 0D_NOT_logical
  void set_full_6D_coupling_calc(bool value);
  bool use_particle_start() const; // 0D_NOT_logical
  void set_use_particle_start(bool value);
  bool use_t_coords() const; // 0D_NOT_logical
  void set_use_t_coords(bool value);
  bool use_z_as_t() const; // 0D_NOT_logical
  void set_use_z_as_t(bool value);
  std::string file_name() const; // 0D_NOT_character
  void set_file_name(const std::string& value);
};

template <>
struct FortranTraits<LatParamProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_lat_param_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_lat_param_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_lat_param_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "lat_param_struct";
  }
};

class LatParamProxy : public FortranProxy<LatParamProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double n_part() const; // 0D_NOT_real
  void set_n_part(double value);
  double total_length() const; // 0D_NOT_real
  void set_total_length(double value);
  double unstable_factor() const; // 0D_NOT_real
  void set_unstable_factor(double value);
  FArray2D<double> t1_with_RF() const; // 2D_NOT_real
  FArray2D<double> t1_no_RF() const; // 2D_NOT_real
  double spin_tune() const; // 0D_NOT_real
  void set_spin_tune(double value);
  int particle() const; // 0D_NOT_integer
  void set_particle(int value);
  int default_tracking_species() const; // 0D_NOT_integer
  void set_default_tracking_species(int value);
  int geometry() const; // 0D_NOT_integer
  void set_geometry(int value);
  int ixx() const; // 0D_NOT_integer
  void set_ixx(int value);
  bool stable() const; // 0D_NOT_logical
  void set_stable(bool value);
  bool live_branch() const; // 0D_NOT_logical
  void set_live_branch(bool value);
  double g1_integral() const; // 0D_NOT_real
  void set_g1_integral(double value);
  double g2_integral() const; // 0D_NOT_real
  void set_g2_integral(double value);
  double g3_integral() const; // 0D_NOT_real
  void set_g3_integral(double value);
  BookkeepingStateProxy bookkeeping_state() const; // 0D_NOT_type
  void set_bookkeeping_state(const BookkeepingStateProxy& src);
  BeamInitProxy beam_init() const; // 0D_NOT_type
  void set_beam_init(const BeamInitProxy& src);
};

template <>
struct FortranTraits<ModeInfoProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_mode_info_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_mode_info_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_mode_info_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "mode_info_struct";
  }
};

class ModeInfoProxy : public FortranProxy<ModeInfoProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool stable() const; // 0D_NOT_logical
  void set_stable(bool value);
  double tune() const; // 0D_NOT_real
  void set_tune(double value);
  double emit() const; // 0D_NOT_real
  void set_emit(double value);
  double chrom() const; // 0D_NOT_real
  void set_chrom(double value);
  double sigma() const; // 0D_NOT_real
  void set_sigma(double value);
  double sigmap() const; // 0D_NOT_real
  void set_sigmap(double value);
};

template <>
struct FortranTraits<PreTrackerProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_pre_tracker_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_pre_tracker_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_pre_tracker_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "pre_tracker_struct";
  }
};

class PreTrackerProxy : public FortranProxy<PreTrackerProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int who() const; // 0D_NOT_integer
  void set_who(int value);
  int ix_ele_start() const; // 0D_NOT_integer
  void set_ix_ele_start(int value);
  int ix_ele_end() const; // 0D_NOT_integer
  void set_ix_ele_end(int value);
  std::string input_file() const; // 0D_NOT_character
  void set_input_file(const std::string& value);
};

template <>
struct FortranTraits<AnormalModeProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_anormal_mode_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_anormal_mode_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_anormal_mode_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "anormal_mode_struct";
  }
};

class AnormalModeProxy : public FortranProxy<AnormalModeProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double emittance() const; // 0D_NOT_real
  void set_emittance(double value);
  double emittance_no_vert() const; // 0D_NOT_real
  void set_emittance_no_vert(double value);
  FArray1D<double> synch_int() const; // 1D_NOT_real
  double j_damp() const; // 0D_NOT_real
  void set_j_damp(double value);
  double alpha_damp() const; // 0D_NOT_real
  void set_alpha_damp(double value);
  double chrom() const; // 0D_NOT_real
  void set_chrom(double value);
  double tune() const; // 0D_NOT_real
  void set_tune(double value);
};

template <>
struct FortranTraits<LinacNormalModeProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_linac_normal_mode_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_linac_normal_mode_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_linac_normal_mode_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "linac_normal_mode_struct";
  }
};

class LinacNormalModeProxy : public FortranProxy<LinacNormalModeProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double i2_E4() const; // 0D_NOT_real
  void set_i2_E4(double value);
  double i3_E7() const; // 0D_NOT_real
  void set_i3_E7(double value);
  double i5a_E6() const; // 0D_NOT_real
  void set_i5a_E6(double value);
  double i5b_E6() const; // 0D_NOT_real
  void set_i5b_E6(double value);
  double sig_E1() const; // 0D_NOT_real
  void set_sig_E1(double value);
  double a_emittance_end() const; // 0D_NOT_real
  void set_a_emittance_end(double value);
  double b_emittance_end() const; // 0D_NOT_real
  void set_b_emittance_end(double value);
};

template <>
struct FortranTraits<NormalModesProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_normal_modes_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_normal_modes_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_normal_modes_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "normal_modes_struct";
  }
};

class NormalModesProxy : public FortranProxy<NormalModesProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> synch_int() const; // 1D_NOT_real
  double sigE_E() const; // 0D_NOT_real
  void set_sigE_E(double value);
  double sig_z() const; // 0D_NOT_real
  void set_sig_z(double value);
  double e_loss() const; // 0D_NOT_real
  void set_e_loss(double value);
  double rf_voltage() const; // 0D_NOT_real
  void set_rf_voltage(double value);
  double pz_aperture() const; // 0D_NOT_real
  void set_pz_aperture(double value);
  double pz_average() const; // 0D_NOT_real
  void set_pz_average(double value);
  double momentum_compaction() const; // 0D_NOT_real
  void set_momentum_compaction(double value);
  double dpz_damp() const; // 0D_NOT_real
  void set_dpz_damp(double value);
  AnormalModeProxy a() const; // 0D_NOT_type
  void set_a(const AnormalModeProxy& src);
  AnormalModeProxy b() const; // 0D_NOT_type
  void set_b(const AnormalModeProxy& src);
  AnormalModeProxy z() const; // 0D_NOT_type
  void set_z(const AnormalModeProxy& src);
  LinacNormalModeProxy lin() const; // 0D_NOT_type
  void set_lin(const LinacNormalModeProxy& src);
};

template <>
struct FortranTraits<EmFieldProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_em_field_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_em_field_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_em_field_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "em_field_struct";
  }
};

class EmFieldProxy : public FortranProxy<EmFieldProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> E() const; // 1D_NOT_real
  FArray1D<double> B() const; // 1D_NOT_real
  FArray2D<double> dE() const; // 2D_NOT_real
  FArray2D<double> dB() const; // 2D_NOT_real
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
  double phi_B() const; // 0D_NOT_real
  void set_phi_B(double value);
  FArray1D<double> A() const; // 1D_NOT_real
};

template <>
struct FortranTraits<StrongBeamProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_strong_beam_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_strong_beam_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_strong_beam_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "strong_beam_struct";
  }
};

class StrongBeamProxy : public FortranProxy<StrongBeamProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int ix_slice() const; // 0D_NOT_integer
  void set_ix_slice(int value);
  double x_center() const; // 0D_NOT_real
  void set_x_center(double value);
  double y_center() const; // 0D_NOT_real
  void set_y_center(double value);
  double x_sigma() const; // 0D_NOT_real
  void set_x_sigma(double value);
  double y_sigma() const; // 0D_NOT_real
  void set_y_sigma(double value);
  double dx() const; // 0D_NOT_real
  void set_dx(double value);
  double dy() const; // 0D_NOT_real
  void set_dy(double value);
};

template <>
struct FortranTraits<TrackPointProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_track_point_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_track_point_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_track_point_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "track_point_struct";
  }
};

class TrackPointProxy : public FortranProxy<TrackPointProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double s_lab() const; // 0D_NOT_real
  void set_s_lab(double value);
  double s_body() const; // 0D_NOT_real
  void set_s_body(double value);
  CoordProxy orb() const; // 0D_NOT_type
  void set_orb(const CoordProxy& src);
  EmFieldProxy field() const; // 0D_NOT_type
  void set_field(const EmFieldProxy& src);
  StrongBeamProxy strong_beam() const; // 0D_NOT_type
  void set_strong_beam(const StrongBeamProxy& src);
  FArray1D<double> vec0() const; // 1D_NOT_real
  FArray2D<double> mat6() const; // 2D_NOT_real
};

template <>
struct FortranTraits<TrackProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_track_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_track_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_track_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "track_struct";
  }
};

class TrackProxy : public FortranProxy<TrackProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TrackPointProxyArray1D pt() const; // 1D_ALLOC_type
  double ds_save() const; // 0D_NOT_real
  void set_ds_save(double value);
  int n_pt() const; // 0D_NOT_integer
  void set_n_pt(int value);
  int n_bad() const; // 0D_NOT_integer
  void set_n_bad(int value);
  int n_ok() const; // 0D_NOT_integer
  void set_n_ok(int value);
};

template <>
struct FortranTraits<SpaceChargeCommonProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_space_charge_common_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_space_charge_common_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_space_charge_common_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "space_charge_common_struct";
  }
};

class SpaceChargeCommonProxy : public FortranProxy<SpaceChargeCommonProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double ds_track_step() const; // 0D_NOT_real
  void set_ds_track_step(double value);
  double dt_track_step() const; // 0D_NOT_real
  void set_dt_track_step(double value);
  double cathode_strength_cutoff() const; // 0D_NOT_real
  void set_cathode_strength_cutoff(double value);
  double rel_tol_tracking() const; // 0D_NOT_real
  void set_rel_tol_tracking(double value);
  double abs_tol_tracking() const; // 0D_NOT_real
  void set_abs_tol_tracking(double value);
  double beam_chamber_height() const; // 0D_NOT_real
  void set_beam_chamber_height(double value);
  double lsc_sigma_cutoff() const; // 0D_NOT_real
  void set_lsc_sigma_cutoff(double value);
  double particle_sigma_cutoff() const; // 0D_NOT_real
  void set_particle_sigma_cutoff(double value);
  FArray1D<int> space_charge_mesh_size() const; // 1D_NOT_integer
  FArray1D<int> csr3d_mesh_size() const; // 1D_NOT_integer
  int n_bin() const; // 0D_NOT_integer
  void set_n_bin(int value);
  int particle_bin_span() const; // 0D_NOT_integer
  void set_particle_bin_span(int value);
  int n_shield_images() const; // 0D_NOT_integer
  void set_n_shield_images(int value);
  int sc_min_in_bin() const; // 0D_NOT_integer
  void set_sc_min_in_bin(int value);
  bool lsc_kick_transverse_dependence() const; // 0D_NOT_logical
  void set_lsc_kick_transverse_dependence(bool value);
  bool debug() const; // 0D_NOT_logical
  void set_debug(bool value);
  std::string diagnostic_output_file() const; // 0D_NOT_character
  void set_diagnostic_output_file(const std::string& value);
};

template <>
struct FortranTraits<BmadCommonProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bmad_common_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bmad_common_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bmad_common_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bmad_common_struct";
  }
};

class BmadCommonProxy : public FortranProxy<BmadCommonProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double max_aperture_limit() const; // 0D_NOT_real
  void set_max_aperture_limit(double value);
  FArray1D<double> d_orb() const; // 1D_NOT_real
  double default_ds_step() const; // 0D_NOT_real
  void set_default_ds_step(double value);
  double significant_length() const; // 0D_NOT_real
  void set_significant_length(double value);
  double rel_tol_tracking() const; // 0D_NOT_real
  void set_rel_tol_tracking(double value);
  double abs_tol_tracking() const; // 0D_NOT_real
  void set_abs_tol_tracking(double value);
  double rel_tol_adaptive_tracking() const; // 0D_NOT_real
  void set_rel_tol_adaptive_tracking(double value);
  double abs_tol_adaptive_tracking() const; // 0D_NOT_real
  void set_abs_tol_adaptive_tracking(double value);
  double init_ds_adaptive_tracking() const; // 0D_NOT_real
  void set_init_ds_adaptive_tracking(double value);
  double min_ds_adaptive_tracking() const; // 0D_NOT_real
  void set_min_ds_adaptive_tracking(double value);
  double fatal_ds_adaptive_tracking() const; // 0D_NOT_real
  void set_fatal_ds_adaptive_tracking(double value);
  double autoscale_amp_abs_tol() const; // 0D_NOT_real
  void set_autoscale_amp_abs_tol(double value);
  double autoscale_amp_rel_tol() const; // 0D_NOT_real
  void set_autoscale_amp_rel_tol(double value);
  double autoscale_phase_tol() const; // 0D_NOT_real
  void set_autoscale_phase_tol(double value);
  double electric_dipole_moment() const; // 0D_NOT_real
  void set_electric_dipole_moment(double value);
  double synch_rad_scale() const; // 0D_NOT_real
  void set_synch_rad_scale(double value);
  double sad_eps_scale() const; // 0D_NOT_real
  void set_sad_eps_scale(double value);
  double sad_amp_max() const; // 0D_NOT_real
  void set_sad_amp_max(double value);
  int sad_n_div_max() const; // 0D_NOT_integer
  void set_sad_n_div_max(int value);
  int taylor_order() const; // 0D_NOT_integer
  void set_taylor_order(int value);
  int runge_kutta_order() const; // 0D_NOT_integer
  void set_runge_kutta_order(int value);
  int default_integ_order() const; // 0D_NOT_integer
  void set_default_integ_order(int value);
  int max_num_runge_kutta_step() const; // 0D_NOT_integer
  void set_max_num_runge_kutta_step(int value);
  bool rf_phase_below_transition_ref() const; // 0D_NOT_logical
  void set_rf_phase_below_transition_ref(bool value);
  bool sr_wakes_on() const; // 0D_NOT_logical
  void set_sr_wakes_on(bool value);
  bool lr_wakes_on() const; // 0D_NOT_logical
  void set_lr_wakes_on(bool value);
  bool auto_bookkeeper() const; // 0D_NOT_logical
  void set_auto_bookkeeper(bool value);
  bool high_energy_space_charge_on() const; // 0D_NOT_logical
  void set_high_energy_space_charge_on(bool value);
  bool csr_and_space_charge_on() const; // 0D_NOT_logical
  void set_csr_and_space_charge_on(bool value);
  bool spin_tracking_on() const; // 0D_NOT_logical
  void set_spin_tracking_on(bool value);
  bool spin_sokolov_ternov_flipping_on() const; // 0D_NOT_logical
  void set_spin_sokolov_ternov_flipping_on(bool value);
  bool radiation_damping_on() const; // 0D_NOT_logical
  void set_radiation_damping_on(bool value);
  bool radiation_zero_average() const; // 0D_NOT_logical
  void set_radiation_zero_average(bool value);
  bool radiation_fluctuations_on() const; // 0D_NOT_logical
  void set_radiation_fluctuations_on(bool value);
  bool conserve_taylor_maps() const; // 0D_NOT_logical
  void set_conserve_taylor_maps(bool value);
  bool absolute_time_tracking() const; // 0D_NOT_logical
  void set_absolute_time_tracking(bool value);
  bool absolute_time_ref_shift() const; // 0D_NOT_logical
  void set_absolute_time_ref_shift(bool value);
  bool convert_to_kinetic_momentum() const; // 0D_NOT_logical
  void set_convert_to_kinetic_momentum(bool value);
  bool normalize_twiss() const; // 0D_NOT_logical
  void set_normalize_twiss(bool value);
  bool aperture_limit_on() const; // 0D_NOT_logical
  void set_aperture_limit_on(bool value);
  bool spin_n0_direction_user_set() const; // 0D_NOT_logical
  void set_spin_n0_direction_user_set(bool value);
  bool debug() const; // 0D_NOT_logical
  void set_debug(bool value);
};

template <>
struct FortranTraits<RadInt1Proxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_rad_int1_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_rad_int1_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_rad_int1_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "rad_int1_struct";
  }
};

class RadInt1Proxy : public FortranProxy<RadInt1Proxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double i0() const; // 0D_NOT_real
  void set_i0(double value);
  double i1() const; // 0D_NOT_real
  void set_i1(double value);
  double i2() const; // 0D_NOT_real
  void set_i2(double value);
  double i3() const; // 0D_NOT_real
  void set_i3(double value);
  double i4a() const; // 0D_NOT_real
  void set_i4a(double value);
  double i4b() const; // 0D_NOT_real
  void set_i4b(double value);
  double i4z() const; // 0D_NOT_real
  void set_i4z(double value);
  double i5a() const; // 0D_NOT_real
  void set_i5a(double value);
  double i5b() const; // 0D_NOT_real
  void set_i5b(double value);
  double i6b() const; // 0D_NOT_real
  void set_i6b(double value);
  double lin_i2_E4() const; // 0D_NOT_real
  void set_lin_i2_E4(double value);
  double lin_i3_E7() const; // 0D_NOT_real
  void set_lin_i3_E7(double value);
  double lin_i5a_E6() const; // 0D_NOT_real
  void set_lin_i5a_E6(double value);
  double lin_i5b_E6() const; // 0D_NOT_real
  void set_lin_i5b_E6(double value);
  double lin_norm_emit_a() const; // 0D_NOT_real
  void set_lin_norm_emit_a(double value);
  double lin_norm_emit_b() const; // 0D_NOT_real
  void set_lin_norm_emit_b(double value);
  double lin_sig_E() const; // 0D_NOT_real
  void set_lin_sig_E(double value);
  double n_steps() const; // 0D_NOT_real
  void set_n_steps(double value);
};

template <>
struct FortranTraits<RadIntBranchProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_rad_int_branch_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_rad_int_branch_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_rad_int_branch_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "rad_int_branch_struct";
  }
};

class RadIntBranchProxy : public FortranProxy<RadIntBranchProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  RadInt1ProxyArray1D ele() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<RadIntAllEleProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_rad_int_all_ele_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_rad_int_all_ele_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_rad_int_all_ele_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "rad_int_all_ele_struct";
  }
};

class RadIntAllEleProxy : public FortranProxy<RadIntAllEleProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  RadIntBranchProxyArray1D branch() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<RfStairStepProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_rf_stair_step_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_rf_stair_step_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_rf_stair_step_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "rf_stair_step_struct";
  }
};

class RfStairStepProxy : public FortranProxy<RfStairStepProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double E_tot0() const; // 0D_NOT_real
  void set_E_tot0(double value);
  double E_tot1() const; // 0D_NOT_real
  void set_E_tot1(double value);
  double p0c() const; // 0D_NOT_real
  void set_p0c(double value);
  double p1c() const; // 0D_NOT_real
  void set_p1c(double value);
  double scale() const; // 0D_NOT_real
  void set_scale(double value);
  double time() const; // 0D_NOT_real
  void set_time(double value);
  double s0() const; // 0D_NOT_real
  void set_s0(double value);
  double s() const; // 0D_NOT_real
  void set_s(double value);
  int ix_step() const; // 0D_NOT_integer
  void set_ix_step(int value);
};

template <>
struct FortranTraits<RfEleProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_rf_ele_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_rf_ele_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_rf_ele_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "rf_ele_struct";
  }
};

class RfEleProxy : public FortranProxy<RfEleProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  RfStairStepProxyArray1D steps() const; // 1D_ALLOC_type
  double ds_step() const; // 0D_NOT_real
  void set_ds_step(double value);
};

template <>
struct FortranTraits<EleProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_ele_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_ele_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_ele_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "ele_struct";
  }
};

class EleProxy : public FortranProxy<EleProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  std::string type() const; // 0D_NOT_character
  void set_type(const std::string& value);
  std::string alias() const; // 0D_NOT_character
  void set_alias(const std::string& value);
  std::string component_name() const; // 0D_NOT_character
  void set_component_name(const std::string& value);
  std::string descrip() const; // 0D_PTR_character
  void set_descrip(const std::string& value);
  TwissProxy a() const; // 0D_NOT_type
  void set_a(const TwissProxy& src);
  TwissProxy b() const; // 0D_NOT_type
  void set_b(const TwissProxy& src);
  TwissProxy z() const; // 0D_NOT_type
  void set_z(const TwissProxy& src);
  XyDispProxy x() const; // 0D_NOT_type
  void set_x(const XyDispProxy& src);
  XyDispProxy y() const; // 0D_NOT_type
  void set_y(const XyDispProxy& src);
  std::optional<AcKickerProxy> ac_kick() const; // 0D_PTR_type
  void set_ac_kick(const AcKickerProxy& src);
  BookkeepingStateProxy bookkeeping_state() const; // 0D_NOT_type
  void set_bookkeeping_state(const BookkeepingStateProxy& src);
  std::optional<BranchProxy> branch() const; // 0D_PTR_type
  void set_branch(const BranchProxy& src);
  std::optional<ControllerProxy> control() const; // 0D_PTR_type
  void set_control(const ControllerProxy& src);
  std::optional<RfEleProxy> rf() const; // 0D_PTR_type
  void set_rf(const RfEleProxy& src);
  std::optional<EleProxy> lord() const; // 0D_PTR_type
  void set_lord(const EleProxy& src);
  FloorPositionProxy floor() const; // 0D_NOT_type
  void set_floor(const FloorPositionProxy& src);
  std::optional<HighEnergySpaceChargeProxy> high_energy_space_charge()
      const; // 0D_PTR_type
  void set_high_energy_space_charge(const HighEnergySpaceChargeProxy& src);
  std::optional<Mode3Proxy> mode3() const; // 0D_PTR_type
  void set_mode3(const Mode3Proxy& src);
  std::optional<PhotonElementProxy> photon() const; // 0D_PTR_type
  void set_photon(const PhotonElementProxy& src);
  std::optional<RadMapEleProxy> rad_map() const; // 0D_PTR_type
  void set_rad_map(const RadMapEleProxy& src);
  TaylorProxyArray1D taylor() const; // 1D_NOT_type
  FArray1D<double> spin_taylor_ref_orb_in() const; // 1D_NOT_real
  TaylorProxyArray1D spin_taylor() const; // 1D_NOT_type
  std::optional<WakeProxy> wake() const; // 0D_PTR_type
  void set_wake(const WakeProxy& src);
  Wall3dProxyArray1D wall3d() const; // 1D_PTR_type
  CartesianMapProxyArray1D cartesian_map() const; // 1D_PTR_type
  CylindricalMapProxyArray1D cylindrical_map() const; // 1D_PTR_type
  GenGradMapProxyArray1D gen_grad_map() const; // 1D_PTR_type
  GridFieldProxyArray1D grid_field() const; // 1D_PTR_type
  CoordProxy map_ref_orb_in() const; // 0D_NOT_type
  void set_map_ref_orb_in(const CoordProxy& src);
  CoordProxy map_ref_orb_out() const; // 0D_NOT_type
  void set_map_ref_orb_out(const CoordProxy& src);
  CoordProxy time_ref_orb_in() const; // 0D_NOT_type
  void set_time_ref_orb_in(const CoordProxy& src);
  CoordProxy time_ref_orb_out() const; // 0D_NOT_type
  void set_time_ref_orb_out(const CoordProxy& src);
  FArray1D<double> value() const; // 1D_NOT_real
  FArray1D<double> old_value() const; // 1D_NOT_real
  FArray2D<double> spin_q() const; // 2D_NOT_real
  FArray1D<double> vec0() const; // 1D_NOT_real
  FArray2D<double> mat6() const; // 2D_NOT_real
  FArray2D<double> c_mat() const; // 2D_NOT_real
  FArray2D<double> dc_mat_dpz() const; // 2D_NOT_real
  double gamma_c() const; // 0D_NOT_real
  void set_gamma_c(double value);
  double s_start() const; // 0D_NOT_real
  void set_s_start(double value);
  double s() const; // 0D_NOT_real
  void set_s(double value);
  double ref_time() const; // 0D_NOT_real
  void set_ref_time(double value);
  FArray1D<double> a_pole() const; // 1D_PTR_real
  FArray1D<double> b_pole() const; // 1D_PTR_real
  FArray1D<double> a_pole_elec() const; // 1D_PTR_real
  FArray1D<double> b_pole_elec() const; // 1D_PTR_real
  FArray1D<double> custom() const; // 1D_PTR_real
  FArray3D<double> r() const; // 3D_PTR_real
  int key() const; // 0D_NOT_integer
  void set_key(int value);
  int sub_key() const; // 0D_NOT_integer
  void set_sub_key(int value);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int lord_status() const; // 0D_NOT_integer
  void set_lord_status(int value);
  int n_slave() const; // 0D_NOT_integer
  void set_n_slave(int value);
  int n_slave_field() const; // 0D_NOT_integer
  void set_n_slave_field(int value);
  int ix1_slave() const; // 0D_NOT_integer
  void set_ix1_slave(int value);
  int slave_status() const; // 0D_NOT_integer
  void set_slave_status(int value);
  int n_lord() const; // 0D_NOT_integer
  void set_n_lord(int value);
  int n_lord_field() const; // 0D_NOT_integer
  void set_n_lord_field(int value);
  int n_lord_ramper() const; // 0D_NOT_integer
  void set_n_lord_ramper(int value);
  int ic1_lord() const; // 0D_NOT_integer
  void set_ic1_lord(int value);
  int ix_pointer() const; // 0D_NOT_integer
  void set_ix_pointer(int value);
  int ixx() const; // 0D_NOT_integer
  void set_ixx(int value);
  int iyy() const; // 0D_NOT_integer
  void set_iyy(int value);
  int izz() const; // 0D_NOT_integer
  void set_izz(int value);
  int mat6_calc_method() const; // 0D_NOT_integer
  void set_mat6_calc_method(int value);
  int tracking_method() const; // 0D_NOT_integer
  void set_tracking_method(int value);
  int spin_tracking_method() const; // 0D_NOT_integer
  void set_spin_tracking_method(int value);
  int csr_method() const; // 0D_NOT_integer
  void set_csr_method(int value);
  int space_charge_method() const; // 0D_NOT_integer
  void set_space_charge_method(int value);
  int ptc_integration_type() const; // 0D_NOT_integer
  void set_ptc_integration_type(int value);
  int field_calc() const; // 0D_NOT_integer
  void set_field_calc(int value);
  int aperture_at() const; // 0D_NOT_integer
  void set_aperture_at(int value);
  int aperture_type() const; // 0D_NOT_integer
  void set_aperture_type(int value);
  int ref_species() const; // 0D_NOT_integer
  void set_ref_species(int value);
  int orientation() const; // 0D_NOT_integer
  void set_orientation(int value);
  bool symplectify() const; // 0D_NOT_logical
  void set_symplectify(bool value);
  bool mode_flip() const; // 0D_NOT_logical
  void set_mode_flip(bool value);
  bool multipoles_on() const; // 0D_NOT_logical
  void set_multipoles_on(bool value);
  bool scale_multipoles() const; // 0D_NOT_logical
  void set_scale_multipoles(bool value);
  bool taylor_map_includes_offsets() const; // 0D_NOT_logical
  void set_taylor_map_includes_offsets(bool value);
  bool field_master() const; // 0D_NOT_logical
  void set_field_master(bool value);
  bool is_on() const; // 0D_NOT_logical
  void set_is_on(bool value);
  bool logic() const; // 0D_NOT_logical
  void set_logic(bool value);
  bool bmad_logic() const; // 0D_NOT_logical
  void set_bmad_logic(bool value);
  bool select() const; // 0D_NOT_logical
  void set_select(bool value);
  bool offset_moves_aperture() const; // 0D_NOT_logical
  void set_offset_moves_aperture(bool value);
};

template <>
struct FortranTraits<ComplexTaylorTermProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_complex_taylor_term_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_complex_taylor_term_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_complex_taylor_term_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "complex_taylor_term_struct";
  }
};

class ComplexTaylorTermProxy : public FortranProxy<ComplexTaylorTermProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::complex<double> coef() const; // 0D_NOT_complex
  void set_coef(std::complex<double> value);
  FArray1D<int> expn() const; // 1D_NOT_integer
};

template <>
struct FortranTraits<ComplexTaylorProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_complex_taylor_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_complex_taylor_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_complex_taylor_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "complex_taylor_struct";
  }
};

class ComplexTaylorProxy : public FortranProxy<ComplexTaylorProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::complex<double> ref() const; // 0D_NOT_complex
  void set_ref(std::complex<double> value);
  ComplexTaylorTermProxyArray1D term() const; // 1D_PTR_type
};

template <>
struct FortranTraits<BranchProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_branch_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_branch_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_branch_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "branch_struct";
  }
};

class BranchProxy : public FortranProxy<BranchProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_from_branch() const; // 0D_NOT_integer
  void set_ix_from_branch(int value);
  int ix_from_ele() const; // 0D_NOT_integer
  void set_ix_from_ele(int value);
  int ix_to_ele() const; // 0D_NOT_integer
  void set_ix_to_ele(int value);
  int ix_fixer() const; // 0D_NOT_integer
  void set_ix_fixer(int value);
  int n_ele_track() const; // 0D_NOT_integer
  void set_n_ele_track(int value);
  int n_ele_max() const; // 0D_NOT_integer
  void set_n_ele_max(int value);
  std::optional<LatProxy> lat() const; // 0D_PTR_type
  void set_lat(const LatProxy& src);
  ModeInfoProxy a() const; // 0D_NOT_type
  void set_a(const ModeInfoProxy& src);
  ModeInfoProxy b() const; // 0D_NOT_type
  void set_b(const ModeInfoProxy& src);
  ModeInfoProxy z() const; // 0D_NOT_type
  void set_z(const ModeInfoProxy& src);
  EleProxyArray1D ele() const; // 1D_PTR_type
  LatParamProxy param() const; // 0D_NOT_type
  void set_param(const LatParamProxy& src);
  CoordProxy particle_start() const; // 0D_NOT_type
  void set_particle_start(const CoordProxy& src);
  Wall3dProxyArray1D wall3d() const; // 1D_PTR_type
};

template <>
struct FortranTraits<LatProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_lat_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_lat_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_lat_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "lat_struct";
  }
};

class LatProxy : public FortranProxy<LatProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string use_name() const; // 0D_NOT_character
  void set_use_name(const std::string& value);
  std::string lattice() const; // 0D_NOT_character
  void set_lattice(const std::string& value);
  std::string machine() const; // 0D_NOT_character
  void set_machine(const std::string& value);
  std::string input_file_name() const; // 0D_NOT_character
  void set_input_file_name(const std::string& value);
  std::string title() const; // 0D_NOT_character
  void set_title(const std::string& value);
  FCharArray1D print_str() const; // 1D_ALLOC_character
  ExpressionAtomProxyArray1D constant() const; // 1D_ALLOC_type
  std::optional<ModeInfoProxy> a() const; // 0D_PTR_type
  void set_a(const ModeInfoProxy& src);
  std::optional<ModeInfoProxy> b() const; // 0D_PTR_type
  void set_b(const ModeInfoProxy& src);
  std::optional<ModeInfoProxy> z() const; // 0D_PTR_type
  void set_z(const ModeInfoProxy& src);
  std::optional<LatParamProxy> param() const; // 0D_PTR_type
  void set_param(const LatParamProxy& src);
  BookkeepingStateProxy lord_state() const; // 0D_NOT_type
  void set_lord_state(const BookkeepingStateProxy& src);
  EleProxy ele_init() const; // 0D_NOT_type
  void set_ele_init(const EleProxy& src);
  EleProxyArray1D ele() const; // 1D_PTR_type
  BranchProxyArray1D branch() const; // 1D_ALLOC_type
  ControlProxyArray1D control() const; // 1D_ALLOC_type
  std::optional<CoordProxy> particle_start() const; // 0D_PTR_type
  void set_particle_start(const CoordProxy& src);
  BeamInitProxy beam_init() const; // 0D_NOT_type
  void set_beam_init(const BeamInitProxy& src);
  PreTrackerProxy pre_tracker() const; // 0D_NOT_type
  void set_pre_tracker(const PreTrackerProxy& src);
  FArray1D<double> custom() const; // 1D_ALLOC_real
  int version() const; // 0D_NOT_integer
  void set_version(int value);
  int* n_ele_track() const; // 0D_PTR_integer
  void set_n_ele_track(int value);
  int* n_ele_max() const; // 0D_PTR_integer
  void set_n_ele_max(int value);
  int n_control_max() const; // 0D_NOT_integer
  void set_n_control_max(int value);
  int n_ic_max() const; // 0D_NOT_integer
  void set_n_ic_max(int value);
  int input_taylor_order() const; // 0D_NOT_integer
  void set_input_taylor_order(int value);
  FArray1D<int> ic() const; // 1D_ALLOC_integer
  int photon_type() const; // 0D_NOT_integer
  void set_photon_type(int value);
  int creation_hash() const; // 0D_NOT_integer
  void set_creation_hash(int value);
  int ramper_slave_bookkeeping() const; // 0D_NOT_integer
  void set_ramper_slave_bookkeeping(int value);
};

template <>
struct FortranTraits<BunchProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bunch_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bunch_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bunch_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bunch_struct";
  }
};

class BunchProxy : public FortranProxy<BunchProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  CoordProxyArray1D particle() const; // 1D_ALLOC_type
  FArray1D<int> ix_z() const; // 1D_ALLOC_integer
  double charge_tot() const; // 0D_NOT_real
  void set_charge_tot(double value);
  double charge_live() const; // 0D_NOT_real
  void set_charge_live(double value);
  double z_center() const; // 0D_NOT_real
  void set_z_center(double value);
  double t_center() const; // 0D_NOT_real
  void set_t_center(double value);
  double t0() const; // 0D_NOT_real
  void set_t0(double value);
  bool drift_between_t_and_s() const; // 0D_NOT_logical
  void set_drift_between_t_and_s(bool value);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_bunch() const; // 0D_NOT_integer
  void set_ix_bunch(int value);
  int ix_turn() const; // 0D_NOT_integer
  void set_ix_turn(int value);
  int n_live() const; // 0D_NOT_integer
  void set_n_live(int value);
  int n_good() const; // 0D_NOT_integer
  void set_n_good(int value);
  int n_bad() const; // 0D_NOT_integer
  void set_n_bad(int value);
};

template <>
struct FortranTraits<BunchParamsProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bunch_params_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bunch_params_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bunch_params_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bunch_params_struct";
  }
};

class BunchParamsProxy : public FortranProxy<BunchParamsProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  CoordProxy centroid() const; // 0D_NOT_type
  void set_centroid(const CoordProxy& src);
  TwissProxy x() const; // 0D_NOT_type
  void set_x(const TwissProxy& src);
  TwissProxy y() const; // 0D_NOT_type
  void set_y(const TwissProxy& src);
  TwissProxy z() const; // 0D_NOT_type
  void set_z(const TwissProxy& src);
  TwissProxy a() const; // 0D_NOT_type
  void set_a(const TwissProxy& src);
  TwissProxy b() const; // 0D_NOT_type
  void set_b(const TwissProxy& src);
  TwissProxy c() const; // 0D_NOT_type
  void set_c(const TwissProxy& src);
  FArray2D<double> sigma() const; // 2D_NOT_real
  FArray1D<double> rel_max() const; // 1D_NOT_real
  FArray1D<double> rel_min() const; // 1D_NOT_real
  double s() const; // 0D_NOT_real
  void set_s(double value);
  double t() const; // 0D_NOT_real
  void set_t(double value);
  double sigma_t() const; // 0D_NOT_real
  void set_sigma_t(double value);
  double charge_live() const; // 0D_NOT_real
  void set_charge_live(double value);
  double charge_tot() const; // 0D_NOT_real
  void set_charge_tot(double value);
  int n_particle_tot() const; // 0D_NOT_integer
  void set_n_particle_tot(int value);
  int n_particle_live() const; // 0D_NOT_integer
  void set_n_particle_live(int value);
  int n_particle_lost_in_ele() const; // 0D_NOT_integer
  void set_n_particle_lost_in_ele(int value);
  int n_good_steps() const; // 0D_NOT_integer
  void set_n_good_steps(int value);
  int n_bad_steps() const; // 0D_NOT_integer
  void set_n_bad_steps(int value);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int location() const; // 0D_NOT_integer
  void set_location(int value);
  bool twiss_valid() const; // 0D_NOT_logical
  void set_twiss_valid(bool value);
};

template <>
struct FortranTraits<BeamProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_beam_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_beam_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_beam_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "beam_struct";
  }
};

class BeamProxy : public FortranProxy<BeamProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  BunchProxyArray1D bunch() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<AperturePointProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_aperture_point_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_aperture_point_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_aperture_point_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "aperture_point_struct";
  }
};

class AperturePointProxy : public FortranProxy<AperturePointProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double x() const; // 0D_NOT_real
  void set_x(double value);
  double y() const; // 0D_NOT_real
  void set_y(double value);
  int plane() const; // 0D_NOT_integer
  void set_plane(int value);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int i_turn() const; // 0D_NOT_integer
  void set_i_turn(int value);
};

template <>
struct FortranTraits<ApertureParamProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_aperture_param_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_aperture_param_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_aperture_param_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "aperture_param_struct";
  }
};

class ApertureParamProxy : public FortranProxy<ApertureParamProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double min_angle() const; // 0D_NOT_real
  void set_min_angle(double value);
  double max_angle() const; // 0D_NOT_real
  void set_max_angle(double value);
  int n_angle() const; // 0D_NOT_integer
  void set_n_angle(int value);
  int n_turn() const; // 0D_NOT_integer
  void set_n_turn(int value);
  double x_init() const; // 0D_NOT_real
  void set_x_init(double value);
  double y_init() const; // 0D_NOT_real
  void set_y_init(double value);
  double rel_accuracy() const; // 0D_NOT_real
  void set_rel_accuracy(double value);
  double abs_accuracy() const; // 0D_NOT_real
  void set_abs_accuracy(double value);
  std::string start_ele() const; // 0D_NOT_character
  void set_start_ele(const std::string& value);
};

template <>
struct FortranTraits<ApertureScanProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_aperture_scan_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_aperture_scan_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_aperture_scan_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "aperture_scan_struct";
  }
};

class ApertureScanProxy : public FortranProxy<ApertureScanProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  AperturePointProxyArray1D point() const; // 1D_ALLOC_type
  CoordProxy ref_orb() const; // 0D_NOT_type
  void set_ref_orb(const CoordProxy& src);
  double pz_start() const; // 0D_NOT_real
  void set_pz_start(double value);
};

template <>
struct FortranTraits<ElePointerProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_ele_pointer_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_ele_pointer_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_ele_pointer_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "ele_pointer_struct";
  }
};

class ElePointerProxy : public FortranProxy<ElePointerProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<EleProxy> ele() const; // 0D_PTR_type
  void set_ele(const EleProxy& src);
  LatEleLocProxy loc() const; // 0D_NOT_type
  void set_loc(const LatEleLocProxy& src);
  int id() const; // 0D_NOT_integer
  void set_id(int value);
};

template <>
struct FortranTraits<ExpressionTreeProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_expression_tree_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_expression_tree_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_expression_tree_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "expression_tree_struct";
  }
};

class ExpressionTreeProxy : public FortranProxy<ExpressionTreeProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  int type() const; // 0D_NOT_integer
  void set_type(int value);
  double value() const; // 0D_NOT_real
  void set_value(double value);
  ExpressionTreeProxyArray1D node() const; // 1D_PTR_type
};

template <>
struct FortranTraits<NametableProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_nametable_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_nametable_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_nametable_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "nametable_struct";
  }
};

class NametableProxy : public FortranProxy<NametableProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FCharArray1D name() const; // 1D_ALLOC_character
  FArray1D<int> index() const; // 1D_ALLOC_integer
  int n_min() const; // 0D_NOT_integer
  void set_n_min(int value);
  int n_max() const; // 0D_NOT_integer
  void set_n_max(int value);
};

template <>
struct FortranTraits<TaoSpinDnDpzProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_spin_dn_dpz_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_spin_dn_dpz_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_spin_dn_dpz_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_spin_dn_dpz_struct";
  }
};

class TaoSpinDnDpzProxy : public FortranProxy<TaoSpinDnDpzProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> vec() const; // 1D_NOT_real
  FArray2D<double> partial() const; // 2D_NOT_real
  FArray2D<double> partial2() const; // 2D_NOT_real
};

template <>
struct FortranTraits<ResonanceHProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_resonance_h_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_resonance_h_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_resonance_h_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "resonance_h_struct";
  }
};

class ResonanceHProxy : public FortranProxy<ResonanceHProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string id() const; // 0D_NOT_character
  void set_id(const std::string& value);
  std::complex<double> c_val() const; // 0D_NOT_complex
  void set_c_val(std::complex<double> value);
};

template <>
struct FortranTraits<SpinOrbitMap1Proxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_spin_orbit_map1_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_spin_orbit_map1_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_spin_orbit_map1_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "spin_orbit_map1_struct";
  }
};

class SpinOrbitMap1Proxy : public FortranProxy<SpinOrbitMap1Proxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray2D<double> orb_mat() const; // 2D_NOT_real
  FArray1D<double> vec0() const; // 1D_NOT_real
  FArray2D<double> spin_q() const; // 2D_NOT_real
};

template <>
struct FortranTraits<SpinAxisProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_spin_axis_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_spin_axis_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_spin_axis_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "spin_axis_struct";
  }
};

class SpinAxisProxy : public FortranProxy<SpinAxisProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> l() const; // 1D_NOT_real
  FArray1D<double> n0() const; // 1D_NOT_real
  FArray1D<double> m() const; // 1D_NOT_real
};

template <>
struct FortranTraits<PtcNormalFormProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_ptc_normal_form_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_ptc_normal_form_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_ptc_normal_form_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "ptc_normal_form_struct";
  }
};

class PtcNormalFormProxy : public FortranProxy<PtcNormalFormProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<EleProxy> ele_origin() const; // 0D_PTR_type
  void set_ele_origin(const EleProxy& src);
  FArray1D<double> orb0() const; // 1D_NOT_real
  bool valid_map() const; // 0D_NOT_logical
  void set_valid_map(bool value);
};

template <>
struct FortranTraits<BmadNormalFormProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bmad_normal_form_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bmad_normal_form_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bmad_normal_form_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bmad_normal_form_struct";
  }
};

class BmadNormalFormProxy : public FortranProxy<BmadNormalFormProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<EleProxy> ele_origin() const; // 0D_PTR_type
  void set_ele_origin(const EleProxy& src);
  TaylorProxyArray1D M() const; // 1D_NOT_type
  TaylorProxyArray1D A() const; // 1D_NOT_type
  TaylorProxyArray1D A_inv() const; // 1D_NOT_type
  TaylorProxyArray1D dhdj() const; // 1D_NOT_type
  ComplexTaylorProxyArray1D F() const; // 1D_NOT_type
  ComplexTaylorProxyArray1D L() const; // 1D_NOT_type
  ResonanceHProxyArray1D h() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<BunchTrackProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bunch_track_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bunch_track_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bunch_track_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bunch_track_struct";
  }
};

class BunchTrackProxy : public FortranProxy<BunchTrackProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  BunchParamsProxyArray1D pt() const; // 1D_ALLOC_type
  double ds_save() const; // 0D_NOT_real
  void set_ds_save(double value);
  int n_pt() const; // 0D_NOT_integer
  void set_n_pt(int value);
};

template <>
struct FortranTraits<SummationRdtProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_summation_rdt_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_summation_rdt_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_summation_rdt_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "summation_rdt_struct";
  }
};

class SummationRdtProxy : public FortranProxy<SummationRdtProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::complex<double> h11001() const; // 0D_NOT_complex
  void set_h11001(std::complex<double> value);
  std::complex<double> h00111() const; // 0D_NOT_complex
  void set_h00111(std::complex<double> value);
  std::complex<double> h20001() const; // 0D_NOT_complex
  void set_h20001(std::complex<double> value);
  std::complex<double> h00201() const; // 0D_NOT_complex
  void set_h00201(std::complex<double> value);
  std::complex<double> h10002() const; // 0D_NOT_complex
  void set_h10002(std::complex<double> value);
  std::complex<double> h21000() const; // 0D_NOT_complex
  void set_h21000(std::complex<double> value);
  std::complex<double> h30000() const; // 0D_NOT_complex
  void set_h30000(std::complex<double> value);
  std::complex<double> h10110() const; // 0D_NOT_complex
  void set_h10110(std::complex<double> value);
  std::complex<double> h10020() const; // 0D_NOT_complex
  void set_h10020(std::complex<double> value);
  std::complex<double> h10200() const; // 0D_NOT_complex
  void set_h10200(std::complex<double> value);
  std::complex<double> h31000() const; // 0D_NOT_complex
  void set_h31000(std::complex<double> value);
  std::complex<double> h40000() const; // 0D_NOT_complex
  void set_h40000(std::complex<double> value);
  std::complex<double> h20110() const; // 0D_NOT_complex
  void set_h20110(std::complex<double> value);
  std::complex<double> h11200() const; // 0D_NOT_complex
  void set_h11200(std::complex<double> value);
  std::complex<double> h20020() const; // 0D_NOT_complex
  void set_h20020(std::complex<double> value);
  std::complex<double> h20200() const; // 0D_NOT_complex
  void set_h20200(std::complex<double> value);
  std::complex<double> h00310() const; // 0D_NOT_complex
  void set_h00310(std::complex<double> value);
  std::complex<double> h00400() const; // 0D_NOT_complex
  void set_h00400(std::complex<double> value);
  std::complex<double> h22000() const; // 0D_NOT_complex
  void set_h22000(std::complex<double> value);
  std::complex<double> h00220() const; // 0D_NOT_complex
  void set_h00220(std::complex<double> value);
  std::complex<double> h11110() const; // 0D_NOT_complex
  void set_h11110(std::complex<double> value);
};

template <>
struct FortranTraits<TaoEleShapeProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_ele_shape_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_ele_shape_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_ele_shape_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_ele_shape_struct";
  }
};

class TaoEleShapeProxy : public FortranProxy<TaoEleShapeProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string ele_id() const; // 0D_NOT_character
  void set_ele_id(const std::string& value);
  std::string shape() const; // 0D_NOT_character
  void set_shape(const std::string& value);
  std::string color() const; // 0D_NOT_character
  void set_color(const std::string& value);
  double size() const; // 0D_NOT_real
  void set_size(double value);
  std::string label() const; // 0D_NOT_character
  void set_label(const std::string& value);
  bool draw() const; // 0D_NOT_logical
  void set_draw(bool value);
  bool multi() const; // 0D_NOT_logical
  void set_multi(bool value);
  int line_width() const; // 0D_NOT_integer
  void set_line_width(int value);
  double offset() const; // 0D_NOT_real
  void set_offset(double value);
  int ix_key() const; // 0D_NOT_integer
  void set_ix_key(int value);
  std::string name_ele() const; // 0D_NOT_character
  void set_name_ele(const std::string& value);
  TaoElePointerProxyArray1D uni() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoElePointerProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_ele_pointer_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_ele_pointer_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_ele_pointer_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_ele_pointer_struct";
  }
};

class TaoElePointerProxy : public FortranProxy<TaoElePointerProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  ElePointerProxyArray1D eles() const; // 1D_ALLOC_type
  int n_loc() const; // 0D_NOT_integer
  void set_n_loc(int value);
};

template <>
struct FortranTraits<TaoCurveProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_curve_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_curve_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_curve_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_curve_struct";
  }
};

class TaoCurveProxy : public FortranProxy<TaoCurveProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  std::string data_source() const; // 0D_NOT_character
  void set_data_source(const std::string& value);
  std::string data_index() const; // 0D_NOT_character
  void set_data_index(const std::string& value);
  std::string data_type_x() const; // 0D_NOT_character
  void set_data_type_x(const std::string& value);
  std::string data_type() const; // 0D_ALLOC_character
  void set_data_type(const std::string& value);
  std::string ele_ref_name() const; // 0D_NOT_character
  void set_ele_ref_name(const std::string& value);
  std::string legend_text() const; // 0D_NOT_character
  void set_legend_text(const std::string& value);
  std::string message_text() const; // 0D_NOT_character
  void set_message_text(const std::string& value);
  std::string component() const; // 0D_NOT_character
  void set_component(const std::string& value);
  std::string why_invalid() const; // 0D_NOT_character
  void set_why_invalid(const std::string& value);
  std::optional<TaoGraphProxy> g() const; // 0D_PTR_type
  void set_g(const TaoGraphProxy& src);
  TaoHistogramProxy hist() const; // 0D_NOT_type
  void set_hist(const TaoHistogramProxy& src);
  TaoCurveColorProxy z_color() const; // 0D_NOT_type
  void set_z_color(const TaoCurveColorProxy& src);
  FArray1D<double> x_line() const; // 1D_ALLOC_real
  FArray1D<double> y_line() const; // 1D_ALLOC_real
  FArray1D<double> y2_line() const; // 1D_ALLOC_real
  FArray1D<int> ix_line() const; // 1D_ALLOC_integer
  FArray1D<double> x_symb() const; // 1D_ALLOC_real
  FArray1D<double> y_symb() const; // 1D_ALLOC_real
  FArray1D<double> z_symb() const; // 1D_ALLOC_real
  FArray1D<double> err_symb() const; // 1D_ALLOC_real
  FArray1D<double> symb_size() const; // 1D_ALLOC_real
  FArray1D<int> ix_symb() const; // 1D_ALLOC_integer
  double y_axis_scale_factor() const; // 0D_NOT_real
  void set_y_axis_scale_factor(double value);
  QpLineProxy line() const; // 0D_NOT_type
  void set_line(const QpLineProxy& src);
  QpSymbolProxy symbol() const; // 0D_NOT_type
  void set_symbol(const QpSymbolProxy& src);
  TaoCurveOrbitProxy orbit() const; // 0D_NOT_type
  void set_orbit(const TaoCurveOrbitProxy& src);
  int ix_universe() const; // 0D_NOT_integer
  void set_ix_universe(int value);
  int symbol_every() const; // 0D_NOT_integer
  void set_symbol_every(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_bunch() const; // 0D_NOT_integer
  void set_ix_bunch(int value);
  int n_turn() const; // 0D_NOT_integer
  void set_n_turn(int value);
  bool use_y2() const; // 0D_NOT_logical
  void set_use_y2(bool value);
  bool draw_line() const; // 0D_NOT_logical
  void set_draw_line(bool value);
  bool draw_symbols() const; // 0D_NOT_logical
  void set_draw_symbols(bool value);
  bool draw_symbol_index() const; // 0D_NOT_logical
  void set_draw_symbol_index(bool value);
  bool draw_error_bars() const; // 0D_NOT_logical
  void set_draw_error_bars(bool value);
  bool smooth_line_calc() const; // 0D_NOT_logical
  void set_smooth_line_calc(bool value);
  bool valid() const; // 0D_NOT_logical
  void set_valid(bool value);
};

template <>
struct FortranTraits<TaoCurveColorProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_curve_color_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_curve_color_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_curve_color_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_curve_color_struct";
  }
};

class TaoCurveColorProxy : public FortranProxy<TaoCurveColorProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string data_type() const; // 0D_NOT_character
  void set_data_type(const std::string& value);
  bool is_on() const; // 0D_NOT_logical
  void set_is_on(bool value);
  double min() const; // 0D_NOT_real
  void set_min(double value);
  double max() const; // 0D_NOT_real
  void set_max(double value);
  bool autoscale() const; // 0D_NOT_logical
  void set_autoscale(bool value);
};

template <>
struct FortranTraits<TaoCurveOrbitProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_curve_orbit_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_curve_orbit_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_curve_orbit_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_curve_orbit_struct";
  }
};

class TaoCurveOrbitProxy : public FortranProxy<TaoCurveOrbitProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double x() const; // 0D_NOT_real
  void set_x(double value);
  double y() const; // 0D_NOT_real
  void set_y(double value);
  double t() const; // 0D_NOT_real
  void set_t(double value);
};

template <>
struct FortranTraits<TaoHistogramProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_histogram_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_histogram_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_histogram_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_histogram_struct";
  }
};

class TaoHistogramProxy : public FortranProxy<TaoHistogramProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool density_normalized() const; // 0D_NOT_logical
  void set_density_normalized(bool value);
  bool weight_by_charge() const; // 0D_NOT_logical
  void set_weight_by_charge(bool value);
  double minimum() const; // 0D_NOT_real
  void set_minimum(double value);
  double maximum() const; // 0D_NOT_real
  void set_maximum(double value);
  double width() const; // 0D_NOT_real
  void set_width(double value);
  double center() const; // 0D_NOT_real
  void set_center(double value);
  int number() const; // 0D_NOT_integer
  void set_number(int value);
};

template <>
struct FortranTraits<LatEleOrder1Proxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_lat_ele_order1_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_lat_ele_order1_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_lat_ele_order1_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "lat_ele_order1_struct";
  }
};

class LatEleOrder1Proxy : public FortranProxy<LatEleOrder1Proxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_order() const; // 0D_NOT_integer
  void set_ix_order(int value);
};

template <>
struct FortranTraits<LatEleOrderArrayProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_lat_ele_order_array_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_lat_ele_order_array_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_lat_ele_order_array_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "lat_ele_order_array_struct";
  }
};

class LatEleOrderArrayProxy : public FortranProxy<LatEleOrderArrayProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  LatEleOrder1ProxyArray1D ele() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoLatSigmaProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_lat_sigma_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_lat_sigma_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_lat_sigma_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_lat_sigma_struct";
  }
};

class TaoLatSigmaProxy : public FortranProxy<TaoLatSigmaProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray2D<double> mat() const; // 2D_NOT_real
};

template <>
struct FortranTraits<TaoSpinEleProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_spin_ele_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_spin_ele_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_spin_ele_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_spin_ele_struct";
  }
};

class TaoSpinEleProxy : public FortranProxy<TaoSpinEleProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoSpinDnDpzProxy dn_dpz() const; // 0D_NOT_type
  void set_dn_dpz(const TaoSpinDnDpzProxy& src);
  FArray1D<double> orb_eigen_val() const; // 1D_NOT_real
  FArray2D<double> orb_eigen_vec() const; // 2D_NOT_real
  FArray2D<double> spin_eigen_vec() const; // 2D_NOT_real
  bool valid() const; // 0D_NOT_logical
  void set_valid(bool value);
};

template <>
struct FortranTraits<TaoPlotCacheProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_plot_cache_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_plot_cache_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_plot_cache_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_plot_cache_struct";
  }
};

class TaoPlotCacheProxy : public FortranProxy<TaoPlotCacheProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  EleProxy ele_to_s() const; // 0D_NOT_type
  void set_ele_to_s(const EleProxy& src);
  CoordProxy orbit() const; // 0D_NOT_type
  void set_orbit(const CoordProxy& src);
  bool err() const; // 0D_NOT_logical
  void set_err(bool value);
};

template <>
struct FortranTraits<TaoSpinPolarizationProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_spin_polarization_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_spin_polarization_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_spin_polarization_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_spin_polarization_struct";
  }
};

class TaoSpinPolarizationProxy : public FortranProxy<TaoSpinPolarizationProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double tune() const; // 0D_NOT_real
  void set_tune(double value);
  double pol_limit_st() const; // 0D_NOT_real
  void set_pol_limit_st(double value);
  double pol_limit_dk() const; // 0D_NOT_real
  void set_pol_limit_dk(double value);
  FArray1D<double> pol_limit_dk_partial() const; // 1D_NOT_real
  FArray1D<double> pol_limit_dk_partial2() const; // 1D_NOT_real
  double pol_rate_bks() const; // 0D_NOT_real
  void set_pol_rate_bks(double value);
  double depol_rate() const; // 0D_NOT_real
  void set_depol_rate(double value);
  FArray1D<double> depol_rate_partial() const; // 1D_NOT_real
  FArray1D<double> depol_rate_partial2() const; // 1D_NOT_real
  double integral_bn() const; // 0D_NOT_real
  void set_integral_bn(double value);
  double integral_bdn() const; // 0D_NOT_real
  void set_integral_bdn(double value);
  double integral_1ns() const; // 0D_NOT_real
  void set_integral_1ns(double value);
  double integral_dn2() const; // 0D_NOT_real
  void set_integral_dn2(double value);
  bool valid() const; // 0D_NOT_logical
  void set_valid(bool value);
  SpinOrbitMap1Proxy q_1turn() const; // 0D_NOT_type
  void set_q_1turn(const SpinOrbitMap1Proxy& src);
  SpinOrbitMap1ProxyArray1D q_ele() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoLatticeBranchProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_lattice_branch_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_lattice_branch_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_lattice_branch_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_lattice_branch_struct";
  }
};

class TaoLatticeBranchProxy : public FortranProxy<TaoLatticeBranchProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<TaoLatticeProxy> tao_lat() const; // 0D_PTR_type
  void set_tao_lat(const TaoLatticeProxy& src);
  TaoLatSigmaProxyArray1D lat_sigma() const; // 1D_ALLOC_type
  TaoSpinEleProxyArray1D spin_ele() const; // 1D_ALLOC_type
  BunchParamsProxyArray1D bunch_params() const; // 1D_ALLOC_type
  BunchTrackProxyArray1D bunch_params_comb() const; // 1D_ALLOC_type
  CoordProxyArray1D orbit() const; // 1D_ALLOC_type
  TaoPlotCacheProxyArray1D plot_cache() const; // 1D_ALLOC_type
  TaoSpinPolarizationProxy spin() const; // 0D_NOT_type
  void set_spin(const TaoSpinPolarizationProxy& src);
  SummationRdtProxy srdt() const; // 0D_NOT_type
  void set_srdt(const SummationRdtProxy& src);
  CoordProxy orb0() const; // 0D_NOT_type
  void set_orb0(const CoordProxy& src);
  NormalModesProxy modes_ri() const; // 0D_NOT_type
  void set_modes_ri(const NormalModesProxy& src);
  NormalModesProxy modes_6d() const; // 0D_NOT_type
  void set_modes_6d(const NormalModesProxy& src);
  PtcNormalFormProxy ptc_normal_form() const; // 0D_NOT_type
  void set_ptc_normal_form(const PtcNormalFormProxy& src);
  BmadNormalFormProxy bmad_normal_form() const; // 0D_NOT_type
  void set_bmad_normal_form(const BmadNormalFormProxy& src);
  CoordProxyArray1D high_E_orb() const; // 1D_ALLOC_type
  CoordProxyArray1D low_E_orb() const; // 1D_ALLOC_type
  TaylorProxyArray1D taylor_save() const; // 1D_NOT_type
  double cache_x_min() const; // 0D_NOT_real
  void set_cache_x_min(double value);
  double cache_x_max() const; // 0D_NOT_real
  void set_cache_x_max(double value);
  double comb_ds_save() const; // 0D_NOT_real
  void set_comb_ds_save(double value);
  int ix_ref_taylor() const; // 0D_NOT_integer
  void set_ix_ref_taylor(int value);
  int ix_ele_taylor() const; // 0D_NOT_integer
  void set_ix_ele_taylor(int value);
  int track_state() const; // 0D_NOT_integer
  void set_track_state(int value);
  int cache_n_pts() const; // 0D_NOT_integer
  void set_cache_n_pts(int value);
  int ix_rad_int_cache() const; // 0D_NOT_integer
  void set_ix_rad_int_cache(int value);
  bool has_open_match_element() const; // 0D_NOT_logical
  void set_has_open_match_element(bool value);
  bool plot_cache_valid() const; // 0D_NOT_logical
  void set_plot_cache_valid(bool value);
  bool spin_map_valid() const; // 0D_NOT_logical
  void set_spin_map_valid(bool value);
  bool twiss_valid() const; // 0D_NOT_logical
  void set_twiss_valid(bool value);
  bool mode_flip_here() const; // 0D_NOT_logical
  void set_mode_flip_here(bool value);
  bool chrom_calc_ok() const; // 0D_NOT_logical
  void set_chrom_calc_ok(bool value);
  bool rad_int_calc_ok() const; // 0D_NOT_logical
  void set_rad_int_calc_ok(bool value);
  bool emit_6d_calc_ok() const; // 0D_NOT_logical
  void set_emit_6d_calc_ok(bool value);
  bool sigma_track_ok() const; // 0D_NOT_logical
  void set_sigma_track_ok(bool value);
};

template <>
struct FortranTraits<TaoModelElementProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_model_element_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_model_element_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_model_element_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_model_element_struct";
  }
};

class TaoModelElementProxy : public FortranProxy<TaoModelElementProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  BeamProxy beam() const; // 0D_NOT_type
  void set_beam(const BeamProxy& src);
  bool save_beam_internally() const; // 0D_NOT_logical
  void set_save_beam_internally(bool value);
  bool save_beam_to_file() const; // 0D_NOT_logical
  void set_save_beam_to_file(bool value);
};

template <>
struct FortranTraits<TaoBeamBranchProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_beam_branch_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_beam_branch_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_beam_branch_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_beam_branch_struct";
  }
};

class TaoBeamBranchProxy : public FortranProxy<TaoBeamBranchProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  BeamProxy beam_at_start() const; // 0D_NOT_type
  void set_beam_at_start(const BeamProxy& src);
  BeamInitProxy beam_init() const; // 0D_NOT_type
  void set_beam_init(const BeamInitProxy& src);
  BeamInitProxy beam_init_used() const; // 0D_NOT_type
  void set_beam_init_used(const BeamInitProxy& src);
  bool init_starting_distribution() const; // 0D_NOT_logical
  void set_init_starting_distribution(bool value);
  std::string track_start() const; // 0D_NOT_character
  void set_track_start(const std::string& value);
  std::string track_end() const; // 0D_NOT_character
  void set_track_end(const std::string& value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_track_start() const; // 0D_NOT_integer
  void set_ix_track_start(int value);
  int ix_track_end() const; // 0D_NOT_integer
  void set_ix_track_end(int value);
};

template <>
struct FortranTraits<TaoD1DataProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_d1_data_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_d1_data_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_d1_data_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_d1_data_struct";
  }
};

class TaoD1DataProxy : public FortranProxy<TaoD1DataProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  std::optional<TaoD2DataProxy> d2() const; // 0D_PTR_type
  void set_d2(const TaoD2DataProxy& src);
  TaoDataProxyArray1D d() const; // 1D_PTR_type
};

template <>
struct FortranTraits<TaoD2DataProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_d2_data_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_d2_data_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_d2_data_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_d2_data_struct";
  }
};

class TaoD2DataProxy : public FortranProxy<TaoD2DataProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  std::string data_file_name() const; // 0D_NOT_character
  void set_data_file_name(const std::string& value);
  std::string ref_file_name() const; // 0D_NOT_character
  void set_ref_file_name(const std::string& value);
  std::string data_date() const; // 0D_NOT_character
  void set_data_date(const std::string& value);
  std::string ref_date() const; // 0D_NOT_character
  void set_ref_date(const std::string& value);
  FCharArray1D descrip() const; // 1D_NOT_character
  TaoD1DataProxyArray1D d1() const; // 1D_ALLOC_type
  int ix_universe() const; // 0D_NOT_integer
  void set_ix_universe(int value);
  int ix_d2_data() const; // 0D_NOT_integer
  void set_ix_d2_data(int value);
  int ix_ref() const; // 0D_NOT_integer
  void set_ix_ref(int value);
  bool data_read_in() const; // 0D_NOT_logical
  void set_data_read_in(bool value);
  bool ref_read_in() const; // 0D_NOT_logical
  void set_ref_read_in(bool value);
};

template <>
struct FortranTraits<TaoDataVarComponentProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_data_var_component_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_data_var_component_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_data_var_component_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_data_var_component_struct";
  }
};

class TaoDataVarComponentProxy : public FortranProxy<TaoDataVarComponentProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  double sign() const; // 0D_NOT_real
  void set_sign(double value);
};

template <>
struct FortranTraits<TaoGraphProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_graph_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_graph_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_graph_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_graph_struct";
  }
};

class TaoGraphProxy : public FortranProxy<TaoGraphProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  std::string type() const; // 0D_NOT_character
  void set_type(const std::string& value);
  std::string title() const; // 0D_NOT_character
  void set_title(const std::string& value);
  std::string title_suffix() const; // 0D_NOT_character
  void set_title_suffix(const std::string& value);
  FCharArray1D text_legend() const; // 1D_NOT_character
  FCharArray1D text_legend_out() const; // 1D_NOT_character
  std::string why_invalid() const; // 0D_NOT_character
  void set_why_invalid(const std::string& value);
  TaoCurveProxyArray1D curve() const; // 1D_ALLOC_type
  std::optional<TaoPlotProxy> p() const; // 0D_PTR_type
  void set_p(const TaoPlotProxy& src);
  TaoFloorPlanProxy floor_plan() const; // 0D_NOT_type
  void set_floor_plan(const TaoFloorPlanProxy& src);
  QpPointProxy text_legend_origin() const; // 0D_NOT_type
  void set_text_legend_origin(const QpPointProxy& src);
  QpPointProxy curve_legend_origin() const; // 0D_NOT_type
  void set_curve_legend_origin(const QpPointProxy& src);
  QpLegendProxy curve_legend() const; // 0D_NOT_type
  void set_curve_legend(const QpLegendProxy& src);
  QpAxisProxy x() const; // 0D_NOT_type
  void set_x(const QpAxisProxy& src);
  QpAxisProxy y() const; // 0D_NOT_type
  void set_y(const QpAxisProxy& src);
  QpAxisProxy x2() const; // 0D_NOT_type
  void set_x2(const QpAxisProxy& src);
  QpAxisProxy y2() const; // 0D_NOT_type
  void set_y2(const QpAxisProxy& src);
  QpRectProxy margin() const; // 0D_NOT_type
  void set_margin(const QpRectProxy& src);
  QpRectProxy scale_margin() const; // 0D_NOT_type
  void set_scale_margin(const QpRectProxy& src);
  double x_axis_scale_factor() const; // 0D_NOT_real
  void set_x_axis_scale_factor(double value);
  double symbol_size_scale() const; // 0D_NOT_real
  void set_symbol_size_scale(double value);
  FArray1D<int> box() const; // 1D_NOT_integer
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_universe() const; // 0D_NOT_integer
  void set_ix_universe(int value);
  bool clip() const; // 0D_NOT_logical
  void set_clip(bool value);
  bool y2_mirrors_y() const; // 0D_NOT_logical
  void set_y2_mirrors_y(bool value);
  bool limited() const; // 0D_NOT_logical
  void set_limited(bool value);
  bool draw_axes() const; // 0D_NOT_logical
  void set_draw_axes(bool value);
  bool draw_curve_legend() const; // 0D_NOT_logical
  void set_draw_curve_legend(bool value);
  bool draw_grid() const; // 0D_NOT_logical
  void set_draw_grid(bool value);
  bool draw_title() const; // 0D_NOT_logical
  void set_draw_title(bool value);
  bool draw_only_good_user_data_or_vars() const; // 0D_NOT_logical
  void set_draw_only_good_user_data_or_vars(bool value);
  bool allow_wrap_around() const; // 0D_NOT_logical
  void set_allow_wrap_around(bool value);
  bool is_valid() const; // 0D_NOT_logical
  void set_is_valid(bool value);
};

template <>
struct FortranTraits<TaoPlotProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_plot_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_plot_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_plot_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_plot_struct";
  }
};

class TaoPlotProxy : public FortranProxy<TaoPlotProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  std::string description() const; // 0D_NOT_character
  void set_description(const std::string& value);
  TaoGraphProxyArray1D graph() const; // 1D_ALLOC_type
  std::optional<TaoPlotRegionProxy> r() const; // 0D_PTR_type
  void set_r(const TaoPlotRegionProxy& src);
  int ix_plot() const; // 0D_NOT_integer
  void set_ix_plot(int value);
  int n_curve_pts() const; // 0D_NOT_integer
  void set_n_curve_pts(int value);
  std::string type() const; // 0D_NOT_character
  void set_type(const std::string& value);
  std::string x_axis_type() const; // 0D_NOT_character
  void set_x_axis_type(const std::string& value);
  bool autoscale_x() const; // 0D_NOT_logical
  void set_autoscale_x(bool value);
  bool autoscale_y() const; // 0D_NOT_logical
  void set_autoscale_y(bool value);
  bool autoscale_gang_x() const; // 0D_NOT_logical
  void set_autoscale_gang_x(bool value);
  bool autoscale_gang_y() const; // 0D_NOT_logical
  void set_autoscale_gang_y(bool value);
  bool list_with_show_plot_command() const; // 0D_NOT_logical
  void set_list_with_show_plot_command(bool value);
  bool phantom() const; // 0D_NOT_logical
  void set_phantom(bool value);
  bool default_plot() const; // 0D_NOT_logical
  void set_default_plot(bool value);
};

template <>
struct FortranTraits<TaoPlotRegionProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_plot_region_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_plot_region_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_plot_region_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_plot_region_struct";
  }
};

class TaoPlotRegionProxy : public FortranProxy<TaoPlotRegionProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  TaoPlotProxy plot() const; // 0D_NOT_type
  void set_plot(const TaoPlotProxy& src);
  FArray1D<double> location() const; // 1D_NOT_real
  bool visible() const; // 0D_NOT_logical
  void set_visible(bool value);
  bool list_with_show_plot_command() const; // 0D_NOT_logical
  void set_list_with_show_plot_command(bool value);
  bool setup_done() const; // 0D_NOT_logical
  void set_setup_done(bool value);
};

template <>
struct FortranTraits<TaoUniversePointerProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_universe_pointer_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_universe_pointer_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_universe_pointer_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_universe_pointer_struct";
  }
};

class TaoUniversePointerProxy : public FortranProxy<TaoUniversePointerProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<TaoUniverseProxy> u() const; // 0D_PTR_type
  void set_u(const TaoUniverseProxy& src);
};

template <>
struct FortranTraits<TaoSuperUniverseProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_super_universe_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_super_universe_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_super_universe_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_super_universe_struct";
  }
};

class TaoSuperUniverseProxy : public FortranProxy<TaoSuperUniverseProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoGlobalProxy global() const; // 0D_NOT_type
  void set_global(const TaoGlobalProxy& src);
  TaoInitProxy init() const; // 0D_NOT_type
  void set_init(const TaoInitProxy& src);
  TaoCommonProxy com() const; // 0D_NOT_type
  void set_com(const TaoCommonProxy& src);
  TaoPlotPageProxy plot_page() const; // 0D_NOT_type
  void set_plot_page(const TaoPlotPageProxy& src);
  TaoV1VarProxyArray1D v1_var() const; // 1D_ALLOC_type
  TaoVarProxyArray1D var() const; // 1D_ALLOC_type
  TaoUniverseProxyArray1D u() const; // 1D_ALLOC_type
  FArray1D<int> key() const; // 1D_ALLOC_integer
  TaoBuildingWallProxy building_wall() const; // 0D_NOT_type
  void set_building_wall(const TaoBuildingWallProxy& src);
  TaoWaveProxy wave() const; // 0D_NOT_type
  void set_wave(const TaoWaveProxy& src);
  int n_var_used() const; // 0D_NOT_integer
  void set_n_var_used(int value);
  int n_v1_var_used() const; // 0D_NOT_integer
  void set_n_v1_var_used(int value);
  TaoCmdHistoryProxyArray1D history() const; // 1D_NOT_type
  bool initialized() const; // 0D_NOT_logical
  void set_initialized(bool value);
};

template <>
struct FortranTraits<TaoVarProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_var_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_var_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_var_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_var_struct";
  }
};

class TaoVarProxy : public FortranProxy<TaoVarProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string ele_name() const; // 0D_NOT_character
  void set_ele_name(const std::string& value);
  std::string attrib_name() const; // 0D_NOT_character
  void set_attrib_name(const std::string& value);
  std::string id() const; // 0D_NOT_character
  void set_id(const std::string& value);
  TaoVarSlaveProxyArray1D slave() const; // 1D_ALLOC_type
  int ix_v1() const; // 0D_NOT_integer
  void set_ix_v1(int value);
  int ix_var() const; // 0D_NOT_integer
  void set_ix_var(int value);
  int ix_dvar() const; // 0D_NOT_integer
  void set_ix_dvar(int value);
  int ix_attrib() const; // 0D_NOT_integer
  void set_ix_attrib(int value);
  int ix_key_table() const; // 0D_NOT_integer
  void set_ix_key_table(int value);
  double* model_value() const; // 0D_PTR_real
  void set_model_value(double value);
  double* base_value() const; // 0D_PTR_real
  void set_base_value(double value);
  double design_value() const; // 0D_NOT_real
  void set_design_value(double value);
  double scratch_value() const; // 0D_NOT_real
  void set_scratch_value(double value);
  double old_value() const; // 0D_NOT_real
  void set_old_value(double value);
  double meas_value() const; // 0D_NOT_real
  void set_meas_value(double value);
  double ref_value() const; // 0D_NOT_real
  void set_ref_value(double value);
  double correction_value() const; // 0D_NOT_real
  void set_correction_value(double value);
  double high_lim() const; // 0D_NOT_real
  void set_high_lim(double value);
  double low_lim() const; // 0D_NOT_real
  void set_low_lim(double value);
  double step() const; // 0D_NOT_real
  void set_step(double value);
  double weight() const; // 0D_NOT_real
  void set_weight(double value);
  double delta_merit() const; // 0D_NOT_real
  void set_delta_merit(double value);
  double merit() const; // 0D_NOT_real
  void set_merit(double value);
  double dMerit_dVar() const; // 0D_NOT_real
  void set_dMerit_dVar(double value);
  double key_val0() const; // 0D_NOT_real
  void set_key_val0(double value);
  double key_delta() const; // 0D_NOT_real
  void set_key_delta(double value);
  double s() const; // 0D_NOT_real
  void set_s(double value);
  double extend_val() const; // 0D_NOT_real
  void set_extend_val(double value);
  std::string merit_type() const; // 0D_NOT_character
  void set_merit_type(const std::string& value);
  bool exists() const; // 0D_NOT_logical
  void set_exists(bool value);
  bool good_var() const; // 0D_NOT_logical
  void set_good_var(bool value);
  bool good_user() const; // 0D_NOT_logical
  void set_good_user(bool value);
  bool good_opt() const; // 0D_NOT_logical
  void set_good_opt(bool value);
  bool good_plot() const; // 0D_NOT_logical
  void set_good_plot(bool value);
  bool useit_opt() const; // 0D_NOT_logical
  void set_useit_opt(bool value);
  bool useit_plot() const; // 0D_NOT_logical
  void set_useit_plot(bool value);
  bool key_bound() const; // 0D_NOT_logical
  void set_key_bound(bool value);
  std::optional<TaoV1VarProxy> v1() const; // 0D_PTR_type
  void set_v1(const TaoV1VarProxy& src);
};

template <>
struct FortranTraits<TaoVarSlaveProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_var_slave_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_var_slave_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_var_slave_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_var_slave_struct";
  }
};

class TaoVarSlaveProxy : public FortranProxy<TaoVarSlaveProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int ix_uni() const; // 0D_NOT_integer
  void set_ix_uni(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  double* model_value() const; // 0D_PTR_real
  void set_model_value(double value);
  double* base_value() const; // 0D_PTR_real
  void set_base_value(double value);
};

template <>
struct FortranTraits<TaoLatticeProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_lattice_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_lattice_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_lattice_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_lattice_struct";
  }
};

class TaoLatticeProxy : public FortranProxy<TaoLatticeProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  LatProxy lat() const; // 0D_NOT_type
  void set_lat(const LatProxy& src);
  LatProxy high_E_lat() const; // 0D_NOT_type
  void set_high_E_lat(const LatProxy& src);
  LatProxy low_E_lat() const; // 0D_NOT_type
  void set_low_E_lat(const LatProxy& src);
  RadIntAllEleProxy rad_int_by_ele_ri() const; // 0D_NOT_type
  void set_rad_int_by_ele_ri(const RadIntAllEleProxy& src);
  RadIntAllEleProxy rad_int_by_ele_6d() const; // 0D_NOT_type
  void set_rad_int_by_ele_6d(const RadIntAllEleProxy& src);
  TaoLatticeBranchProxyArray1D tao_branch() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoBeamUniProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_beam_uni_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_beam_uni_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_beam_uni_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_beam_uni_struct";
  }
};

class TaoBeamUniProxy : public FortranProxy<TaoBeamUniProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string saved_at() const; // 0D_NOT_character
  void set_saved_at(const std::string& value);
  std::string dump_file() const; // 0D_NOT_character
  void set_dump_file(const std::string& value);
  std::string dump_at() const; // 0D_NOT_character
  void set_dump_at(const std::string& value);
  bool track_beam_in_universe() const; // 0D_NOT_logical
  void set_track_beam_in_universe(bool value);
  bool always_reinit() const; // 0D_NOT_logical
  void set_always_reinit(bool value);
};

template <>
struct FortranTraits<TaoDynamicApertureProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_dynamic_aperture_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_dynamic_aperture_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_dynamic_aperture_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_dynamic_aperture_struct";
  }
};

class TaoDynamicApertureProxy : public FortranProxy<TaoDynamicApertureProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  ApertureParamProxy param() const; // 0D_NOT_type
  void set_param(const ApertureParamProxy& src);
  ApertureScanProxyArray1D scan() const; // 1D_ALLOC_type
  FArray1D<double> pz() const; // 1D_ALLOC_real
  double ellipse_scale() const; // 0D_NOT_real
  void set_ellipse_scale(double value);
  double a_emit() const; // 0D_NOT_real
  void set_a_emit(double value);
  double b_emit() const; // 0D_NOT_real
  void set_b_emit(double value);
};

template <>
struct FortranTraits<TaoModelBranchProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_model_branch_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_model_branch_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_model_branch_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_model_branch_struct";
  }
};

class TaoModelBranchProxy : public FortranProxy<TaoModelBranchProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoModelElementProxyArray1D ele() const; // 1D_ALLOC_type
  TaoBeamBranchProxy beam() const; // 0D_NOT_type
  void set_beam(const TaoBeamBranchProxy& src);
};

template <>
struct FortranTraits<TaoSpinMapProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_spin_map_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_spin_map_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_spin_map_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_spin_map_struct";
  }
};

class TaoSpinMapProxy : public FortranProxy<TaoSpinMapProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool valid() const; // 0D_NOT_logical
  void set_valid(bool value);
  SpinOrbitMap1Proxy map1() const; // 0D_NOT_type
  void set_map1(const SpinOrbitMap1Proxy& src);
  SpinAxisProxy axis_input() const; // 0D_NOT_type
  void set_axis_input(const SpinAxisProxy& src);
  SpinAxisProxy axis0() const; // 0D_NOT_type
  void set_axis0(const SpinAxisProxy& src);
  SpinAxisProxy axis1() const; // 0D_NOT_type
  void set_axis1(const SpinAxisProxy& src);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_ref() const; // 0D_NOT_integer
  void set_ix_ref(int value);
  int ix_uni() const; // 0D_NOT_integer
  void set_ix_uni(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  FArray2D<double> mat8() const; // 2D_NOT_real
};

template <>
struct FortranTraits<TaoDataProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_data_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_data_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_data_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_data_struct";
  }
};

class TaoDataProxy : public FortranProxy<TaoDataProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string ele_name() const; // 0D_NOT_character
  void set_ele_name(const std::string& value);
  std::string ele_start_name() const; // 0D_NOT_character
  void set_ele_start_name(const std::string& value);
  std::string ele_ref_name() const; // 0D_NOT_character
  void set_ele_ref_name(const std::string& value);
  std::string data_type() const; // 0D_ALLOC_character
  void set_data_type(const std::string& value);
  std::string merit_type() const; // 0D_NOT_character
  void set_merit_type(const std::string& value);
  std::string id() const; // 0D_NOT_character
  void set_id(const std::string& value);
  std::string data_source() const; // 0D_NOT_character
  void set_data_source(const std::string& value);
  std::string why_invalid() const; // 0D_NOT_character
  void set_why_invalid(const std::string& value);
  int ix_uni() const; // 0D_NOT_integer
  void set_ix_uni(int value);
  int ix_bunch() const; // 0D_NOT_integer
  void set_ix_bunch(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_ele_start() const; // 0D_NOT_integer
  void set_ix_ele_start(int value);
  int ix_ele_ref() const; // 0D_NOT_integer
  void set_ix_ele_ref(int value);
  int ix_ele_merit() const; // 0D_NOT_integer
  void set_ix_ele_merit(int value);
  int ix_d1() const; // 0D_NOT_integer
  void set_ix_d1(int value);
  int ix_data() const; // 0D_NOT_integer
  void set_ix_data(int value);
  int ix_dModel() const; // 0D_NOT_integer
  void set_ix_dModel(int value);
  int eval_point() const; // 0D_NOT_integer
  void set_eval_point(int value);
  double meas_value() const; // 0D_NOT_real
  void set_meas_value(double value);
  double ref_value() const; // 0D_NOT_real
  void set_ref_value(double value);
  double model_value() const; // 0D_NOT_real
  void set_model_value(double value);
  double design_value() const; // 0D_NOT_real
  void set_design_value(double value);
  double old_value() const; // 0D_NOT_real
  void set_old_value(double value);
  double base_value() const; // 0D_NOT_real
  void set_base_value(double value);
  double error_rms() const; // 0D_NOT_real
  void set_error_rms(double value);
  double delta_merit() const; // 0D_NOT_real
  void set_delta_merit(double value);
  double weight() const; // 0D_NOT_real
  void set_weight(double value);
  double invalid_value() const; // 0D_NOT_real
  void set_invalid_value(double value);
  double merit() const; // 0D_NOT_real
  void set_merit(double value);
  double s() const; // 0D_NOT_real
  void set_s(double value);
  double s_offset() const; // 0D_NOT_real
  void set_s_offset(double value);
  double ref_s_offset() const; // 0D_NOT_real
  void set_ref_s_offset(double value);
  bool err_message_printed() const; // 0D_NOT_logical
  void set_err_message_printed(bool value);
  bool exists() const; // 0D_NOT_logical
  void set_exists(bool value);
  bool good_model() const; // 0D_NOT_logical
  void set_good_model(bool value);
  bool good_base() const; // 0D_NOT_logical
  void set_good_base(bool value);
  bool good_design() const; // 0D_NOT_logical
  void set_good_design(bool value);
  bool good_meas() const; // 0D_NOT_logical
  void set_good_meas(bool value);
  bool good_ref() const; // 0D_NOT_logical
  void set_good_ref(bool value);
  bool good_user() const; // 0D_NOT_logical
  void set_good_user(bool value);
  bool good_opt() const; // 0D_NOT_logical
  void set_good_opt(bool value);
  bool good_plot() const; // 0D_NOT_logical
  void set_good_plot(bool value);
  bool useit_plot() const; // 0D_NOT_logical
  void set_useit_plot(bool value);
  bool useit_opt() const; // 0D_NOT_logical
  void set_useit_opt(bool value);
  TaoSpinMapProxy spin_map() const; // 0D_NOT_type
  void set_spin_map(const TaoSpinMapProxy& src);
  std::optional<TaoD1DataProxy> d1() const; // 0D_PTR_type
  void set_d1(const TaoD1DataProxy& src);
};

template <>
struct FortranTraits<TaoPingScaleProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_ping_scale_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_ping_scale_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_ping_scale_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_ping_scale_struct";
  }
};

class TaoPingScaleProxy : public FortranProxy<TaoPingScaleProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double a_mode_meas() const; // 0D_NOT_real
  void set_a_mode_meas(double value);
  double a_mode_ref() const; // 0D_NOT_real
  void set_a_mode_ref(double value);
  double b_mode_meas() const; // 0D_NOT_real
  void set_b_mode_meas(double value);
  double b_mode_ref() const; // 0D_NOT_real
  void set_b_mode_ref(double value);
};

template <>
struct FortranTraits<TaoUniverseCalcProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_universe_calc_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_universe_calc_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_universe_calc_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_universe_calc_struct";
  }
};

class TaoUniverseCalcProxy : public FortranProxy<TaoUniverseCalcProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int srdt_for_data() const; // 0D_NOT_integer
  void set_srdt_for_data(int value);
  bool rad_int_for_data() const; // 0D_NOT_logical
  void set_rad_int_for_data(bool value);
  bool rad_int_for_plotting() const; // 0D_NOT_logical
  void set_rad_int_for_plotting(bool value);
  bool chrom_for_data() const; // 0D_NOT_logical
  void set_chrom_for_data(bool value);
  bool chrom_for_plotting() const; // 0D_NOT_logical
  void set_chrom_for_plotting(bool value);
  bool lat_sigma_for_data() const; // 0D_NOT_logical
  void set_lat_sigma_for_data(bool value);
  bool lat_sigma_for_plotting() const; // 0D_NOT_logical
  void set_lat_sigma_for_plotting(bool value);
  bool dynamic_aperture() const; // 0D_NOT_logical
  void set_dynamic_aperture(bool value);
  bool one_turn_map() const; // 0D_NOT_logical
  void set_one_turn_map(bool value);
  bool lattice() const; // 0D_NOT_logical
  void set_lattice(bool value);
  bool twiss() const; // 0D_NOT_logical
  void set_twiss(bool value);
  bool track() const; // 0D_NOT_logical
  void set_track(bool value);
  bool spin_matrices() const; // 0D_NOT_logical
  void set_spin_matrices(bool value);
};

template <>
struct FortranTraits<LatEleOrderProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_lat_ele_order_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_lat_ele_order_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_lat_ele_order_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "lat_ele_order_struct";
  }
};

class LatEleOrderProxy : public FortranProxy<LatEleOrderProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  LatEleOrderArrayProxyArray1D branch() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoTitleProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_title_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_title_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_title_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_title_struct";
  }
};

class TaoTitleProxy : public FortranProxy<TaoTitleProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string string() const; // 0D_NOT_character
  void set_string(const std::string& value);
  double x() const; // 0D_NOT_real
  void set_x(double value);
  double y() const; // 0D_NOT_real
  void set_y(double value);
  std::string units() const; // 0D_NOT_character
  void set_units(const std::string& value);
  std::string justify() const; // 0D_NOT_character
  void set_justify(const std::string& value);
  bool draw_it() const; // 0D_NOT_logical
  void set_draw_it(bool value);
};

template <>
struct FortranTraits<QpRectProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_qp_rect_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_qp_rect_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_qp_rect_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "qp_rect_struct";
  }
};

class QpRectProxy : public FortranProxy<QpRectProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double x1() const; // 0D_NOT_real
  void set_x1(double value);
  double x2() const; // 0D_NOT_real
  void set_x2(double value);
  double y1() const; // 0D_NOT_real
  void set_y1(double value);
  double y2() const; // 0D_NOT_real
  void set_y2(double value);
  std::string units() const; // 0D_NOT_character
  void set_units(const std::string& value);
};

template <>
struct FortranTraits<TaoDrawingProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_drawing_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_drawing_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_drawing_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_drawing_struct";
  }
};

class TaoDrawingProxy : public FortranProxy<TaoDrawingProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoEleShapeProxyArray1D ele_shape() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoShapePatternProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_shape_pattern_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_shape_pattern_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_shape_pattern_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_shape_pattern_struct";
  }
};

class TaoShapePatternProxy : public FortranProxy<TaoShapePatternProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  QpLineProxy line() const; // 0D_NOT_type
  void set_line(const QpLineProxy& src);
  TaoShapePatternPointProxyArray1D pt() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoShapePatternPointProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_shape_pattern_point_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_shape_pattern_point_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_shape_pattern_point_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_shape_pattern_point_struct";
  }
};

class TaoShapePatternPointProxy
    : public FortranProxy<TaoShapePatternPointProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double s() const; // 0D_NOT_real
  void set_s(double value);
  double y() const; // 0D_NOT_real
  void set_y(double value);
  double radius() const; // 0D_NOT_real
  void set_radius(double value);
};

template <>
struct FortranTraits<QpAxisProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_qp_axis_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_qp_axis_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_qp_axis_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "qp_axis_struct";
  }
};

class QpAxisProxy : public FortranProxy<QpAxisProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string label() const; // 0D_NOT_character
  void set_label(const std::string& value);
  double min() const; // 0D_NOT_real
  void set_min(double value);
  double max() const; // 0D_NOT_real
  void set_max(double value);
  double tick_min() const; // 0D_NOT_real
  void set_tick_min(double value);
  double tick_max() const; // 0D_NOT_real
  void set_tick_max(double value);
  double eval_min() const; // 0D_NOT_real
  void set_eval_min(double value);
  double eval_max() const; // 0D_NOT_real
  void set_eval_max(double value);
  double dtick() const; // 0D_NOT_real
  void set_dtick(double value);
  double number_offset() const; // 0D_NOT_real
  void set_number_offset(double value);
  double label_offset() const; // 0D_NOT_real
  void set_label_offset(double value);
  double major_tick_len() const; // 0D_NOT_real
  void set_major_tick_len(double value);
  double minor_tick_len() const; // 0D_NOT_real
  void set_minor_tick_len(double value);
  std::string label_color() const; // 0D_NOT_character
  void set_label_color(const std::string& value);
  int major_div() const; // 0D_NOT_integer
  void set_major_div(int value);
  int major_div_nominal() const; // 0D_NOT_integer
  void set_major_div_nominal(int value);
  int minor_div() const; // 0D_NOT_integer
  void set_minor_div(int value);
  int minor_div_max() const; // 0D_NOT_integer
  void set_minor_div_max(int value);
  int places() const; // 0D_NOT_integer
  void set_places(int value);
  std::string type() const; // 0D_NOT_character
  void set_type(const std::string& value);
  std::string bounds() const; // 0D_NOT_character
  void set_bounds(const std::string& value);
  int tick_side() const; // 0D_NOT_integer
  void set_tick_side(int value);
  int number_side() const; // 0D_NOT_integer
  void set_number_side(int value);
  bool draw_label() const; // 0D_NOT_logical
  void set_draw_label(bool value);
  bool draw_numbers() const; // 0D_NOT_logical
  void set_draw_numbers(bool value);
};

template <>
struct FortranTraits<QpLegendProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_qp_legend_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_qp_legend_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_qp_legend_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "qp_legend_struct";
  }
};

class QpLegendProxy : public FortranProxy<QpLegendProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double row_spacing() const; // 0D_NOT_real
  void set_row_spacing(double value);
  double line_length() const; // 0D_NOT_real
  void set_line_length(double value);
  double text_offset() const; // 0D_NOT_real
  void set_text_offset(double value);
  bool draw_line() const; // 0D_NOT_logical
  void set_draw_line(bool value);
  bool draw_symbol() const; // 0D_NOT_logical
  void set_draw_symbol(bool value);
  bool draw_text() const; // 0D_NOT_logical
  void set_draw_text(bool value);
};

template <>
struct FortranTraits<QpPointProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_qp_point_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_qp_point_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_qp_point_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "qp_point_struct";
  }
};

class QpPointProxy : public FortranProxy<QpPointProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double x() const; // 0D_NOT_real
  void set_x(double value);
  double y() const; // 0D_NOT_real
  void set_y(double value);
  std::string units() const; // 0D_NOT_character
  void set_units(const std::string& value);
};

template <>
struct FortranTraits<QpLineProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_qp_line_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_qp_line_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_qp_line_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "qp_line_struct";
  }
};

class QpLineProxy : public FortranProxy<QpLineProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int width() const; // 0D_NOT_integer
  void set_width(int value);
  std::string color() const; // 0D_NOT_character
  void set_color(const std::string& value);
  std::string pattern() const; // 0D_NOT_character
  void set_pattern(const std::string& value);
};

template <>
struct FortranTraits<QpSymbolProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_qp_symbol_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_qp_symbol_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_qp_symbol_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "qp_symbol_struct";
  }
};

class QpSymbolProxy : public FortranProxy<QpSymbolProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string type() const; // 0D_NOT_character
  void set_type(const std::string& value);
  double height() const; // 0D_NOT_real
  void set_height(double value);
  std::string color() const; // 0D_NOT_character
  void set_color(const std::string& value);
  std::string fill_pattern() const; // 0D_NOT_character
  void set_fill_pattern(const std::string& value);
  int line_width() const; // 0D_NOT_integer
  void set_line_width(int value);
};

template <>
struct FortranTraits<TaoFloorPlanProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_floor_plan_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_floor_plan_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_floor_plan_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_floor_plan_struct";
  }
};

class TaoFloorPlanProxy : public FortranProxy<TaoFloorPlanProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string view() const; // 0D_NOT_character
  void set_view(const std::string& value);
  double rotation() const; // 0D_NOT_real
  void set_rotation(double value);
  bool correct_distortion() const; // 0D_NOT_logical
  void set_correct_distortion(bool value);
  bool flip_label_side() const; // 0D_NOT_logical
  void set_flip_label_side(bool value);
  bool size_is_absolute() const; // 0D_NOT_logical
  void set_size_is_absolute(bool value);
  bool draw_only_first_pass() const; // 0D_NOT_logical
  void set_draw_only_first_pass(bool value);
  bool draw_building_wall() const; // 0D_NOT_logical
  void set_draw_building_wall(bool value);
  double orbit_scale() const; // 0D_NOT_real
  void set_orbit_scale(double value);
  std::string orbit_color() const; // 0D_NOT_character
  void set_orbit_color(const std::string& value);
  std::string orbit_pattern() const; // 0D_NOT_character
  void set_orbit_pattern(const std::string& value);
  std::string orbit_lattice() const; // 0D_NOT_character
  void set_orbit_lattice(const std::string& value);
  int orbit_width() const; // 0D_NOT_integer
  void set_orbit_width(int value);
};

template <>
struct FortranTraits<TaoV1VarProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_v1_var_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_v1_var_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_v1_var_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_v1_var_struct";
  }
};

class TaoV1VarProxy : public FortranProxy<TaoV1VarProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  int ix_v1_var() const; // 0D_NOT_integer
  void set_ix_v1_var(int value);
  TaoVarProxyArray1D v() const; // 1D_PTR_type
};

template <>
struct FortranTraits<TaoGlobalProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_global_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_global_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_global_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_global_struct";
  }
};

class TaoGlobalProxy : public FortranProxy<TaoGlobalProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double beam_dead_cutoff() const; // 0D_NOT_real
  void set_beam_dead_cutoff(double value);
  double lm_opt_deriv_reinit() const; // 0D_NOT_real
  void set_lm_opt_deriv_reinit(double value);
  double de_lm_step_ratio() const; // 0D_NOT_real
  void set_de_lm_step_ratio(double value);
  double de_var_to_population_factor() const; // 0D_NOT_real
  void set_de_var_to_population_factor(double value);
  double lmdif_eps() const; // 0D_NOT_real
  void set_lmdif_eps(double value);
  double lmdif_negligible_merit() const; // 0D_NOT_real
  void set_lmdif_negligible_merit(double value);
  double svd_cutoff() const; // 0D_NOT_real
  void set_svd_cutoff(double value);
  double unstable_penalty() const; // 0D_NOT_real
  void set_unstable_penalty(double value);
  double merit_stop_value() const; // 0D_NOT_real
  void set_merit_stop_value(double value);
  double dmerit_stop_value() const; // 0D_NOT_real
  void set_dmerit_stop_value(double value);
  double random_sigma_cutoff() const; // 0D_NOT_real
  void set_random_sigma_cutoff(double value);
  double delta_e_chrom() const; // 0D_NOT_real
  void set_delta_e_chrom(double value);
  double max_plot_time() const; // 0D_NOT_real
  void set_max_plot_time(double value);
  int default_universe() const; // 0D_NOT_integer
  void set_default_universe(int value);
  int default_branch() const; // 0D_NOT_integer
  void set_default_branch(int value);
  int n_opti_cycles() const; // 0D_NOT_integer
  void set_n_opti_cycles(int value);
  int n_opti_loops() const; // 0D_NOT_integer
  void set_n_opti_loops(int value);
  int n_threads() const; // 0D_NOT_integer
  void set_n_threads(int value);
  int phase_units() const; // 0D_NOT_integer
  void set_phase_units(int value);
  int bunch_to_plot() const; // 0D_NOT_integer
  void set_bunch_to_plot(int value);
  int random_seed() const; // 0D_NOT_integer
  void set_random_seed(int value);
  int n_top10_merit() const; // 0D_NOT_integer
  void set_n_top10_merit(int value);
  int srdt_gen_n_slices() const; // 0D_NOT_integer
  void set_srdt_gen_n_slices(int value);
  int datum_err_messages_max() const; // 0D_NOT_integer
  void set_datum_err_messages_max(int value);
  int srdt_sxt_n_slices() const; // 0D_NOT_integer
  void set_srdt_sxt_n_slices(int value);
  bool srdt_use_cache() const; // 0D_NOT_logical
  void set_srdt_use_cache(bool value);
  std::string quiet() const; // 0D_NOT_character
  void set_quiet(const std::string& value);
  std::string random_engine() const; // 0D_NOT_character
  void set_random_engine(const std::string& value);
  std::string random_gauss_converter() const; // 0D_NOT_character
  void set_random_gauss_converter(const std::string& value);
  std::string track_type() const; // 0D_NOT_character
  void set_track_type(const std::string& value);
  std::string lat_sigma_calc_uses_emit_from() const; // 0D_NOT_character
  void set_lat_sigma_calc_uses_emit_from(const std::string& value);
  std::string prompt_string() const; // 0D_NOT_character
  void set_prompt_string(const std::string& value);
  std::string prompt_color() const; // 0D_NOT_character
  void set_prompt_color(const std::string& value);
  std::string optimizer() const; // 0D_NOT_character
  void set_optimizer(const std::string& value);
  std::string print_command() const; // 0D_NOT_character
  void set_print_command(const std::string& value);
  std::string var_out_file() const; // 0D_NOT_character
  void set_var_out_file(const std::string& value);
  std::string history_file() const; // 0D_NOT_character
  void set_history_file(const std::string& value);
  bool beam_timer_on() const; // 0D_NOT_logical
  void set_beam_timer_on(bool value);
  bool box_plots() const; // 0D_NOT_logical
  void set_box_plots(bool value);
  bool blank_line_between_commands() const; // 0D_NOT_logical
  void set_blank_line_between_commands(bool value);
  bool cmd_file_abort_on_error() const; // 0D_NOT_logical
  void set_cmd_file_abort_on_error(bool value);
  bool concatenate_maps() const; // 0D_NOT_logical
  void set_concatenate_maps(bool value);
  bool derivative_recalc() const; // 0D_NOT_logical
  void set_derivative_recalc(bool value);
  bool derivative_uses_design() const; // 0D_NOT_logical
  void set_derivative_uses_design(bool value);
  bool disable_smooth_line_calc() const; // 0D_NOT_logical
  void set_disable_smooth_line_calc(bool value);
  bool draw_curve_off_scale_warn() const; // 0D_NOT_logical
  void set_draw_curve_off_scale_warn(bool value);
  bool external_plotting() const; // 0D_NOT_logical
  void set_external_plotting(bool value);
  bool label_lattice_elements() const; // 0D_NOT_logical
  void set_label_lattice_elements(bool value);
  bool label_keys() const; // 0D_NOT_logical
  void set_label_keys(bool value);
  bool lattice_calc_on() const; // 0D_NOT_logical
  void set_lattice_calc_on(bool value);
  bool only_limit_opt_vars() const; // 0D_NOT_logical
  void set_only_limit_opt_vars(bool value);
  bool opt_with_ref() const; // 0D_NOT_logical
  void set_opt_with_ref(bool value);
  bool opt_with_base() const; // 0D_NOT_logical
  void set_opt_with_base(bool value);
  bool opt_match_auto_recalc() const; // 0D_NOT_logical
  void set_opt_match_auto_recalc(bool value);
  bool opti_write_var_file() const; // 0D_NOT_logical
  void set_opti_write_var_file(bool value);
  bool optimizer_allow_user_abort() const; // 0D_NOT_logical
  void set_optimizer_allow_user_abort(bool value);
  bool optimizer_var_limit_warn() const; // 0D_NOT_logical
  void set_optimizer_var_limit_warn(bool value);
  bool plot_on() const; // 0D_NOT_logical
  void set_plot_on(bool value);
  bool rad_int_user_calc_on() const; // 0D_NOT_logical
  void set_rad_int_user_calc_on(bool value);
  bool rf_on() const; // 0D_NOT_logical
  void set_rf_on(bool value);
  bool single_step() const; // 0D_NOT_logical
  void set_single_step(bool value);
  bool stop_on_error() const; // 0D_NOT_logical
  void set_stop_on_error(bool value);
  bool svd_retreat_on_merit_increase() const; // 0D_NOT_logical
  void set_svd_retreat_on_merit_increase(bool value);
  bool var_limits_on() const; // 0D_NOT_logical
  void set_var_limits_on(bool value);
  bool wait_for_CR_in_single_mode() const; // 0D_NOT_logical
  void set_wait_for_CR_in_single_mode(bool value);
  bool symbol_import() const; // 0D_NOT_logical
  void set_symbol_import(bool value);
  bool debug_on() const; // 0D_NOT_logical
  void set_debug_on(bool value);
  bool expression_tree_on() const; // 0D_NOT_logical
  void set_expression_tree_on(bool value);
  bool verbose_on() const; // 0D_NOT_logical
  void set_verbose_on(bool value);
};

template <>
struct FortranTraits<TaoInitProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_init_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_init_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_init_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_init_struct";
  }
};

class TaoInitProxy : public FortranProxy<TaoInitProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool parse_cmd_args() const; // 0D_NOT_logical
  void set_parse_cmd_args(bool value);
  bool debug_switch() const; // 0D_NOT_logical
  void set_debug_switch(bool value);
  bool external_plotting_switch() const; // 0D_NOT_logical
  void set_external_plotting_switch(bool value);
  std::string init_name() const; // 0D_NOT_character
  void set_init_name(const std::string& value);
  std::string hook_init_file() const; // 0D_NOT_character
  void set_hook_init_file(const std::string& value);
  std::string hook_lat_file() const; // 0D_NOT_character
  void set_hook_lat_file(const std::string& value);
  std::string hook_beam_file() const; // 0D_NOT_character
  void set_hook_beam_file(const std::string& value);
  std::string hook_data_file() const; // 0D_NOT_character
  void set_hook_data_file(const std::string& value);
  std::string hook_plot_file() const; // 0D_NOT_character
  void set_hook_plot_file(const std::string& value);
  std::string hook_startup_file() const; // 0D_NOT_character
  void set_hook_startup_file(const std::string& value);
  std::string hook_var_file() const; // 0D_NOT_character
  void set_hook_var_file(const std::string& value);
  std::string hook_building_wall_file() const; // 0D_NOT_character
  void set_hook_building_wall_file(const std::string& value);
  std::string init_file_arg_path() const; // 0D_NOT_character
  void set_init_file_arg_path(const std::string& value);
  std::string lattice_file_arg() const; // 0D_NOT_character
  void set_lattice_file_arg(const std::string& value);
  std::string hook_init_file_arg() const; // 0D_NOT_character
  void set_hook_init_file_arg(const std::string& value);
  std::string init_file_arg() const; // 0D_NOT_character
  void set_init_file_arg(const std::string& value);
  std::string beam_file_arg() const; // 0D_NOT_character
  void set_beam_file_arg(const std::string& value);
  std::string beam_init_position_file_arg() const; // 0D_NOT_character
  void set_beam_init_position_file_arg(const std::string& value);
  std::string command_arg() const; // 0D_NOT_character
  void set_command_arg(const std::string& value);
  std::string data_file_arg() const; // 0D_NOT_character
  void set_data_file_arg(const std::string& value);
  std::string plot_file_arg() const; // 0D_NOT_character
  void set_plot_file_arg(const std::string& value);
  std::string startup_file_arg() const; // 0D_NOT_character
  void set_startup_file_arg(const std::string& value);
  std::string var_file_arg() const; // 0D_NOT_character
  void set_var_file_arg(const std::string& value);
  std::string building_wall_file_arg() const; // 0D_NOT_character
  void set_building_wall_file_arg(const std::string& value);
  std::string geometry_arg() const; // 0D_NOT_character
  void set_geometry_arg(const std::string& value);
  std::string slice_lattice_arg() const; // 0D_NOT_character
  void set_slice_lattice_arg(const std::string& value);
  std::string start_branch_at_arg() const; // 0D_NOT_character
  void set_start_branch_at_arg(const std::string& value);
  std::string log_startup_arg() const; // 0D_NOT_character
  void set_log_startup_arg(const std::string& value);
  std::string no_stopping_arg() const; // 0D_NOT_character
  void set_no_stopping_arg(const std::string& value);
  std::string noplot_arg() const; // 0D_NOT_character
  void set_noplot_arg(const std::string& value);
  std::string no_rad_int_arg() const; // 0D_NOT_character
  void set_no_rad_int_arg(const std::string& value);
  std::string reverse_arg() const; // 0D_NOT_character
  void set_reverse_arg(const std::string& value);
  std::string debug_arg() const; // 0D_NOT_character
  void set_debug_arg(const std::string& value);
  std::string disable_smooth_line_calc_arg() const; // 0D_NOT_character
  void set_disable_smooth_line_calc_arg(const std::string& value);
  std::string rf_on_arg() const; // 0D_NOT_character
  void set_rf_on_arg(const std::string& value);
  std::string prompt_color_arg() const; // 0D_NOT_character
  void set_prompt_color_arg(const std::string& value);
  std::string quiet_arg() const; // 0D_NOT_character
  void set_quiet_arg(const std::string& value);
  std::string noinit_arg() const; // 0D_NOT_character
  void set_noinit_arg(const std::string& value);
  std::string nostartup_arg() const; // 0D_NOT_character
  void set_nostartup_arg(const std::string& value);
  std::string symbol_import_arg() const; // 0D_NOT_character
  void set_symbol_import_arg(const std::string& value);
  std::string unique_name_suffix() const; // 0D_NOT_character
  void set_unique_name_suffix(const std::string& value);
};

template <>
struct FortranTraits<TaoCommonProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_common_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_common_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_common_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_common_struct";
  }
};

class TaoCommonProxy : public FortranProxy<TaoCommonProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoPlotRegionProxyArray1D plot_place_buffer() const; // 1D_ALLOC_type
  FArray2D<double> covar() const; // 2D_ALLOC_real
  FArray2D<double> alpha() const; // 2D_ALLOC_real
  double dummy_target() const; // 0D_NOT_real
  void set_dummy_target(double value);
  int n_alias() const; // 0D_NOT_integer
  void set_n_alias(int value);
  int cmd_file_level() const; // 0D_NOT_integer
  void set_cmd_file_level(int value);
  int ix_key_bank() const; // 0D_NOT_integer
  void set_ix_key_bank(int value);
  int ix_history() const; // 0D_NOT_integer
  void set_ix_history(int value);
  int n_history() const; // 0D_NOT_integer
  void set_n_history(int value);
  int lev_loop() const; // 0D_NOT_integer
  void set_lev_loop(int value);
  int n_err_messages_printed() const; // 0D_NOT_integer
  void set_n_err_messages_printed(int value);
  int n_universes() const; // 0D_NOT_integer
  void set_n_universes(int value);
  int ix_beam_track_active_element() const; // 0D_NOT_integer
  void set_ix_beam_track_active_element(int value);
  bool cmd_file_paused() const; // 0D_NOT_logical
  void set_cmd_file_paused(bool value);
  bool use_cmd_here() const; // 0D_NOT_logical
  void set_use_cmd_here(bool value);
  bool cmd_from_cmd_file() const; // 0D_NOT_logical
  void set_cmd_from_cmd_file(bool value);
  bool use_saved_beam_in_tracking() const; // 0D_NOT_logical
  void set_use_saved_beam_in_tracking(bool value);
  bool single_mode() const; // 0D_NOT_logical
  void set_single_mode(bool value);
  bool combine_consecutive_elements_of_like_name() const; // 0D_NOT_logical
  void set_combine_consecutive_elements_of_like_name(bool value);
  bool have_tracked_beam() const; // 0D_NOT_logical
  void set_have_tracked_beam(bool value);
  bool init_plot_needed() const; // 0D_NOT_logical
  void set_init_plot_needed(bool value);
  bool init_beam() const; // 0D_NOT_logical
  void set_init_beam(bool value);
  bool init_var() const; // 0D_NOT_logical
  void set_init_var(bool value);
  bool init_read_lat_info() const; // 0D_NOT_logical
  void set_init_read_lat_info(bool value);
  bool optimizer_running() const; // 0D_NOT_logical
  void set_optimizer_running(bool value);
  bool have_datums_using_expressions() const; // 0D_NOT_logical
  void set_have_datums_using_expressions(bool value);
  bool print_to_terminal() const; // 0D_NOT_logical
  void set_print_to_terminal(bool value);
  bool lattice_calc_done() const; // 0D_NOT_logical
  void set_lattice_calc_done(bool value);
  bool add_measurement_noise() const; // 0D_NOT_logical
  void set_add_measurement_noise(bool value);
  bool command_arg_has_been_executed() const; // 0D_NOT_logical
  void set_command_arg_has_been_executed(bool value);
  bool all_merit_weights_positive() const; // 0D_NOT_logical
  void set_all_merit_weights_positive(bool value);
  bool multi_turn_orbit_is_plotted() const; // 0D_NOT_logical
  void set_multi_turn_orbit_is_plotted(bool value);
  bool force_chrom_calc() const; // 0D_NOT_logical
  void set_force_chrom_calc(bool value);
  bool force_rad_int_calc() const; // 0D_NOT_logical
  void set_force_rad_int_calc(bool value);
  bool rad_int_ri_calc_on() const; // 0D_NOT_logical
  void set_rad_int_ri_calc_on(bool value);
  bool rad_int_6d_calc_on() const; // 0D_NOT_logical
  void set_rad_int_6d_calc_on(bool value);
  FCharArray1D valid_plot_who() const; // 1D_NOT_character
  std::string single_mode_buffer() const; // 0D_NOT_character
  void set_single_mode_buffer(const std::string& value);
  std::string cmd() const; // 0D_NOT_character
  void set_cmd(const std::string& value);
  std::string saved_cmd_line() const; // 0D_NOT_character
  void set_saved_cmd_line(const std::string& value);
};

template <>
struct FortranTraits<TaoPlotPageProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_plot_page_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_plot_page_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_plot_page_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_plot_page_struct";
  }
};

class TaoPlotPageProxy : public FortranProxy<TaoPlotPageProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoTitleProxy title() const; // 0D_NOT_type
  void set_title(const TaoTitleProxy& src);
  TaoTitleProxy subtitle() const; // 0D_NOT_type
  void set_subtitle(const TaoTitleProxy& src);
  QpRectProxy border() const; // 0D_NOT_type
  void set_border(const QpRectProxy& src);
  TaoDrawingProxy floor_plan() const; // 0D_NOT_type
  void set_floor_plan(const TaoDrawingProxy& src);
  TaoDrawingProxy lat_layout() const; // 0D_NOT_type
  void set_lat_layout(const TaoDrawingProxy& src);
  TaoShapePatternProxyArray1D pattern() const; // 1D_ALLOC_type
  TaoPlotProxyArray1D template_() const; // 1D_ALLOC_type
  TaoPlotRegionProxyArray1D region() const; // 1D_ALLOC_type
  std::string plot_display_type() const; // 0D_NOT_character
  void set_plot_display_type(const std::string& value);
  FArray1D<double> size() const; // 1D_NOT_real
  double text_height() const; // 0D_NOT_real
  void set_text_height(double value);
  double main_title_text_scale() const; // 0D_NOT_real
  void set_main_title_text_scale(double value);
  double graph_title_text_scale() const; // 0D_NOT_real
  void set_graph_title_text_scale(double value);
  double axis_number_text_scale() const; // 0D_NOT_real
  void set_axis_number_text_scale(double value);
  double axis_label_text_scale() const; // 0D_NOT_real
  void set_axis_label_text_scale(double value);
  double legend_text_scale() const; // 0D_NOT_real
  void set_legend_text_scale(double value);
  double key_table_text_scale() const; // 0D_NOT_real
  void set_key_table_text_scale(double value);
  double floor_plan_shape_scale() const; // 0D_NOT_real
  void set_floor_plan_shape_scale(double value);
  double floor_plan_text_scale() const; // 0D_NOT_real
  void set_floor_plan_text_scale(double value);
  double lat_layout_shape_scale() const; // 0D_NOT_real
  void set_lat_layout_shape_scale(double value);
  double lat_layout_text_scale() const; // 0D_NOT_real
  void set_lat_layout_text_scale(double value);
  int n_curve_pts() const; // 0D_NOT_integer
  void set_n_curve_pts(int value);
  int id_window() const; // 0D_NOT_integer
  void set_id_window(int value);
  bool delete_overlapping_plots() const; // 0D_NOT_logical
  void set_delete_overlapping_plots(bool value);
  bool draw_graph_title_suffix() const; // 0D_NOT_logical
  void set_draw_graph_title_suffix(bool value);
};

template <>
struct FortranTraits<TaoBuildingWallProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_building_wall_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_building_wall_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_building_wall_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_building_wall_struct";
  }
};

class TaoBuildingWallProxy : public FortranProxy<TaoBuildingWallProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoBuildingWallOrientationProxy orientation() const; // 0D_NOT_type
  void set_orientation(const TaoBuildingWallOrientationProxy& src);
  TaoBuildingWallSectionProxyArray1D section() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoBuildingWallOrientationProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_building_wall_orientation_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_building_wall_orientation_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_building_wall_orientation_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_building_wall_orientation_struct";
  }
};

class TaoBuildingWallOrientationProxy
    : public FortranProxy<TaoBuildingWallOrientationProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double theta() const; // 0D_NOT_real
  void set_theta(double value);
  double x_offset() const; // 0D_NOT_real
  void set_x_offset(double value);
  double z_offset() const; // 0D_NOT_real
  void set_z_offset(double value);
};

template <>
struct FortranTraits<TaoBuildingWallSectionProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_building_wall_section_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_building_wall_section_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_building_wall_section_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_building_wall_section_struct";
  }
};

class TaoBuildingWallSectionProxy
    : public FortranProxy<TaoBuildingWallSectionProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string& value);
  std::string constraint() const; // 0D_NOT_character
  void set_constraint(const std::string& value);
  TaoBuildingWallPointProxyArray1D point() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoBuildingWallPointProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_building_wall_point_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_building_wall_point_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_building_wall_point_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_building_wall_point_struct";
  }
};

class TaoBuildingWallPointProxy
    : public FortranProxy<TaoBuildingWallPointProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double z() const; // 0D_NOT_real
  void set_z(double value);
  double x() const; // 0D_NOT_real
  void set_x(double value);
  double radius() const; // 0D_NOT_real
  void set_radius(double value);
  double z_center() const; // 0D_NOT_real
  void set_z_center(double value);
  double x_center() const; // 0D_NOT_real
  void set_x_center(double value);
};

template <>
struct FortranTraits<TaoWaveProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_wave_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_wave_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_wave_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_wave_struct";
  }
};

class TaoWaveProxy : public FortranProxy<TaoWaveProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string data_type() const; // 0D_NOT_character
  void set_data_type(const std::string& value);
  double rms_rel_a() const; // 0D_NOT_real
  void set_rms_rel_a(double value);
  double rms_rel_b() const; // 0D_NOT_real
  void set_rms_rel_b(double value);
  double rms_rel_as() const; // 0D_NOT_real
  void set_rms_rel_as(double value);
  double rms_rel_bs() const; // 0D_NOT_real
  void set_rms_rel_bs(double value);
  double rms_rel_ar() const; // 0D_NOT_real
  void set_rms_rel_ar(double value);
  double rms_rel_br() const; // 0D_NOT_real
  void set_rms_rel_br(double value);
  double rms_rel_k() const; // 0D_NOT_real
  void set_rms_rel_k(double value);
  double rms_rel_ks() const; // 0D_NOT_real
  void set_rms_rel_ks(double value);
  double rms_rel_kr() const; // 0D_NOT_real
  void set_rms_rel_kr(double value);
  double rms_phi() const; // 0D_NOT_real
  void set_rms_phi(double value);
  double rms_phi_s() const; // 0D_NOT_real
  void set_rms_phi_s(double value);
  double rms_phi_r() const; // 0D_NOT_real
  void set_rms_phi_r(double value);
  double amp_ba_s() const; // 0D_NOT_real
  void set_amp_ba_s(double value);
  double amp_ba_r() const; // 0D_NOT_real
  void set_amp_ba_r(double value);
  double chi_a() const; // 0D_NOT_real
  void set_chi_a(double value);
  double chi_c() const; // 0D_NOT_real
  void set_chi_c(double value);
  double chi_ba() const; // 0D_NOT_real
  void set_chi_ba(double value);
  FArray1D<double> amp_a() const; // 1D_NOT_real
  FArray1D<double> amp_b() const; // 1D_NOT_real
  FArray1D<double> amp_ba() const; // 1D_NOT_real
  FArray1D<double> coef_a() const; // 1D_NOT_real
  FArray1D<double> coef_b() const; // 1D_NOT_real
  FArray1D<double> coef_ba() const; // 1D_NOT_real
  int n_func() const; // 0D_NOT_integer
  void set_n_func(int value);
  int ix_a1() const; // 0D_NOT_integer
  void set_ix_a1(int value);
  int ix_a2() const; // 0D_NOT_integer
  void set_ix_a2(int value);
  int ix_b1() const; // 0D_NOT_integer
  void set_ix_b1(int value);
  int ix_b2() const; // 0D_NOT_integer
  void set_ix_b2(int value);
  int i_a1() const; // 0D_NOT_integer
  void set_i_a1(int value);
  int i_a2() const; // 0D_NOT_integer
  void set_i_a2(int value);
  int i_b1() const; // 0D_NOT_integer
  void set_i_b1(int value);
  int i_b2() const; // 0D_NOT_integer
  void set_i_b2(int value);
  int n_a() const; // 0D_NOT_integer
  void set_n_a(int value);
  int n_b() const; // 0D_NOT_integer
  void set_n_b(int value);
  int i_curve_wrap_pt() const; // 0D_NOT_integer
  void set_i_curve_wrap_pt(int value);
  FArray1D<int> ix_data() const; // 1D_ALLOC_integer
  int n_kick() const; // 0D_NOT_integer
  void set_n_kick(int value);
  TaoWaveKickPtProxyArray1D kick() const; // 1D_ALLOC_type
  TaoGraphProxy base_graph() const; // 0D_NOT_type
  void set_base_graph(const TaoGraphProxy& src);
  std::optional<TaoPlotRegionProxy> region() const; // 0D_PTR_type
  void set_region(const TaoPlotRegionProxy& src);
  std::optional<TaoD1DataProxy> d1_dat() const; // 0D_PTR_type
  void set_d1_dat(const TaoD1DataProxy& src);
};

template <>
struct FortranTraits<TaoWaveKickPtProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_wave_kick_pt_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_wave_kick_pt_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_wave_kick_pt_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_wave_kick_pt_struct";
  }
};

class TaoWaveKickPtProxy : public FortranProxy<TaoWaveKickPtProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double phi_s() const; // 0D_NOT_real
  void set_phi_s(double value);
  double phi_r() const; // 0D_NOT_real
  void set_phi_r(double value);
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
  double amp() const; // 0D_NOT_real
  void set_amp(double value);
  double s() const; // 0D_NOT_real
  void set_s(double value);
  int ix_dat_before_kick() const; // 0D_NOT_integer
  void set_ix_dat_before_kick(int value);
  std::optional<EleProxy> ele() const; // 0D_PTR_type
  void set_ele(const EleProxy& src);
};

template <>
struct FortranTraits<TaoCmdHistoryProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_cmd_history_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_cmd_history_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_cmd_history_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_cmd_history_struct";
  }
};

class TaoCmdHistoryProxy : public FortranProxy<TaoCmdHistoryProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string cmd() const; // 0D_ALLOC_character
  void set_cmd(const std::string& value);
  int ix() const; // 0D_NOT_integer
  void set_ix(int value);
};

template <>
struct FortranTraits<TaoUniverseProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_tao_universe_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_tao_universe_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_tao_universe_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "tao_universe_struct";
  }
};

class TaoUniverseProxy : public FortranProxy<TaoUniverseProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<TaoLatticeProxy> model() const; // 0D_PTR_type
  void set_model(const TaoLatticeProxy& src);
  std::optional<TaoLatticeProxy> design() const; // 0D_PTR_type
  void set_design(const TaoLatticeProxy& src);
  std::optional<TaoLatticeProxy> base() const; // 0D_PTR_type
  void set_base(const TaoLatticeProxy& src);
  TaoBeamUniProxy beam() const; // 0D_NOT_type
  void set_beam(const TaoBeamUniProxy& src);
  TaoDynamicApertureProxy dynamic_aperture() const; // 0D_NOT_type
  void set_dynamic_aperture(const TaoDynamicApertureProxy& src);
  TaoModelBranchProxyArray1D model_branch() const; // 1D_PTR_type
  TaoD2DataProxyArray1D d2_data() const; // 1D_ALLOC_type
  TaoDataProxyArray1D data() const; // 1D_ALLOC_type
  TaoPingScaleProxy ping_scale() const; // 0D_NOT_type
  void set_ping_scale(const TaoPingScaleProxy& src);
  LatProxy scratch_lat() const; // 0D_NOT_type
  void set_scratch_lat(const LatProxy& src);
  TaoUniverseCalcProxy calc() const; // 0D_NOT_type
  void set_calc(const TaoUniverseCalcProxy& src);
  LatEleOrderProxy ele_order() const; // 0D_NOT_type
  void set_ele_order(const LatEleOrderProxy& src);
  TaoSpinMapProxy spin_map() const; // 0D_NOT_type
  void set_spin_map(const TaoSpinMapProxy& src);
  FArray2D<double> dModel_dVar() const; // 2D_ALLOC_real
  int ix_uni() const; // 0D_NOT_integer
  void set_ix_uni(int value);
  int n_d2_data_used() const; // 0D_NOT_integer
  void set_n_d2_data_used(int value);
  int n_data_used() const; // 0D_NOT_integer
  void set_n_data_used(int value);
  bool is_on() const; // 0D_NOT_logical
  void set_is_on(bool value);
  bool design_same_as_previous() const; // 0D_NOT_logical
  void set_design_same_as_previous(bool value);
  bool picked_uni() const; // 0D_NOT_logical
  void set_picked_uni(bool value);
};

template <>
struct FortranTraits<MadEnergyProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_mad_energy_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_mad_energy_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_mad_energy_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "mad_energy_struct";
  }
};

class MadEnergyProxy : public FortranProxy<MadEnergyProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double total() const; // 0D_NOT_real
  void set_total(double value);
  double beta() const; // 0D_NOT_real
  void set_beta(double value);
  double gamma() const; // 0D_NOT_real
  void set_gamma(double value);
  double kinetic() const; // 0D_NOT_real
  void set_kinetic(double value);
  double p0c() const; // 0D_NOT_real
  void set_p0c(double value);
  int particle() const; // 0D_NOT_integer
  void set_particle(int value);
};

template <>
struct FortranTraits<MadMapProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_mad_map_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_mad_map_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_mad_map_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "mad_map_struct";
  }
};

class MadMapProxy : public FortranProxy<MadMapProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> k() const; // 1D_NOT_real
  FArray2D<double> r() const; // 2D_NOT_real
  FArray3D<double> t() const; // 3D_NOT_real
};

template <>
struct FortranTraits<RandomStateProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_random_state_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_random_state_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_random_state_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "random_state_struct";
  }
};

class RandomStateProxy : public FortranProxy<RandomStateProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int64_t ix() const; // 0D_NOT_integer8
  void set_ix(int64_t value);
  int64_t iy() const; // 0D_NOT_integer8
  void set_iy(int64_t value);
  bool number_stored() const; // 0D_NOT_logical
  void set_number_stored(bool value);
  double h_saved() const; // 0D_NOT_real
  void set_h_saved(double value);
  int engine() const; // 0D_NOT_integer
  void set_engine(int value);
  int seed() const; // 0D_NOT_integer
  void set_seed(int value);
  double am() const; // 0D_NOT_real
  void set_am(double value);
  int gauss_converter() const; // 0D_NOT_integer
  void set_gauss_converter(int value);
  double gauss_sigma_cut() const; // 0D_NOT_real
  void set_gauss_sigma_cut(double value);
  int64_t in_sobseq() const; // 0D_NOT_integer8
  void set_in_sobseq(int64_t value);
  FArray1D<double> x_sobseq() const; // 1D_NOT_real
};

template <>
struct FortranTraits<BbuStageProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bbu_stage_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bbu_stage_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bbu_stage_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bbu_stage_struct";
  }
};

class BbuStageProxy : public FortranProxy<BbuStageProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int ix_ele_lr_wake() const; // 0D_NOT_integer
  void set_ix_ele_lr_wake(int value);
  int ix_ele_stage_end() const; // 0D_NOT_integer
  void set_ix_ele_stage_end(int value);
  int ix_pass() const; // 0D_NOT_integer
  void set_ix_pass(int value);
  int ix_stage_pass1() const; // 0D_NOT_integer
  void set_ix_stage_pass1(int value);
  int ix_head_bunch() const; // 0D_NOT_integer
  void set_ix_head_bunch(int value);
  int ix_hom_max() const; // 0D_NOT_integer
  void set_ix_hom_max(int value);
  double hom_voltage_max() const; // 0D_NOT_real
  void set_hom_voltage_max(double value);
  double time_at_wake_ele() const; // 0D_NOT_real
  void set_time_at_wake_ele(double value);
  FArray1D<double> ave_orb() const; // 1D_NOT_real
  FArray1D<double> rms_orb() const; // 1D_NOT_real
  FArray1D<double> min_orb() const; // 1D_NOT_real
  FArray1D<double> max_orb() const; // 1D_NOT_real
  int n_orb() const; // 0D_NOT_integer
  void set_n_orb(int value);
};

template <>
struct FortranTraits<BbuBeamProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bbu_beam_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bbu_beam_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bbu_beam_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bbu_beam_struct";
  }
};

class BbuBeamProxy : public FortranProxy<BbuBeamProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  BunchProxyArray1D bunch() const; // 1D_ALLOC_type
  BbuStageProxyArray1D stage() const; // 1D_ALLOC_type
  FArray1D<int> ix_ele_bunch() const; // 1D_ALLOC_integer
  int ix_bunch_head() const; // 0D_NOT_integer
  void set_ix_bunch_head(int value);
  int ix_bunch_end() const; // 0D_NOT_integer
  void set_ix_bunch_end(int value);
  int n_bunch_in_lat() const; // 0D_NOT_integer
  void set_n_bunch_in_lat(int value);
  int ix_stage_voltage_max() const; // 0D_NOT_integer
  void set_ix_stage_voltage_max(int value);
  double hom_voltage_max() const; // 0D_NOT_real
  void set_hom_voltage_max(double value);
  double time_now() const; // 0D_NOT_real
  void set_time_now(double value);
  double one_turn_time() const; // 0D_NOT_real
  void set_one_turn_time(double value);
  double rf_wavelength_max() const; // 0D_NOT_real
  void set_rf_wavelength_max(double value);
};

template <>
struct FortranTraits<BbuParamProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_bbu_param_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_bbu_param_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_bbu_param_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "bbu_param_struct";
  }
};

class BbuParamProxy : public FortranProxy<BbuParamProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string lat_filename() const; // 0D_NOT_character
  void set_lat_filename(const std::string& value);
  std::string lat2_filename() const; // 0D_NOT_character
  void set_lat2_filename(const std::string& value);
  std::string bunch_by_bunch_info_file() const; // 0D_NOT_character
  void set_bunch_by_bunch_info_file(const std::string& value);
  bool hybridize() const; // 0D_NOT_logical
  void set_hybridize(bool value);
  bool write_digested_hybrid_lat() const; // 0D_NOT_logical
  void set_write_digested_hybrid_lat(bool value);
  bool write_voltage_vs_time_dat() const; // 0D_NOT_logical
  void set_write_voltage_vs_time_dat(bool value);
  bool keep_overlays_and_groups() const; // 0D_NOT_logical
  void set_keep_overlays_and_groups(bool value);
  bool keep_all_lcavities() const; // 0D_NOT_logical
  void set_keep_all_lcavities(bool value);
  bool use_taylor_for_hybrids() const; // 0D_NOT_logical
  void set_use_taylor_for_hybrids(bool value);
  bool stable_orbit_anal() const; // 0D_NOT_logical
  void set_stable_orbit_anal(bool value);
  double limit_factor() const; // 0D_NOT_real
  void set_limit_factor(double value);
  double simulation_turns_max() const; // 0D_NOT_real
  void set_simulation_turns_max(double value);
  double bunch_freq() const; // 0D_NOT_real
  void set_bunch_freq(double value);
  double init_particle_offset() const; // 0D_NOT_real
  void set_init_particle_offset(double value);
  double current() const; // 0D_NOT_real
  void set_current(double value);
  double rel_tol() const; // 0D_NOT_real
  void set_rel_tol(double value);
  bool drscan() const; // 0D_NOT_logical
  void set_drscan(bool value);
  bool use_interpolated_threshold() const; // 0D_NOT_logical
  void set_use_interpolated_threshold(bool value);
  bool write_hom_info() const; // 0D_NOT_logical
  void set_write_hom_info(bool value);
  int elindex() const; // 0D_NOT_integer
  void set_elindex(int value);
  std::string elname() const; // 0D_NOT_character
  void set_elname(const std::string& value);
  int nstep() const; // 0D_NOT_integer
  void set_nstep(int value);
  double begdr() const; // 0D_NOT_real
  void set_begdr(double value);
  double enddr() const; // 0D_NOT_real
  void set_enddr(double value);
  int nrep() const; // 0D_NOT_integer
  void set_nrep(int value);
  int ran_seed() const; // 0D_NOT_integer
  void set_ran_seed(int value);
  int hom_order_cutoff() const; // 0D_NOT_integer
  void set_hom_order_cutoff(int value);
  double ran_gauss_sigma_cut() const; // 0D_NOT_real
  void set_ran_gauss_sigma_cut(double value);
  std::string ele_track_end() const; // 0D_NOT_character
  void set_ele_track_end(const std::string& value);
  int ix_ele_track_end() const; // 0D_NOT_integer
  void set_ix_ele_track_end(int value);
  bool regression() const; // 0D_NOT_logical
  void set_regression(bool value);
  bool normalize_z_to_rf() const; // 0D_NOT_logical
  void set_normalize_z_to_rf(bool value);
  bool ramp_on() const; // 0D_NOT_logical
  void set_ramp_on(bool value);
  FArray1D<double> ramp_pattern() const; // 1D_NOT_real
  int ramp_n_start() const; // 0D_NOT_integer
  void set_ramp_n_start(int value);
  int n_ramp_pattern() const; // 0D_NOT_integer
  void set_n_ramp_pattern(int value);
};

template <>
struct FortranTraits<AllEncompassingProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_all_encompassing_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_all_encompassing_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_all_encompassing_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "all_encompassing_struct";
  }
};

class AllEncompassingProxy : public FortranProxy<AllEncompassingProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double real_rp_0d() const; // 0D_NOT_real
  void set_real_rp_0d(double value);
  FArray1D<double> real_rp_1d() const; // 1D_NOT_real
  FArray2D<double> real_rp_2d() const; // 2D_NOT_real
  FArray3D<double> real_rp_3d() const; // 3D_NOT_real
  double* real_rp_0d_ptr() const; // 0D_PTR_real
  void set_real_rp_0d_ptr(double value);
  FArray1D<double> real_rp_1d_ptr() const; // 1D_PTR_real
  FArray2D<double> real_rp_2d_ptr() const; // 2D_PTR_real
  FArray3D<double> real_rp_3d_ptr() const; // 3D_PTR_real
  FArray1D<double> real_rp_1d_alloc() const; // 1D_ALLOC_real
  FArray2D<double> real_rp_2d_alloc() const; // 2D_ALLOC_real
  FArray3D<double> real_rp_3d_alloc() const; // 3D_ALLOC_real
  double real_dp_0d() const; // 0D_NOT_real
  void set_real_dp_0d(double value);
  FArray1D<double> real_dp_1d() const; // 1D_NOT_real
  FArray2D<double> real_dp_2d() const; // 2D_NOT_real
  FArray3D<double> real_dp_3d() const; // 3D_NOT_real
  double* real_dp_0d_ptr() const; // 0D_PTR_real
  void set_real_dp_0d_ptr(double value);
  FArray1D<double> real_dp_1d_ptr() const; // 1D_PTR_real
  FArray2D<double> real_dp_2d_ptr() const; // 2D_PTR_real
  FArray3D<double> real_dp_3d_ptr() const; // 3D_PTR_real
  FArray1D<double> real_dp_1d_alloc() const; // 1D_ALLOC_real
  FArray2D<double> real_dp_2d_alloc() const; // 2D_ALLOC_real
  FArray3D<double> real_dp_3d_alloc() const; // 3D_ALLOC_real
  std::complex<double> complex_dp_0d() const; // 0D_NOT_complex
  void set_complex_dp_0d(std::complex<double> value);
  FArray1D<std::complex<double>> complex_dp_1d() const; // 1D_NOT_complex
  FArray2D<std::complex<double>> complex_dp_2d() const; // 2D_NOT_complex
  FArray3D<std::complex<double>> complex_dp_3d() const; // 3D_NOT_complex
  FArray1D<std::complex<double>> complex_dp_1d_ptr() const; // 1D_PTR_complex
  FArray2D<std::complex<double>> complex_dp_2d_ptr() const; // 2D_PTR_complex
  FArray3D<std::complex<double>> complex_dp_3d_ptr() const; // 3D_PTR_complex
  FArray1D<std::complex<double>> complex_dp_1d_alloc()
      const; // 1D_ALLOC_complex
  FArray2D<std::complex<double>> complex_dp_2d_alloc()
      const; // 2D_ALLOC_complex
  FArray3D<std::complex<double>> complex_dp_3d_alloc()
      const; // 3D_ALLOC_complex
  int int_0d() const; // 0D_NOT_integer
  void set_int_0d(int value);
  FArray1D<int> int_1d() const; // 1D_NOT_integer
  FArray2D<int> int_2d() const; // 2D_NOT_integer
  FArray3D<int> int_3d() const; // 3D_NOT_integer
  int* int_0d_ptr() const; // 0D_PTR_integer
  void set_int_0d_ptr(int value);
  FArray1D<int> int_1d_ptr() const; // 1D_PTR_integer
  FArray2D<int> int_2d_ptr() const; // 2D_PTR_integer
  FArray3D<int> int_3d_ptr() const; // 3D_PTR_integer
  FArray1D<int> int_1d_alloc() const; // 1D_ALLOC_integer
  FArray2D<int> int_2d_alloc() const; // 2D_ALLOC_integer
  FArray3D<int> int_3d_alloc() const; // 3D_ALLOC_integer
  int64_t int8_0d() const; // 0D_NOT_integer8
  void set_int8_0d(int64_t value);
  int64_t* int8_0d_ptr() const; // 0D_PTR_integer8
  void set_int8_0d_ptr(int64_t value);
  bool logical_0d() const; // 0D_NOT_logical
  void set_logical_0d(bool value);
  bool* logical_0d_ptr() const; // 0D_PTR_logical
  void set_logical_0d_ptr(bool value);
  TestSubProxy type_0d() const; // 0D_NOT_type
  void set_type_0d(const TestSubProxy& src);
  TestSubProxyArray1D type_1d() const; // 1D_NOT_type
  TestSubProxyArray2D type_2d() const; // 2D_NOT_type
  TestSubProxyArray3D type_3d() const; // 3D_NOT_type
  std::optional<TestSubProxy> type_0d_ptr() const; // 0D_PTR_type
  void set_type_0d_ptr(const TestSubProxy& src);
  TestSubProxyArray1D type_1d_ptr() const; // 1D_PTR_type
  TestSubProxyArray2D type_2d_ptr() const; // 2D_PTR_type
  TestSubProxyArray3D type_3d_ptr() const; // 3D_PTR_type
  TestSubProxyArray1D type_1d_alloc() const; // 1D_ALLOC_type
  TestSubProxyArray2D type_2d_alloc() const; // 2D_ALLOC_type
  TestSubProxyArray3D type_3d_alloc() const; // 3D_ALLOC_type
};

template <>
struct FortranTraits<TestSubProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_test_sub_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_test_sub_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_test_sub_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "test_sub_struct";
  }
};

class TestSubProxy : public FortranProxy<TestSubProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TestSubSubProxy sr() const; // 0D_NOT_type
  void set_sr(const TestSubSubProxy& src);
};

template <>
struct FortranTraits<TestSubSubProxy> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_test_sub_sub_struct(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_test_sub_sub_struct(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_test_sub_sub_struct(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "test_sub_sub_struct";
  }
};

class TestSubSubProxy : public FortranProxy<TestSubSubProxy> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int64_t aaa() const; // 0D_NOT_integer8
  void set_aaa(int64_t value);
  int bbb() const; // 0D_NOT_integer
  void set_bbb(int value);
  std::string file() const; // 0D_NOT_character
  void set_file(const std::string& value);
  double t_ref() const; // 0D_NOT_real
  void set_t_ref(double value);
  double freq_spread() const; // 0D_NOT_real
  void set_freq_spread(double value);
};

} // namespace Bmad
