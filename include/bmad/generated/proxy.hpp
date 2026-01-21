#pragma once

#include <complex>
#include <memory>
#include <string>

#include "bmad/fortran_arrays.hpp"
#include "bmad/proxy_base.hpp"

extern "C" {
// Forward declarations for Fortran interface
void spline_struct_get_x0(const void *struct_obj, double *value_out);
void spline_struct_set_x0(void *struct_obj, double value_in);
void spline_struct_get_y0(const void *struct_obj, double *value_out);
void spline_struct_set_y0(void *struct_obj, double value_in);
void spline_struct_get_x1(const void *struct_obj, double *value_out);
void spline_struct_set_x1(void *struct_obj, double value_in);
void spline_struct_get_coef_info(const void *s, double **d, int *bounds, bool *is_alloc);
void spline_struct_set_coef(void *s, const void *d, const int *shape);
void spin_polar_struct_get_polarization(const void *struct_obj, double *value_out);
void spin_polar_struct_set_polarization(void *struct_obj, double value_in);
void spin_polar_struct_get_theta(const void *struct_obj, double *value_out);
void spin_polar_struct_set_theta(void *struct_obj, double value_in);
void spin_polar_struct_get_phi(const void *struct_obj, double *value_out);
void spin_polar_struct_set_phi(void *struct_obj, double value_in);
void spin_polar_struct_get_xi(const void *struct_obj, double *value_out);
void spin_polar_struct_set_xi(void *struct_obj, double value_in);
void ac_kicker_time_struct_get_amp(const void *struct_obj, double *value_out);
void ac_kicker_time_struct_set_amp(void *struct_obj, double value_in);
void ac_kicker_time_struct_get_time(const void *struct_obj, double *value_out);
void ac_kicker_time_struct_set_time(void *struct_obj, double value_in);
void ac_kicker_time_struct_get_spline(const void *struct_obj, void **ptr_out);
void ac_kicker_time_struct_set_spline(void *struct_obj, const void *src_ptr);
void ac_kicker_freq_struct_get_f(const void *struct_obj, double *value_out);
void ac_kicker_freq_struct_set_f(void *struct_obj, double value_in);
void ac_kicker_freq_struct_get_amp(const void *struct_obj, double *value_out);
void ac_kicker_freq_struct_set_amp(void *struct_obj, double value_in);
void ac_kicker_freq_struct_get_phi(const void *struct_obj, double *value_out);
void ac_kicker_freq_struct_set_phi(void *struct_obj, double value_in);

void ac_kicker_struct_get_amp_vs_time_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void ac_kicker_struct_get_frequency_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void interval1_coef_struct_get_c0(const void *struct_obj, double *value_out);
void interval1_coef_struct_set_c0(void *struct_obj, double value_in);
void interval1_coef_struct_get_c1(const void *struct_obj, double *value_out);
void interval1_coef_struct_set_c1(void *struct_obj, double value_in);
void interval1_coef_struct_get_n_exp(const void *struct_obj, double *value_out);
void interval1_coef_struct_set_n_exp(void *struct_obj, double value_in);
void photon_reflect_table_struct_get_angle_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void photon_reflect_table_struct_set_angle(void *s, const void *d, const int *shape);
void photon_reflect_table_struct_get_energy_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void photon_reflect_table_struct_set_energy(void *s, const void *d, const int *shape);

void photon_reflect_table_struct_get_int1_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void photon_reflect_table_struct_get_p_reflect_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void photon_reflect_table_struct_set_p_reflect(void *s, const void *d, const int *shape);
void photon_reflect_table_struct_get_max_energy(const void *struct_obj, double *value_out);
void photon_reflect_table_struct_set_max_energy(void *struct_obj, double value_in);
void photon_reflect_table_struct_get_p_reflect_scratch_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void photon_reflect_table_struct_set_p_reflect_scratch(void *s, const void *d, const int *shape);
void photon_reflect_table_struct_get_bragg_angle_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void photon_reflect_table_struct_set_bragg_angle(void *s, const void *d, const int *shape);
void photon_reflect_surface_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void photon_reflect_surface_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void photon_reflect_surface_struct_get_description_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void photon_reflect_surface_struct_set_description(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);
void photon_reflect_surface_struct_get_reflectivity_file_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void photon_reflect_surface_struct_set_reflectivity_file(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);

void photon_reflect_surface_struct_get_table_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void photon_reflect_surface_struct_get_surface_roughness_rms(
    const void *struct_obj,
    double *value_out
);
void photon_reflect_surface_struct_set_surface_roughness_rms(void *struct_obj, double value_in);
void photon_reflect_surface_struct_get_roughness_correlation_len(
    const void *struct_obj,
    double *value_out
);
void photon_reflect_surface_struct_set_roughness_correlation_len(void *struct_obj, double value_in);
void photon_reflect_surface_struct_get_ix_surface(const void *struct_obj, int *value_out);
void photon_reflect_surface_struct_set_ix_surface(void *struct_obj, int value_in);
void coord_struct_get_vec_info(const void *s, double **d, int *bounds, bool *is_alloc);
void coord_struct_set_vec(void *s, const void *d, const int *shape);
void coord_struct_get_s(const void *struct_obj, double *value_out);
void coord_struct_set_s(void *struct_obj, double value_in);
void coord_struct_get_t(const void *struct_obj, long double *value_out);
void coord_struct_set_t(void *struct_obj, long double value_in);
void coord_struct_get_spin_info(const void *s, double **d, int *bounds, bool *is_alloc);
void coord_struct_set_spin(void *s, const void *d, const int *shape);
void coord_struct_get_field_info(const void *s, double **d, int *bounds, bool *is_alloc);
void coord_struct_set_field(void *s, const void *d, const int *shape);
void coord_struct_get_phase_info(const void *s, double **d, int *bounds, bool *is_alloc);
void coord_struct_set_phase(void *s, const void *d, const int *shape);
void coord_struct_get_charge(const void *struct_obj, double *value_out);
void coord_struct_set_charge(void *struct_obj, double value_in);
void coord_struct_get_dt_ref(const void *struct_obj, double *value_out);
void coord_struct_set_dt_ref(void *struct_obj, double value_in);
void coord_struct_get_r(const void *struct_obj, double *value_out);
void coord_struct_set_r(void *struct_obj, double value_in);
void coord_struct_get_p0c(const void *struct_obj, double *value_out);
void coord_struct_set_p0c(void *struct_obj, double value_in);
void coord_struct_get_E_potential(const void *struct_obj, double *value_out);
void coord_struct_set_E_potential(void *struct_obj, double value_in);
void coord_struct_get_beta(const void *struct_obj, double *value_out);
void coord_struct_set_beta(void *struct_obj, double value_in);
void coord_struct_get_ix_ele(const void *struct_obj, int *value_out);
void coord_struct_set_ix_ele(void *struct_obj, int value_in);
void coord_struct_get_ix_branch(const void *struct_obj, int *value_out);
void coord_struct_set_ix_branch(void *struct_obj, int value_in);
void coord_struct_get_ix_turn(const void *struct_obj, int *value_out);
void coord_struct_set_ix_turn(void *struct_obj, int value_in);
void coord_struct_get_ix_user(const void *struct_obj, int *value_out);
void coord_struct_set_ix_user(void *struct_obj, int value_in);
void coord_struct_get_state(const void *struct_obj, int *value_out);
void coord_struct_set_state(void *struct_obj, int value_in);
void coord_struct_get_direction(const void *struct_obj, int *value_out);
void coord_struct_set_direction(void *struct_obj, int value_in);
void coord_struct_get_time_dir(const void *struct_obj, int *value_out);
void coord_struct_set_time_dir(void *struct_obj, int value_in);
void coord_struct_get_species(const void *struct_obj, int *value_out);
void coord_struct_set_species(void *struct_obj, int value_in);
void coord_struct_get_location(const void *struct_obj, int *value_out);
void coord_struct_set_location(void *struct_obj, int value_in);

void coord_array_struct_get_orbit_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bpm_phase_coupling_struct_get_K_22a(const void *struct_obj, double *value_out);
void bpm_phase_coupling_struct_set_K_22a(void *struct_obj, double value_in);
void bpm_phase_coupling_struct_get_K_12a(const void *struct_obj, double *value_out);
void bpm_phase_coupling_struct_set_K_12a(void *struct_obj, double value_in);
void bpm_phase_coupling_struct_get_K_11b(const void *struct_obj, double *value_out);
void bpm_phase_coupling_struct_set_K_11b(void *struct_obj, double value_in);
void bpm_phase_coupling_struct_get_K_12b(const void *struct_obj, double *value_out);
void bpm_phase_coupling_struct_set_K_12b(void *struct_obj, double value_in);
void bpm_phase_coupling_struct_get_Cbar22_a(const void *struct_obj, double *value_out);
void bpm_phase_coupling_struct_set_Cbar22_a(void *struct_obj, double value_in);
void bpm_phase_coupling_struct_get_Cbar12_a(const void *struct_obj, double *value_out);
void bpm_phase_coupling_struct_set_Cbar12_a(void *struct_obj, double value_in);
void bpm_phase_coupling_struct_get_Cbar11_b(const void *struct_obj, double *value_out);
void bpm_phase_coupling_struct_set_Cbar11_b(void *struct_obj, double value_in);
void bpm_phase_coupling_struct_get_Cbar12_b(const void *struct_obj, double *value_out);
void bpm_phase_coupling_struct_set_Cbar12_b(void *struct_obj, double value_in);
void bpm_phase_coupling_struct_get_phi_a(const void *struct_obj, double *value_out);
void bpm_phase_coupling_struct_set_phi_a(void *struct_obj, double value_in);
void bpm_phase_coupling_struct_get_phi_b(const void *struct_obj, double *value_out);
void bpm_phase_coupling_struct_set_phi_b(void *struct_obj, double value_in);
void expression_atom_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void expression_atom_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void expression_atom_struct_get_type(const void *struct_obj, int *value_out);
void expression_atom_struct_set_type(void *struct_obj, int value_in);
void expression_atom_struct_get_value(const void *struct_obj, double *value_out);
void expression_atom_struct_set_value(void *struct_obj, double value_in);
void wake_sr_z_long_struct_get_w_info(const void *s, double **d, int *bounds, bool *is_alloc);
void wake_sr_z_long_struct_set_w(void *s, const void *d, const int *shape);
void wake_sr_z_long_struct_get_fw_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    bool *is_alloc
);
void wake_sr_z_long_struct_set_fw(void *s, const void *d, const int *shape);
void wake_sr_z_long_struct_get_fbunch_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    bool *is_alloc
);
void wake_sr_z_long_struct_set_fbunch(void *s, const void *d, const int *shape);
void wake_sr_z_long_struct_get_w_out_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    bool *is_alloc
);
void wake_sr_z_long_struct_set_w_out(void *s, const void *d, const int *shape);
void wake_sr_z_long_struct_get_dz(const void *struct_obj, double *value_out);
void wake_sr_z_long_struct_set_dz(void *struct_obj, double value_in);
void wake_sr_z_long_struct_get_z0(const void *struct_obj, double *value_out);
void wake_sr_z_long_struct_set_z0(void *struct_obj, double value_in);
void wake_sr_z_long_struct_get_smoothing_sigma(const void *struct_obj, double *value_out);
void wake_sr_z_long_struct_set_smoothing_sigma(void *struct_obj, double value_in);
void wake_sr_z_long_struct_get_position_dependence(const void *struct_obj, int *value_out);
void wake_sr_z_long_struct_set_position_dependence(void *struct_obj, int value_in);
void wake_sr_z_long_struct_get_time_based(const void *struct_obj, bool *value_out);
void wake_sr_z_long_struct_set_time_based(void *struct_obj, bool value_in);
void wake_sr_mode_struct_get_amp(const void *struct_obj, double *value_out);
void wake_sr_mode_struct_set_amp(void *struct_obj, double value_in);
void wake_sr_mode_struct_get_damp(const void *struct_obj, double *value_out);
void wake_sr_mode_struct_set_damp(void *struct_obj, double value_in);
void wake_sr_mode_struct_get_k(const void *struct_obj, double *value_out);
void wake_sr_mode_struct_set_k(void *struct_obj, double value_in);
void wake_sr_mode_struct_get_phi(const void *struct_obj, double *value_out);
void wake_sr_mode_struct_set_phi(void *struct_obj, double value_in);
void wake_sr_mode_struct_get_b_sin(const void *struct_obj, double *value_out);
void wake_sr_mode_struct_set_b_sin(void *struct_obj, double value_in);
void wake_sr_mode_struct_get_b_cos(const void *struct_obj, double *value_out);
void wake_sr_mode_struct_set_b_cos(void *struct_obj, double value_in);
void wake_sr_mode_struct_get_a_sin(const void *struct_obj, double *value_out);
void wake_sr_mode_struct_set_a_sin(void *struct_obj, double value_in);
void wake_sr_mode_struct_get_a_cos(const void *struct_obj, double *value_out);
void wake_sr_mode_struct_set_a_cos(void *struct_obj, double value_in);
void wake_sr_mode_struct_get_polarization(const void *struct_obj, int *value_out);
void wake_sr_mode_struct_set_polarization(void *struct_obj, int value_in);
void wake_sr_mode_struct_get_position_dependence(const void *struct_obj, int *value_out);
void wake_sr_mode_struct_set_position_dependence(void *struct_obj, int value_in);
void wake_sr_struct_get_file_info(const void *s, char **d, int *bounds, bool *a);
void wake_sr_struct_set_file(void *struct_obj, const char *str_ptr, int str_len);
void wake_sr_struct_get_z_long(const void *struct_obj, void **ptr_out);
void wake_sr_struct_set_z_long(void *struct_obj, const void *src_ptr);

void wake_sr_struct_get_long_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void wake_sr_struct_get_trans_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void wake_sr_struct_get_z_ref_long(const void *struct_obj, double *value_out);
void wake_sr_struct_set_z_ref_long(void *struct_obj, double value_in);
void wake_sr_struct_get_z_ref_trans(const void *struct_obj, double *value_out);
void wake_sr_struct_set_z_ref_trans(void *struct_obj, double value_in);
void wake_sr_struct_get_z_max(const void *struct_obj, double *value_out);
void wake_sr_struct_set_z_max(void *struct_obj, double value_in);
void wake_sr_struct_get_amp_scale(const void *struct_obj, double *value_out);
void wake_sr_struct_set_amp_scale(void *struct_obj, double value_in);
void wake_sr_struct_get_z_scale(const void *struct_obj, double *value_out);
void wake_sr_struct_set_z_scale(void *struct_obj, double value_in);
void wake_sr_struct_get_scale_with_length(const void *struct_obj, bool *value_out);
void wake_sr_struct_set_scale_with_length(void *struct_obj, bool value_in);
void wake_lr_mode_struct_get_freq(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_freq(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_freq_in(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_freq_in(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_R_over_Q(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_R_over_Q(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_Q(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_Q(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_damp(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_damp(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_phi(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_phi(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_angle(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_angle(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_b_sin(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_b_sin(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_b_cos(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_b_cos(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_a_sin(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_a_sin(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_a_cos(const void *struct_obj, double *value_out);
void wake_lr_mode_struct_set_a_cos(void *struct_obj, double value_in);
void wake_lr_mode_struct_get_m(const void *struct_obj, int *value_out);
void wake_lr_mode_struct_set_m(void *struct_obj, int value_in);
void wake_lr_mode_struct_get_polarized(const void *struct_obj, bool *value_out);
void wake_lr_mode_struct_set_polarized(void *struct_obj, bool value_in);
void wake_lr_struct_get_file_info(const void *s, char **d, int *bounds, bool *a);
void wake_lr_struct_set_file(void *struct_obj, const char *str_ptr, int str_len);

void wake_lr_struct_get_mode_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void wake_lr_struct_get_t_ref(const void *struct_obj, double *value_out);
void wake_lr_struct_set_t_ref(void *struct_obj, double value_in);
void wake_lr_struct_get_freq_spread(const void *struct_obj, double *value_out);
void wake_lr_struct_set_freq_spread(void *struct_obj, double value_in);
void wake_lr_struct_get_amp_scale(const void *struct_obj, double *value_out);
void wake_lr_struct_set_amp_scale(void *struct_obj, double value_in);
void wake_lr_struct_get_time_scale(const void *struct_obj, double *value_out);
void wake_lr_struct_set_time_scale(void *struct_obj, double value_in);
void wake_lr_struct_get_self_wake_on(const void *struct_obj, bool *value_out);
void wake_lr_struct_set_self_wake_on(void *struct_obj, bool value_in);
void lat_ele_loc_struct_get_ix_ele(const void *struct_obj, int *value_out);
void lat_ele_loc_struct_set_ix_ele(void *struct_obj, int value_in);
void lat_ele_loc_struct_get_ix_branch(const void *struct_obj, int *value_out);
void lat_ele_loc_struct_set_ix_branch(void *struct_obj, int value_in);
void wake_struct_get_sr(const void *struct_obj, void **ptr_out);
void wake_struct_set_sr(void *struct_obj, const void *src_ptr);
void wake_struct_get_lr(const void *struct_obj, void **ptr_out);
void wake_struct_set_lr(void *struct_obj, const void *src_ptr);
void taylor_term_struct_get_coef(const void *struct_obj, double *value_out);
void taylor_term_struct_set_coef(void *struct_obj, double value_in);
void taylor_term_struct_get_expn_info(const void *s, int **d, int *bounds, bool *is_alloc);
void taylor_term_struct_set_expn(void *s, const void *d, const int *shape);
void taylor_struct_get_ref(const void *struct_obj, double *value_out);
void taylor_struct_set_ref(void *struct_obj, double value_in);

void taylor_struct_get_term_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void em_taylor_term_struct_get_coef(const void *struct_obj, double *value_out);
void em_taylor_term_struct_set_coef(void *struct_obj, double value_in);
void em_taylor_term_struct_get_expn_info(const void *s, int **d, int *bounds, bool *is_alloc);
void em_taylor_term_struct_set_expn(void *s, const void *d, const int *shape);
void em_taylor_struct_get_ref(const void *struct_obj, double *value_out);
void em_taylor_struct_set_ref(void *struct_obj, double value_in);

void em_taylor_struct_get_term_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void cartesian_map_term1_struct_get_coef(const void *struct_obj, double *value_out);
void cartesian_map_term1_struct_set_coef(void *struct_obj, double value_in);
void cartesian_map_term1_struct_get_kx(const void *struct_obj, double *value_out);
void cartesian_map_term1_struct_set_kx(void *struct_obj, double value_in);
void cartesian_map_term1_struct_get_ky(const void *struct_obj, double *value_out);
void cartesian_map_term1_struct_set_ky(void *struct_obj, double value_in);
void cartesian_map_term1_struct_get_kz(const void *struct_obj, double *value_out);
void cartesian_map_term1_struct_set_kz(void *struct_obj, double value_in);
void cartesian_map_term1_struct_get_x0(const void *struct_obj, double *value_out);
void cartesian_map_term1_struct_set_x0(void *struct_obj, double value_in);
void cartesian_map_term1_struct_get_y0(const void *struct_obj, double *value_out);
void cartesian_map_term1_struct_set_y0(void *struct_obj, double value_in);
void cartesian_map_term1_struct_get_phi_z(const void *struct_obj, double *value_out);
void cartesian_map_term1_struct_set_phi_z(void *struct_obj, double value_in);
void cartesian_map_term1_struct_get_family(const void *struct_obj, int *value_out);
void cartesian_map_term1_struct_set_family(void *struct_obj, int value_in);
void cartesian_map_term1_struct_get_form(const void *struct_obj, int *value_out);
void cartesian_map_term1_struct_set_form(void *struct_obj, int value_in);
void cartesian_map_term_struct_get_file_info(const void *s, char **d, int *bounds, bool *a);
void cartesian_map_term_struct_set_file(void *struct_obj, const char *str_ptr, int str_len);
void cartesian_map_term_struct_get_n_link(const void *struct_obj, int *value_out);
void cartesian_map_term_struct_set_n_link(void *struct_obj, int value_in);

void cartesian_map_term_struct_get_term_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void cartesian_map_struct_get_field_scale(const void *struct_obj, double *value_out);
void cartesian_map_struct_set_field_scale(void *struct_obj, double value_in);
void cartesian_map_struct_get_r0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void cartesian_map_struct_set_r0(void *s, const void *d, const int *shape);
void cartesian_map_struct_get_master_parameter(const void *struct_obj, int *value_out);
void cartesian_map_struct_set_master_parameter(void *struct_obj, int value_in);
void cartesian_map_struct_get_ele_anchor_pt(const void *struct_obj, int *value_out);
void cartesian_map_struct_set_ele_anchor_pt(void *struct_obj, int value_in);
void cartesian_map_struct_get_field_type(const void *struct_obj, int *value_out);
void cartesian_map_struct_set_field_type(void *struct_obj, int value_in);
void cartesian_map_struct_get_ptr(const void *struct_obj, void **ptr_out);
void cartesian_map_struct_set_ptr(void *struct_obj, const void *src_ptr);
void cylindrical_map_term1_struct_get_e_coef(
    const void *struct_obj,
    std::complex<double> *value_out
);
void cylindrical_map_term1_struct_set_e_coef(void *struct_obj, std::complex<double> value_in);
void cylindrical_map_term1_struct_get_b_coef(
    const void *struct_obj,
    std::complex<double> *value_out
);
void cylindrical_map_term1_struct_set_b_coef(void *struct_obj, std::complex<double> value_in);
void cylindrical_map_term_struct_get_file_info(const void *s, char **d, int *bounds, bool *a);
void cylindrical_map_term_struct_set_file(void *struct_obj, const char *str_ptr, int str_len);
void cylindrical_map_term_struct_get_n_link(const void *struct_obj, int *value_out);
void cylindrical_map_term_struct_set_n_link(void *struct_obj, int value_in);

void cylindrical_map_term_struct_get_term_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void cylindrical_map_struct_get_m(const void *struct_obj, int *value_out);
void cylindrical_map_struct_set_m(void *struct_obj, int value_in);
void cylindrical_map_struct_get_harmonic(const void *struct_obj, int *value_out);
void cylindrical_map_struct_set_harmonic(void *struct_obj, int value_in);
void cylindrical_map_struct_get_phi0_fieldmap(const void *struct_obj, double *value_out);
void cylindrical_map_struct_set_phi0_fieldmap(void *struct_obj, double value_in);
void cylindrical_map_struct_get_theta0_azimuth(const void *struct_obj, double *value_out);
void cylindrical_map_struct_set_theta0_azimuth(void *struct_obj, double value_in);
void cylindrical_map_struct_get_field_scale(const void *struct_obj, double *value_out);
void cylindrical_map_struct_set_field_scale(void *struct_obj, double value_in);
void cylindrical_map_struct_get_master_parameter(const void *struct_obj, int *value_out);
void cylindrical_map_struct_set_master_parameter(void *struct_obj, int value_in);
void cylindrical_map_struct_get_ele_anchor_pt(const void *struct_obj, int *value_out);
void cylindrical_map_struct_set_ele_anchor_pt(void *struct_obj, int value_in);
void cylindrical_map_struct_get_dz(const void *struct_obj, double *value_out);
void cylindrical_map_struct_set_dz(void *struct_obj, double value_in);
void cylindrical_map_struct_get_r0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void cylindrical_map_struct_set_r0(void *s, const void *d, const int *shape);
void cylindrical_map_struct_get_ptr(const void *struct_obj, void **ptr_out);
void cylindrical_map_struct_set_ptr(void *struct_obj, const void *src_ptr);
void bicubic_cmplx_coef_struct_get_coef_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void bicubic_cmplx_coef_struct_set_coef(void *s, const void *d, const int *shape);
void bicubic_cmplx_coef_struct_get_i_box_info(const void *s, int **d, int *bounds, bool *is_alloc);
void bicubic_cmplx_coef_struct_set_i_box(void *s, const void *d, const int *shape);
void tricubic_cmplx_coef_struct_get_coef_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void tricubic_cmplx_coef_struct_set_coef(void *s, const void *d, const int *shape);
void tricubic_cmplx_coef_struct_get_i_box_info(const void *s, int **d, int *bounds, bool *is_alloc);
void tricubic_cmplx_coef_struct_set_i_box(void *s, const void *d, const int *shape);
void grid_field_pt1_struct_get_E_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    bool *is_alloc
);
void grid_field_pt1_struct_set_E(void *s, const void *d, const int *shape);
void grid_field_pt1_struct_get_B_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    bool *is_alloc
);
void grid_field_pt1_struct_set_B(void *s, const void *d, const int *shape);
void grid_field_pt_struct_get_file_info(const void *s, char **d, int *bounds, bool *a);
void grid_field_pt_struct_set_file(void *struct_obj, const char *str_ptr, int str_len);
void grid_field_pt_struct_get_n_link(const void *struct_obj, int *value_out);
void grid_field_pt_struct_set_n_link(void *struct_obj, int value_in);

void grid_field_pt_struct_get_pt_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void grid_field_struct_get_geometry(const void *struct_obj, int *value_out);
void grid_field_struct_set_geometry(void *struct_obj, int value_in);
void grid_field_struct_get_harmonic(const void *struct_obj, int *value_out);
void grid_field_struct_set_harmonic(void *struct_obj, int value_in);
void grid_field_struct_get_phi0_fieldmap(const void *struct_obj, double *value_out);
void grid_field_struct_set_phi0_fieldmap(void *struct_obj, double value_in);
void grid_field_struct_get_field_scale(const void *struct_obj, double *value_out);
void grid_field_struct_set_field_scale(void *struct_obj, double value_in);
void grid_field_struct_get_field_type(const void *struct_obj, int *value_out);
void grid_field_struct_set_field_type(void *struct_obj, int value_in);
void grid_field_struct_get_master_parameter(const void *struct_obj, int *value_out);
void grid_field_struct_set_master_parameter(void *struct_obj, int value_in);
void grid_field_struct_get_ele_anchor_pt(const void *struct_obj, int *value_out);
void grid_field_struct_set_ele_anchor_pt(void *struct_obj, int value_in);
void grid_field_struct_get_interpolation_order(const void *struct_obj, int *value_out);
void grid_field_struct_set_interpolation_order(void *struct_obj, int value_in);
void grid_field_struct_get_dr_info(const void *s, double **d, int *bounds, bool *is_alloc);
void grid_field_struct_set_dr(void *s, const void *d, const int *shape);
void grid_field_struct_get_r0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void grid_field_struct_set_r0(void *s, const void *d, const int *shape);
void grid_field_struct_get_curved_ref_frame(const void *struct_obj, bool *value_out);
void grid_field_struct_set_curved_ref_frame(void *struct_obj, bool value_in);
void grid_field_struct_get_ptr(const void *struct_obj, void **ptr_out);
void grid_field_struct_set_ptr(void *struct_obj, const void *src_ptr);

void grid_field_struct_get_bi_coef_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void grid_field_struct_get_tri_coef_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void floor_position_struct_get_r_info(const void *s, double **d, int *bounds, bool *is_alloc);
void floor_position_struct_set_r(void *s, const void *d, const int *shape);
void floor_position_struct_get_w_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void floor_position_struct_set_w(void *s, const void *d, const int *shape);
void floor_position_struct_get_theta(const void *struct_obj, double *value_out);
void floor_position_struct_set_theta(void *struct_obj, double value_in);
void floor_position_struct_get_phi(const void *struct_obj, double *value_out);
void floor_position_struct_set_phi(void *struct_obj, double value_in);
void floor_position_struct_get_psi(const void *struct_obj, double *value_out);
void floor_position_struct_set_psi(void *struct_obj, double value_in);
void high_energy_space_charge_struct_get_closed_orb(const void *struct_obj, void **ptr_out);
void high_energy_space_charge_struct_set_closed_orb(void *struct_obj, const void *src_ptr);
void high_energy_space_charge_struct_get_kick_const(const void *struct_obj, double *value_out);
void high_energy_space_charge_struct_set_kick_const(void *struct_obj, double value_in);
void high_energy_space_charge_struct_get_sig_x(const void *struct_obj, double *value_out);
void high_energy_space_charge_struct_set_sig_x(void *struct_obj, double value_in);
void high_energy_space_charge_struct_get_sig_y(const void *struct_obj, double *value_out);
void high_energy_space_charge_struct_set_sig_y(void *struct_obj, double value_in);
void high_energy_space_charge_struct_get_phi(const void *struct_obj, double *value_out);
void high_energy_space_charge_struct_set_phi(void *struct_obj, double value_in);
void high_energy_space_charge_struct_get_sin_phi(const void *struct_obj, double *value_out);
void high_energy_space_charge_struct_set_sin_phi(void *struct_obj, double value_in);
void high_energy_space_charge_struct_get_cos_phi(const void *struct_obj, double *value_out);
void high_energy_space_charge_struct_set_cos_phi(void *struct_obj, double value_in);
void high_energy_space_charge_struct_get_sig_z(const void *struct_obj, double *value_out);
void high_energy_space_charge_struct_set_sig_z(void *struct_obj, double value_in);
void xy_disp_struct_get_eta(const void *struct_obj, double *value_out);
void xy_disp_struct_set_eta(void *struct_obj, double value_in);
void xy_disp_struct_get_etap(const void *struct_obj, double *value_out);
void xy_disp_struct_set_etap(void *struct_obj, double value_in);
void xy_disp_struct_get_deta_ds(const void *struct_obj, double *value_out);
void xy_disp_struct_set_deta_ds(void *struct_obj, double value_in);
void xy_disp_struct_get_sigma(const void *struct_obj, double *value_out);
void xy_disp_struct_set_sigma(void *struct_obj, double value_in);
void xy_disp_struct_get_deta_dpz(const void *struct_obj, double *value_out);
void xy_disp_struct_set_deta_dpz(void *struct_obj, double value_in);
void xy_disp_struct_get_detap_dpz(const void *struct_obj, double *value_out);
void xy_disp_struct_set_detap_dpz(void *struct_obj, double value_in);
void twiss_struct_get_beta(const void *struct_obj, double *value_out);
void twiss_struct_set_beta(void *struct_obj, double value_in);
void twiss_struct_get_alpha(const void *struct_obj, double *value_out);
void twiss_struct_set_alpha(void *struct_obj, double value_in);
void twiss_struct_get_gamma(const void *struct_obj, double *value_out);
void twiss_struct_set_gamma(void *struct_obj, double value_in);
void twiss_struct_get_phi(const void *struct_obj, double *value_out);
void twiss_struct_set_phi(void *struct_obj, double value_in);
void twiss_struct_get_eta(const void *struct_obj, double *value_out);
void twiss_struct_set_eta(void *struct_obj, double value_in);
void twiss_struct_get_etap(const void *struct_obj, double *value_out);
void twiss_struct_set_etap(void *struct_obj, double value_in);
void twiss_struct_get_deta_ds(const void *struct_obj, double *value_out);
void twiss_struct_set_deta_ds(void *struct_obj, double value_in);
void twiss_struct_get_sigma(const void *struct_obj, double *value_out);
void twiss_struct_set_sigma(void *struct_obj, double value_in);
void twiss_struct_get_sigma_p(const void *struct_obj, double *value_out);
void twiss_struct_set_sigma_p(void *struct_obj, double value_in);
void twiss_struct_get_emit(const void *struct_obj, double *value_out);
void twiss_struct_set_emit(void *struct_obj, double value_in);
void twiss_struct_get_norm_emit(const void *struct_obj, double *value_out);
void twiss_struct_set_norm_emit(void *struct_obj, double value_in);
void twiss_struct_get_chrom(const void *struct_obj, double *value_out);
void twiss_struct_set_chrom(void *struct_obj, double value_in);
void twiss_struct_get_dbeta_dpz(const void *struct_obj, double *value_out);
void twiss_struct_set_dbeta_dpz(void *struct_obj, double value_in);
void twiss_struct_get_dalpha_dpz(const void *struct_obj, double *value_out);
void twiss_struct_set_dalpha_dpz(void *struct_obj, double value_in);
void twiss_struct_get_deta_dpz(const void *struct_obj, double *value_out);
void twiss_struct_set_deta_dpz(void *struct_obj, double value_in);
void twiss_struct_get_detap_dpz(const void *struct_obj, double *value_out);
void twiss_struct_set_detap_dpz(void *struct_obj, double value_in);
void mode3_struct_get_v_info(const void *s, double **d, int *bounds, int *strides, bool *is_alloc);
void mode3_struct_set_v(void *s, const void *d, const int *shape);
void mode3_struct_get_a(const void *struct_obj, void **ptr_out);
void mode3_struct_set_a(void *struct_obj, const void *src_ptr);
void mode3_struct_get_b(const void *struct_obj, void **ptr_out);
void mode3_struct_set_b(void *struct_obj, const void *src_ptr);
void mode3_struct_get_c(const void *struct_obj, void **ptr_out);
void mode3_struct_set_c(void *struct_obj, const void *src_ptr);
void mode3_struct_get_x(const void *struct_obj, void **ptr_out);
void mode3_struct_set_x(void *struct_obj, const void *src_ptr);
void mode3_struct_get_y(const void *struct_obj, void **ptr_out);
void mode3_struct_set_y(void *struct_obj, const void *src_ptr);
void bookkeeping_state_struct_get_attributes(const void *struct_obj, int *value_out);
void bookkeeping_state_struct_set_attributes(void *struct_obj, int value_in);
void bookkeeping_state_struct_get_control(const void *struct_obj, int *value_out);
void bookkeeping_state_struct_set_control(void *struct_obj, int value_in);
void bookkeeping_state_struct_get_floor_position(const void *struct_obj, int *value_out);
void bookkeeping_state_struct_set_floor_position(void *struct_obj, int value_in);
void bookkeeping_state_struct_get_s_position(const void *struct_obj, int *value_out);
void bookkeeping_state_struct_set_s_position(void *struct_obj, int value_in);
void bookkeeping_state_struct_get_ref_energy(const void *struct_obj, int *value_out);
void bookkeeping_state_struct_set_ref_energy(void *struct_obj, int value_in);
void bookkeeping_state_struct_get_mat6(const void *struct_obj, int *value_out);
void bookkeeping_state_struct_set_mat6(void *struct_obj, int value_in);
void bookkeeping_state_struct_get_rad_int(const void *struct_obj, int *value_out);
void bookkeeping_state_struct_set_rad_int(void *struct_obj, int value_in);
void bookkeeping_state_struct_get_ptc(const void *struct_obj, int *value_out);
void bookkeeping_state_struct_set_ptc(void *struct_obj, int value_in);
void bookkeeping_state_struct_get_has_misalign(const void *struct_obj, bool *value_out);
void bookkeeping_state_struct_set_has_misalign(void *struct_obj, bool value_in);
void rad_map_struct_get_ref_orb_info(const void *s, double **d, int *bounds, bool *is_alloc);
void rad_map_struct_set_ref_orb(void *s, const void *d, const int *shape);
void rad_map_struct_get_damp_dmat_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void rad_map_struct_set_damp_dmat(void *s, const void *d, const int *shape);
void rad_map_struct_get_xfer_damp_vec_info(const void *s, double **d, int *bounds, bool *is_alloc);
void rad_map_struct_set_xfer_damp_vec(void *s, const void *d, const int *shape);
void rad_map_struct_get_xfer_damp_mat_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void rad_map_struct_set_xfer_damp_mat(void *s, const void *d, const int *shape);
void rad_map_struct_get_stoc_mat_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void rad_map_struct_set_stoc_mat(void *s, const void *d, const int *shape);
void rad_map_ele_struct_get_rm0(const void *struct_obj, void **ptr_out);
void rad_map_ele_struct_set_rm0(void *struct_obj, const void *src_ptr);
void rad_map_ele_struct_get_rm1(const void *struct_obj, void **ptr_out);
void rad_map_ele_struct_set_rm1(void *struct_obj, const void *src_ptr);
void rad_map_ele_struct_get_stale(const void *struct_obj, bool *value_out);
void rad_map_ele_struct_set_stale(void *struct_obj, bool value_in);
void gen_grad1_struct_get_m(const void *struct_obj, int *value_out);
void gen_grad1_struct_set_m(void *struct_obj, int value_in);
void gen_grad1_struct_get_sincos(const void *struct_obj, int *value_out);
void gen_grad1_struct_set_sincos(void *struct_obj, int value_in);
void gen_grad1_struct_get_n_deriv_max(const void *struct_obj, int *value_out);
void gen_grad1_struct_set_n_deriv_max(void *struct_obj, int value_in);
void gen_grad1_struct_get_deriv_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void gen_grad1_struct_set_deriv(void *s, const void *d, const int *shape);
void gen_grad_map_struct_get_file_info(const void *s, char **d, int *bounds, bool *a);
void gen_grad_map_struct_set_file(void *struct_obj, const char *str_ptr, int str_len);

void gen_grad_map_struct_get_gg_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void gen_grad_map_struct_get_ele_anchor_pt(const void *struct_obj, int *value_out);
void gen_grad_map_struct_set_ele_anchor_pt(void *struct_obj, int value_in);
void gen_grad_map_struct_get_field_type(const void *struct_obj, int *value_out);
void gen_grad_map_struct_set_field_type(void *struct_obj, int value_in);
void gen_grad_map_struct_get_iz0(const void *struct_obj, int *value_out);
void gen_grad_map_struct_set_iz0(void *struct_obj, int value_in);
void gen_grad_map_struct_get_iz1(const void *struct_obj, int *value_out);
void gen_grad_map_struct_set_iz1(void *struct_obj, int value_in);
void gen_grad_map_struct_get_dz(const void *struct_obj, double *value_out);
void gen_grad_map_struct_set_dz(void *struct_obj, double value_in);
void gen_grad_map_struct_get_r0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void gen_grad_map_struct_set_r0(void *s, const void *d, const int *shape);
void gen_grad_map_struct_get_field_scale(const void *struct_obj, double *value_out);
void gen_grad_map_struct_set_field_scale(void *struct_obj, double value_in);
void gen_grad_map_struct_get_master_parameter(const void *struct_obj, int *value_out);
void gen_grad_map_struct_set_master_parameter(void *struct_obj, int value_in);
void gen_grad_map_struct_get_curved_ref_frame(const void *struct_obj, bool *value_out);
void gen_grad_map_struct_set_curved_ref_frame(void *struct_obj, bool value_in);
void surface_segmented_pt_struct_get_x0(const void *struct_obj, double *value_out);
void surface_segmented_pt_struct_set_x0(void *struct_obj, double value_in);
void surface_segmented_pt_struct_get_y0(const void *struct_obj, double *value_out);
void surface_segmented_pt_struct_set_y0(void *struct_obj, double value_in);
void surface_segmented_pt_struct_get_z0(const void *struct_obj, double *value_out);
void surface_segmented_pt_struct_set_z0(void *struct_obj, double value_in);
void surface_segmented_pt_struct_get_dz_dx(const void *struct_obj, double *value_out);
void surface_segmented_pt_struct_set_dz_dx(void *struct_obj, double value_in);
void surface_segmented_pt_struct_get_dz_dy(const void *struct_obj, double *value_out);
void surface_segmented_pt_struct_set_dz_dy(void *struct_obj, double value_in);
void surface_segmented_struct_get_active(const void *struct_obj, bool *value_out);
void surface_segmented_struct_set_active(void *struct_obj, bool value_in);
void surface_segmented_struct_get_dr_info(const void *s, double **d, int *bounds, bool *is_alloc);
void surface_segmented_struct_set_dr(void *s, const void *d, const int *shape);
void surface_segmented_struct_get_r0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void surface_segmented_struct_set_r0(void *s, const void *d, const int *shape);

void surface_segmented_struct_get_pt_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void surface_h_misalign_pt_struct_get_x0(const void *struct_obj, double *value_out);
void surface_h_misalign_pt_struct_set_x0(void *struct_obj, double value_in);
void surface_h_misalign_pt_struct_get_y0(const void *struct_obj, double *value_out);
void surface_h_misalign_pt_struct_set_y0(void *struct_obj, double value_in);
void surface_h_misalign_pt_struct_get_rot_y(const void *struct_obj, double *value_out);
void surface_h_misalign_pt_struct_set_rot_y(void *struct_obj, double value_in);
void surface_h_misalign_pt_struct_get_rot_t(const void *struct_obj, double *value_out);
void surface_h_misalign_pt_struct_set_rot_t(void *struct_obj, double value_in);
void surface_h_misalign_pt_struct_get_rot_y_rms(const void *struct_obj, double *value_out);
void surface_h_misalign_pt_struct_set_rot_y_rms(void *struct_obj, double value_in);
void surface_h_misalign_pt_struct_get_rot_t_rms(const void *struct_obj, double *value_out);
void surface_h_misalign_pt_struct_set_rot_t_rms(void *struct_obj, double value_in);
void surface_h_misalign_struct_get_active(const void *struct_obj, bool *value_out);
void surface_h_misalign_struct_set_active(void *struct_obj, bool value_in);
void surface_h_misalign_struct_get_dr_info(const void *s, double **d, int *bounds, bool *is_alloc);
void surface_h_misalign_struct_set_dr(void *s, const void *d, const int *shape);
void surface_h_misalign_struct_get_r0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void surface_h_misalign_struct_set_r0(void *s, const void *d, const int *shape);

void surface_h_misalign_struct_get_pt_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void surface_displacement_pt_struct_get_x0(const void *struct_obj, double *value_out);
void surface_displacement_pt_struct_set_x0(void *struct_obj, double value_in);
void surface_displacement_pt_struct_get_y0(const void *struct_obj, double *value_out);
void surface_displacement_pt_struct_set_y0(void *struct_obj, double value_in);
void surface_displacement_pt_struct_get_z0(const void *struct_obj, double *value_out);
void surface_displacement_pt_struct_set_z0(void *struct_obj, double value_in);
void surface_displacement_pt_struct_get_dz_dx(const void *struct_obj, double *value_out);
void surface_displacement_pt_struct_set_dz_dx(void *struct_obj, double value_in);
void surface_displacement_pt_struct_get_dz_dy(const void *struct_obj, double *value_out);
void surface_displacement_pt_struct_set_dz_dy(void *struct_obj, double value_in);
void surface_displacement_pt_struct_get_d2z_dxdy(const void *struct_obj, double *value_out);
void surface_displacement_pt_struct_set_d2z_dxdy(void *struct_obj, double value_in);
void surface_displacement_struct_get_active(const void *struct_obj, bool *value_out);
void surface_displacement_struct_set_active(void *struct_obj, bool value_in);
void surface_displacement_struct_get_dr_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void surface_displacement_struct_set_dr(void *s, const void *d, const int *shape);
void surface_displacement_struct_get_r0_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void surface_displacement_struct_set_r0(void *s, const void *d, const int *shape);

void surface_displacement_struct_get_pt_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void target_point_struct_get_r_info(const void *s, double **d, int *bounds, bool *is_alloc);
void target_point_struct_set_r(void *s, const void *d, const int *shape);
void surface_curvature_struct_get_xy_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void surface_curvature_struct_set_xy(void *s, const void *d, const int *shape);
void surface_curvature_struct_get_spherical(const void *struct_obj, double *value_out);
void surface_curvature_struct_set_spherical(void *struct_obj, double value_in);
void surface_curvature_struct_get_elliptical_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void surface_curvature_struct_set_elliptical(void *s, const void *d, const int *shape);
void surface_curvature_struct_get_has_curvature(const void *struct_obj, bool *value_out);
void surface_curvature_struct_set_has_curvature(void *struct_obj, bool value_in);
void photon_target_struct_get_type(const void *struct_obj, int *value_out);
void photon_target_struct_set_type(void *struct_obj, int value_in);
void photon_target_struct_get_n_corner(const void *struct_obj, int *value_out);
void photon_target_struct_set_n_corner(void *struct_obj, int value_in);
void photon_target_struct_get_ele_loc(const void *struct_obj, void **ptr_out);
void photon_target_struct_set_ele_loc(void *struct_obj, const void *src_ptr);

void photon_target_struct_get_corner_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void photon_target_struct_get_center(const void *struct_obj, void **ptr_out);
void photon_target_struct_set_center(void *struct_obj, const void *src_ptr);
void photon_material_struct_get_f0_m1(const void *struct_obj, std::complex<double> *value_out);
void photon_material_struct_set_f0_m1(void *struct_obj, std::complex<double> value_in);
void photon_material_struct_get_f0_m2(const void *struct_obj, std::complex<double> *value_out);
void photon_material_struct_set_f0_m2(void *struct_obj, std::complex<double> value_in);
void photon_material_struct_get_f_0(const void *struct_obj, std::complex<double> *value_out);
void photon_material_struct_set_f_0(void *struct_obj, std::complex<double> value_in);
void photon_material_struct_get_f_h(const void *struct_obj, std::complex<double> *value_out);
void photon_material_struct_set_f_h(void *struct_obj, std::complex<double> value_in);
void photon_material_struct_get_f_hbar(const void *struct_obj, std::complex<double> *value_out);
void photon_material_struct_set_f_hbar(void *struct_obj, std::complex<double> value_in);
void photon_material_struct_get_f_hkl(const void *struct_obj, std::complex<double> *value_out);
void photon_material_struct_set_f_hkl(void *struct_obj, std::complex<double> value_in);
void photon_material_struct_get_h_norm_info(const void *s, double **d, int *bounds, bool *is_alloc);
void photon_material_struct_set_h_norm(void *s, const void *d, const int *shape);
void photon_material_struct_get_l_ref_info(const void *s, double **d, int *bounds, bool *is_alloc);
void photon_material_struct_set_l_ref(void *s, const void *d, const int *shape);
void pixel_pt_struct_get_n_photon(const void *struct_obj, int64_t *value_out);
void pixel_pt_struct_set_n_photon(void *struct_obj, int64_t value_in);
void pixel_pt_struct_get_E_x(const void *struct_obj, std::complex<double> *value_out);
void pixel_pt_struct_set_E_x(void *struct_obj, std::complex<double> value_in);
void pixel_pt_struct_get_E_y(const void *struct_obj, std::complex<double> *value_out);
void pixel_pt_struct_set_E_y(void *struct_obj, std::complex<double> value_in);
void pixel_pt_struct_get_intensity_x(const void *struct_obj, double *value_out);
void pixel_pt_struct_set_intensity_x(void *struct_obj, double value_in);
void pixel_pt_struct_get_intensity_y(const void *struct_obj, double *value_out);
void pixel_pt_struct_set_intensity_y(void *struct_obj, double value_in);
void pixel_pt_struct_get_intensity(const void *struct_obj, double *value_out);
void pixel_pt_struct_set_intensity(void *struct_obj, double value_in);
void pixel_pt_struct_get_orbit_info(const void *s, double **d, int *bounds, bool *is_alloc);
void pixel_pt_struct_set_orbit(void *s, const void *d, const int *shape);
void pixel_pt_struct_get_orbit_rms_info(const void *s, double **d, int *bounds, bool *is_alloc);
void pixel_pt_struct_set_orbit_rms(void *s, const void *d, const int *shape);
void pixel_pt_struct_get_init_orbit_info(const void *s, double **d, int *bounds, bool *is_alloc);
void pixel_pt_struct_set_init_orbit(void *s, const void *d, const int *shape);
void pixel_pt_struct_get_init_orbit_rms_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void pixel_pt_struct_set_init_orbit_rms(void *s, const void *d, const int *shape);
void pixel_detec_struct_get_dr_info(const void *s, double **d, int *bounds, bool *is_alloc);
void pixel_detec_struct_set_dr(void *s, const void *d, const int *shape);
void pixel_detec_struct_get_r0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void pixel_detec_struct_set_r0(void *s, const void *d, const int *shape);
void pixel_detec_struct_get_n_track_tot(const void *struct_obj, int64_t *value_out);
void pixel_detec_struct_set_n_track_tot(void *struct_obj, int64_t value_in);
void pixel_detec_struct_get_n_hit_detec(const void *struct_obj, int64_t *value_out);
void pixel_detec_struct_set_n_hit_detec(void *struct_obj, int64_t value_in);
void pixel_detec_struct_get_n_hit_pixel(const void *struct_obj, int64_t *value_out);
void pixel_detec_struct_set_n_hit_pixel(void *struct_obj, int64_t value_in);

void pixel_detec_struct_get_pt_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void photon_element_struct_get_curvature(const void *struct_obj, void **ptr_out);
void photon_element_struct_set_curvature(void *struct_obj, const void *src_ptr);
void photon_element_struct_get_target(const void *struct_obj, void **ptr_out);
void photon_element_struct_set_target(void *struct_obj, const void *src_ptr);
void photon_element_struct_get_material(const void *struct_obj, void **ptr_out);
void photon_element_struct_set_material(void *struct_obj, const void *src_ptr);
void photon_element_struct_get_segmented(const void *struct_obj, void **ptr_out);
void photon_element_struct_set_segmented(void *struct_obj, const void *src_ptr);
void photon_element_struct_get_h_misalign(const void *struct_obj, void **ptr_out);
void photon_element_struct_set_h_misalign(void *struct_obj, const void *src_ptr);
void photon_element_struct_get_displacement(const void *struct_obj, void **ptr_out);
void photon_element_struct_set_displacement(void *struct_obj, const void *src_ptr);
void photon_element_struct_get_pixel(const void *struct_obj, void **ptr_out);
void photon_element_struct_set_pixel(void *struct_obj, const void *src_ptr);
void photon_element_struct_get_reflectivity_table_type(const void *struct_obj, int *value_out);
void photon_element_struct_set_reflectivity_table_type(void *struct_obj, int value_in);
void photon_element_struct_get_reflectivity_table_sigma(const void *struct_obj, void **ptr_out);
void photon_element_struct_set_reflectivity_table_sigma(void *struct_obj, const void *src_ptr);
void photon_element_struct_get_reflectivity_table_pi(const void *struct_obj, void **ptr_out);
void photon_element_struct_set_reflectivity_table_pi(void *struct_obj, const void *src_ptr);

void photon_element_struct_get_init_energy_prob_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void photon_element_struct_get_integrated_init_energy_prob_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void photon_element_struct_set_integrated_init_energy_prob(
    void *s,
    const void *d,
    const int *shape
);
void wall3d_vertex_struct_get_x(const void *struct_obj, double *value_out);
void wall3d_vertex_struct_set_x(void *struct_obj, double value_in);
void wall3d_vertex_struct_get_y(const void *struct_obj, double *value_out);
void wall3d_vertex_struct_set_y(void *struct_obj, double value_in);
void wall3d_vertex_struct_get_radius_x(const void *struct_obj, double *value_out);
void wall3d_vertex_struct_set_radius_x(void *struct_obj, double value_in);
void wall3d_vertex_struct_get_radius_y(const void *struct_obj, double *value_out);
void wall3d_vertex_struct_set_radius_y(void *struct_obj, double value_in);
void wall3d_vertex_struct_get_tilt(const void *struct_obj, double *value_out);
void wall3d_vertex_struct_set_tilt(void *struct_obj, double value_in);
void wall3d_vertex_struct_get_angle(const void *struct_obj, double *value_out);
void wall3d_vertex_struct_set_angle(void *struct_obj, double value_in);
void wall3d_vertex_struct_get_x0(const void *struct_obj, double *value_out);
void wall3d_vertex_struct_set_x0(void *struct_obj, double value_in);
void wall3d_vertex_struct_get_y0(const void *struct_obj, double *value_out);
void wall3d_vertex_struct_set_y0(void *struct_obj, double value_in);
void wall3d_vertex_struct_get_type(const void *struct_obj, int *value_out);
void wall3d_vertex_struct_set_type(void *struct_obj, int value_in);
void wall3d_section_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void wall3d_section_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void wall3d_section_struct_get_material_info(const void *s, char **d, int *bounds, bool *a);
void wall3d_section_struct_set_material(void *struct_obj, const char *str_ptr, int str_len);

void wall3d_section_struct_get_v_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void wall3d_section_struct_get_surface(const void *struct_obj, void **ptr_out);
void wall3d_section_struct_set_surface(void *struct_obj, const void *src_ptr);
void wall3d_section_struct_get_type(const void *struct_obj, int *value_out);
void wall3d_section_struct_set_type(void *struct_obj, int value_in);
void wall3d_section_struct_get_n_vertex_input(const void *struct_obj, int *value_out);
void wall3d_section_struct_set_n_vertex_input(void *struct_obj, int value_in);
void wall3d_section_struct_get_ix_ele(const void *struct_obj, int *value_out);
void wall3d_section_struct_set_ix_ele(void *struct_obj, int value_in);
void wall3d_section_struct_get_ix_branch(const void *struct_obj, int *value_out);
void wall3d_section_struct_set_ix_branch(void *struct_obj, int value_in);
void wall3d_section_struct_get_vertices_state(const void *struct_obj, int *value_out);
void wall3d_section_struct_set_vertices_state(void *struct_obj, int value_in);
void wall3d_section_struct_get_patch_in_region(const void *struct_obj, bool *value_out);
void wall3d_section_struct_set_patch_in_region(void *struct_obj, bool value_in);
void wall3d_section_struct_get_thickness(const void *struct_obj, double *value_out);
void wall3d_section_struct_set_thickness(void *struct_obj, double value_in);
void wall3d_section_struct_get_s(const void *struct_obj, double *value_out);
void wall3d_section_struct_set_s(void *struct_obj, double value_in);
void wall3d_section_struct_get_r0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void wall3d_section_struct_set_r0(void *s, const void *d, const int *shape);
void wall3d_section_struct_get_dx0_ds(const void *struct_obj, double *value_out);
void wall3d_section_struct_set_dx0_ds(void *struct_obj, double value_in);
void wall3d_section_struct_get_dy0_ds(const void *struct_obj, double *value_out);
void wall3d_section_struct_set_dy0_ds(void *struct_obj, double value_in);
void wall3d_section_struct_get_x0_coef_info(const void *s, double **d, int *bounds, bool *is_alloc);
void wall3d_section_struct_set_x0_coef(void *s, const void *d, const int *shape);
void wall3d_section_struct_get_y0_coef_info(const void *s, double **d, int *bounds, bool *is_alloc);
void wall3d_section_struct_set_y0_coef(void *s, const void *d, const int *shape);
void wall3d_section_struct_get_dr_ds(const void *struct_obj, double *value_out);
void wall3d_section_struct_set_dr_ds(void *struct_obj, double value_in);
void wall3d_section_struct_get_p1_coef_info(const void *s, double **d, int *bounds, bool *is_alloc);
void wall3d_section_struct_set_p1_coef(void *s, const void *d, const int *shape);
void wall3d_section_struct_get_p2_coef_info(const void *s, double **d, int *bounds, bool *is_alloc);
void wall3d_section_struct_set_p2_coef(void *s, const void *d, const int *shape);
void wall3d_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void wall3d_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void wall3d_struct_get_type(const void *struct_obj, int *value_out);
void wall3d_struct_set_type(void *struct_obj, int value_in);
void wall3d_struct_get_ix_wall3d(const void *struct_obj, int *value_out);
void wall3d_struct_set_ix_wall3d(void *struct_obj, int value_in);
void wall3d_struct_get_n_link(const void *struct_obj, int *value_out);
void wall3d_struct_set_n_link(void *struct_obj, int value_in);
void wall3d_struct_get_thickness(const void *struct_obj, double *value_out);
void wall3d_struct_set_thickness(void *struct_obj, double value_in);
void wall3d_struct_get_clear_material_info(const void *s, char **d, int *bounds, bool *a);
void wall3d_struct_set_clear_material(void *struct_obj, const char *str_ptr, int str_len);
void wall3d_struct_get_opaque_material_info(const void *s, char **d, int *bounds, bool *a);
void wall3d_struct_set_opaque_material(void *struct_obj, const char *str_ptr, int str_len);
void wall3d_struct_get_superimpose(const void *struct_obj, bool *value_out);
void wall3d_struct_set_superimpose(void *struct_obj, bool value_in);
void wall3d_struct_get_ele_anchor_pt(const void *struct_obj, int *value_out);
void wall3d_struct_set_ele_anchor_pt(void *struct_obj, int value_in);

void wall3d_struct_get_section_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void ramper_lord_struct_get_ix_ele(const void *struct_obj, int *value_out);
void ramper_lord_struct_set_ix_ele(void *struct_obj, int value_in);
void ramper_lord_struct_get_ix_con(const void *struct_obj, int *value_out);
void ramper_lord_struct_set_ix_con(void *struct_obj, int value_in);
void ramper_lord_struct_get_attrib_ptr(const void *struct_obj, double **ptr_out);
void ramper_lord_struct_set_attrib_ptr(void *struct_obj, double value_in);
void control_struct_get_value(const void *struct_obj, double *value_out);
void control_struct_set_value(void *struct_obj, double value_in);
void control_struct_get_y_knot_info(const void *s, double **d, int *bounds, bool *is_alloc);
void control_struct_set_y_knot(void *s, const void *d, const int *shape);

void control_struct_get_stack_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void control_struct_get_slave(const void *struct_obj, void **ptr_out);
void control_struct_set_slave(void *struct_obj, const void *src_ptr);
void control_struct_get_lord(const void *struct_obj, void **ptr_out);
void control_struct_set_lord(void *struct_obj, const void *src_ptr);
void control_struct_get_slave_name_info(const void *s, char **d, int *bounds, bool *a);
void control_struct_set_slave_name(void *struct_obj, const char *str_ptr, int str_len);
void control_struct_get_attribute_info(const void *s, char **d, int *bounds, bool *a);
void control_struct_set_attribute(void *struct_obj, const char *str_ptr, int str_len);
void control_struct_get_ix_attrib(const void *struct_obj, int *value_out);
void control_struct_set_ix_attrib(void *struct_obj, int value_in);
void control_var1_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void control_var1_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void control_var1_struct_get_value(const void *struct_obj, double *value_out);
void control_var1_struct_set_value(void *struct_obj, double value_in);
void control_var1_struct_get_old_value(const void *struct_obj, double *value_out);
void control_var1_struct_set_old_value(void *struct_obj, double value_in);
void control_ramp1_struct_get_y_knot_info(const void *s, double **d, int *bounds, bool *is_alloc);
void control_ramp1_struct_set_y_knot(void *s, const void *d, const int *shape);

void control_ramp1_struct_get_stack_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void control_ramp1_struct_get_attribute_info(const void *s, char **d, int *bounds, bool *a);
void control_ramp1_struct_set_attribute(void *struct_obj, const char *str_ptr, int str_len);
void control_ramp1_struct_get_slave_name_info(const void *s, char **d, int *bounds, bool *a);
void control_ramp1_struct_set_slave_name(void *struct_obj, const char *str_ptr, int str_len);
void control_ramp1_struct_get_is_controller(const void *struct_obj, bool *value_out);
void control_ramp1_struct_set_is_controller(void *struct_obj, bool value_in);

void controller_struct_get_var_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void controller_struct_get_ramp_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void controller_struct_get_ramper_lord_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void controller_struct_get_x_knot_info(const void *s, double **d, int *bounds, bool *is_alloc);
void controller_struct_set_x_knot(void *s, const void *d, const int *shape);
void ellipse_beam_init_struct_get_part_per_ellipse(const void *struct_obj, int *value_out);
void ellipse_beam_init_struct_set_part_per_ellipse(void *struct_obj, int value_in);
void ellipse_beam_init_struct_get_n_ellipse(const void *struct_obj, int *value_out);
void ellipse_beam_init_struct_set_n_ellipse(void *struct_obj, int value_in);
void ellipse_beam_init_struct_get_sigma_cutoff(const void *struct_obj, double *value_out);
void ellipse_beam_init_struct_set_sigma_cutoff(void *struct_obj, double value_in);
void kv_beam_init_struct_get_part_per_phi_info(const void *s, int **d, int *bounds, bool *is_alloc);
void kv_beam_init_struct_set_part_per_phi(void *s, const void *d, const int *shape);
void kv_beam_init_struct_get_n_I2(const void *struct_obj, int *value_out);
void kv_beam_init_struct_set_n_I2(void *struct_obj, int value_in);
void kv_beam_init_struct_get_A(const void *struct_obj, double *value_out);
void kv_beam_init_struct_set_A(void *struct_obj, double value_in);
void grid_beam_init_struct_get_n_x(const void *struct_obj, int *value_out);
void grid_beam_init_struct_set_n_x(void *struct_obj, int value_in);
void grid_beam_init_struct_get_n_px(const void *struct_obj, int *value_out);
void grid_beam_init_struct_set_n_px(void *struct_obj, int value_in);
void grid_beam_init_struct_get_x_min(const void *struct_obj, double *value_out);
void grid_beam_init_struct_set_x_min(void *struct_obj, double value_in);
void grid_beam_init_struct_get_x_max(const void *struct_obj, double *value_out);
void grid_beam_init_struct_set_x_max(void *struct_obj, double value_in);
void grid_beam_init_struct_get_px_min(const void *struct_obj, double *value_out);
void grid_beam_init_struct_set_px_min(void *struct_obj, double value_in);
void grid_beam_init_struct_get_px_max(const void *struct_obj, double *value_out);
void grid_beam_init_struct_set_px_max(void *struct_obj, double value_in);
void beam_init_struct_get_position_file_info(const void *s, char **d, int *bounds, bool *a);
void beam_init_struct_set_position_file(void *struct_obj, const char *str_ptr, int str_len);

void beam_init_struct_get_distribution_type_info(
    const void *s,
    char **d,
    int *bounds, // [lower, upper]
    int *str_len,
    bool *is_alloc
);

void beam_init_struct_get_spin_info(const void *s, double **d, int *bounds, bool *is_alloc);
void beam_init_struct_set_spin(void *s, const void *d, const int *shape);

void beam_init_struct_get_ellipse_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void beam_init_struct_get_KV(const void *struct_obj, void **ptr_out);
void beam_init_struct_set_KV(void *struct_obj, const void *src_ptr);

void beam_init_struct_get_grid_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void beam_init_struct_get_center_jitter_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void beam_init_struct_set_center_jitter(void *s, const void *d, const int *shape);
void beam_init_struct_get_emit_jitter_info(const void *s, double **d, int *bounds, bool *is_alloc);
void beam_init_struct_set_emit_jitter(void *s, const void *d, const int *shape);
void beam_init_struct_get_sig_z_jitter(const void *struct_obj, double *value_out);
void beam_init_struct_set_sig_z_jitter(void *struct_obj, double value_in);
void beam_init_struct_get_sig_pz_jitter(const void *struct_obj, double *value_out);
void beam_init_struct_set_sig_pz_jitter(void *struct_obj, double value_in);
void beam_init_struct_get_n_particle(const void *struct_obj, int *value_out);
void beam_init_struct_set_n_particle(void *struct_obj, int value_in);
void beam_init_struct_get_renorm_center(const void *struct_obj, bool *value_out);
void beam_init_struct_set_renorm_center(void *struct_obj, bool value_in);
void beam_init_struct_get_renorm_sigma(const void *struct_obj, bool *value_out);
void beam_init_struct_set_renorm_sigma(void *struct_obj, bool value_in);
void beam_init_struct_get_random_engine_info(const void *s, char **d, int *bounds, bool *a);
void beam_init_struct_set_random_engine(void *struct_obj, const char *str_ptr, int str_len);
void beam_init_struct_get_random_gauss_converter_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void beam_init_struct_set_random_gauss_converter(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);
void beam_init_struct_get_random_sigma_cutoff(const void *struct_obj, double *value_out);
void beam_init_struct_set_random_sigma_cutoff(void *struct_obj, double value_in);
void beam_init_struct_get_a_norm_emit(const void *struct_obj, double *value_out);
void beam_init_struct_set_a_norm_emit(void *struct_obj, double value_in);
void beam_init_struct_get_b_norm_emit(const void *struct_obj, double *value_out);
void beam_init_struct_set_b_norm_emit(void *struct_obj, double value_in);
void beam_init_struct_get_a_emit(const void *struct_obj, double *value_out);
void beam_init_struct_set_a_emit(void *struct_obj, double value_in);
void beam_init_struct_get_b_emit(const void *struct_obj, double *value_out);
void beam_init_struct_set_b_emit(void *struct_obj, double value_in);
void beam_init_struct_get_dPz_dz(const void *struct_obj, double *value_out);
void beam_init_struct_set_dPz_dz(void *struct_obj, double value_in);
void beam_init_struct_get_center_info(const void *s, double **d, int *bounds, bool *is_alloc);
void beam_init_struct_set_center(void *s, const void *d, const int *shape);
void beam_init_struct_get_t_offset(const void *struct_obj, double *value_out);
void beam_init_struct_set_t_offset(void *struct_obj, double value_in);
void beam_init_struct_get_dt_bunch(const void *struct_obj, double *value_out);
void beam_init_struct_set_dt_bunch(void *struct_obj, double value_in);
void beam_init_struct_get_sig_z(const void *struct_obj, double *value_out);
void beam_init_struct_set_sig_z(void *struct_obj, double value_in);
void beam_init_struct_get_sig_pz(const void *struct_obj, double *value_out);
void beam_init_struct_set_sig_pz(void *struct_obj, double value_in);
void beam_init_struct_get_bunch_charge(const void *struct_obj, double *value_out);
void beam_init_struct_set_bunch_charge(void *struct_obj, double value_in);
void beam_init_struct_get_n_bunch(const void *struct_obj, int *value_out);
void beam_init_struct_set_n_bunch(void *struct_obj, int value_in);
void beam_init_struct_get_ix_turn(const void *struct_obj, int *value_out);
void beam_init_struct_set_ix_turn(void *struct_obj, int value_in);
void beam_init_struct_get_species_info(const void *s, char **d, int *bounds, bool *a);
void beam_init_struct_set_species(void *struct_obj, const char *str_ptr, int str_len);
void beam_init_struct_get_full_6D_coupling_calc(const void *struct_obj, bool *value_out);
void beam_init_struct_set_full_6D_coupling_calc(void *struct_obj, bool value_in);
void beam_init_struct_get_use_particle_start(const void *struct_obj, bool *value_out);
void beam_init_struct_set_use_particle_start(void *struct_obj, bool value_in);
void beam_init_struct_get_use_t_coords(const void *struct_obj, bool *value_out);
void beam_init_struct_set_use_t_coords(void *struct_obj, bool value_in);
void beam_init_struct_get_use_z_as_t(const void *struct_obj, bool *value_out);
void beam_init_struct_set_use_z_as_t(void *struct_obj, bool value_in);
void beam_init_struct_get_file_name_info(const void *s, char **d, int *bounds, bool *a);
void beam_init_struct_set_file_name(void *struct_obj, const char *str_ptr, int str_len);
void lat_param_struct_get_n_part(const void *struct_obj, double *value_out);
void lat_param_struct_set_n_part(void *struct_obj, double value_in);
void lat_param_struct_get_total_length(const void *struct_obj, double *value_out);
void lat_param_struct_set_total_length(void *struct_obj, double value_in);
void lat_param_struct_get_unstable_factor(const void *struct_obj, double *value_out);
void lat_param_struct_set_unstable_factor(void *struct_obj, double value_in);
void lat_param_struct_get_t1_with_RF_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void lat_param_struct_set_t1_with_RF(void *s, const void *d, const int *shape);
void lat_param_struct_get_t1_no_RF_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void lat_param_struct_set_t1_no_RF(void *s, const void *d, const int *shape);
void lat_param_struct_get_spin_tune(const void *struct_obj, double *value_out);
void lat_param_struct_set_spin_tune(void *struct_obj, double value_in);
void lat_param_struct_get_particle(const void *struct_obj, int *value_out);
void lat_param_struct_set_particle(void *struct_obj, int value_in);
void lat_param_struct_get_default_tracking_species(const void *struct_obj, int *value_out);
void lat_param_struct_set_default_tracking_species(void *struct_obj, int value_in);
void lat_param_struct_get_geometry(const void *struct_obj, int *value_out);
void lat_param_struct_set_geometry(void *struct_obj, int value_in);
void lat_param_struct_get_ixx(const void *struct_obj, int *value_out);
void lat_param_struct_set_ixx(void *struct_obj, int value_in);
void lat_param_struct_get_stable(const void *struct_obj, bool *value_out);
void lat_param_struct_set_stable(void *struct_obj, bool value_in);
void lat_param_struct_get_live_branch(const void *struct_obj, bool *value_out);
void lat_param_struct_set_live_branch(void *struct_obj, bool value_in);
void lat_param_struct_get_g1_integral(const void *struct_obj, double *value_out);
void lat_param_struct_set_g1_integral(void *struct_obj, double value_in);
void lat_param_struct_get_g2_integral(const void *struct_obj, double *value_out);
void lat_param_struct_set_g2_integral(void *struct_obj, double value_in);
void lat_param_struct_get_g3_integral(const void *struct_obj, double *value_out);
void lat_param_struct_set_g3_integral(void *struct_obj, double value_in);
void lat_param_struct_get_bookkeeping_state(const void *struct_obj, void **ptr_out);
void lat_param_struct_set_bookkeeping_state(void *struct_obj, const void *src_ptr);
void lat_param_struct_get_beam_init(const void *struct_obj, void **ptr_out);
void lat_param_struct_set_beam_init(void *struct_obj, const void *src_ptr);
void mode_info_struct_get_stable(const void *struct_obj, bool *value_out);
void mode_info_struct_set_stable(void *struct_obj, bool value_in);
void mode_info_struct_get_tune(const void *struct_obj, double *value_out);
void mode_info_struct_set_tune(void *struct_obj, double value_in);
void mode_info_struct_get_emit(const void *struct_obj, double *value_out);
void mode_info_struct_set_emit(void *struct_obj, double value_in);
void mode_info_struct_get_chrom(const void *struct_obj, double *value_out);
void mode_info_struct_set_chrom(void *struct_obj, double value_in);
void mode_info_struct_get_sigma(const void *struct_obj, double *value_out);
void mode_info_struct_set_sigma(void *struct_obj, double value_in);
void mode_info_struct_get_sigmap(const void *struct_obj, double *value_out);
void mode_info_struct_set_sigmap(void *struct_obj, double value_in);
void pre_tracker_struct_get_who(const void *struct_obj, int *value_out);
void pre_tracker_struct_set_who(void *struct_obj, int value_in);
void pre_tracker_struct_get_ix_ele_start(const void *struct_obj, int *value_out);
void pre_tracker_struct_set_ix_ele_start(void *struct_obj, int value_in);
void pre_tracker_struct_get_ix_ele_end(const void *struct_obj, int *value_out);
void pre_tracker_struct_set_ix_ele_end(void *struct_obj, int value_in);
void pre_tracker_struct_get_input_file_info(const void *s, char **d, int *bounds, bool *a);
void pre_tracker_struct_set_input_file(void *struct_obj, const char *str_ptr, int str_len);
void anormal_mode_struct_get_emittance(const void *struct_obj, double *value_out);
void anormal_mode_struct_set_emittance(void *struct_obj, double value_in);
void anormal_mode_struct_get_emittance_no_vert(const void *struct_obj, double *value_out);
void anormal_mode_struct_set_emittance_no_vert(void *struct_obj, double value_in);
void anormal_mode_struct_get_synch_int_info(const void *s, double **d, int *bounds, bool *is_alloc);
void anormal_mode_struct_set_synch_int(void *s, const void *d, const int *shape);
void anormal_mode_struct_get_j_damp(const void *struct_obj, double *value_out);
void anormal_mode_struct_set_j_damp(void *struct_obj, double value_in);
void anormal_mode_struct_get_alpha_damp(const void *struct_obj, double *value_out);
void anormal_mode_struct_set_alpha_damp(void *struct_obj, double value_in);
void anormal_mode_struct_get_chrom(const void *struct_obj, double *value_out);
void anormal_mode_struct_set_chrom(void *struct_obj, double value_in);
void anormal_mode_struct_get_tune(const void *struct_obj, double *value_out);
void anormal_mode_struct_set_tune(void *struct_obj, double value_in);
void linac_normal_mode_struct_get_i2_E4(const void *struct_obj, double *value_out);
void linac_normal_mode_struct_set_i2_E4(void *struct_obj, double value_in);
void linac_normal_mode_struct_get_i3_E7(const void *struct_obj, double *value_out);
void linac_normal_mode_struct_set_i3_E7(void *struct_obj, double value_in);
void linac_normal_mode_struct_get_i5a_E6(const void *struct_obj, double *value_out);
void linac_normal_mode_struct_set_i5a_E6(void *struct_obj, double value_in);
void linac_normal_mode_struct_get_i5b_E6(const void *struct_obj, double *value_out);
void linac_normal_mode_struct_set_i5b_E6(void *struct_obj, double value_in);
void linac_normal_mode_struct_get_sig_E1(const void *struct_obj, double *value_out);
void linac_normal_mode_struct_set_sig_E1(void *struct_obj, double value_in);
void linac_normal_mode_struct_get_a_emittance_end(const void *struct_obj, double *value_out);
void linac_normal_mode_struct_set_a_emittance_end(void *struct_obj, double value_in);
void linac_normal_mode_struct_get_b_emittance_end(const void *struct_obj, double *value_out);
void linac_normal_mode_struct_set_b_emittance_end(void *struct_obj, double value_in);
void normal_modes_struct_get_synch_int_info(const void *s, double **d, int *bounds, bool *is_alloc);
void normal_modes_struct_set_synch_int(void *s, const void *d, const int *shape);
void normal_modes_struct_get_sigE_E(const void *struct_obj, double *value_out);
void normal_modes_struct_set_sigE_E(void *struct_obj, double value_in);
void normal_modes_struct_get_sig_z(const void *struct_obj, double *value_out);
void normal_modes_struct_set_sig_z(void *struct_obj, double value_in);
void normal_modes_struct_get_e_loss(const void *struct_obj, double *value_out);
void normal_modes_struct_set_e_loss(void *struct_obj, double value_in);
void normal_modes_struct_get_rf_voltage(const void *struct_obj, double *value_out);
void normal_modes_struct_set_rf_voltage(void *struct_obj, double value_in);
void normal_modes_struct_get_pz_aperture(const void *struct_obj, double *value_out);
void normal_modes_struct_set_pz_aperture(void *struct_obj, double value_in);
void normal_modes_struct_get_pz_average(const void *struct_obj, double *value_out);
void normal_modes_struct_set_pz_average(void *struct_obj, double value_in);
void normal_modes_struct_get_momentum_compaction(const void *struct_obj, double *value_out);
void normal_modes_struct_set_momentum_compaction(void *struct_obj, double value_in);
void normal_modes_struct_get_dpz_damp(const void *struct_obj, double *value_out);
void normal_modes_struct_set_dpz_damp(void *struct_obj, double value_in);
void normal_modes_struct_get_a(const void *struct_obj, void **ptr_out);
void normal_modes_struct_set_a(void *struct_obj, const void *src_ptr);
void normal_modes_struct_get_b(const void *struct_obj, void **ptr_out);
void normal_modes_struct_set_b(void *struct_obj, const void *src_ptr);
void normal_modes_struct_get_z(const void *struct_obj, void **ptr_out);
void normal_modes_struct_set_z(void *struct_obj, const void *src_ptr);
void normal_modes_struct_get_lin(const void *struct_obj, void **ptr_out);
void normal_modes_struct_set_lin(void *struct_obj, const void *src_ptr);
void em_field_struct_get_E_info(const void *s, double **d, int *bounds, bool *is_alloc);
void em_field_struct_set_E(void *s, const void *d, const int *shape);
void em_field_struct_get_B_info(const void *s, double **d, int *bounds, bool *is_alloc);
void em_field_struct_set_B(void *s, const void *d, const int *shape);
void em_field_struct_get_dE_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void em_field_struct_set_dE(void *s, const void *d, const int *shape);
void em_field_struct_get_dB_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void em_field_struct_set_dB(void *s, const void *d, const int *shape);
void em_field_struct_get_phi(const void *struct_obj, double *value_out);
void em_field_struct_set_phi(void *struct_obj, double value_in);
void em_field_struct_get_phi_B(const void *struct_obj, double *value_out);
void em_field_struct_set_phi_B(void *struct_obj, double value_in);
void em_field_struct_get_A_info(const void *s, double **d, int *bounds, bool *is_alloc);
void em_field_struct_set_A(void *s, const void *d, const int *shape);
void strong_beam_struct_get_ix_slice(const void *struct_obj, int *value_out);
void strong_beam_struct_set_ix_slice(void *struct_obj, int value_in);
void strong_beam_struct_get_x_center(const void *struct_obj, double *value_out);
void strong_beam_struct_set_x_center(void *struct_obj, double value_in);
void strong_beam_struct_get_y_center(const void *struct_obj, double *value_out);
void strong_beam_struct_set_y_center(void *struct_obj, double value_in);
void strong_beam_struct_get_x_sigma(const void *struct_obj, double *value_out);
void strong_beam_struct_set_x_sigma(void *struct_obj, double value_in);
void strong_beam_struct_get_y_sigma(const void *struct_obj, double *value_out);
void strong_beam_struct_set_y_sigma(void *struct_obj, double value_in);
void strong_beam_struct_get_dx(const void *struct_obj, double *value_out);
void strong_beam_struct_set_dx(void *struct_obj, double value_in);
void strong_beam_struct_get_dy(const void *struct_obj, double *value_out);
void strong_beam_struct_set_dy(void *struct_obj, double value_in);
void track_point_struct_get_s_lab(const void *struct_obj, double *value_out);
void track_point_struct_set_s_lab(void *struct_obj, double value_in);
void track_point_struct_get_s_body(const void *struct_obj, double *value_out);
void track_point_struct_set_s_body(void *struct_obj, double value_in);
void track_point_struct_get_orb(const void *struct_obj, void **ptr_out);
void track_point_struct_set_orb(void *struct_obj, const void *src_ptr);
void track_point_struct_get_field(const void *struct_obj, void **ptr_out);
void track_point_struct_set_field(void *struct_obj, const void *src_ptr);
void track_point_struct_get_strong_beam(const void *struct_obj, void **ptr_out);
void track_point_struct_set_strong_beam(void *struct_obj, const void *src_ptr);
void track_point_struct_get_vec0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void track_point_struct_set_vec0(void *s, const void *d, const int *shape);
void track_point_struct_get_mat6_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void track_point_struct_set_mat6(void *s, const void *d, const int *shape);

void track_struct_get_pt_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void track_struct_get_ds_save(const void *struct_obj, double *value_out);
void track_struct_set_ds_save(void *struct_obj, double value_in);
void track_struct_get_n_pt(const void *struct_obj, int *value_out);
void track_struct_set_n_pt(void *struct_obj, int value_in);
void track_struct_get_n_bad(const void *struct_obj, int *value_out);
void track_struct_set_n_bad(void *struct_obj, int value_in);
void track_struct_get_n_ok(const void *struct_obj, int *value_out);
void track_struct_set_n_ok(void *struct_obj, int value_in);
void space_charge_common_struct_get_ds_track_step(const void *struct_obj, double *value_out);
void space_charge_common_struct_set_ds_track_step(void *struct_obj, double value_in);
void space_charge_common_struct_get_dt_track_step(const void *struct_obj, double *value_out);
void space_charge_common_struct_set_dt_track_step(void *struct_obj, double value_in);
void space_charge_common_struct_get_cathode_strength_cutoff(
    const void *struct_obj,
    double *value_out
);
void space_charge_common_struct_set_cathode_strength_cutoff(void *struct_obj, double value_in);
void space_charge_common_struct_get_rel_tol_tracking(const void *struct_obj, double *value_out);
void space_charge_common_struct_set_rel_tol_tracking(void *struct_obj, double value_in);
void space_charge_common_struct_get_abs_tol_tracking(const void *struct_obj, double *value_out);
void space_charge_common_struct_set_abs_tol_tracking(void *struct_obj, double value_in);
void space_charge_common_struct_get_beam_chamber_height(const void *struct_obj, double *value_out);
void space_charge_common_struct_set_beam_chamber_height(void *struct_obj, double value_in);
void space_charge_common_struct_get_lsc_sigma_cutoff(const void *struct_obj, double *value_out);
void space_charge_common_struct_set_lsc_sigma_cutoff(void *struct_obj, double value_in);
void space_charge_common_struct_get_particle_sigma_cutoff(
    const void *struct_obj,
    double *value_out
);
void space_charge_common_struct_set_particle_sigma_cutoff(void *struct_obj, double value_in);
void space_charge_common_struct_get_space_charge_mesh_size_info(
    const void *s,
    int **d,
    int *bounds,
    bool *is_alloc
);
void space_charge_common_struct_set_space_charge_mesh_size(
    void *s,
    const void *d,
    const int *shape
);
void space_charge_common_struct_get_csr3d_mesh_size_info(
    const void *s,
    int **d,
    int *bounds,
    bool *is_alloc
);
void space_charge_common_struct_set_csr3d_mesh_size(void *s, const void *d, const int *shape);
void space_charge_common_struct_get_n_bin(const void *struct_obj, int *value_out);
void space_charge_common_struct_set_n_bin(void *struct_obj, int value_in);
void space_charge_common_struct_get_particle_bin_span(const void *struct_obj, int *value_out);
void space_charge_common_struct_set_particle_bin_span(void *struct_obj, int value_in);
void space_charge_common_struct_get_n_shield_images(const void *struct_obj, int *value_out);
void space_charge_common_struct_set_n_shield_images(void *struct_obj, int value_in);
void space_charge_common_struct_get_sc_min_in_bin(const void *struct_obj, int *value_out);
void space_charge_common_struct_set_sc_min_in_bin(void *struct_obj, int value_in);
void space_charge_common_struct_get_lsc_kick_transverse_dependence(
    const void *struct_obj,
    bool *value_out
);
void space_charge_common_struct_set_lsc_kick_transverse_dependence(void *struct_obj, bool value_in);
void space_charge_common_struct_get_debug(const void *struct_obj, bool *value_out);
void space_charge_common_struct_set_debug(void *struct_obj, bool value_in);
void space_charge_common_struct_get_diagnostic_output_file_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void space_charge_common_struct_set_diagnostic_output_file(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);
void bmad_common_struct_get_max_aperture_limit(const void *struct_obj, double *value_out);
void bmad_common_struct_set_max_aperture_limit(void *struct_obj, double value_in);
void bmad_common_struct_get_d_orb_info(const void *s, double **d, int *bounds, bool *is_alloc);
void bmad_common_struct_set_d_orb(void *s, const void *d, const int *shape);
void bmad_common_struct_get_default_ds_step(const void *struct_obj, double *value_out);
void bmad_common_struct_set_default_ds_step(void *struct_obj, double value_in);
void bmad_common_struct_get_significant_length(const void *struct_obj, double *value_out);
void bmad_common_struct_set_significant_length(void *struct_obj, double value_in);
void bmad_common_struct_get_rel_tol_tracking(const void *struct_obj, double *value_out);
void bmad_common_struct_set_rel_tol_tracking(void *struct_obj, double value_in);
void bmad_common_struct_get_abs_tol_tracking(const void *struct_obj, double *value_out);
void bmad_common_struct_set_abs_tol_tracking(void *struct_obj, double value_in);
void bmad_common_struct_get_rel_tol_adaptive_tracking(const void *struct_obj, double *value_out);
void bmad_common_struct_set_rel_tol_adaptive_tracking(void *struct_obj, double value_in);
void bmad_common_struct_get_abs_tol_adaptive_tracking(const void *struct_obj, double *value_out);
void bmad_common_struct_set_abs_tol_adaptive_tracking(void *struct_obj, double value_in);
void bmad_common_struct_get_init_ds_adaptive_tracking(const void *struct_obj, double *value_out);
void bmad_common_struct_set_init_ds_adaptive_tracking(void *struct_obj, double value_in);
void bmad_common_struct_get_min_ds_adaptive_tracking(const void *struct_obj, double *value_out);
void bmad_common_struct_set_min_ds_adaptive_tracking(void *struct_obj, double value_in);
void bmad_common_struct_get_fatal_ds_adaptive_tracking(const void *struct_obj, double *value_out);
void bmad_common_struct_set_fatal_ds_adaptive_tracking(void *struct_obj, double value_in);
void bmad_common_struct_get_autoscale_amp_abs_tol(const void *struct_obj, double *value_out);
void bmad_common_struct_set_autoscale_amp_abs_tol(void *struct_obj, double value_in);
void bmad_common_struct_get_autoscale_amp_rel_tol(const void *struct_obj, double *value_out);
void bmad_common_struct_set_autoscale_amp_rel_tol(void *struct_obj, double value_in);
void bmad_common_struct_get_autoscale_phase_tol(const void *struct_obj, double *value_out);
void bmad_common_struct_set_autoscale_phase_tol(void *struct_obj, double value_in);
void bmad_common_struct_get_electric_dipole_moment(const void *struct_obj, double *value_out);
void bmad_common_struct_set_electric_dipole_moment(void *struct_obj, double value_in);
void bmad_common_struct_get_synch_rad_scale(const void *struct_obj, double *value_out);
void bmad_common_struct_set_synch_rad_scale(void *struct_obj, double value_in);
void bmad_common_struct_get_sad_eps_scale(const void *struct_obj, double *value_out);
void bmad_common_struct_set_sad_eps_scale(void *struct_obj, double value_in);
void bmad_common_struct_get_sad_amp_max(const void *struct_obj, double *value_out);
void bmad_common_struct_set_sad_amp_max(void *struct_obj, double value_in);
void bmad_common_struct_get_sad_n_div_max(const void *struct_obj, int *value_out);
void bmad_common_struct_set_sad_n_div_max(void *struct_obj, int value_in);
void bmad_common_struct_get_taylor_order(const void *struct_obj, int *value_out);
void bmad_common_struct_set_taylor_order(void *struct_obj, int value_in);
void bmad_common_struct_get_runge_kutta_order(const void *struct_obj, int *value_out);
void bmad_common_struct_set_runge_kutta_order(void *struct_obj, int value_in);
void bmad_common_struct_get_default_integ_order(const void *struct_obj, int *value_out);
void bmad_common_struct_set_default_integ_order(void *struct_obj, int value_in);
void bmad_common_struct_get_max_num_runge_kutta_step(const void *struct_obj, int *value_out);
void bmad_common_struct_set_max_num_runge_kutta_step(void *struct_obj, int value_in);
void bmad_common_struct_get_rf_phase_below_transition_ref(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_rf_phase_below_transition_ref(void *struct_obj, bool value_in);
void bmad_common_struct_get_sr_wakes_on(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_sr_wakes_on(void *struct_obj, bool value_in);
void bmad_common_struct_get_lr_wakes_on(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_lr_wakes_on(void *struct_obj, bool value_in);
void bmad_common_struct_get_auto_bookkeeper(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_auto_bookkeeper(void *struct_obj, bool value_in);
void bmad_common_struct_get_high_energy_space_charge_on(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_high_energy_space_charge_on(void *struct_obj, bool value_in);
void bmad_common_struct_get_csr_and_space_charge_on(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_csr_and_space_charge_on(void *struct_obj, bool value_in);
void bmad_common_struct_get_spin_tracking_on(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_spin_tracking_on(void *struct_obj, bool value_in);
void bmad_common_struct_get_spin_sokolov_ternov_flipping_on(
    const void *struct_obj,
    bool *value_out
);
void bmad_common_struct_set_spin_sokolov_ternov_flipping_on(void *struct_obj, bool value_in);
void bmad_common_struct_get_radiation_damping_on(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_radiation_damping_on(void *struct_obj, bool value_in);
void bmad_common_struct_get_radiation_zero_average(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_radiation_zero_average(void *struct_obj, bool value_in);
void bmad_common_struct_get_radiation_fluctuations_on(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_radiation_fluctuations_on(void *struct_obj, bool value_in);
void bmad_common_struct_get_conserve_taylor_maps(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_conserve_taylor_maps(void *struct_obj, bool value_in);
void bmad_common_struct_get_absolute_time_tracking(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_absolute_time_tracking(void *struct_obj, bool value_in);
void bmad_common_struct_get_absolute_time_ref_shift(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_absolute_time_ref_shift(void *struct_obj, bool value_in);
void bmad_common_struct_get_convert_to_kinetic_momentum(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_convert_to_kinetic_momentum(void *struct_obj, bool value_in);
void bmad_common_struct_get_normalize_twiss(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_normalize_twiss(void *struct_obj, bool value_in);
void bmad_common_struct_get_aperture_limit_on(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_aperture_limit_on(void *struct_obj, bool value_in);
void bmad_common_struct_get_spin_n0_direction_user_set(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_spin_n0_direction_user_set(void *struct_obj, bool value_in);
void bmad_common_struct_get_debug(const void *struct_obj, bool *value_out);
void bmad_common_struct_set_debug(void *struct_obj, bool value_in);
void rad_int1_struct_get_i0(const void *struct_obj, double *value_out);
void rad_int1_struct_set_i0(void *struct_obj, double value_in);
void rad_int1_struct_get_i1(const void *struct_obj, double *value_out);
void rad_int1_struct_set_i1(void *struct_obj, double value_in);
void rad_int1_struct_get_i2(const void *struct_obj, double *value_out);
void rad_int1_struct_set_i2(void *struct_obj, double value_in);
void rad_int1_struct_get_i3(const void *struct_obj, double *value_out);
void rad_int1_struct_set_i3(void *struct_obj, double value_in);
void rad_int1_struct_get_i4a(const void *struct_obj, double *value_out);
void rad_int1_struct_set_i4a(void *struct_obj, double value_in);
void rad_int1_struct_get_i4b(const void *struct_obj, double *value_out);
void rad_int1_struct_set_i4b(void *struct_obj, double value_in);
void rad_int1_struct_get_i4z(const void *struct_obj, double *value_out);
void rad_int1_struct_set_i4z(void *struct_obj, double value_in);
void rad_int1_struct_get_i5a(const void *struct_obj, double *value_out);
void rad_int1_struct_set_i5a(void *struct_obj, double value_in);
void rad_int1_struct_get_i5b(const void *struct_obj, double *value_out);
void rad_int1_struct_set_i5b(void *struct_obj, double value_in);
void rad_int1_struct_get_i6b(const void *struct_obj, double *value_out);
void rad_int1_struct_set_i6b(void *struct_obj, double value_in);
void rad_int1_struct_get_lin_i2_E4(const void *struct_obj, double *value_out);
void rad_int1_struct_set_lin_i2_E4(void *struct_obj, double value_in);
void rad_int1_struct_get_lin_i3_E7(const void *struct_obj, double *value_out);
void rad_int1_struct_set_lin_i3_E7(void *struct_obj, double value_in);
void rad_int1_struct_get_lin_i5a_E6(const void *struct_obj, double *value_out);
void rad_int1_struct_set_lin_i5a_E6(void *struct_obj, double value_in);
void rad_int1_struct_get_lin_i5b_E6(const void *struct_obj, double *value_out);
void rad_int1_struct_set_lin_i5b_E6(void *struct_obj, double value_in);
void rad_int1_struct_get_lin_norm_emit_a(const void *struct_obj, double *value_out);
void rad_int1_struct_set_lin_norm_emit_a(void *struct_obj, double value_in);
void rad_int1_struct_get_lin_norm_emit_b(const void *struct_obj, double *value_out);
void rad_int1_struct_set_lin_norm_emit_b(void *struct_obj, double value_in);
void rad_int1_struct_get_lin_sig_E(const void *struct_obj, double *value_out);
void rad_int1_struct_set_lin_sig_E(void *struct_obj, double value_in);
void rad_int1_struct_get_n_steps(const void *struct_obj, double *value_out);
void rad_int1_struct_set_n_steps(void *struct_obj, double value_in);

void rad_int_branch_struct_get_ele_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void rad_int_all_ele_struct_get_branch_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void rf_stair_step_struct_get_E_tot0(const void *struct_obj, double *value_out);
void rf_stair_step_struct_set_E_tot0(void *struct_obj, double value_in);
void rf_stair_step_struct_get_E_tot1(const void *struct_obj, double *value_out);
void rf_stair_step_struct_set_E_tot1(void *struct_obj, double value_in);
void rf_stair_step_struct_get_p0c(const void *struct_obj, double *value_out);
void rf_stair_step_struct_set_p0c(void *struct_obj, double value_in);
void rf_stair_step_struct_get_p1c(const void *struct_obj, double *value_out);
void rf_stair_step_struct_set_p1c(void *struct_obj, double value_in);
void rf_stair_step_struct_get_scale(const void *struct_obj, double *value_out);
void rf_stair_step_struct_set_scale(void *struct_obj, double value_in);
void rf_stair_step_struct_get_time(const void *struct_obj, double *value_out);
void rf_stair_step_struct_set_time(void *struct_obj, double value_in);
void rf_stair_step_struct_get_s0(const void *struct_obj, double *value_out);
void rf_stair_step_struct_set_s0(void *struct_obj, double value_in);
void rf_stair_step_struct_get_s(const void *struct_obj, double *value_out);
void rf_stair_step_struct_set_s(void *struct_obj, double value_in);
void rf_stair_step_struct_get_ix_step(const void *struct_obj, int *value_out);
void rf_stair_step_struct_set_ix_step(void *struct_obj, int value_in);

void rf_ele_struct_get_steps_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void rf_ele_struct_get_ds_step(const void *struct_obj, double *value_out);
void rf_ele_struct_set_ds_step(void *struct_obj, double value_in);
void ele_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void ele_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void ele_struct_get_type_info(const void *s, char **d, int *bounds, bool *a);
void ele_struct_set_type(void *struct_obj, const char *str_ptr, int str_len);
void ele_struct_get_alias_info(const void *s, char **d, int *bounds, bool *a);
void ele_struct_set_alias(void *struct_obj, const char *str_ptr, int str_len);
void ele_struct_get_component_name_info(const void *s, char **d, int *bounds, bool *a);
void ele_struct_set_component_name(void *struct_obj, const char *str_ptr, int str_len);

void ele_struct_get_descrip_info(const void *s, char **d, int *len, bool *is_alloc);

void ele_struct_set_descrip(void *struct_obj, const char *str_ptr, int str_len);
void ele_struct_get_a(const void *struct_obj, void **ptr_out);
void ele_struct_set_a(void *struct_obj, const void *src_ptr);
void ele_struct_get_b(const void *struct_obj, void **ptr_out);
void ele_struct_set_b(void *struct_obj, const void *src_ptr);
void ele_struct_get_z(const void *struct_obj, void **ptr_out);
void ele_struct_set_z(void *struct_obj, const void *src_ptr);
void ele_struct_get_x(const void *struct_obj, void **ptr_out);
void ele_struct_set_x(void *struct_obj, const void *src_ptr);
void ele_struct_get_y(const void *struct_obj, void **ptr_out);
void ele_struct_set_y(void *struct_obj, const void *src_ptr);
void ele_struct_get_ac_kick(const void *struct_obj, void **ptr_out);
void ele_struct_set_ac_kick(void *struct_obj, const void *src_ptr);
void ele_struct_get_bookkeeping_state(const void *struct_obj, void **ptr_out);
void ele_struct_set_bookkeeping_state(void *struct_obj, const void *src_ptr);
void ele_struct_get_branch(const void *struct_obj, void **ptr_out);
void ele_struct_set_branch(void *struct_obj, const void *src_ptr);
void ele_struct_get_control(const void *struct_obj, void **ptr_out);
void ele_struct_set_control(void *struct_obj, const void *src_ptr);
void ele_struct_get_rf(const void *struct_obj, void **ptr_out);
void ele_struct_set_rf(void *struct_obj, const void *src_ptr);
void ele_struct_get_lord(const void *struct_obj, void **ptr_out);
void ele_struct_set_lord(void *struct_obj, const void *src_ptr);
void ele_struct_get_ptc_fibre(const void *struct_obj, void **ptr_out);
void ele_struct_set_ptc_fibre(void *struct_obj, const void *src_ptr);
void ele_struct_get_floor(const void *struct_obj, void **ptr_out);
void ele_struct_set_floor(void *struct_obj, const void *src_ptr);
void ele_struct_get_high_energy_space_charge(const void *struct_obj, void **ptr_out);
void ele_struct_set_high_energy_space_charge(void *struct_obj, const void *src_ptr);
void ele_struct_get_mode3(const void *struct_obj, void **ptr_out);
void ele_struct_set_mode3(void *struct_obj, const void *src_ptr);
void ele_struct_get_photon(const void *struct_obj, void **ptr_out);
void ele_struct_set_photon(void *struct_obj, const void *src_ptr);
void ele_struct_get_rad_map(const void *struct_obj, void **ptr_out);
void ele_struct_set_rad_map(void *struct_obj, const void *src_ptr);

void ele_struct_get_taylor_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void ele_struct_get_spin_taylor_ref_orb_in_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void ele_struct_set_spin_taylor_ref_orb_in(void *s, const void *d, const int *shape);

void ele_struct_get_spin_taylor_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void ele_struct_get_wake(const void *struct_obj, void **ptr_out);
void ele_struct_set_wake(void *struct_obj, const void *src_ptr);

void ele_struct_get_wall3d_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void ele_struct_get_cartesian_map_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void ele_struct_get_cylindrical_map_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void ele_struct_get_gen_grad_map_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void ele_struct_get_grid_field_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void ele_struct_get_map_ref_orb_in(const void *struct_obj, void **ptr_out);
void ele_struct_set_map_ref_orb_in(void *struct_obj, const void *src_ptr);
void ele_struct_get_map_ref_orb_out(const void *struct_obj, void **ptr_out);
void ele_struct_set_map_ref_orb_out(void *struct_obj, const void *src_ptr);
void ele_struct_get_time_ref_orb_in(const void *struct_obj, void **ptr_out);
void ele_struct_set_time_ref_orb_in(void *struct_obj, const void *src_ptr);
void ele_struct_get_time_ref_orb_out(const void *struct_obj, void **ptr_out);
void ele_struct_set_time_ref_orb_out(void *struct_obj, const void *src_ptr);
void ele_struct_get_value_info(const void *s, double **d, int *bounds, bool *is_alloc);
void ele_struct_set_value(void *s, const void *d, const int *shape);
void ele_struct_get_old_value_info(const void *s, double **d, int *bounds, bool *is_alloc);
void ele_struct_set_old_value(void *s, const void *d, const int *shape);
void ele_struct_get_spin_q_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void ele_struct_set_spin_q(void *s, const void *d, const int *shape);
void ele_struct_get_vec0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void ele_struct_set_vec0(void *s, const void *d, const int *shape);
void ele_struct_get_mat6_info(const void *s, double **d, int *bounds, int *strides, bool *is_alloc);
void ele_struct_set_mat6(void *s, const void *d, const int *shape);
void ele_struct_get_c_mat_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void ele_struct_set_c_mat(void *s, const void *d, const int *shape);
void ele_struct_get_dc_mat_dpz_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void ele_struct_set_dc_mat_dpz(void *s, const void *d, const int *shape);
void ele_struct_get_gamma_c(const void *struct_obj, double *value_out);
void ele_struct_set_gamma_c(void *struct_obj, double value_in);
void ele_struct_get_s_start(const void *struct_obj, double *value_out);
void ele_struct_set_s_start(void *struct_obj, double value_in);
void ele_struct_get_s(const void *struct_obj, double *value_out);
void ele_struct_set_s(void *struct_obj, double value_in);
void ele_struct_get_ref_time(const void *struct_obj, double *value_out);
void ele_struct_set_ref_time(void *struct_obj, double value_in);
void ele_struct_get_a_pole_info(const void *s, double **d, int *bounds, bool *is_alloc);
void ele_struct_set_a_pole(void *s, const void *d, const int *shape);
void ele_struct_get_b_pole_info(const void *s, double **d, int *bounds, bool *is_alloc);
void ele_struct_set_b_pole(void *s, const void *d, const int *shape);
void ele_struct_get_a_pole_elec_info(const void *s, double **d, int *bounds, bool *is_alloc);
void ele_struct_set_a_pole_elec(void *s, const void *d, const int *shape);
void ele_struct_get_b_pole_elec_info(const void *s, double **d, int *bounds, bool *is_alloc);
void ele_struct_set_b_pole_elec(void *s, const void *d, const int *shape);
void ele_struct_get_custom_info(const void *s, double **d, int *bounds, bool *is_alloc);
void ele_struct_set_custom(void *s, const void *d, const int *shape);
void ele_struct_get_r_info(const void *s, double **d, int *bounds, int *strides, bool *is_alloc);
void ele_struct_set_r(void *s, const void *d, const int *shape);
void ele_struct_get_key(const void *struct_obj, int *value_out);
void ele_struct_set_key(void *struct_obj, int value_in);
void ele_struct_get_sub_key(const void *struct_obj, int *value_out);
void ele_struct_set_sub_key(void *struct_obj, int value_in);
void ele_struct_get_ix_ele(const void *struct_obj, int *value_out);
void ele_struct_set_ix_ele(void *struct_obj, int value_in);
void ele_struct_get_ix_branch(const void *struct_obj, int *value_out);
void ele_struct_set_ix_branch(void *struct_obj, int value_in);
void ele_struct_get_lord_status(const void *struct_obj, int *value_out);
void ele_struct_set_lord_status(void *struct_obj, int value_in);
void ele_struct_get_n_slave(const void *struct_obj, int *value_out);
void ele_struct_set_n_slave(void *struct_obj, int value_in);
void ele_struct_get_n_slave_field(const void *struct_obj, int *value_out);
void ele_struct_set_n_slave_field(void *struct_obj, int value_in);
void ele_struct_get_ix1_slave(const void *struct_obj, int *value_out);
void ele_struct_set_ix1_slave(void *struct_obj, int value_in);
void ele_struct_get_slave_status(const void *struct_obj, int *value_out);
void ele_struct_set_slave_status(void *struct_obj, int value_in);
void ele_struct_get_n_lord(const void *struct_obj, int *value_out);
void ele_struct_set_n_lord(void *struct_obj, int value_in);
void ele_struct_get_n_lord_field(const void *struct_obj, int *value_out);
void ele_struct_set_n_lord_field(void *struct_obj, int value_in);
void ele_struct_get_n_lord_ramper(const void *struct_obj, int *value_out);
void ele_struct_set_n_lord_ramper(void *struct_obj, int value_in);
void ele_struct_get_ic1_lord(const void *struct_obj, int *value_out);
void ele_struct_set_ic1_lord(void *struct_obj, int value_in);
void ele_struct_get_ix_pointer(const void *struct_obj, int *value_out);
void ele_struct_set_ix_pointer(void *struct_obj, int value_in);
void ele_struct_get_ixx(const void *struct_obj, int *value_out);
void ele_struct_set_ixx(void *struct_obj, int value_in);
void ele_struct_get_iyy(const void *struct_obj, int *value_out);
void ele_struct_set_iyy(void *struct_obj, int value_in);
void ele_struct_get_izz(const void *struct_obj, int *value_out);
void ele_struct_set_izz(void *struct_obj, int value_in);
void ele_struct_get_mat6_calc_method(const void *struct_obj, int *value_out);
void ele_struct_set_mat6_calc_method(void *struct_obj, int value_in);
void ele_struct_get_tracking_method(const void *struct_obj, int *value_out);
void ele_struct_set_tracking_method(void *struct_obj, int value_in);
void ele_struct_get_spin_tracking_method(const void *struct_obj, int *value_out);
void ele_struct_set_spin_tracking_method(void *struct_obj, int value_in);
void ele_struct_get_csr_method(const void *struct_obj, int *value_out);
void ele_struct_set_csr_method(void *struct_obj, int value_in);
void ele_struct_get_space_charge_method(const void *struct_obj, int *value_out);
void ele_struct_set_space_charge_method(void *struct_obj, int value_in);
void ele_struct_get_ptc_integration_type(const void *struct_obj, int *value_out);
void ele_struct_set_ptc_integration_type(void *struct_obj, int value_in);
void ele_struct_get_field_calc(const void *struct_obj, int *value_out);
void ele_struct_set_field_calc(void *struct_obj, int value_in);
void ele_struct_get_aperture_at(const void *struct_obj, int *value_out);
void ele_struct_set_aperture_at(void *struct_obj, int value_in);
void ele_struct_get_aperture_type(const void *struct_obj, int *value_out);
void ele_struct_set_aperture_type(void *struct_obj, int value_in);
void ele_struct_get_ref_species(const void *struct_obj, int *value_out);
void ele_struct_set_ref_species(void *struct_obj, int value_in);
void ele_struct_get_orientation(const void *struct_obj, int *value_out);
void ele_struct_set_orientation(void *struct_obj, int value_in);
void ele_struct_get_symplectify(const void *struct_obj, bool *value_out);
void ele_struct_set_symplectify(void *struct_obj, bool value_in);
void ele_struct_get_mode_flip(const void *struct_obj, bool *value_out);
void ele_struct_set_mode_flip(void *struct_obj, bool value_in);
void ele_struct_get_multipoles_on(const void *struct_obj, bool *value_out);
void ele_struct_set_multipoles_on(void *struct_obj, bool value_in);
void ele_struct_get_scale_multipoles(const void *struct_obj, bool *value_out);
void ele_struct_set_scale_multipoles(void *struct_obj, bool value_in);
void ele_struct_get_taylor_map_includes_offsets(const void *struct_obj, bool *value_out);
void ele_struct_set_taylor_map_includes_offsets(void *struct_obj, bool value_in);
void ele_struct_get_field_master(const void *struct_obj, bool *value_out);
void ele_struct_set_field_master(void *struct_obj, bool value_in);
void ele_struct_get_is_on(const void *struct_obj, bool *value_out);
void ele_struct_set_is_on(void *struct_obj, bool value_in);
void ele_struct_get_logic(const void *struct_obj, bool *value_out);
void ele_struct_set_logic(void *struct_obj, bool value_in);
void ele_struct_get_bmad_logic(const void *struct_obj, bool *value_out);
void ele_struct_set_bmad_logic(void *struct_obj, bool value_in);
void ele_struct_get_select(const void *struct_obj, bool *value_out);
void ele_struct_set_select(void *struct_obj, bool value_in);
void ele_struct_get_offset_moves_aperture(const void *struct_obj, bool *value_out);
void ele_struct_set_offset_moves_aperture(void *struct_obj, bool value_in);
void complex_taylor_term_struct_get_coef(const void *struct_obj, std::complex<double> *value_out);
void complex_taylor_term_struct_set_coef(void *struct_obj, std::complex<double> value_in);
void complex_taylor_term_struct_get_expn_info(const void *s, int **d, int *bounds, bool *is_alloc);
void complex_taylor_term_struct_set_expn(void *s, const void *d, const int *shape);
void complex_taylor_struct_get_ref(const void *struct_obj, std::complex<double> *value_out);
void complex_taylor_struct_set_ref(void *struct_obj, std::complex<double> value_in);

void complex_taylor_struct_get_term_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void branch_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void branch_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void branch_struct_get_ix_branch(const void *struct_obj, int *value_out);
void branch_struct_set_ix_branch(void *struct_obj, int value_in);
void branch_struct_get_ix_from_branch(const void *struct_obj, int *value_out);
void branch_struct_set_ix_from_branch(void *struct_obj, int value_in);
void branch_struct_get_ix_from_ele(const void *struct_obj, int *value_out);
void branch_struct_set_ix_from_ele(void *struct_obj, int value_in);
void branch_struct_get_ix_to_ele(const void *struct_obj, int *value_out);
void branch_struct_set_ix_to_ele(void *struct_obj, int value_in);
void branch_struct_get_ix_fixer(const void *struct_obj, int *value_out);
void branch_struct_set_ix_fixer(void *struct_obj, int value_in);
void branch_struct_get_n_ele_track(const void *struct_obj, int *value_out);
void branch_struct_set_n_ele_track(void *struct_obj, int value_in);
void branch_struct_get_n_ele_max(const void *struct_obj, int *value_out);
void branch_struct_set_n_ele_max(void *struct_obj, int value_in);
void branch_struct_get_lat(const void *struct_obj, void **ptr_out);
void branch_struct_set_lat(void *struct_obj, const void *src_ptr);
void branch_struct_get_a(const void *struct_obj, void **ptr_out);
void branch_struct_set_a(void *struct_obj, const void *src_ptr);
void branch_struct_get_b(const void *struct_obj, void **ptr_out);
void branch_struct_set_b(void *struct_obj, const void *src_ptr);
void branch_struct_get_z(const void *struct_obj, void **ptr_out);
void branch_struct_set_z(void *struct_obj, const void *src_ptr);

void branch_struct_get_ele_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void branch_struct_get_param(const void *struct_obj, void **ptr_out);
void branch_struct_set_param(void *struct_obj, const void *src_ptr);
void branch_struct_get_particle_start(const void *struct_obj, void **ptr_out);
void branch_struct_set_particle_start(void *struct_obj, const void *src_ptr);

void branch_struct_get_wall3d_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void lat_struct_get_use_name_info(const void *s, char **d, int *bounds, bool *a);
void lat_struct_set_use_name(void *struct_obj, const char *str_ptr, int str_len);
void lat_struct_get_lattice_info(const void *s, char **d, int *bounds, bool *a);
void lat_struct_set_lattice(void *struct_obj, const char *str_ptr, int str_len);
void lat_struct_get_machine_info(const void *s, char **d, int *bounds, bool *a);
void lat_struct_set_machine(void *struct_obj, const char *str_ptr, int str_len);
void lat_struct_get_input_file_name_info(const void *s, char **d, int *bounds, bool *a);
void lat_struct_set_input_file_name(void *struct_obj, const char *str_ptr, int str_len);
void lat_struct_get_title_info(const void *s, char **d, int *bounds, bool *a);
void lat_struct_set_title(void *struct_obj, const char *str_ptr, int str_len);

void lat_struct_get_print_str_info(
    const void *s,
    char **d,
    int *bounds, // [lower, upper]
    int *str_len,
    bool *is_alloc
);

void lat_struct_get_constant_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void lat_struct_get_a(const void *struct_obj, void **ptr_out);
void lat_struct_set_a(void *struct_obj, const void *src_ptr);
void lat_struct_get_b(const void *struct_obj, void **ptr_out);
void lat_struct_set_b(void *struct_obj, const void *src_ptr);
void lat_struct_get_z(const void *struct_obj, void **ptr_out);
void lat_struct_set_z(void *struct_obj, const void *src_ptr);
void lat_struct_get_param(const void *struct_obj, void **ptr_out);
void lat_struct_set_param(void *struct_obj, const void *src_ptr);
void lat_struct_get_lord_state(const void *struct_obj, void **ptr_out);
void lat_struct_set_lord_state(void *struct_obj, const void *src_ptr);
void lat_struct_get_ele_init(const void *struct_obj, void **ptr_out);
void lat_struct_set_ele_init(void *struct_obj, const void *src_ptr);

void lat_struct_get_ele_info(const void *s, void **d, int *bounds, bool *is_alloc, size_t *el_size);

void lat_struct_get_branch_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void lat_struct_get_control_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void lat_struct_get_particle_start(const void *struct_obj, void **ptr_out);
void lat_struct_set_particle_start(void *struct_obj, const void *src_ptr);
void lat_struct_get_beam_init(const void *struct_obj, void **ptr_out);
void lat_struct_set_beam_init(void *struct_obj, const void *src_ptr);
void lat_struct_get_pre_tracker(const void *struct_obj, void **ptr_out);
void lat_struct_set_pre_tracker(void *struct_obj, const void *src_ptr);
void lat_struct_get_custom_info(const void *s, double **d, int *bounds, bool *is_alloc);
void lat_struct_set_custom(void *s, const void *d, const int *shape);
void lat_struct_get_version(const void *struct_obj, int *value_out);
void lat_struct_set_version(void *struct_obj, int value_in);
void lat_struct_get_n_ele_track(const void *struct_obj, int **ptr_out);
void lat_struct_set_n_ele_track(void *struct_obj, int value_in);
void lat_struct_get_n_ele_max(const void *struct_obj, int **ptr_out);
void lat_struct_set_n_ele_max(void *struct_obj, int value_in);
void lat_struct_get_n_control_max(const void *struct_obj, int *value_out);
void lat_struct_set_n_control_max(void *struct_obj, int value_in);
void lat_struct_get_n_ic_max(const void *struct_obj, int *value_out);
void lat_struct_set_n_ic_max(void *struct_obj, int value_in);
void lat_struct_get_input_taylor_order(const void *struct_obj, int *value_out);
void lat_struct_set_input_taylor_order(void *struct_obj, int value_in);
void lat_struct_get_ic_info(const void *s, int **d, int *bounds, bool *is_alloc);
void lat_struct_set_ic(void *s, const void *d, const int *shape);
void lat_struct_get_photon_type(const void *struct_obj, int *value_out);
void lat_struct_set_photon_type(void *struct_obj, int value_in);
void lat_struct_get_creation_hash(const void *struct_obj, int *value_out);
void lat_struct_set_creation_hash(void *struct_obj, int value_in);
void lat_struct_get_ramper_slave_bookkeeping(const void *struct_obj, int *value_out);
void lat_struct_set_ramper_slave_bookkeeping(void *struct_obj, int value_in);

void bunch_struct_get_particle_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bunch_struct_get_ix_z_info(const void *s, int **d, int *bounds, bool *is_alloc);
void bunch_struct_set_ix_z(void *s, const void *d, const int *shape);
void bunch_struct_get_charge_tot(const void *struct_obj, double *value_out);
void bunch_struct_set_charge_tot(void *struct_obj, double value_in);
void bunch_struct_get_charge_live(const void *struct_obj, double *value_out);
void bunch_struct_set_charge_live(void *struct_obj, double value_in);
void bunch_struct_get_z_center(const void *struct_obj, double *value_out);
void bunch_struct_set_z_center(void *struct_obj, double value_in);
void bunch_struct_get_t_center(const void *struct_obj, double *value_out);
void bunch_struct_set_t_center(void *struct_obj, double value_in);
void bunch_struct_get_t0(const void *struct_obj, double *value_out);
void bunch_struct_set_t0(void *struct_obj, double value_in);
void bunch_struct_get_drift_between_t_and_s(const void *struct_obj, bool *value_out);
void bunch_struct_set_drift_between_t_and_s(void *struct_obj, bool value_in);
void bunch_struct_get_ix_ele(const void *struct_obj, int *value_out);
void bunch_struct_set_ix_ele(void *struct_obj, int value_in);
void bunch_struct_get_ix_bunch(const void *struct_obj, int *value_out);
void bunch_struct_set_ix_bunch(void *struct_obj, int value_in);
void bunch_struct_get_ix_turn(const void *struct_obj, int *value_out);
void bunch_struct_set_ix_turn(void *struct_obj, int value_in);
void bunch_struct_get_n_live(const void *struct_obj, int *value_out);
void bunch_struct_set_n_live(void *struct_obj, int value_in);
void bunch_struct_get_n_good(const void *struct_obj, int *value_out);
void bunch_struct_set_n_good(void *struct_obj, int value_in);
void bunch_struct_get_n_bad(const void *struct_obj, int *value_out);
void bunch_struct_set_n_bad(void *struct_obj, int value_in);
void bunch_params_struct_get_centroid(const void *struct_obj, void **ptr_out);
void bunch_params_struct_set_centroid(void *struct_obj, const void *src_ptr);
void bunch_params_struct_get_x(const void *struct_obj, void **ptr_out);
void bunch_params_struct_set_x(void *struct_obj, const void *src_ptr);
void bunch_params_struct_get_y(const void *struct_obj, void **ptr_out);
void bunch_params_struct_set_y(void *struct_obj, const void *src_ptr);
void bunch_params_struct_get_z(const void *struct_obj, void **ptr_out);
void bunch_params_struct_set_z(void *struct_obj, const void *src_ptr);
void bunch_params_struct_get_a(const void *struct_obj, void **ptr_out);
void bunch_params_struct_set_a(void *struct_obj, const void *src_ptr);
void bunch_params_struct_get_b(const void *struct_obj, void **ptr_out);
void bunch_params_struct_set_b(void *struct_obj, const void *src_ptr);
void bunch_params_struct_get_c(const void *struct_obj, void **ptr_out);
void bunch_params_struct_set_c(void *struct_obj, const void *src_ptr);
void bunch_params_struct_get_sigma_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void bunch_params_struct_set_sigma(void *s, const void *d, const int *shape);
void bunch_params_struct_get_rel_max_info(const void *s, double **d, int *bounds, bool *is_alloc);
void bunch_params_struct_set_rel_max(void *s, const void *d, const int *shape);
void bunch_params_struct_get_rel_min_info(const void *s, double **d, int *bounds, bool *is_alloc);
void bunch_params_struct_set_rel_min(void *s, const void *d, const int *shape);
void bunch_params_struct_get_s(const void *struct_obj, double *value_out);
void bunch_params_struct_set_s(void *struct_obj, double value_in);
void bunch_params_struct_get_t(const void *struct_obj, double *value_out);
void bunch_params_struct_set_t(void *struct_obj, double value_in);
void bunch_params_struct_get_sigma_t(const void *struct_obj, double *value_out);
void bunch_params_struct_set_sigma_t(void *struct_obj, double value_in);
void bunch_params_struct_get_charge_live(const void *struct_obj, double *value_out);
void bunch_params_struct_set_charge_live(void *struct_obj, double value_in);
void bunch_params_struct_get_charge_tot(const void *struct_obj, double *value_out);
void bunch_params_struct_set_charge_tot(void *struct_obj, double value_in);
void bunch_params_struct_get_n_particle_tot(const void *struct_obj, int *value_out);
void bunch_params_struct_set_n_particle_tot(void *struct_obj, int value_in);
void bunch_params_struct_get_n_particle_live(const void *struct_obj, int *value_out);
void bunch_params_struct_set_n_particle_live(void *struct_obj, int value_in);
void bunch_params_struct_get_n_particle_lost_in_ele(const void *struct_obj, int *value_out);
void bunch_params_struct_set_n_particle_lost_in_ele(void *struct_obj, int value_in);
void bunch_params_struct_get_n_good_steps(const void *struct_obj, int *value_out);
void bunch_params_struct_set_n_good_steps(void *struct_obj, int value_in);
void bunch_params_struct_get_n_bad_steps(const void *struct_obj, int *value_out);
void bunch_params_struct_set_n_bad_steps(void *struct_obj, int value_in);
void bunch_params_struct_get_ix_ele(const void *struct_obj, int *value_out);
void bunch_params_struct_set_ix_ele(void *struct_obj, int value_in);
void bunch_params_struct_get_location(const void *struct_obj, int *value_out);
void bunch_params_struct_set_location(void *struct_obj, int value_in);
void bunch_params_struct_get_twiss_valid(const void *struct_obj, bool *value_out);
void bunch_params_struct_set_twiss_valid(void *struct_obj, bool value_in);

void beam_struct_get_bunch_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void aperture_point_struct_get_x(const void *struct_obj, double *value_out);
void aperture_point_struct_set_x(void *struct_obj, double value_in);
void aperture_point_struct_get_y(const void *struct_obj, double *value_out);
void aperture_point_struct_set_y(void *struct_obj, double value_in);
void aperture_point_struct_get_plane(const void *struct_obj, int *value_out);
void aperture_point_struct_set_plane(void *struct_obj, int value_in);
void aperture_point_struct_get_ix_ele(const void *struct_obj, int *value_out);
void aperture_point_struct_set_ix_ele(void *struct_obj, int value_in);
void aperture_point_struct_get_i_turn(const void *struct_obj, int *value_out);
void aperture_point_struct_set_i_turn(void *struct_obj, int value_in);
void aperture_param_struct_get_min_angle(const void *struct_obj, double *value_out);
void aperture_param_struct_set_min_angle(void *struct_obj, double value_in);
void aperture_param_struct_get_max_angle(const void *struct_obj, double *value_out);
void aperture_param_struct_set_max_angle(void *struct_obj, double value_in);
void aperture_param_struct_get_n_angle(const void *struct_obj, int *value_out);
void aperture_param_struct_set_n_angle(void *struct_obj, int value_in);
void aperture_param_struct_get_n_turn(const void *struct_obj, int *value_out);
void aperture_param_struct_set_n_turn(void *struct_obj, int value_in);
void aperture_param_struct_get_x_init(const void *struct_obj, double *value_out);
void aperture_param_struct_set_x_init(void *struct_obj, double value_in);
void aperture_param_struct_get_y_init(const void *struct_obj, double *value_out);
void aperture_param_struct_set_y_init(void *struct_obj, double value_in);
void aperture_param_struct_get_rel_accuracy(const void *struct_obj, double *value_out);
void aperture_param_struct_set_rel_accuracy(void *struct_obj, double value_in);
void aperture_param_struct_get_abs_accuracy(const void *struct_obj, double *value_out);
void aperture_param_struct_set_abs_accuracy(void *struct_obj, double value_in);
void aperture_param_struct_get_start_ele_info(const void *s, char **d, int *bounds, bool *a);
void aperture_param_struct_set_start_ele(void *struct_obj, const char *str_ptr, int str_len);

void aperture_scan_struct_get_point_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void aperture_scan_struct_get_ref_orb(const void *struct_obj, void **ptr_out);
void aperture_scan_struct_set_ref_orb(void *struct_obj, const void *src_ptr);
void aperture_scan_struct_get_pz_start(const void *struct_obj, double *value_out);
void aperture_scan_struct_set_pz_start(void *struct_obj, double value_in);
void ele_pointer_struct_get_ele(const void *struct_obj, void **ptr_out);
void ele_pointer_struct_set_ele(void *struct_obj, const void *src_ptr);
void ele_pointer_struct_get_loc(const void *struct_obj, void **ptr_out);
void ele_pointer_struct_set_loc(void *struct_obj, const void *src_ptr);
void ele_pointer_struct_get_id(const void *struct_obj, int *value_out);
void ele_pointer_struct_set_id(void *struct_obj, int value_in);
void expression_tree_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void expression_tree_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void expression_tree_struct_get_type(const void *struct_obj, int *value_out);
void expression_tree_struct_set_type(void *struct_obj, int value_in);
void expression_tree_struct_get_value(const void *struct_obj, double *value_out);
void expression_tree_struct_set_value(void *struct_obj, double value_in);

void expression_tree_struct_get_node_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void nametable_struct_get_name_info(
    const void *s,
    char **d,
    int *bounds, // [lower, upper]
    int *str_len,
    bool *is_alloc
);

void nametable_struct_get_index_info(const void *s, int **d, int *bounds, bool *is_alloc);
void nametable_struct_set_index(void *s, const void *d, const int *shape);
void nametable_struct_get_n_min(const void *struct_obj, int *value_out);
void nametable_struct_set_n_min(void *struct_obj, int value_in);
void nametable_struct_get_n_max(const void *struct_obj, int *value_out);
void nametable_struct_set_n_max(void *struct_obj, int value_in);
void tao_spin_dn_dpz_struct_get_vec_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_spin_dn_dpz_struct_set_vec(void *s, const void *d, const int *shape);
void tao_spin_dn_dpz_struct_get_partial_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void tao_spin_dn_dpz_struct_set_partial(void *s, const void *d, const int *shape);
void tao_spin_dn_dpz_struct_get_partial2_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void tao_spin_dn_dpz_struct_set_partial2(void *s, const void *d, const int *shape);
void resonance_h_struct_get_id_info(const void *s, char **d, int *bounds, bool *a);
void resonance_h_struct_set_id(void *struct_obj, const char *str_ptr, int str_len);
void resonance_h_struct_get_c_val(const void *struct_obj, std::complex<double> *value_out);
void resonance_h_struct_set_c_val(void *struct_obj, std::complex<double> value_in);
void spin_orbit_map1_struct_get_orb_mat_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void spin_orbit_map1_struct_set_orb_mat(void *s, const void *d, const int *shape);
void spin_orbit_map1_struct_get_vec0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void spin_orbit_map1_struct_set_vec0(void *s, const void *d, const int *shape);
void spin_orbit_map1_struct_get_spin_q_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void spin_orbit_map1_struct_set_spin_q(void *s, const void *d, const int *shape);
void spin_axis_struct_get_l_info(const void *s, double **d, int *bounds, bool *is_alloc);
void spin_axis_struct_set_l(void *s, const void *d, const int *shape);
void spin_axis_struct_get_n0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void spin_axis_struct_set_n0(void *s, const void *d, const int *shape);
void spin_axis_struct_get_m_info(const void *s, double **d, int *bounds, bool *is_alloc);
void spin_axis_struct_set_m(void *s, const void *d, const int *shape);
void ptc_normal_form_struct_get_ele_origin(const void *struct_obj, void **ptr_out);
void ptc_normal_form_struct_set_ele_origin(void *struct_obj, const void *src_ptr);
void ptc_normal_form_struct_get_orb0_info(const void *s, double **d, int *bounds, bool *is_alloc);
void ptc_normal_form_struct_set_orb0(void *s, const void *d, const int *shape);
void ptc_normal_form_struct_get_valid_map(const void *struct_obj, bool *value_out);
void ptc_normal_form_struct_set_valid_map(void *struct_obj, bool value_in);
void bmad_normal_form_struct_get_ele_origin(const void *struct_obj, void **ptr_out);
void bmad_normal_form_struct_set_ele_origin(void *struct_obj, const void *src_ptr);

void bmad_normal_form_struct_get_M_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bmad_normal_form_struct_get_A_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bmad_normal_form_struct_get_A_inv_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bmad_normal_form_struct_get_dhdj_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bmad_normal_form_struct_get_F_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bmad_normal_form_struct_get_L_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bmad_normal_form_struct_get_h_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bunch_track_struct_get_pt_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bunch_track_struct_get_ds_save(const void *struct_obj, double *value_out);
void bunch_track_struct_set_ds_save(void *struct_obj, double value_in);
void bunch_track_struct_get_n_pt(const void *struct_obj, int *value_out);
void bunch_track_struct_set_n_pt(void *struct_obj, int value_in);
void summation_rdt_struct_get_h11001(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h11001(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h00111(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h00111(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h20001(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h20001(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h00201(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h00201(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h10002(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h10002(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h21000(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h21000(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h30000(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h30000(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h10110(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h10110(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h10020(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h10020(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h10200(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h10200(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h31000(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h31000(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h40000(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h40000(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h20110(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h20110(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h11200(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h11200(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h20020(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h20020(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h20200(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h20200(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h00310(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h00310(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h00400(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h00400(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h22000(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h22000(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h00220(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h00220(void *struct_obj, std::complex<double> value_in);
void summation_rdt_struct_get_h11110(const void *struct_obj, std::complex<double> *value_out);
void summation_rdt_struct_set_h11110(void *struct_obj, std::complex<double> value_in);
void tao_ele_shape_struct_get_ele_id_info(const void *s, char **d, int *bounds, bool *a);
void tao_ele_shape_struct_set_ele_id(void *struct_obj, const char *str_ptr, int str_len);
void tao_ele_shape_struct_get_shape_info(const void *s, char **d, int *bounds, bool *a);
void tao_ele_shape_struct_set_shape(void *struct_obj, const char *str_ptr, int str_len);
void tao_ele_shape_struct_get_color_info(const void *s, char **d, int *bounds, bool *a);
void tao_ele_shape_struct_set_color(void *struct_obj, const char *str_ptr, int str_len);
void tao_ele_shape_struct_get_size(const void *struct_obj, double *value_out);
void tao_ele_shape_struct_set_size(void *struct_obj, double value_in);
void tao_ele_shape_struct_get_label_info(const void *s, char **d, int *bounds, bool *a);
void tao_ele_shape_struct_set_label(void *struct_obj, const char *str_ptr, int str_len);
void tao_ele_shape_struct_get_draw(const void *struct_obj, bool *value_out);
void tao_ele_shape_struct_set_draw(void *struct_obj, bool value_in);
void tao_ele_shape_struct_get_multi(const void *struct_obj, bool *value_out);
void tao_ele_shape_struct_set_multi(void *struct_obj, bool value_in);
void tao_ele_shape_struct_get_line_width(const void *struct_obj, int *value_out);
void tao_ele_shape_struct_set_line_width(void *struct_obj, int value_in);
void tao_ele_shape_struct_get_offset(const void *struct_obj, double *value_out);
void tao_ele_shape_struct_set_offset(void *struct_obj, double value_in);
void tao_ele_shape_struct_get_ix_key(const void *struct_obj, int *value_out);
void tao_ele_shape_struct_set_ix_key(void *struct_obj, int value_in);
void tao_ele_shape_struct_get_name_ele_info(const void *s, char **d, int *bounds, bool *a);
void tao_ele_shape_struct_set_name_ele(void *struct_obj, const char *str_ptr, int str_len);

void tao_ele_shape_struct_get_uni_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_ele_pointer_struct_get_eles_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_ele_pointer_struct_get_n_loc(const void *struct_obj, int *value_out);
void tao_ele_pointer_struct_set_n_loc(void *struct_obj, int value_in);
void tao_curve_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_curve_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_curve_struct_get_data_source_info(const void *s, char **d, int *bounds, bool *a);
void tao_curve_struct_set_data_source(void *struct_obj, const char *str_ptr, int str_len);
void tao_curve_struct_get_data_index_info(const void *s, char **d, int *bounds, bool *a);
void tao_curve_struct_set_data_index(void *struct_obj, const char *str_ptr, int str_len);
void tao_curve_struct_get_data_type_x_info(const void *s, char **d, int *bounds, bool *a);
void tao_curve_struct_set_data_type_x(void *struct_obj, const char *str_ptr, int str_len);

void tao_curve_struct_get_data_type_info(const void *s, char **d, int *len, bool *is_alloc);

void tao_curve_struct_set_data_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_curve_struct_get_ele_ref_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_curve_struct_set_ele_ref_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_curve_struct_get_legend_text_info(const void *s, char **d, int *bounds, bool *a);
void tao_curve_struct_set_legend_text(void *struct_obj, const char *str_ptr, int str_len);
void tao_curve_struct_get_message_text_info(const void *s, char **d, int *bounds, bool *a);
void tao_curve_struct_set_message_text(void *struct_obj, const char *str_ptr, int str_len);
void tao_curve_struct_get_component_info(const void *s, char **d, int *bounds, bool *a);
void tao_curve_struct_set_component(void *struct_obj, const char *str_ptr, int str_len);
void tao_curve_struct_get_why_invalid_info(const void *s, char **d, int *bounds, bool *a);
void tao_curve_struct_set_why_invalid(void *struct_obj, const char *str_ptr, int str_len);
void tao_curve_struct_get_g(const void *struct_obj, void **ptr_out);
void tao_curve_struct_set_g(void *struct_obj, const void *src_ptr);
void tao_curve_struct_get_hist(const void *struct_obj, void **ptr_out);
void tao_curve_struct_set_hist(void *struct_obj, const void *src_ptr);
void tao_curve_struct_get_z_color(const void *struct_obj, void **ptr_out);
void tao_curve_struct_set_z_color(void *struct_obj, const void *src_ptr);
void tao_curve_struct_get_x_line_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_curve_struct_set_x_line(void *s, const void *d, const int *shape);
void tao_curve_struct_get_y_line_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_curve_struct_set_y_line(void *s, const void *d, const int *shape);
void tao_curve_struct_get_y2_line_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_curve_struct_set_y2_line(void *s, const void *d, const int *shape);
void tao_curve_struct_get_ix_line_info(const void *s, int **d, int *bounds, bool *is_alloc);
void tao_curve_struct_set_ix_line(void *s, const void *d, const int *shape);
void tao_curve_struct_get_x_symb_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_curve_struct_set_x_symb(void *s, const void *d, const int *shape);
void tao_curve_struct_get_y_symb_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_curve_struct_set_y_symb(void *s, const void *d, const int *shape);
void tao_curve_struct_get_z_symb_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_curve_struct_set_z_symb(void *s, const void *d, const int *shape);
void tao_curve_struct_get_err_symb_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_curve_struct_set_err_symb(void *s, const void *d, const int *shape);
void tao_curve_struct_get_symb_size_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_curve_struct_set_symb_size(void *s, const void *d, const int *shape);
void tao_curve_struct_get_ix_symb_info(const void *s, int **d, int *bounds, bool *is_alloc);
void tao_curve_struct_set_ix_symb(void *s, const void *d, const int *shape);
void tao_curve_struct_get_y_axis_scale_factor(const void *struct_obj, double *value_out);
void tao_curve_struct_set_y_axis_scale_factor(void *struct_obj, double value_in);
void tao_curve_struct_get_line(const void *struct_obj, void **ptr_out);
void tao_curve_struct_set_line(void *struct_obj, const void *src_ptr);
void tao_curve_struct_get_symbol(const void *struct_obj, void **ptr_out);
void tao_curve_struct_set_symbol(void *struct_obj, const void *src_ptr);
void tao_curve_struct_get_orbit(const void *struct_obj, void **ptr_out);
void tao_curve_struct_set_orbit(void *struct_obj, const void *src_ptr);
void tao_curve_struct_get_ix_universe(const void *struct_obj, int *value_out);
void tao_curve_struct_set_ix_universe(void *struct_obj, int value_in);
void tao_curve_struct_get_symbol_every(const void *struct_obj, int *value_out);
void tao_curve_struct_set_symbol_every(void *struct_obj, int value_in);
void tao_curve_struct_get_ix_branch(const void *struct_obj, int *value_out);
void tao_curve_struct_set_ix_branch(void *struct_obj, int value_in);
void tao_curve_struct_get_ix_bunch(const void *struct_obj, int *value_out);
void tao_curve_struct_set_ix_bunch(void *struct_obj, int value_in);
void tao_curve_struct_get_n_turn(const void *struct_obj, int *value_out);
void tao_curve_struct_set_n_turn(void *struct_obj, int value_in);
void tao_curve_struct_get_use_y2(const void *struct_obj, bool *value_out);
void tao_curve_struct_set_use_y2(void *struct_obj, bool value_in);
void tao_curve_struct_get_draw_line(const void *struct_obj, bool *value_out);
void tao_curve_struct_set_draw_line(void *struct_obj, bool value_in);
void tao_curve_struct_get_draw_symbols(const void *struct_obj, bool *value_out);
void tao_curve_struct_set_draw_symbols(void *struct_obj, bool value_in);
void tao_curve_struct_get_draw_symbol_index(const void *struct_obj, bool *value_out);
void tao_curve_struct_set_draw_symbol_index(void *struct_obj, bool value_in);
void tao_curve_struct_get_draw_error_bars(const void *struct_obj, bool *value_out);
void tao_curve_struct_set_draw_error_bars(void *struct_obj, bool value_in);
void tao_curve_struct_get_smooth_line_calc(const void *struct_obj, bool *value_out);
void tao_curve_struct_set_smooth_line_calc(void *struct_obj, bool value_in);
void tao_curve_struct_get_valid(const void *struct_obj, bool *value_out);
void tao_curve_struct_set_valid(void *struct_obj, bool value_in);
void tao_curve_color_struct_get_data_type_info(const void *s, char **d, int *bounds, bool *a);
void tao_curve_color_struct_set_data_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_curve_color_struct_get_is_on(const void *struct_obj, bool *value_out);
void tao_curve_color_struct_set_is_on(void *struct_obj, bool value_in);
void tao_curve_color_struct_get_min(const void *struct_obj, double *value_out);
void tao_curve_color_struct_set_min(void *struct_obj, double value_in);
void tao_curve_color_struct_get_max(const void *struct_obj, double *value_out);
void tao_curve_color_struct_set_max(void *struct_obj, double value_in);
void tao_curve_color_struct_get_autoscale(const void *struct_obj, bool *value_out);
void tao_curve_color_struct_set_autoscale(void *struct_obj, bool value_in);
void tao_curve_orbit_struct_get_x(const void *struct_obj, double *value_out);
void tao_curve_orbit_struct_set_x(void *struct_obj, double value_in);
void tao_curve_orbit_struct_get_y(const void *struct_obj, double *value_out);
void tao_curve_orbit_struct_set_y(void *struct_obj, double value_in);
void tao_curve_orbit_struct_get_t(const void *struct_obj, double *value_out);
void tao_curve_orbit_struct_set_t(void *struct_obj, double value_in);
void tao_histogram_struct_get_density_normalized(const void *struct_obj, bool *value_out);
void tao_histogram_struct_set_density_normalized(void *struct_obj, bool value_in);
void tao_histogram_struct_get_weight_by_charge(const void *struct_obj, bool *value_out);
void tao_histogram_struct_set_weight_by_charge(void *struct_obj, bool value_in);
void tao_histogram_struct_get_minimum(const void *struct_obj, double *value_out);
void tao_histogram_struct_set_minimum(void *struct_obj, double value_in);
void tao_histogram_struct_get_maximum(const void *struct_obj, double *value_out);
void tao_histogram_struct_set_maximum(void *struct_obj, double value_in);
void tao_histogram_struct_get_width(const void *struct_obj, double *value_out);
void tao_histogram_struct_set_width(void *struct_obj, double value_in);
void tao_histogram_struct_get_center(const void *struct_obj, double *value_out);
void tao_histogram_struct_set_center(void *struct_obj, double value_in);
void tao_histogram_struct_get_number(const void *struct_obj, int *value_out);
void tao_histogram_struct_set_number(void *struct_obj, int value_in);
void lat_ele_order1_struct_get_ix_branch(const void *struct_obj, int *value_out);
void lat_ele_order1_struct_set_ix_branch(void *struct_obj, int value_in);
void lat_ele_order1_struct_get_ix_order(const void *struct_obj, int *value_out);
void lat_ele_order1_struct_set_ix_order(void *struct_obj, int value_in);

void lat_ele_order_array_struct_get_ele_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lat_sigma_struct_get_mat_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void tao_lat_sigma_struct_set_mat(void *s, const void *d, const int *shape);
void tao_spin_ele_struct_get_dn_dpz(const void *struct_obj, void **ptr_out);
void tao_spin_ele_struct_set_dn_dpz(void *struct_obj, const void *src_ptr);
void tao_spin_ele_struct_get_orb_eigen_val_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void tao_spin_ele_struct_set_orb_eigen_val(void *s, const void *d, const int *shape);
void tao_spin_ele_struct_get_orb_eigen_vec_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void tao_spin_ele_struct_set_orb_eigen_vec(void *s, const void *d, const int *shape);
void tao_spin_ele_struct_get_spin_eigen_vec_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void tao_spin_ele_struct_set_spin_eigen_vec(void *s, const void *d, const int *shape);
void tao_spin_ele_struct_get_valid(const void *struct_obj, bool *value_out);
void tao_spin_ele_struct_set_valid(void *struct_obj, bool value_in);
void tao_plot_cache_struct_get_ele_to_s(const void *struct_obj, void **ptr_out);
void tao_plot_cache_struct_set_ele_to_s(void *struct_obj, const void *src_ptr);
void tao_plot_cache_struct_get_orbit(const void *struct_obj, void **ptr_out);
void tao_plot_cache_struct_set_orbit(void *struct_obj, const void *src_ptr);
void tao_plot_cache_struct_get_err(const void *struct_obj, bool *value_out);
void tao_plot_cache_struct_set_err(void *struct_obj, bool value_in);
void tao_spin_polarization_struct_get_tune(const void *struct_obj, double *value_out);
void tao_spin_polarization_struct_set_tune(void *struct_obj, double value_in);
void tao_spin_polarization_struct_get_pol_limit_st(const void *struct_obj, double *value_out);
void tao_spin_polarization_struct_set_pol_limit_st(void *struct_obj, double value_in);
void tao_spin_polarization_struct_get_pol_limit_dk(const void *struct_obj, double *value_out);
void tao_spin_polarization_struct_set_pol_limit_dk(void *struct_obj, double value_in);
void tao_spin_polarization_struct_get_pol_limit_dk_partial_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void tao_spin_polarization_struct_set_pol_limit_dk_partial(
    void *s,
    const void *d,
    const int *shape
);
void tao_spin_polarization_struct_get_pol_limit_dk_partial2_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void tao_spin_polarization_struct_set_pol_limit_dk_partial2(
    void *s,
    const void *d,
    const int *shape
);
void tao_spin_polarization_struct_get_pol_rate_bks(const void *struct_obj, double *value_out);
void tao_spin_polarization_struct_set_pol_rate_bks(void *struct_obj, double value_in);
void tao_spin_polarization_struct_get_depol_rate(const void *struct_obj, double *value_out);
void tao_spin_polarization_struct_set_depol_rate(void *struct_obj, double value_in);
void tao_spin_polarization_struct_get_depol_rate_partial_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void tao_spin_polarization_struct_set_depol_rate_partial(void *s, const void *d, const int *shape);
void tao_spin_polarization_struct_get_depol_rate_partial2_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void tao_spin_polarization_struct_set_depol_rate_partial2(void *s, const void *d, const int *shape);
void tao_spin_polarization_struct_get_integral_bn(const void *struct_obj, double *value_out);
void tao_spin_polarization_struct_set_integral_bn(void *struct_obj, double value_in);
void tao_spin_polarization_struct_get_integral_bdn(const void *struct_obj, double *value_out);
void tao_spin_polarization_struct_set_integral_bdn(void *struct_obj, double value_in);
void tao_spin_polarization_struct_get_integral_1ns(const void *struct_obj, double *value_out);
void tao_spin_polarization_struct_set_integral_1ns(void *struct_obj, double value_in);
void tao_spin_polarization_struct_get_integral_dn2(const void *struct_obj, double *value_out);
void tao_spin_polarization_struct_set_integral_dn2(void *struct_obj, double value_in);
void tao_spin_polarization_struct_get_valid(const void *struct_obj, bool *value_out);
void tao_spin_polarization_struct_set_valid(void *struct_obj, bool value_in);
void tao_spin_polarization_struct_get_q_1turn(const void *struct_obj, void **ptr_out);
void tao_spin_polarization_struct_set_q_1turn(void *struct_obj, const void *src_ptr);

void tao_spin_polarization_struct_get_q_ele_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lattice_branch_struct_get_tao_lat(const void *struct_obj, void **ptr_out);
void tao_lattice_branch_struct_set_tao_lat(void *struct_obj, const void *src_ptr);

void tao_lattice_branch_struct_get_lat_sigma_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lattice_branch_struct_get_spin_ele_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lattice_branch_struct_get_bunch_params_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lattice_branch_struct_get_bunch_params_comb_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lattice_branch_struct_get_orbit_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lattice_branch_struct_get_plot_cache_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lattice_branch_struct_get_spin(const void *struct_obj, void **ptr_out);
void tao_lattice_branch_struct_set_spin(void *struct_obj, const void *src_ptr);
void tao_lattice_branch_struct_get_srdt(const void *struct_obj, void **ptr_out);
void tao_lattice_branch_struct_set_srdt(void *struct_obj, const void *src_ptr);
void tao_lattice_branch_struct_get_orb0(const void *struct_obj, void **ptr_out);
void tao_lattice_branch_struct_set_orb0(void *struct_obj, const void *src_ptr);
void tao_lattice_branch_struct_get_modes_ri(const void *struct_obj, void **ptr_out);
void tao_lattice_branch_struct_set_modes_ri(void *struct_obj, const void *src_ptr);
void tao_lattice_branch_struct_get_modes_6d(const void *struct_obj, void **ptr_out);
void tao_lattice_branch_struct_set_modes_6d(void *struct_obj, const void *src_ptr);
void tao_lattice_branch_struct_get_ptc_normal_form(const void *struct_obj, void **ptr_out);
void tao_lattice_branch_struct_set_ptc_normal_form(void *struct_obj, const void *src_ptr);
void tao_lattice_branch_struct_get_bmad_normal_form(const void *struct_obj, void **ptr_out);
void tao_lattice_branch_struct_set_bmad_normal_form(void *struct_obj, const void *src_ptr);

void tao_lattice_branch_struct_get_high_E_orb_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lattice_branch_struct_get_low_E_orb_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lattice_branch_struct_get_taylor_save_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_lattice_branch_struct_get_cache_x_min(const void *struct_obj, double *value_out);
void tao_lattice_branch_struct_set_cache_x_min(void *struct_obj, double value_in);
void tao_lattice_branch_struct_get_cache_x_max(const void *struct_obj, double *value_out);
void tao_lattice_branch_struct_set_cache_x_max(void *struct_obj, double value_in);
void tao_lattice_branch_struct_get_comb_ds_save(const void *struct_obj, double *value_out);
void tao_lattice_branch_struct_set_comb_ds_save(void *struct_obj, double value_in);
void tao_lattice_branch_struct_get_ix_ref_taylor(const void *struct_obj, int *value_out);
void tao_lattice_branch_struct_set_ix_ref_taylor(void *struct_obj, int value_in);
void tao_lattice_branch_struct_get_ix_ele_taylor(const void *struct_obj, int *value_out);
void tao_lattice_branch_struct_set_ix_ele_taylor(void *struct_obj, int value_in);
void tao_lattice_branch_struct_get_track_state(const void *struct_obj, int *value_out);
void tao_lattice_branch_struct_set_track_state(void *struct_obj, int value_in);
void tao_lattice_branch_struct_get_cache_n_pts(const void *struct_obj, int *value_out);
void tao_lattice_branch_struct_set_cache_n_pts(void *struct_obj, int value_in);
void tao_lattice_branch_struct_get_ix_rad_int_cache(const void *struct_obj, int *value_out);
void tao_lattice_branch_struct_set_ix_rad_int_cache(void *struct_obj, int value_in);
void tao_lattice_branch_struct_get_has_open_match_element(const void *struct_obj, bool *value_out);
void tao_lattice_branch_struct_set_has_open_match_element(void *struct_obj, bool value_in);
void tao_lattice_branch_struct_get_plot_cache_valid(const void *struct_obj, bool *value_out);
void tao_lattice_branch_struct_set_plot_cache_valid(void *struct_obj, bool value_in);
void tao_lattice_branch_struct_get_spin_map_valid(const void *struct_obj, bool *value_out);
void tao_lattice_branch_struct_set_spin_map_valid(void *struct_obj, bool value_in);
void tao_lattice_branch_struct_get_twiss_valid(const void *struct_obj, bool *value_out);
void tao_lattice_branch_struct_set_twiss_valid(void *struct_obj, bool value_in);
void tao_lattice_branch_struct_get_mode_flip_here(const void *struct_obj, bool *value_out);
void tao_lattice_branch_struct_set_mode_flip_here(void *struct_obj, bool value_in);
void tao_lattice_branch_struct_get_chrom_calc_ok(const void *struct_obj, bool *value_out);
void tao_lattice_branch_struct_set_chrom_calc_ok(void *struct_obj, bool value_in);
void tao_lattice_branch_struct_get_rad_int_calc_ok(const void *struct_obj, bool *value_out);
void tao_lattice_branch_struct_set_rad_int_calc_ok(void *struct_obj, bool value_in);
void tao_lattice_branch_struct_get_emit_6d_calc_ok(const void *struct_obj, bool *value_out);
void tao_lattice_branch_struct_set_emit_6d_calc_ok(void *struct_obj, bool value_in);
void tao_lattice_branch_struct_get_sigma_track_ok(const void *struct_obj, bool *value_out);
void tao_lattice_branch_struct_set_sigma_track_ok(void *struct_obj, bool value_in);
void tao_model_element_struct_get_beam(const void *struct_obj, void **ptr_out);
void tao_model_element_struct_set_beam(void *struct_obj, const void *src_ptr);
void tao_model_element_struct_get_save_beam_internally(const void *struct_obj, bool *value_out);
void tao_model_element_struct_set_save_beam_internally(void *struct_obj, bool value_in);
void tao_model_element_struct_get_save_beam_to_file(const void *struct_obj, bool *value_out);
void tao_model_element_struct_set_save_beam_to_file(void *struct_obj, bool value_in);
void tao_beam_branch_struct_get_beam_at_start(const void *struct_obj, void **ptr_out);
void tao_beam_branch_struct_set_beam_at_start(void *struct_obj, const void *src_ptr);
void tao_beam_branch_struct_get_beam_init(const void *struct_obj, void **ptr_out);
void tao_beam_branch_struct_set_beam_init(void *struct_obj, const void *src_ptr);
void tao_beam_branch_struct_get_beam_init_used(const void *struct_obj, void **ptr_out);
void tao_beam_branch_struct_set_beam_init_used(void *struct_obj, const void *src_ptr);
void tao_beam_branch_struct_get_init_starting_distribution(const void *struct_obj, bool *value_out);
void tao_beam_branch_struct_set_init_starting_distribution(void *struct_obj, bool value_in);
void tao_beam_branch_struct_get_track_start_info(const void *s, char **d, int *bounds, bool *a);
void tao_beam_branch_struct_set_track_start(void *struct_obj, const char *str_ptr, int str_len);
void tao_beam_branch_struct_get_track_end_info(const void *s, char **d, int *bounds, bool *a);
void tao_beam_branch_struct_set_track_end(void *struct_obj, const char *str_ptr, int str_len);
void tao_beam_branch_struct_get_ix_branch(const void *struct_obj, int *value_out);
void tao_beam_branch_struct_set_ix_branch(void *struct_obj, int value_in);
void tao_beam_branch_struct_get_ix_track_start(const void *struct_obj, int *value_out);
void tao_beam_branch_struct_set_ix_track_start(void *struct_obj, int value_in);
void tao_beam_branch_struct_get_ix_track_end(const void *struct_obj, int *value_out);
void tao_beam_branch_struct_set_ix_track_end(void *struct_obj, int value_in);
void tao_d1_data_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_d1_data_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_d1_data_struct_get_d2(const void *struct_obj, void **ptr_out);
void tao_d1_data_struct_set_d2(void *struct_obj, const void *src_ptr);

void tao_d1_data_struct_get_d_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_d2_data_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_d2_data_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_d2_data_struct_get_data_file_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_d2_data_struct_set_data_file_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_d2_data_struct_get_ref_file_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_d2_data_struct_set_ref_file_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_d2_data_struct_get_data_date_info(const void *s, char **d, int *bounds, bool *a);
void tao_d2_data_struct_set_data_date(void *struct_obj, const char *str_ptr, int str_len);
void tao_d2_data_struct_get_ref_date_info(const void *s, char **d, int *bounds, bool *a);
void tao_d2_data_struct_set_ref_date(void *struct_obj, const char *str_ptr, int str_len);

void tao_d2_data_struct_get_descrip_info(
    const void *s,
    char **d,
    int *bounds, // [lower, upper]
    int *str_len,
    bool *is_alloc
);

void tao_d2_data_struct_get_d1_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_d2_data_struct_get_ix_universe(const void *struct_obj, int *value_out);
void tao_d2_data_struct_set_ix_universe(void *struct_obj, int value_in);
void tao_d2_data_struct_get_ix_d2_data(const void *struct_obj, int *value_out);
void tao_d2_data_struct_set_ix_d2_data(void *struct_obj, int value_in);
void tao_d2_data_struct_get_ix_ref(const void *struct_obj, int *value_out);
void tao_d2_data_struct_set_ix_ref(void *struct_obj, int value_in);
void tao_d2_data_struct_get_data_read_in(const void *struct_obj, bool *value_out);
void tao_d2_data_struct_set_data_read_in(void *struct_obj, bool value_in);
void tao_d2_data_struct_get_ref_read_in(const void *struct_obj, bool *value_out);
void tao_d2_data_struct_set_ref_read_in(void *struct_obj, bool value_in);
void tao_data_var_component_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_data_var_component_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_data_var_component_struct_get_sign(const void *struct_obj, double *value_out);
void tao_data_var_component_struct_set_sign(void *struct_obj, double value_in);
void tao_graph_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_graph_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_graph_struct_get_type_info(const void *s, char **d, int *bounds, bool *a);
void tao_graph_struct_set_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_graph_struct_get_title_info(const void *s, char **d, int *bounds, bool *a);
void tao_graph_struct_set_title(void *struct_obj, const char *str_ptr, int str_len);
void tao_graph_struct_get_title_suffix_info(const void *s, char **d, int *bounds, bool *a);
void tao_graph_struct_set_title_suffix(void *struct_obj, const char *str_ptr, int str_len);

void tao_graph_struct_get_text_legend_info(
    const void *s,
    char **d,
    int *bounds, // [lower, upper]
    int *str_len,
    bool *is_alloc
);

void tao_graph_struct_get_text_legend_out_info(
    const void *s,
    char **d,
    int *bounds, // [lower, upper]
    int *str_len,
    bool *is_alloc
);

void tao_graph_struct_get_why_invalid_info(const void *s, char **d, int *bounds, bool *a);
void tao_graph_struct_set_why_invalid(void *struct_obj, const char *str_ptr, int str_len);

void tao_graph_struct_get_curve_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_graph_struct_get_p(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_p(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_floor_plan(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_floor_plan(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_text_legend_origin(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_text_legend_origin(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_curve_legend_origin(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_curve_legend_origin(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_curve_legend(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_curve_legend(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_x(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_x(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_y(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_y(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_x2(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_x2(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_y2(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_y2(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_margin(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_margin(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_scale_margin(const void *struct_obj, void **ptr_out);
void tao_graph_struct_set_scale_margin(void *struct_obj, const void *src_ptr);
void tao_graph_struct_get_x_axis_scale_factor(const void *struct_obj, double *value_out);
void tao_graph_struct_set_x_axis_scale_factor(void *struct_obj, double value_in);
void tao_graph_struct_get_symbol_size_scale(const void *struct_obj, double *value_out);
void tao_graph_struct_set_symbol_size_scale(void *struct_obj, double value_in);
void tao_graph_struct_get_box_info(const void *s, int **d, int *bounds, bool *is_alloc);
void tao_graph_struct_set_box(void *s, const void *d, const int *shape);
void tao_graph_struct_get_ix_branch(const void *struct_obj, int *value_out);
void tao_graph_struct_set_ix_branch(void *struct_obj, int value_in);
void tao_graph_struct_get_ix_universe(const void *struct_obj, int *value_out);
void tao_graph_struct_set_ix_universe(void *struct_obj, int value_in);
void tao_graph_struct_get_clip(const void *struct_obj, bool *value_out);
void tao_graph_struct_set_clip(void *struct_obj, bool value_in);
void tao_graph_struct_get_y2_mirrors_y(const void *struct_obj, bool *value_out);
void tao_graph_struct_set_y2_mirrors_y(void *struct_obj, bool value_in);
void tao_graph_struct_get_limited(const void *struct_obj, bool *value_out);
void tao_graph_struct_set_limited(void *struct_obj, bool value_in);
void tao_graph_struct_get_draw_axes(const void *struct_obj, bool *value_out);
void tao_graph_struct_set_draw_axes(void *struct_obj, bool value_in);
void tao_graph_struct_get_draw_curve_legend(const void *struct_obj, bool *value_out);
void tao_graph_struct_set_draw_curve_legend(void *struct_obj, bool value_in);
void tao_graph_struct_get_draw_grid(const void *struct_obj, bool *value_out);
void tao_graph_struct_set_draw_grid(void *struct_obj, bool value_in);
void tao_graph_struct_get_draw_title(const void *struct_obj, bool *value_out);
void tao_graph_struct_set_draw_title(void *struct_obj, bool value_in);
void tao_graph_struct_get_draw_only_good_user_data_or_vars(const void *struct_obj, bool *value_out);
void tao_graph_struct_set_draw_only_good_user_data_or_vars(void *struct_obj, bool value_in);
void tao_graph_struct_get_allow_wrap_around(const void *struct_obj, bool *value_out);
void tao_graph_struct_set_allow_wrap_around(void *struct_obj, bool value_in);
void tao_graph_struct_get_is_valid(const void *struct_obj, bool *value_out);
void tao_graph_struct_set_is_valid(void *struct_obj, bool value_in);
void tao_plot_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_plot_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_plot_struct_get_description_info(const void *s, char **d, int *bounds, bool *a);
void tao_plot_struct_set_description(void *struct_obj, const char *str_ptr, int str_len);

void tao_plot_struct_get_graph_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_plot_struct_get_r(const void *struct_obj, void **ptr_out);
void tao_plot_struct_set_r(void *struct_obj, const void *src_ptr);
void tao_plot_struct_get_ix_plot(const void *struct_obj, int *value_out);
void tao_plot_struct_set_ix_plot(void *struct_obj, int value_in);
void tao_plot_struct_get_n_curve_pts(const void *struct_obj, int *value_out);
void tao_plot_struct_set_n_curve_pts(void *struct_obj, int value_in);
void tao_plot_struct_get_type_info(const void *s, char **d, int *bounds, bool *a);
void tao_plot_struct_set_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_plot_struct_get_x_axis_type_info(const void *s, char **d, int *bounds, bool *a);
void tao_plot_struct_set_x_axis_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_plot_struct_get_autoscale_x(const void *struct_obj, bool *value_out);
void tao_plot_struct_set_autoscale_x(void *struct_obj, bool value_in);
void tao_plot_struct_get_autoscale_y(const void *struct_obj, bool *value_out);
void tao_plot_struct_set_autoscale_y(void *struct_obj, bool value_in);
void tao_plot_struct_get_autoscale_gang_x(const void *struct_obj, bool *value_out);
void tao_plot_struct_set_autoscale_gang_x(void *struct_obj, bool value_in);
void tao_plot_struct_get_autoscale_gang_y(const void *struct_obj, bool *value_out);
void tao_plot_struct_set_autoscale_gang_y(void *struct_obj, bool value_in);
void tao_plot_struct_get_list_with_show_plot_command(const void *struct_obj, bool *value_out);
void tao_plot_struct_set_list_with_show_plot_command(void *struct_obj, bool value_in);
void tao_plot_struct_get_phantom(const void *struct_obj, bool *value_out);
void tao_plot_struct_set_phantom(void *struct_obj, bool value_in);
void tao_plot_struct_get_default_plot(const void *struct_obj, bool *value_out);
void tao_plot_struct_set_default_plot(void *struct_obj, bool value_in);
void tao_plot_region_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_plot_region_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_plot_region_struct_get_plot(const void *struct_obj, void **ptr_out);
void tao_plot_region_struct_set_plot(void *struct_obj, const void *src_ptr);
void tao_plot_region_struct_get_location_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void tao_plot_region_struct_set_location(void *s, const void *d, const int *shape);
void tao_plot_region_struct_get_visible(const void *struct_obj, bool *value_out);
void tao_plot_region_struct_set_visible(void *struct_obj, bool value_in);
void tao_plot_region_struct_get_list_with_show_plot_command(
    const void *struct_obj,
    bool *value_out
);
void tao_plot_region_struct_set_list_with_show_plot_command(void *struct_obj, bool value_in);
void tao_plot_region_struct_get_setup_done(const void *struct_obj, bool *value_out);
void tao_plot_region_struct_set_setup_done(void *struct_obj, bool value_in);
void tao_universe_pointer_struct_get_u(const void *struct_obj, void **ptr_out);
void tao_universe_pointer_struct_set_u(void *struct_obj, const void *src_ptr);
void tao_super_universe_struct_get_global(const void *struct_obj, void **ptr_out);
void tao_super_universe_struct_set_global(void *struct_obj, const void *src_ptr);
void tao_super_universe_struct_get_init(const void *struct_obj, void **ptr_out);
void tao_super_universe_struct_set_init(void *struct_obj, const void *src_ptr);
void tao_super_universe_struct_get_com(const void *struct_obj, void **ptr_out);
void tao_super_universe_struct_set_com(void *struct_obj, const void *src_ptr);
void tao_super_universe_struct_get_plot_page(const void *struct_obj, void **ptr_out);
void tao_super_universe_struct_set_plot_page(void *struct_obj, const void *src_ptr);

void tao_super_universe_struct_get_v1_var_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_super_universe_struct_get_var_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_super_universe_struct_get_u_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_super_universe_struct_get_key_info(const void *s, int **d, int *bounds, bool *is_alloc);
void tao_super_universe_struct_set_key(void *s, const void *d, const int *shape);
void tao_super_universe_struct_get_building_wall(const void *struct_obj, void **ptr_out);
void tao_super_universe_struct_set_building_wall(void *struct_obj, const void *src_ptr);
void tao_super_universe_struct_get_wave(const void *struct_obj, void **ptr_out);
void tao_super_universe_struct_set_wave(void *struct_obj, const void *src_ptr);
void tao_super_universe_struct_get_n_var_used(const void *struct_obj, int *value_out);
void tao_super_universe_struct_set_n_var_used(void *struct_obj, int value_in);
void tao_super_universe_struct_get_n_v1_var_used(const void *struct_obj, int *value_out);
void tao_super_universe_struct_set_n_v1_var_used(void *struct_obj, int value_in);

void tao_super_universe_struct_get_history_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_super_universe_struct_get_initialized(const void *struct_obj, bool *value_out);
void tao_super_universe_struct_set_initialized(void *struct_obj, bool value_in);
void tao_var_struct_get_ele_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_var_struct_set_ele_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_var_struct_get_attrib_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_var_struct_set_attrib_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_var_struct_get_id_info(const void *s, char **d, int *bounds, bool *a);
void tao_var_struct_set_id(void *struct_obj, const char *str_ptr, int str_len);

void tao_var_struct_get_slave_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_var_struct_get_ix_v1(const void *struct_obj, int *value_out);
void tao_var_struct_set_ix_v1(void *struct_obj, int value_in);
void tao_var_struct_get_ix_var(const void *struct_obj, int *value_out);
void tao_var_struct_set_ix_var(void *struct_obj, int value_in);
void tao_var_struct_get_ix_dvar(const void *struct_obj, int *value_out);
void tao_var_struct_set_ix_dvar(void *struct_obj, int value_in);
void tao_var_struct_get_ix_attrib(const void *struct_obj, int *value_out);
void tao_var_struct_set_ix_attrib(void *struct_obj, int value_in);
void tao_var_struct_get_ix_key_table(const void *struct_obj, int *value_out);
void tao_var_struct_set_ix_key_table(void *struct_obj, int value_in);
void tao_var_struct_get_model_value(const void *struct_obj, double **ptr_out);
void tao_var_struct_set_model_value(void *struct_obj, double value_in);
void tao_var_struct_get_base_value(const void *struct_obj, double **ptr_out);
void tao_var_struct_set_base_value(void *struct_obj, double value_in);
void tao_var_struct_get_design_value(const void *struct_obj, double *value_out);
void tao_var_struct_set_design_value(void *struct_obj, double value_in);
void tao_var_struct_get_scratch_value(const void *struct_obj, double *value_out);
void tao_var_struct_set_scratch_value(void *struct_obj, double value_in);
void tao_var_struct_get_old_value(const void *struct_obj, double *value_out);
void tao_var_struct_set_old_value(void *struct_obj, double value_in);
void tao_var_struct_get_meas_value(const void *struct_obj, double *value_out);
void tao_var_struct_set_meas_value(void *struct_obj, double value_in);
void tao_var_struct_get_ref_value(const void *struct_obj, double *value_out);
void tao_var_struct_set_ref_value(void *struct_obj, double value_in);
void tao_var_struct_get_correction_value(const void *struct_obj, double *value_out);
void tao_var_struct_set_correction_value(void *struct_obj, double value_in);
void tao_var_struct_get_high_lim(const void *struct_obj, double *value_out);
void tao_var_struct_set_high_lim(void *struct_obj, double value_in);
void tao_var_struct_get_low_lim(const void *struct_obj, double *value_out);
void tao_var_struct_set_low_lim(void *struct_obj, double value_in);
void tao_var_struct_get_step(const void *struct_obj, double *value_out);
void tao_var_struct_set_step(void *struct_obj, double value_in);
void tao_var_struct_get_weight(const void *struct_obj, double *value_out);
void tao_var_struct_set_weight(void *struct_obj, double value_in);
void tao_var_struct_get_delta_merit(const void *struct_obj, double *value_out);
void tao_var_struct_set_delta_merit(void *struct_obj, double value_in);
void tao_var_struct_get_merit(const void *struct_obj, double *value_out);
void tao_var_struct_set_merit(void *struct_obj, double value_in);
void tao_var_struct_get_dMerit_dVar(const void *struct_obj, double *value_out);
void tao_var_struct_set_dMerit_dVar(void *struct_obj, double value_in);
void tao_var_struct_get_key_val0(const void *struct_obj, double *value_out);
void tao_var_struct_set_key_val0(void *struct_obj, double value_in);
void tao_var_struct_get_key_delta(const void *struct_obj, double *value_out);
void tao_var_struct_set_key_delta(void *struct_obj, double value_in);
void tao_var_struct_get_s(const void *struct_obj, double *value_out);
void tao_var_struct_set_s(void *struct_obj, double value_in);
void tao_var_struct_get_extend_val(const void *struct_obj, double *value_out);
void tao_var_struct_set_extend_val(void *struct_obj, double value_in);
void tao_var_struct_get_merit_type_info(const void *s, char **d, int *bounds, bool *a);
void tao_var_struct_set_merit_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_var_struct_get_exists(const void *struct_obj, bool *value_out);
void tao_var_struct_set_exists(void *struct_obj, bool value_in);
void tao_var_struct_get_good_var(const void *struct_obj, bool *value_out);
void tao_var_struct_set_good_var(void *struct_obj, bool value_in);
void tao_var_struct_get_good_user(const void *struct_obj, bool *value_out);
void tao_var_struct_set_good_user(void *struct_obj, bool value_in);
void tao_var_struct_get_good_opt(const void *struct_obj, bool *value_out);
void tao_var_struct_set_good_opt(void *struct_obj, bool value_in);
void tao_var_struct_get_good_plot(const void *struct_obj, bool *value_out);
void tao_var_struct_set_good_plot(void *struct_obj, bool value_in);
void tao_var_struct_get_useit_opt(const void *struct_obj, bool *value_out);
void tao_var_struct_set_useit_opt(void *struct_obj, bool value_in);
void tao_var_struct_get_useit_plot(const void *struct_obj, bool *value_out);
void tao_var_struct_set_useit_plot(void *struct_obj, bool value_in);
void tao_var_struct_get_key_bound(const void *struct_obj, bool *value_out);
void tao_var_struct_set_key_bound(void *struct_obj, bool value_in);
void tao_var_struct_get_v1(const void *struct_obj, void **ptr_out);
void tao_var_struct_set_v1(void *struct_obj, const void *src_ptr);
void tao_var_slave_struct_get_ix_uni(const void *struct_obj, int *value_out);
void tao_var_slave_struct_set_ix_uni(void *struct_obj, int value_in);
void tao_var_slave_struct_get_ix_branch(const void *struct_obj, int *value_out);
void tao_var_slave_struct_set_ix_branch(void *struct_obj, int value_in);
void tao_var_slave_struct_get_ix_ele(const void *struct_obj, int *value_out);
void tao_var_slave_struct_set_ix_ele(void *struct_obj, int value_in);
void tao_var_slave_struct_get_model_value(const void *struct_obj, double **ptr_out);
void tao_var_slave_struct_set_model_value(void *struct_obj, double value_in);
void tao_var_slave_struct_get_base_value(const void *struct_obj, double **ptr_out);
void tao_var_slave_struct_set_base_value(void *struct_obj, double value_in);
void tao_lattice_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_lattice_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_lattice_struct_get_lat(const void *struct_obj, void **ptr_out);
void tao_lattice_struct_set_lat(void *struct_obj, const void *src_ptr);
void tao_lattice_struct_get_high_E_lat(const void *struct_obj, void **ptr_out);
void tao_lattice_struct_set_high_E_lat(void *struct_obj, const void *src_ptr);
void tao_lattice_struct_get_low_E_lat(const void *struct_obj, void **ptr_out);
void tao_lattice_struct_set_low_E_lat(void *struct_obj, const void *src_ptr);
void tao_lattice_struct_get_rad_int_by_ele_ri(const void *struct_obj, void **ptr_out);
void tao_lattice_struct_set_rad_int_by_ele_ri(void *struct_obj, const void *src_ptr);
void tao_lattice_struct_get_rad_int_by_ele_6d(const void *struct_obj, void **ptr_out);
void tao_lattice_struct_set_rad_int_by_ele_6d(void *struct_obj, const void *src_ptr);

void tao_lattice_struct_get_tao_branch_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_beam_uni_struct_get_saved_at_info(const void *s, char **d, int *bounds, bool *a);
void tao_beam_uni_struct_set_saved_at(void *struct_obj, const char *str_ptr, int str_len);
void tao_beam_uni_struct_get_dump_file_info(const void *s, char **d, int *bounds, bool *a);
void tao_beam_uni_struct_set_dump_file(void *struct_obj, const char *str_ptr, int str_len);
void tao_beam_uni_struct_get_dump_at_info(const void *s, char **d, int *bounds, bool *a);
void tao_beam_uni_struct_set_dump_at(void *struct_obj, const char *str_ptr, int str_len);
void tao_beam_uni_struct_get_track_beam_in_universe(const void *struct_obj, bool *value_out);
void tao_beam_uni_struct_set_track_beam_in_universe(void *struct_obj, bool value_in);
void tao_beam_uni_struct_get_always_reinit(const void *struct_obj, bool *value_out);
void tao_beam_uni_struct_set_always_reinit(void *struct_obj, bool value_in);
void tao_dynamic_aperture_struct_get_param(const void *struct_obj, void **ptr_out);
void tao_dynamic_aperture_struct_set_param(void *struct_obj, const void *src_ptr);

void tao_dynamic_aperture_struct_get_scan_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_dynamic_aperture_struct_get_pz_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void tao_dynamic_aperture_struct_set_pz(void *s, const void *d, const int *shape);
void tao_dynamic_aperture_struct_get_ellipse_scale(const void *struct_obj, double *value_out);
void tao_dynamic_aperture_struct_set_ellipse_scale(void *struct_obj, double value_in);
void tao_dynamic_aperture_struct_get_a_emit(const void *struct_obj, double *value_out);
void tao_dynamic_aperture_struct_set_a_emit(void *struct_obj, double value_in);
void tao_dynamic_aperture_struct_get_b_emit(const void *struct_obj, double *value_out);
void tao_dynamic_aperture_struct_set_b_emit(void *struct_obj, double value_in);

void tao_model_branch_struct_get_ele_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_model_branch_struct_get_beam(const void *struct_obj, void **ptr_out);
void tao_model_branch_struct_set_beam(void *struct_obj, const void *src_ptr);
void tao_spin_map_struct_get_valid(const void *struct_obj, bool *value_out);
void tao_spin_map_struct_set_valid(void *struct_obj, bool value_in);
void tao_spin_map_struct_get_map1(const void *struct_obj, void **ptr_out);
void tao_spin_map_struct_set_map1(void *struct_obj, const void *src_ptr);
void tao_spin_map_struct_get_axis_input(const void *struct_obj, void **ptr_out);
void tao_spin_map_struct_set_axis_input(void *struct_obj, const void *src_ptr);
void tao_spin_map_struct_get_axis0(const void *struct_obj, void **ptr_out);
void tao_spin_map_struct_set_axis0(void *struct_obj, const void *src_ptr);
void tao_spin_map_struct_get_axis1(const void *struct_obj, void **ptr_out);
void tao_spin_map_struct_set_axis1(void *struct_obj, const void *src_ptr);
void tao_spin_map_struct_get_ix_ele(const void *struct_obj, int *value_out);
void tao_spin_map_struct_set_ix_ele(void *struct_obj, int value_in);
void tao_spin_map_struct_get_ix_ref(const void *struct_obj, int *value_out);
void tao_spin_map_struct_set_ix_ref(void *struct_obj, int value_in);
void tao_spin_map_struct_get_ix_uni(const void *struct_obj, int *value_out);
void tao_spin_map_struct_set_ix_uni(void *struct_obj, int value_in);
void tao_spin_map_struct_get_ix_branch(const void *struct_obj, int *value_out);
void tao_spin_map_struct_set_ix_branch(void *struct_obj, int value_in);
void tao_spin_map_struct_get_mat8_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void tao_spin_map_struct_set_mat8(void *s, const void *d, const int *shape);
void tao_data_struct_get_ele_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_data_struct_set_ele_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_data_struct_get_ele_start_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_data_struct_set_ele_start_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_data_struct_get_ele_ref_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_data_struct_set_ele_ref_name(void *struct_obj, const char *str_ptr, int str_len);

void tao_data_struct_get_data_type_info(const void *s, char **d, int *len, bool *is_alloc);

void tao_data_struct_set_data_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_data_struct_get_merit_type_info(const void *s, char **d, int *bounds, bool *a);
void tao_data_struct_set_merit_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_data_struct_get_id_info(const void *s, char **d, int *bounds, bool *a);
void tao_data_struct_set_id(void *struct_obj, const char *str_ptr, int str_len);
void tao_data_struct_get_data_source_info(const void *s, char **d, int *bounds, bool *a);
void tao_data_struct_set_data_source(void *struct_obj, const char *str_ptr, int str_len);
void tao_data_struct_get_why_invalid_info(const void *s, char **d, int *bounds, bool *a);
void tao_data_struct_set_why_invalid(void *struct_obj, const char *str_ptr, int str_len);
void tao_data_struct_get_ix_uni(const void *struct_obj, int *value_out);
void tao_data_struct_set_ix_uni(void *struct_obj, int value_in);
void tao_data_struct_get_ix_bunch(const void *struct_obj, int *value_out);
void tao_data_struct_set_ix_bunch(void *struct_obj, int value_in);
void tao_data_struct_get_ix_branch(const void *struct_obj, int *value_out);
void tao_data_struct_set_ix_branch(void *struct_obj, int value_in);
void tao_data_struct_get_ix_ele(const void *struct_obj, int *value_out);
void tao_data_struct_set_ix_ele(void *struct_obj, int value_in);
void tao_data_struct_get_ix_ele_start(const void *struct_obj, int *value_out);
void tao_data_struct_set_ix_ele_start(void *struct_obj, int value_in);
void tao_data_struct_get_ix_ele_ref(const void *struct_obj, int *value_out);
void tao_data_struct_set_ix_ele_ref(void *struct_obj, int value_in);
void tao_data_struct_get_ix_ele_merit(const void *struct_obj, int *value_out);
void tao_data_struct_set_ix_ele_merit(void *struct_obj, int value_in);
void tao_data_struct_get_ix_d1(const void *struct_obj, int *value_out);
void tao_data_struct_set_ix_d1(void *struct_obj, int value_in);
void tao_data_struct_get_ix_data(const void *struct_obj, int *value_out);
void tao_data_struct_set_ix_data(void *struct_obj, int value_in);
void tao_data_struct_get_ix_dModel(const void *struct_obj, int *value_out);
void tao_data_struct_set_ix_dModel(void *struct_obj, int value_in);
void tao_data_struct_get_eval_point(const void *struct_obj, int *value_out);
void tao_data_struct_set_eval_point(void *struct_obj, int value_in);
void tao_data_struct_get_meas_value(const void *struct_obj, double *value_out);
void tao_data_struct_set_meas_value(void *struct_obj, double value_in);
void tao_data_struct_get_ref_value(const void *struct_obj, double *value_out);
void tao_data_struct_set_ref_value(void *struct_obj, double value_in);
void tao_data_struct_get_model_value(const void *struct_obj, double *value_out);
void tao_data_struct_set_model_value(void *struct_obj, double value_in);
void tao_data_struct_get_design_value(const void *struct_obj, double *value_out);
void tao_data_struct_set_design_value(void *struct_obj, double value_in);
void tao_data_struct_get_old_value(const void *struct_obj, double *value_out);
void tao_data_struct_set_old_value(void *struct_obj, double value_in);
void tao_data_struct_get_base_value(const void *struct_obj, double *value_out);
void tao_data_struct_set_base_value(void *struct_obj, double value_in);
void tao_data_struct_get_error_rms(const void *struct_obj, double *value_out);
void tao_data_struct_set_error_rms(void *struct_obj, double value_in);
void tao_data_struct_get_delta_merit(const void *struct_obj, double *value_out);
void tao_data_struct_set_delta_merit(void *struct_obj, double value_in);
void tao_data_struct_get_weight(const void *struct_obj, double *value_out);
void tao_data_struct_set_weight(void *struct_obj, double value_in);
void tao_data_struct_get_invalid_value(const void *struct_obj, double *value_out);
void tao_data_struct_set_invalid_value(void *struct_obj, double value_in);
void tao_data_struct_get_merit(const void *struct_obj, double *value_out);
void tao_data_struct_set_merit(void *struct_obj, double value_in);
void tao_data_struct_get_s(const void *struct_obj, double *value_out);
void tao_data_struct_set_s(void *struct_obj, double value_in);
void tao_data_struct_get_s_offset(const void *struct_obj, double *value_out);
void tao_data_struct_set_s_offset(void *struct_obj, double value_in);
void tao_data_struct_get_ref_s_offset(const void *struct_obj, double *value_out);
void tao_data_struct_set_ref_s_offset(void *struct_obj, double value_in);
void tao_data_struct_get_err_message_printed(const void *struct_obj, bool *value_out);
void tao_data_struct_set_err_message_printed(void *struct_obj, bool value_in);
void tao_data_struct_get_exists(const void *struct_obj, bool *value_out);
void tao_data_struct_set_exists(void *struct_obj, bool value_in);
void tao_data_struct_get_good_model(const void *struct_obj, bool *value_out);
void tao_data_struct_set_good_model(void *struct_obj, bool value_in);
void tao_data_struct_get_good_base(const void *struct_obj, bool *value_out);
void tao_data_struct_set_good_base(void *struct_obj, bool value_in);
void tao_data_struct_get_good_design(const void *struct_obj, bool *value_out);
void tao_data_struct_set_good_design(void *struct_obj, bool value_in);
void tao_data_struct_get_good_meas(const void *struct_obj, bool *value_out);
void tao_data_struct_set_good_meas(void *struct_obj, bool value_in);
void tao_data_struct_get_good_ref(const void *struct_obj, bool *value_out);
void tao_data_struct_set_good_ref(void *struct_obj, bool value_in);
void tao_data_struct_get_good_user(const void *struct_obj, bool *value_out);
void tao_data_struct_set_good_user(void *struct_obj, bool value_in);
void tao_data_struct_get_good_opt(const void *struct_obj, bool *value_out);
void tao_data_struct_set_good_opt(void *struct_obj, bool value_in);
void tao_data_struct_get_good_plot(const void *struct_obj, bool *value_out);
void tao_data_struct_set_good_plot(void *struct_obj, bool value_in);
void tao_data_struct_get_useit_plot(const void *struct_obj, bool *value_out);
void tao_data_struct_set_useit_plot(void *struct_obj, bool value_in);
void tao_data_struct_get_useit_opt(const void *struct_obj, bool *value_out);
void tao_data_struct_set_useit_opt(void *struct_obj, bool value_in);
void tao_data_struct_get_spin_map(const void *struct_obj, void **ptr_out);
void tao_data_struct_set_spin_map(void *struct_obj, const void *src_ptr);
void tao_data_struct_get_d1(const void *struct_obj, void **ptr_out);
void tao_data_struct_set_d1(void *struct_obj, const void *src_ptr);
void tao_ping_scale_struct_get_a_mode_meas(const void *struct_obj, double *value_out);
void tao_ping_scale_struct_set_a_mode_meas(void *struct_obj, double value_in);
void tao_ping_scale_struct_get_a_mode_ref(const void *struct_obj, double *value_out);
void tao_ping_scale_struct_set_a_mode_ref(void *struct_obj, double value_in);
void tao_ping_scale_struct_get_b_mode_meas(const void *struct_obj, double *value_out);
void tao_ping_scale_struct_set_b_mode_meas(void *struct_obj, double value_in);
void tao_ping_scale_struct_get_b_mode_ref(const void *struct_obj, double *value_out);
void tao_ping_scale_struct_set_b_mode_ref(void *struct_obj, double value_in);
void tao_universe_calc_struct_get_srdt_for_data(const void *struct_obj, int *value_out);
void tao_universe_calc_struct_set_srdt_for_data(void *struct_obj, int value_in);
void tao_universe_calc_struct_get_rad_int_for_data(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_rad_int_for_data(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_rad_int_for_plotting(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_rad_int_for_plotting(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_chrom_for_data(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_chrom_for_data(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_chrom_for_plotting(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_chrom_for_plotting(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_lat_sigma_for_data(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_lat_sigma_for_data(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_lat_sigma_for_plotting(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_lat_sigma_for_plotting(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_dynamic_aperture(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_dynamic_aperture(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_one_turn_map(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_one_turn_map(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_lattice(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_lattice(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_twiss(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_twiss(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_track(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_track(void *struct_obj, bool value_in);
void tao_universe_calc_struct_get_spin_matrices(const void *struct_obj, bool *value_out);
void tao_universe_calc_struct_set_spin_matrices(void *struct_obj, bool value_in);

void lat_ele_order_struct_get_branch_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_expression_info_struct_get_good(const void *struct_obj, bool *value_out);
void tao_expression_info_struct_set_good(void *struct_obj, bool value_in);
void tao_expression_info_struct_get_ele(const void *struct_obj, void **ptr_out);
void tao_expression_info_struct_set_ele(void *struct_obj, const void *src_ptr);
void tao_expression_info_struct_get_s(const void *struct_obj, double *value_out);
void tao_expression_info_struct_set_s(void *struct_obj, double value_in);
void tao_eval_node_struct_get_type(const void *struct_obj, int *value_out);
void tao_eval_node_struct_set_type(void *struct_obj, int value_in);
void tao_eval_node_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_eval_node_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_eval_node_struct_get_scale(const void *struct_obj, double *value_out);
void tao_eval_node_struct_set_scale(void *struct_obj, double value_in);
void tao_eval_node_struct_get_value_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_eval_node_struct_set_value(void *s, const void *d, const int *shape);

void tao_eval_node_struct_get_info_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_eval_node_struct_get_node_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_title_struct_get_string_info(const void *s, char **d, int *bounds, bool *a);
void tao_title_struct_set_string(void *struct_obj, const char *str_ptr, int str_len);
void tao_title_struct_get_x(const void *struct_obj, double *value_out);
void tao_title_struct_set_x(void *struct_obj, double value_in);
void tao_title_struct_get_y(const void *struct_obj, double *value_out);
void tao_title_struct_set_y(void *struct_obj, double value_in);
void tao_title_struct_get_units_info(const void *s, char **d, int *bounds, bool *a);
void tao_title_struct_set_units(void *struct_obj, const char *str_ptr, int str_len);
void tao_title_struct_get_justify_info(const void *s, char **d, int *bounds, bool *a);
void tao_title_struct_set_justify(void *struct_obj, const char *str_ptr, int str_len);
void tao_title_struct_get_draw_it(const void *struct_obj, bool *value_out);
void tao_title_struct_set_draw_it(void *struct_obj, bool value_in);
void qp_rect_struct_get_x1(const void *struct_obj, double *value_out);
void qp_rect_struct_set_x1(void *struct_obj, double value_in);
void qp_rect_struct_get_x2(const void *struct_obj, double *value_out);
void qp_rect_struct_set_x2(void *struct_obj, double value_in);
void qp_rect_struct_get_y1(const void *struct_obj, double *value_out);
void qp_rect_struct_set_y1(void *struct_obj, double value_in);
void qp_rect_struct_get_y2(const void *struct_obj, double *value_out);
void qp_rect_struct_set_y2(void *struct_obj, double value_in);
void qp_rect_struct_get_units_info(const void *s, char **d, int *bounds, bool *a);
void qp_rect_struct_set_units(void *struct_obj, const char *str_ptr, int str_len);

void tao_drawing_struct_get_ele_shape_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_shape_pattern_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_shape_pattern_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_shape_pattern_struct_get_line(const void *struct_obj, void **ptr_out);
void tao_shape_pattern_struct_set_line(void *struct_obj, const void *src_ptr);

void tao_shape_pattern_struct_get_pt_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_shape_pattern_point_struct_get_s(const void *struct_obj, double *value_out);
void tao_shape_pattern_point_struct_set_s(void *struct_obj, double value_in);
void tao_shape_pattern_point_struct_get_y(const void *struct_obj, double *value_out);
void tao_shape_pattern_point_struct_set_y(void *struct_obj, double value_in);
void tao_shape_pattern_point_struct_get_radius(const void *struct_obj, double *value_out);
void tao_shape_pattern_point_struct_set_radius(void *struct_obj, double value_in);
void qp_axis_struct_get_label_info(const void *s, char **d, int *bounds, bool *a);
void qp_axis_struct_set_label(void *struct_obj, const char *str_ptr, int str_len);
void qp_axis_struct_get_min(const void *struct_obj, double *value_out);
void qp_axis_struct_set_min(void *struct_obj, double value_in);
void qp_axis_struct_get_max(const void *struct_obj, double *value_out);
void qp_axis_struct_set_max(void *struct_obj, double value_in);
void qp_axis_struct_get_tick_min(const void *struct_obj, double *value_out);
void qp_axis_struct_set_tick_min(void *struct_obj, double value_in);
void qp_axis_struct_get_tick_max(const void *struct_obj, double *value_out);
void qp_axis_struct_set_tick_max(void *struct_obj, double value_in);
void qp_axis_struct_get_eval_min(const void *struct_obj, double *value_out);
void qp_axis_struct_set_eval_min(void *struct_obj, double value_in);
void qp_axis_struct_get_eval_max(const void *struct_obj, double *value_out);
void qp_axis_struct_set_eval_max(void *struct_obj, double value_in);
void qp_axis_struct_get_dtick(const void *struct_obj, double *value_out);
void qp_axis_struct_set_dtick(void *struct_obj, double value_in);
void qp_axis_struct_get_number_offset(const void *struct_obj, double *value_out);
void qp_axis_struct_set_number_offset(void *struct_obj, double value_in);
void qp_axis_struct_get_label_offset(const void *struct_obj, double *value_out);
void qp_axis_struct_set_label_offset(void *struct_obj, double value_in);
void qp_axis_struct_get_major_tick_len(const void *struct_obj, double *value_out);
void qp_axis_struct_set_major_tick_len(void *struct_obj, double value_in);
void qp_axis_struct_get_minor_tick_len(const void *struct_obj, double *value_out);
void qp_axis_struct_set_minor_tick_len(void *struct_obj, double value_in);
void qp_axis_struct_get_label_color_info(const void *s, char **d, int *bounds, bool *a);
void qp_axis_struct_set_label_color(void *struct_obj, const char *str_ptr, int str_len);
void qp_axis_struct_get_major_div(const void *struct_obj, int *value_out);
void qp_axis_struct_set_major_div(void *struct_obj, int value_in);
void qp_axis_struct_get_major_div_nominal(const void *struct_obj, int *value_out);
void qp_axis_struct_set_major_div_nominal(void *struct_obj, int value_in);
void qp_axis_struct_get_minor_div(const void *struct_obj, int *value_out);
void qp_axis_struct_set_minor_div(void *struct_obj, int value_in);
void qp_axis_struct_get_minor_div_max(const void *struct_obj, int *value_out);
void qp_axis_struct_set_minor_div_max(void *struct_obj, int value_in);
void qp_axis_struct_get_places(const void *struct_obj, int *value_out);
void qp_axis_struct_set_places(void *struct_obj, int value_in);
void qp_axis_struct_get_type_info(const void *s, char **d, int *bounds, bool *a);
void qp_axis_struct_set_type(void *struct_obj, const char *str_ptr, int str_len);
void qp_axis_struct_get_bounds_info(const void *s, char **d, int *bounds, bool *a);
void qp_axis_struct_set_bounds(void *struct_obj, const char *str_ptr, int str_len);
void qp_axis_struct_get_tick_side(const void *struct_obj, int *value_out);
void qp_axis_struct_set_tick_side(void *struct_obj, int value_in);
void qp_axis_struct_get_number_side(const void *struct_obj, int *value_out);
void qp_axis_struct_set_number_side(void *struct_obj, int value_in);
void qp_axis_struct_get_draw_label(const void *struct_obj, bool *value_out);
void qp_axis_struct_set_draw_label(void *struct_obj, bool value_in);
void qp_axis_struct_get_draw_numbers(const void *struct_obj, bool *value_out);
void qp_axis_struct_set_draw_numbers(void *struct_obj, bool value_in);
void qp_legend_struct_get_row_spacing(const void *struct_obj, double *value_out);
void qp_legend_struct_set_row_spacing(void *struct_obj, double value_in);
void qp_legend_struct_get_line_length(const void *struct_obj, double *value_out);
void qp_legend_struct_set_line_length(void *struct_obj, double value_in);
void qp_legend_struct_get_text_offset(const void *struct_obj, double *value_out);
void qp_legend_struct_set_text_offset(void *struct_obj, double value_in);
void qp_legend_struct_get_draw_line(const void *struct_obj, bool *value_out);
void qp_legend_struct_set_draw_line(void *struct_obj, bool value_in);
void qp_legend_struct_get_draw_symbol(const void *struct_obj, bool *value_out);
void qp_legend_struct_set_draw_symbol(void *struct_obj, bool value_in);
void qp_legend_struct_get_draw_text(const void *struct_obj, bool *value_out);
void qp_legend_struct_set_draw_text(void *struct_obj, bool value_in);
void qp_point_struct_get_x(const void *struct_obj, double *value_out);
void qp_point_struct_set_x(void *struct_obj, double value_in);
void qp_point_struct_get_y(const void *struct_obj, double *value_out);
void qp_point_struct_set_y(void *struct_obj, double value_in);
void qp_point_struct_get_units_info(const void *s, char **d, int *bounds, bool *a);
void qp_point_struct_set_units(void *struct_obj, const char *str_ptr, int str_len);
void qp_line_struct_get_width(const void *struct_obj, int *value_out);
void qp_line_struct_set_width(void *struct_obj, int value_in);
void qp_line_struct_get_color_info(const void *s, char **d, int *bounds, bool *a);
void qp_line_struct_set_color(void *struct_obj, const char *str_ptr, int str_len);
void qp_line_struct_get_pattern_info(const void *s, char **d, int *bounds, bool *a);
void qp_line_struct_set_pattern(void *struct_obj, const char *str_ptr, int str_len);
void qp_symbol_struct_get_type_info(const void *s, char **d, int *bounds, bool *a);
void qp_symbol_struct_set_type(void *struct_obj, const char *str_ptr, int str_len);
void qp_symbol_struct_get_height(const void *struct_obj, double *value_out);
void qp_symbol_struct_set_height(void *struct_obj, double value_in);
void qp_symbol_struct_get_color_info(const void *s, char **d, int *bounds, bool *a);
void qp_symbol_struct_set_color(void *struct_obj, const char *str_ptr, int str_len);
void qp_symbol_struct_get_fill_pattern_info(const void *s, char **d, int *bounds, bool *a);
void qp_symbol_struct_set_fill_pattern(void *struct_obj, const char *str_ptr, int str_len);
void qp_symbol_struct_get_line_width(const void *struct_obj, int *value_out);
void qp_symbol_struct_set_line_width(void *struct_obj, int value_in);
void tao_floor_plan_struct_get_view_info(const void *s, char **d, int *bounds, bool *a);
void tao_floor_plan_struct_set_view(void *struct_obj, const char *str_ptr, int str_len);
void tao_floor_plan_struct_get_rotation(const void *struct_obj, double *value_out);
void tao_floor_plan_struct_set_rotation(void *struct_obj, double value_in);
void tao_floor_plan_struct_get_correct_distortion(const void *struct_obj, bool *value_out);
void tao_floor_plan_struct_set_correct_distortion(void *struct_obj, bool value_in);
void tao_floor_plan_struct_get_flip_label_side(const void *struct_obj, bool *value_out);
void tao_floor_plan_struct_set_flip_label_side(void *struct_obj, bool value_in);
void tao_floor_plan_struct_get_size_is_absolute(const void *struct_obj, bool *value_out);
void tao_floor_plan_struct_set_size_is_absolute(void *struct_obj, bool value_in);
void tao_floor_plan_struct_get_draw_only_first_pass(const void *struct_obj, bool *value_out);
void tao_floor_plan_struct_set_draw_only_first_pass(void *struct_obj, bool value_in);
void tao_floor_plan_struct_get_draw_building_wall(const void *struct_obj, bool *value_out);
void tao_floor_plan_struct_set_draw_building_wall(void *struct_obj, bool value_in);
void tao_floor_plan_struct_get_orbit_scale(const void *struct_obj, double *value_out);
void tao_floor_plan_struct_set_orbit_scale(void *struct_obj, double value_in);
void tao_floor_plan_struct_get_orbit_color_info(const void *s, char **d, int *bounds, bool *a);
void tao_floor_plan_struct_set_orbit_color(void *struct_obj, const char *str_ptr, int str_len);
void tao_floor_plan_struct_get_orbit_pattern_info(const void *s, char **d, int *bounds, bool *a);
void tao_floor_plan_struct_set_orbit_pattern(void *struct_obj, const char *str_ptr, int str_len);
void tao_floor_plan_struct_get_orbit_lattice_info(const void *s, char **d, int *bounds, bool *a);
void tao_floor_plan_struct_set_orbit_lattice(void *struct_obj, const char *str_ptr, int str_len);
void tao_floor_plan_struct_get_orbit_width(const void *struct_obj, int *value_out);
void tao_floor_plan_struct_set_orbit_width(void *struct_obj, int value_in);
void tao_v1_var_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_v1_var_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_v1_var_struct_get_ix_v1_var(const void *struct_obj, int *value_out);
void tao_v1_var_struct_set_ix_v1_var(void *struct_obj, int value_in);

void tao_v1_var_struct_get_v_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_global_struct_get_beam_dead_cutoff(const void *struct_obj, double *value_out);
void tao_global_struct_set_beam_dead_cutoff(void *struct_obj, double value_in);
void tao_global_struct_get_lm_opt_deriv_reinit(const void *struct_obj, double *value_out);
void tao_global_struct_set_lm_opt_deriv_reinit(void *struct_obj, double value_in);
void tao_global_struct_get_de_lm_step_ratio(const void *struct_obj, double *value_out);
void tao_global_struct_set_de_lm_step_ratio(void *struct_obj, double value_in);
void tao_global_struct_get_de_var_to_population_factor(const void *struct_obj, double *value_out);
void tao_global_struct_set_de_var_to_population_factor(void *struct_obj, double value_in);
void tao_global_struct_get_lmdif_eps(const void *struct_obj, double *value_out);
void tao_global_struct_set_lmdif_eps(void *struct_obj, double value_in);
void tao_global_struct_get_lmdif_negligible_merit(const void *struct_obj, double *value_out);
void tao_global_struct_set_lmdif_negligible_merit(void *struct_obj, double value_in);
void tao_global_struct_get_svd_cutoff(const void *struct_obj, double *value_out);
void tao_global_struct_set_svd_cutoff(void *struct_obj, double value_in);
void tao_global_struct_get_unstable_penalty(const void *struct_obj, double *value_out);
void tao_global_struct_set_unstable_penalty(void *struct_obj, double value_in);
void tao_global_struct_get_merit_stop_value(const void *struct_obj, double *value_out);
void tao_global_struct_set_merit_stop_value(void *struct_obj, double value_in);
void tao_global_struct_get_dmerit_stop_value(const void *struct_obj, double *value_out);
void tao_global_struct_set_dmerit_stop_value(void *struct_obj, double value_in);
void tao_global_struct_get_random_sigma_cutoff(const void *struct_obj, double *value_out);
void tao_global_struct_set_random_sigma_cutoff(void *struct_obj, double value_in);
void tao_global_struct_get_delta_e_chrom(const void *struct_obj, double *value_out);
void tao_global_struct_set_delta_e_chrom(void *struct_obj, double value_in);
void tao_global_struct_get_max_plot_time(const void *struct_obj, double *value_out);
void tao_global_struct_set_max_plot_time(void *struct_obj, double value_in);
void tao_global_struct_get_default_universe(const void *struct_obj, int *value_out);
void tao_global_struct_set_default_universe(void *struct_obj, int value_in);
void tao_global_struct_get_default_branch(const void *struct_obj, int *value_out);
void tao_global_struct_set_default_branch(void *struct_obj, int value_in);
void tao_global_struct_get_n_opti_cycles(const void *struct_obj, int *value_out);
void tao_global_struct_set_n_opti_cycles(void *struct_obj, int value_in);
void tao_global_struct_get_n_opti_loops(const void *struct_obj, int *value_out);
void tao_global_struct_set_n_opti_loops(void *struct_obj, int value_in);
void tao_global_struct_get_n_threads(const void *struct_obj, int *value_out);
void tao_global_struct_set_n_threads(void *struct_obj, int value_in);
void tao_global_struct_get_phase_units(const void *struct_obj, int *value_out);
void tao_global_struct_set_phase_units(void *struct_obj, int value_in);
void tao_global_struct_get_bunch_to_plot(const void *struct_obj, int *value_out);
void tao_global_struct_set_bunch_to_plot(void *struct_obj, int value_in);
void tao_global_struct_get_random_seed(const void *struct_obj, int *value_out);
void tao_global_struct_set_random_seed(void *struct_obj, int value_in);
void tao_global_struct_get_n_top10_merit(const void *struct_obj, int *value_out);
void tao_global_struct_set_n_top10_merit(void *struct_obj, int value_in);
void tao_global_struct_get_srdt_gen_n_slices(const void *struct_obj, int *value_out);
void tao_global_struct_set_srdt_gen_n_slices(void *struct_obj, int value_in);
void tao_global_struct_get_datum_err_messages_max(const void *struct_obj, int *value_out);
void tao_global_struct_set_datum_err_messages_max(void *struct_obj, int value_in);
void tao_global_struct_get_srdt_sxt_n_slices(const void *struct_obj, int *value_out);
void tao_global_struct_set_srdt_sxt_n_slices(void *struct_obj, int value_in);
void tao_global_struct_get_srdt_use_cache(const void *struct_obj, bool *value_out);
void tao_global_struct_set_srdt_use_cache(void *struct_obj, bool value_in);
void tao_global_struct_get_quiet_info(const void *s, char **d, int *bounds, bool *a);
void tao_global_struct_set_quiet(void *struct_obj, const char *str_ptr, int str_len);
void tao_global_struct_get_random_engine_info(const void *s, char **d, int *bounds, bool *a);
void tao_global_struct_set_random_engine(void *struct_obj, const char *str_ptr, int str_len);
void tao_global_struct_get_random_gauss_converter_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void tao_global_struct_set_random_gauss_converter(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);
void tao_global_struct_get_track_type_info(const void *s, char **d, int *bounds, bool *a);
void tao_global_struct_set_track_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_global_struct_get_lat_sigma_calc_uses_emit_from_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void tao_global_struct_set_lat_sigma_calc_uses_emit_from(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);
void tao_global_struct_get_prompt_string_info(const void *s, char **d, int *bounds, bool *a);
void tao_global_struct_set_prompt_string(void *struct_obj, const char *str_ptr, int str_len);
void tao_global_struct_get_prompt_color_info(const void *s, char **d, int *bounds, bool *a);
void tao_global_struct_set_prompt_color(void *struct_obj, const char *str_ptr, int str_len);
void tao_global_struct_get_optimizer_info(const void *s, char **d, int *bounds, bool *a);
void tao_global_struct_set_optimizer(void *struct_obj, const char *str_ptr, int str_len);
void tao_global_struct_get_print_command_info(const void *s, char **d, int *bounds, bool *a);
void tao_global_struct_set_print_command(void *struct_obj, const char *str_ptr, int str_len);
void tao_global_struct_get_var_out_file_info(const void *s, char **d, int *bounds, bool *a);
void tao_global_struct_set_var_out_file(void *struct_obj, const char *str_ptr, int str_len);
void tao_global_struct_get_history_file_info(const void *s, char **d, int *bounds, bool *a);
void tao_global_struct_set_history_file(void *struct_obj, const char *str_ptr, int str_len);
void tao_global_struct_get_beam_timer_on(const void *struct_obj, bool *value_out);
void tao_global_struct_set_beam_timer_on(void *struct_obj, bool value_in);
void tao_global_struct_get_box_plots(const void *struct_obj, bool *value_out);
void tao_global_struct_set_box_plots(void *struct_obj, bool value_in);
void tao_global_struct_get_blank_line_between_commands(const void *struct_obj, bool *value_out);
void tao_global_struct_set_blank_line_between_commands(void *struct_obj, bool value_in);
void tao_global_struct_get_cmd_file_abort_on_error(const void *struct_obj, bool *value_out);
void tao_global_struct_set_cmd_file_abort_on_error(void *struct_obj, bool value_in);
void tao_global_struct_get_concatenate_maps(const void *struct_obj, bool *value_out);
void tao_global_struct_set_concatenate_maps(void *struct_obj, bool value_in);
void tao_global_struct_get_derivative_recalc(const void *struct_obj, bool *value_out);
void tao_global_struct_set_derivative_recalc(void *struct_obj, bool value_in);
void tao_global_struct_get_derivative_uses_design(const void *struct_obj, bool *value_out);
void tao_global_struct_set_derivative_uses_design(void *struct_obj, bool value_in);
void tao_global_struct_get_disable_smooth_line_calc(const void *struct_obj, bool *value_out);
void tao_global_struct_set_disable_smooth_line_calc(void *struct_obj, bool value_in);
void tao_global_struct_get_draw_curve_off_scale_warn(const void *struct_obj, bool *value_out);
void tao_global_struct_set_draw_curve_off_scale_warn(void *struct_obj, bool value_in);
void tao_global_struct_get_external_plotting(const void *struct_obj, bool *value_out);
void tao_global_struct_set_external_plotting(void *struct_obj, bool value_in);
void tao_global_struct_get_label_lattice_elements(const void *struct_obj, bool *value_out);
void tao_global_struct_set_label_lattice_elements(void *struct_obj, bool value_in);
void tao_global_struct_get_label_keys(const void *struct_obj, bool *value_out);
void tao_global_struct_set_label_keys(void *struct_obj, bool value_in);
void tao_global_struct_get_lattice_calc_on(const void *struct_obj, bool *value_out);
void tao_global_struct_set_lattice_calc_on(void *struct_obj, bool value_in);
void tao_global_struct_get_only_limit_opt_vars(const void *struct_obj, bool *value_out);
void tao_global_struct_set_only_limit_opt_vars(void *struct_obj, bool value_in);
void tao_global_struct_get_opt_with_ref(const void *struct_obj, bool *value_out);
void tao_global_struct_set_opt_with_ref(void *struct_obj, bool value_in);
void tao_global_struct_get_opt_with_base(const void *struct_obj, bool *value_out);
void tao_global_struct_set_opt_with_base(void *struct_obj, bool value_in);
void tao_global_struct_get_opt_match_auto_recalc(const void *struct_obj, bool *value_out);
void tao_global_struct_set_opt_match_auto_recalc(void *struct_obj, bool value_in);
void tao_global_struct_get_opti_write_var_file(const void *struct_obj, bool *value_out);
void tao_global_struct_set_opti_write_var_file(void *struct_obj, bool value_in);
void tao_global_struct_get_optimizer_allow_user_abort(const void *struct_obj, bool *value_out);
void tao_global_struct_set_optimizer_allow_user_abort(void *struct_obj, bool value_in);
void tao_global_struct_get_optimizer_var_limit_warn(const void *struct_obj, bool *value_out);
void tao_global_struct_set_optimizer_var_limit_warn(void *struct_obj, bool value_in);
void tao_global_struct_get_plot_on(const void *struct_obj, bool *value_out);
void tao_global_struct_set_plot_on(void *struct_obj, bool value_in);
void tao_global_struct_get_rad_int_user_calc_on(const void *struct_obj, bool *value_out);
void tao_global_struct_set_rad_int_user_calc_on(void *struct_obj, bool value_in);
void tao_global_struct_get_rf_on(const void *struct_obj, bool *value_out);
void tao_global_struct_set_rf_on(void *struct_obj, bool value_in);
void tao_global_struct_get_single_step(const void *struct_obj, bool *value_out);
void tao_global_struct_set_single_step(void *struct_obj, bool value_in);
void tao_global_struct_get_stop_on_error(const void *struct_obj, bool *value_out);
void tao_global_struct_set_stop_on_error(void *struct_obj, bool value_in);
void tao_global_struct_get_svd_retreat_on_merit_increase(const void *struct_obj, bool *value_out);
void tao_global_struct_set_svd_retreat_on_merit_increase(void *struct_obj, bool value_in);
void tao_global_struct_get_var_limits_on(const void *struct_obj, bool *value_out);
void tao_global_struct_set_var_limits_on(void *struct_obj, bool value_in);
void tao_global_struct_get_wait_for_CR_in_single_mode(const void *struct_obj, bool *value_out);
void tao_global_struct_set_wait_for_CR_in_single_mode(void *struct_obj, bool value_in);
void tao_global_struct_get_symbol_import(const void *struct_obj, bool *value_out);
void tao_global_struct_set_symbol_import(void *struct_obj, bool value_in);
void tao_global_struct_get_debug_on(const void *struct_obj, bool *value_out);
void tao_global_struct_set_debug_on(void *struct_obj, bool value_in);
void tao_global_struct_get_expression_tree_on(const void *struct_obj, bool *value_out);
void tao_global_struct_set_expression_tree_on(void *struct_obj, bool value_in);
void tao_global_struct_get_verbose_on(const void *struct_obj, bool *value_out);
void tao_global_struct_set_verbose_on(void *struct_obj, bool value_in);
void tao_init_struct_get_parse_cmd_args(const void *struct_obj, bool *value_out);
void tao_init_struct_set_parse_cmd_args(void *struct_obj, bool value_in);
void tao_init_struct_get_debug_switch(const void *struct_obj, bool *value_out);
void tao_init_struct_set_debug_switch(void *struct_obj, bool value_in);
void tao_init_struct_get_external_plotting_switch(const void *struct_obj, bool *value_out);
void tao_init_struct_set_external_plotting_switch(void *struct_obj, bool value_in);
void tao_init_struct_get_init_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_init_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_hook_init_file_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_hook_init_file(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_hook_lat_file_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_hook_lat_file(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_hook_beam_file_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_hook_beam_file(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_hook_data_file_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_hook_data_file(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_hook_plot_file_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_hook_plot_file(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_hook_startup_file_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_hook_startup_file(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_hook_var_file_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_hook_var_file(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_hook_building_wall_file_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void tao_init_struct_set_hook_building_wall_file(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);
void tao_init_struct_get_init_file_arg_path_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_init_file_arg_path(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_lattice_file_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_lattice_file_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_hook_init_file_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_hook_init_file_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_init_file_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_init_file_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_beam_file_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_beam_file_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_beam_init_position_file_arg_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void tao_init_struct_set_beam_init_position_file_arg(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);
void tao_init_struct_get_command_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_command_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_data_file_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_data_file_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_plot_file_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_plot_file_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_startup_file_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_startup_file_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_var_file_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_var_file_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_building_wall_file_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_building_wall_file_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_geometry_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_geometry_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_slice_lattice_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_slice_lattice_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_start_branch_at_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_start_branch_at_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_log_startup_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_log_startup_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_no_stopping_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_no_stopping_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_noplot_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_noplot_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_no_rad_int_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_no_rad_int_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_reverse_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_reverse_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_debug_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_debug_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_disable_smooth_line_calc_arg_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void tao_init_struct_set_disable_smooth_line_calc_arg(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);
void tao_init_struct_get_rf_on_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_rf_on_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_prompt_color_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_prompt_color_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_quiet_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_quiet_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_noinit_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_noinit_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_nostartup_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_nostartup_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_symbol_import_arg_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_symbol_import_arg(void *struct_obj, const char *str_ptr, int str_len);
void tao_init_struct_get_unique_name_suffix_info(const void *s, char **d, int *bounds, bool *a);
void tao_init_struct_set_unique_name_suffix(void *struct_obj, const char *str_ptr, int str_len);

void tao_common_struct_get_plot_place_buffer_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_common_struct_get_covar_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void tao_common_struct_set_covar(void *s, const void *d, const int *shape);
void tao_common_struct_get_alpha_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void tao_common_struct_set_alpha(void *s, const void *d, const int *shape);
void tao_common_struct_get_dummy_target(const void *struct_obj, double *value_out);
void tao_common_struct_set_dummy_target(void *struct_obj, double value_in);
void tao_common_struct_get_n_alias(const void *struct_obj, int *value_out);
void tao_common_struct_set_n_alias(void *struct_obj, int value_in);
void tao_common_struct_get_cmd_file_level(const void *struct_obj, int *value_out);
void tao_common_struct_set_cmd_file_level(void *struct_obj, int value_in);
void tao_common_struct_get_ix_key_bank(const void *struct_obj, int *value_out);
void tao_common_struct_set_ix_key_bank(void *struct_obj, int value_in);
void tao_common_struct_get_ix_history(const void *struct_obj, int *value_out);
void tao_common_struct_set_ix_history(void *struct_obj, int value_in);
void tao_common_struct_get_n_history(const void *struct_obj, int *value_out);
void tao_common_struct_set_n_history(void *struct_obj, int value_in);
void tao_common_struct_get_lev_loop(const void *struct_obj, int *value_out);
void tao_common_struct_set_lev_loop(void *struct_obj, int value_in);
void tao_common_struct_get_n_err_messages_printed(const void *struct_obj, int *value_out);
void tao_common_struct_set_n_err_messages_printed(void *struct_obj, int value_in);
void tao_common_struct_get_n_universes(const void *struct_obj, int *value_out);
void tao_common_struct_set_n_universes(void *struct_obj, int value_in);
void tao_common_struct_get_ix_beam_track_active_element(const void *struct_obj, int *value_out);
void tao_common_struct_set_ix_beam_track_active_element(void *struct_obj, int value_in);
void tao_common_struct_get_cmd_file_paused(const void *struct_obj, bool *value_out);
void tao_common_struct_set_cmd_file_paused(void *struct_obj, bool value_in);
void tao_common_struct_get_use_cmd_here(const void *struct_obj, bool *value_out);
void tao_common_struct_set_use_cmd_here(void *struct_obj, bool value_in);
void tao_common_struct_get_cmd_from_cmd_file(const void *struct_obj, bool *value_out);
void tao_common_struct_set_cmd_from_cmd_file(void *struct_obj, bool value_in);
void tao_common_struct_get_use_saved_beam_in_tracking(const void *struct_obj, bool *value_out);
void tao_common_struct_set_use_saved_beam_in_tracking(void *struct_obj, bool value_in);
void tao_common_struct_get_single_mode(const void *struct_obj, bool *value_out);
void tao_common_struct_set_single_mode(void *struct_obj, bool value_in);
void tao_common_struct_get_combine_consecutive_elements_of_like_name(
    const void *struct_obj,
    bool *value_out
);
void tao_common_struct_set_combine_consecutive_elements_of_like_name(
    void *struct_obj,
    bool value_in
);
void tao_common_struct_get_have_tracked_beam(const void *struct_obj, bool *value_out);
void tao_common_struct_set_have_tracked_beam(void *struct_obj, bool value_in);
void tao_common_struct_get_init_plot_needed(const void *struct_obj, bool *value_out);
void tao_common_struct_set_init_plot_needed(void *struct_obj, bool value_in);
void tao_common_struct_get_init_beam(const void *struct_obj, bool *value_out);
void tao_common_struct_set_init_beam(void *struct_obj, bool value_in);
void tao_common_struct_get_init_var(const void *struct_obj, bool *value_out);
void tao_common_struct_set_init_var(void *struct_obj, bool value_in);
void tao_common_struct_get_init_read_lat_info(const void *struct_obj, bool *value_out);
void tao_common_struct_set_init_read_lat_info(void *struct_obj, bool value_in);
void tao_common_struct_get_optimizer_running(const void *struct_obj, bool *value_out);
void tao_common_struct_set_optimizer_running(void *struct_obj, bool value_in);
void tao_common_struct_get_have_datums_using_expressions(const void *struct_obj, bool *value_out);
void tao_common_struct_set_have_datums_using_expressions(void *struct_obj, bool value_in);
void tao_common_struct_get_print_to_terminal(const void *struct_obj, bool *value_out);
void tao_common_struct_set_print_to_terminal(void *struct_obj, bool value_in);
void tao_common_struct_get_lattice_calc_done(const void *struct_obj, bool *value_out);
void tao_common_struct_set_lattice_calc_done(void *struct_obj, bool value_in);
void tao_common_struct_get_add_measurement_noise(const void *struct_obj, bool *value_out);
void tao_common_struct_set_add_measurement_noise(void *struct_obj, bool value_in);
void tao_common_struct_get_is_err_message_printed_info(
    const void *s,
    bool **d,
    int *bounds,
    bool *is_alloc
);
void tao_common_struct_set_is_err_message_printed(void *s, const void *d, const int *shape);
void tao_common_struct_get_command_arg_has_been_executed(const void *struct_obj, bool *value_out);
void tao_common_struct_set_command_arg_has_been_executed(void *struct_obj, bool value_in);
void tao_common_struct_get_all_merit_weights_positive(const void *struct_obj, bool *value_out);
void tao_common_struct_set_all_merit_weights_positive(void *struct_obj, bool value_in);
void tao_common_struct_get_multi_turn_orbit_is_plotted(const void *struct_obj, bool *value_out);
void tao_common_struct_set_multi_turn_orbit_is_plotted(void *struct_obj, bool value_in);
void tao_common_struct_get_force_chrom_calc(const void *struct_obj, bool *value_out);
void tao_common_struct_set_force_chrom_calc(void *struct_obj, bool value_in);
void tao_common_struct_get_force_rad_int_calc(const void *struct_obj, bool *value_out);
void tao_common_struct_set_force_rad_int_calc(void *struct_obj, bool value_in);
void tao_common_struct_get_rad_int_ri_calc_on(const void *struct_obj, bool *value_out);
void tao_common_struct_set_rad_int_ri_calc_on(void *struct_obj, bool value_in);
void tao_common_struct_get_rad_int_6d_calc_on(const void *struct_obj, bool *value_out);
void tao_common_struct_set_rad_int_6d_calc_on(void *struct_obj, bool value_in);

void tao_common_struct_get_valid_plot_who_info(
    const void *s,
    char **d,
    int *bounds, // [lower, upper]
    int *str_len,
    bool *is_alloc
);

void tao_common_struct_get_single_mode_buffer_info(const void *s, char **d, int *bounds, bool *a);
void tao_common_struct_set_single_mode_buffer(void *struct_obj, const char *str_ptr, int str_len);
void tao_common_struct_get_cmd_info(const void *s, char **d, int *bounds, bool *a);
void tao_common_struct_set_cmd(void *struct_obj, const char *str_ptr, int str_len);
void tao_common_struct_get_saved_cmd_line_info(const void *s, char **d, int *bounds, bool *a);
void tao_common_struct_set_saved_cmd_line(void *struct_obj, const char *str_ptr, int str_len);
void tao_plot_page_struct_get_title(const void *struct_obj, void **ptr_out);
void tao_plot_page_struct_set_title(void *struct_obj, const void *src_ptr);
void tao_plot_page_struct_get_subtitle(const void *struct_obj, void **ptr_out);
void tao_plot_page_struct_set_subtitle(void *struct_obj, const void *src_ptr);
void tao_plot_page_struct_get_border(const void *struct_obj, void **ptr_out);
void tao_plot_page_struct_set_border(void *struct_obj, const void *src_ptr);
void tao_plot_page_struct_get_floor_plan(const void *struct_obj, void **ptr_out);
void tao_plot_page_struct_set_floor_plan(void *struct_obj, const void *src_ptr);
void tao_plot_page_struct_get_lat_layout(const void *struct_obj, void **ptr_out);
void tao_plot_page_struct_set_lat_layout(void *struct_obj, const void *src_ptr);

void tao_plot_page_struct_get_pattern_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_plot_page_struct_get_template_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_plot_page_struct_get_region_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_plot_page_struct_get_plot_display_type_info(const void *s, char **d, int *bounds, bool *a);
void tao_plot_page_struct_set_plot_display_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_plot_page_struct_get_size_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_plot_page_struct_set_size(void *s, const void *d, const int *shape);
void tao_plot_page_struct_get_text_height(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_text_height(void *struct_obj, double value_in);
void tao_plot_page_struct_get_main_title_text_scale(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_main_title_text_scale(void *struct_obj, double value_in);
void tao_plot_page_struct_get_graph_title_text_scale(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_graph_title_text_scale(void *struct_obj, double value_in);
void tao_plot_page_struct_get_axis_number_text_scale(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_axis_number_text_scale(void *struct_obj, double value_in);
void tao_plot_page_struct_get_axis_label_text_scale(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_axis_label_text_scale(void *struct_obj, double value_in);
void tao_plot_page_struct_get_legend_text_scale(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_legend_text_scale(void *struct_obj, double value_in);
void tao_plot_page_struct_get_key_table_text_scale(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_key_table_text_scale(void *struct_obj, double value_in);
void tao_plot_page_struct_get_floor_plan_shape_scale(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_floor_plan_shape_scale(void *struct_obj, double value_in);
void tao_plot_page_struct_get_floor_plan_text_scale(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_floor_plan_text_scale(void *struct_obj, double value_in);
void tao_plot_page_struct_get_lat_layout_shape_scale(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_lat_layout_shape_scale(void *struct_obj, double value_in);
void tao_plot_page_struct_get_lat_layout_text_scale(const void *struct_obj, double *value_out);
void tao_plot_page_struct_set_lat_layout_text_scale(void *struct_obj, double value_in);
void tao_plot_page_struct_get_n_curve_pts(const void *struct_obj, int *value_out);
void tao_plot_page_struct_set_n_curve_pts(void *struct_obj, int value_in);
void tao_plot_page_struct_get_id_window(const void *struct_obj, int *value_out);
void tao_plot_page_struct_set_id_window(void *struct_obj, int value_in);
void tao_plot_page_struct_get_delete_overlapping_plots(const void *struct_obj, bool *value_out);
void tao_plot_page_struct_set_delete_overlapping_plots(void *struct_obj, bool value_in);
void tao_plot_page_struct_get_draw_graph_title_suffix(const void *struct_obj, bool *value_out);
void tao_plot_page_struct_set_draw_graph_title_suffix(void *struct_obj, bool value_in);
void tao_building_wall_struct_get_orientation(const void *struct_obj, void **ptr_out);
void tao_building_wall_struct_set_orientation(void *struct_obj, const void *src_ptr);

void tao_building_wall_struct_get_section_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_building_wall_orientation_struct_get_theta(const void *struct_obj, double *value_out);
void tao_building_wall_orientation_struct_set_theta(void *struct_obj, double value_in);
void tao_building_wall_orientation_struct_get_x_offset(const void *struct_obj, double *value_out);
void tao_building_wall_orientation_struct_set_x_offset(void *struct_obj, double value_in);
void tao_building_wall_orientation_struct_get_z_offset(const void *struct_obj, double *value_out);
void tao_building_wall_orientation_struct_set_z_offset(void *struct_obj, double value_in);
void tao_building_wall_section_struct_get_name_info(const void *s, char **d, int *bounds, bool *a);
void tao_building_wall_section_struct_set_name(void *struct_obj, const char *str_ptr, int str_len);
void tao_building_wall_section_struct_get_constraint_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void tao_building_wall_section_struct_set_constraint(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);

void tao_building_wall_section_struct_get_point_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_building_wall_point_struct_get_z(const void *struct_obj, double *value_out);
void tao_building_wall_point_struct_set_z(void *struct_obj, double value_in);
void tao_building_wall_point_struct_get_x(const void *struct_obj, double *value_out);
void tao_building_wall_point_struct_set_x(void *struct_obj, double value_in);
void tao_building_wall_point_struct_get_radius(const void *struct_obj, double *value_out);
void tao_building_wall_point_struct_set_radius(void *struct_obj, double value_in);
void tao_building_wall_point_struct_get_z_center(const void *struct_obj, double *value_out);
void tao_building_wall_point_struct_set_z_center(void *struct_obj, double value_in);
void tao_building_wall_point_struct_get_x_center(const void *struct_obj, double *value_out);
void tao_building_wall_point_struct_set_x_center(void *struct_obj, double value_in);
void tao_wave_struct_get_data_type_info(const void *s, char **d, int *bounds, bool *a);
void tao_wave_struct_set_data_type(void *struct_obj, const char *str_ptr, int str_len);
void tao_wave_struct_get_rms_rel_a(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_rel_a(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_b(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_rel_b(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_as(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_rel_as(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_bs(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_rel_bs(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_ar(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_rel_ar(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_br(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_rel_br(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_k(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_rel_k(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_ks(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_rel_ks(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_rel_kr(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_rel_kr(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_phi(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_phi(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_phi_s(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_phi_s(void *struct_obj, double value_in);
void tao_wave_struct_get_rms_phi_r(const void *struct_obj, double *value_out);
void tao_wave_struct_set_rms_phi_r(void *struct_obj, double value_in);
void tao_wave_struct_get_amp_ba_s(const void *struct_obj, double *value_out);
void tao_wave_struct_set_amp_ba_s(void *struct_obj, double value_in);
void tao_wave_struct_get_amp_ba_r(const void *struct_obj, double *value_out);
void tao_wave_struct_set_amp_ba_r(void *struct_obj, double value_in);
void tao_wave_struct_get_chi_a(const void *struct_obj, double *value_out);
void tao_wave_struct_set_chi_a(void *struct_obj, double value_in);
void tao_wave_struct_get_chi_c(const void *struct_obj, double *value_out);
void tao_wave_struct_set_chi_c(void *struct_obj, double value_in);
void tao_wave_struct_get_chi_ba(const void *struct_obj, double *value_out);
void tao_wave_struct_set_chi_ba(void *struct_obj, double value_in);
void tao_wave_struct_get_amp_a_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_wave_struct_set_amp_a(void *s, const void *d, const int *shape);
void tao_wave_struct_get_amp_b_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_wave_struct_set_amp_b(void *s, const void *d, const int *shape);
void tao_wave_struct_get_amp_ba_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_wave_struct_set_amp_ba(void *s, const void *d, const int *shape);
void tao_wave_struct_get_coef_a_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_wave_struct_set_coef_a(void *s, const void *d, const int *shape);
void tao_wave_struct_get_coef_b_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_wave_struct_set_coef_b(void *s, const void *d, const int *shape);
void tao_wave_struct_get_coef_ba_info(const void *s, double **d, int *bounds, bool *is_alloc);
void tao_wave_struct_set_coef_ba(void *s, const void *d, const int *shape);
void tao_wave_struct_get_n_func(const void *struct_obj, int *value_out);
void tao_wave_struct_set_n_func(void *struct_obj, int value_in);
void tao_wave_struct_get_ix_a1(const void *struct_obj, int *value_out);
void tao_wave_struct_set_ix_a1(void *struct_obj, int value_in);
void tao_wave_struct_get_ix_a2(const void *struct_obj, int *value_out);
void tao_wave_struct_set_ix_a2(void *struct_obj, int value_in);
void tao_wave_struct_get_ix_b1(const void *struct_obj, int *value_out);
void tao_wave_struct_set_ix_b1(void *struct_obj, int value_in);
void tao_wave_struct_get_ix_b2(const void *struct_obj, int *value_out);
void tao_wave_struct_set_ix_b2(void *struct_obj, int value_in);
void tao_wave_struct_get_i_a1(const void *struct_obj, int *value_out);
void tao_wave_struct_set_i_a1(void *struct_obj, int value_in);
void tao_wave_struct_get_i_a2(const void *struct_obj, int *value_out);
void tao_wave_struct_set_i_a2(void *struct_obj, int value_in);
void tao_wave_struct_get_i_b1(const void *struct_obj, int *value_out);
void tao_wave_struct_set_i_b1(void *struct_obj, int value_in);
void tao_wave_struct_get_i_b2(const void *struct_obj, int *value_out);
void tao_wave_struct_set_i_b2(void *struct_obj, int value_in);
void tao_wave_struct_get_n_a(const void *struct_obj, int *value_out);
void tao_wave_struct_set_n_a(void *struct_obj, int value_in);
void tao_wave_struct_get_n_b(const void *struct_obj, int *value_out);
void tao_wave_struct_set_n_b(void *struct_obj, int value_in);
void tao_wave_struct_get_i_curve_wrap_pt(const void *struct_obj, int *value_out);
void tao_wave_struct_set_i_curve_wrap_pt(void *struct_obj, int value_in);
void tao_wave_struct_get_ix_data_info(const void *s, int **d, int *bounds, bool *is_alloc);
void tao_wave_struct_set_ix_data(void *s, const void *d, const int *shape);
void tao_wave_struct_get_n_kick(const void *struct_obj, int *value_out);
void tao_wave_struct_set_n_kick(void *struct_obj, int value_in);

void tao_wave_struct_get_kick_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_wave_struct_get_base_graph(const void *struct_obj, void **ptr_out);
void tao_wave_struct_set_base_graph(void *struct_obj, const void *src_ptr);
void tao_wave_struct_get_region(const void *struct_obj, void **ptr_out);
void tao_wave_struct_set_region(void *struct_obj, const void *src_ptr);
void tao_wave_struct_get_d1_dat(const void *struct_obj, void **ptr_out);
void tao_wave_struct_set_d1_dat(void *struct_obj, const void *src_ptr);
void tao_wave_kick_pt_struct_get_phi_s(const void *struct_obj, double *value_out);
void tao_wave_kick_pt_struct_set_phi_s(void *struct_obj, double value_in);
void tao_wave_kick_pt_struct_get_phi_r(const void *struct_obj, double *value_out);
void tao_wave_kick_pt_struct_set_phi_r(void *struct_obj, double value_in);
void tao_wave_kick_pt_struct_get_phi(const void *struct_obj, double *value_out);
void tao_wave_kick_pt_struct_set_phi(void *struct_obj, double value_in);
void tao_wave_kick_pt_struct_get_amp(const void *struct_obj, double *value_out);
void tao_wave_kick_pt_struct_set_amp(void *struct_obj, double value_in);
void tao_wave_kick_pt_struct_get_s(const void *struct_obj, double *value_out);
void tao_wave_kick_pt_struct_set_s(void *struct_obj, double value_in);
void tao_wave_kick_pt_struct_get_ix_dat_before_kick(const void *struct_obj, int *value_out);
void tao_wave_kick_pt_struct_set_ix_dat_before_kick(void *struct_obj, int value_in);
void tao_wave_kick_pt_struct_get_ele(const void *struct_obj, void **ptr_out);
void tao_wave_kick_pt_struct_set_ele(void *struct_obj, const void *src_ptr);

void tao_cmd_history_struct_get_cmd_info(const void *s, char **d, int *len, bool *is_alloc);

void tao_cmd_history_struct_set_cmd(void *struct_obj, const char *str_ptr, int str_len);
void tao_cmd_history_struct_get_ix(const void *struct_obj, int *value_out);
void tao_cmd_history_struct_set_ix(void *struct_obj, int value_in);
void tao_universe_struct_get_model(const void *struct_obj, void **ptr_out);
void tao_universe_struct_set_model(void *struct_obj, const void *src_ptr);
void tao_universe_struct_get_design(const void *struct_obj, void **ptr_out);
void tao_universe_struct_set_design(void *struct_obj, const void *src_ptr);
void tao_universe_struct_get_base(const void *struct_obj, void **ptr_out);
void tao_universe_struct_set_base(void *struct_obj, const void *src_ptr);
void tao_universe_struct_get_beam(const void *struct_obj, void **ptr_out);
void tao_universe_struct_set_beam(void *struct_obj, const void *src_ptr);
void tao_universe_struct_get_dynamic_aperture(const void *struct_obj, void **ptr_out);
void tao_universe_struct_set_dynamic_aperture(void *struct_obj, const void *src_ptr);

void tao_universe_struct_get_model_branch_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_universe_struct_get_d2_data_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_universe_struct_get_data_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void tao_universe_struct_get_ping_scale(const void *struct_obj, void **ptr_out);
void tao_universe_struct_set_ping_scale(void *struct_obj, const void *src_ptr);
void tao_universe_struct_get_scratch_lat(const void *struct_obj, void **ptr_out);
void tao_universe_struct_set_scratch_lat(void *struct_obj, const void *src_ptr);
void tao_universe_struct_get_calc(const void *struct_obj, void **ptr_out);
void tao_universe_struct_set_calc(void *struct_obj, const void *src_ptr);
void tao_universe_struct_get_ele_order(const void *struct_obj, void **ptr_out);
void tao_universe_struct_set_ele_order(void *struct_obj, const void *src_ptr);
void tao_universe_struct_get_spin_map(const void *struct_obj, void **ptr_out);
void tao_universe_struct_set_spin_map(void *struct_obj, const void *src_ptr);
void tao_universe_struct_get_dModel_dVar_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void tao_universe_struct_set_dModel_dVar(void *s, const void *d, const int *shape);
void tao_universe_struct_get_ix_uni(const void *struct_obj, int *value_out);
void tao_universe_struct_set_ix_uni(void *struct_obj, int value_in);
void tao_universe_struct_get_n_d2_data_used(const void *struct_obj, int *value_out);
void tao_universe_struct_set_n_d2_data_used(void *struct_obj, int value_in);
void tao_universe_struct_get_n_data_used(const void *struct_obj, int *value_out);
void tao_universe_struct_set_n_data_used(void *struct_obj, int value_in);
void tao_universe_struct_get_is_on(const void *struct_obj, bool *value_out);
void tao_universe_struct_set_is_on(void *struct_obj, bool value_in);
void tao_universe_struct_get_design_same_as_previous(const void *struct_obj, bool *value_out);
void tao_universe_struct_set_design_same_as_previous(void *struct_obj, bool value_in);
void tao_universe_struct_get_picked_uni(const void *struct_obj, bool *value_out);
void tao_universe_struct_set_picked_uni(void *struct_obj, bool value_in);
void mad_energy_struct_get_total(const void *struct_obj, double *value_out);
void mad_energy_struct_set_total(void *struct_obj, double value_in);
void mad_energy_struct_get_beta(const void *struct_obj, double *value_out);
void mad_energy_struct_set_beta(void *struct_obj, double value_in);
void mad_energy_struct_get_gamma(const void *struct_obj, double *value_out);
void mad_energy_struct_set_gamma(void *struct_obj, double value_in);
void mad_energy_struct_get_kinetic(const void *struct_obj, double *value_out);
void mad_energy_struct_set_kinetic(void *struct_obj, double value_in);
void mad_energy_struct_get_p0c(const void *struct_obj, double *value_out);
void mad_energy_struct_set_p0c(void *struct_obj, double value_in);
void mad_energy_struct_get_particle(const void *struct_obj, int *value_out);
void mad_energy_struct_set_particle(void *struct_obj, int value_in);
void mad_map_struct_get_k_info(const void *s, double **d, int *bounds, bool *is_alloc);
void mad_map_struct_set_k(void *s, const void *d, const int *shape);
void mad_map_struct_get_r_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void mad_map_struct_set_r(void *s, const void *d, const int *shape);
void mad_map_struct_get_t_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void mad_map_struct_set_t(void *s, const void *d, const int *shape);
void random_state_struct_get_ix(const void *struct_obj, int64_t *value_out);
void random_state_struct_set_ix(void *struct_obj, int64_t value_in);
void random_state_struct_get_iy(const void *struct_obj, int64_t *value_out);
void random_state_struct_set_iy(void *struct_obj, int64_t value_in);
void random_state_struct_get_number_stored(const void *struct_obj, bool *value_out);
void random_state_struct_set_number_stored(void *struct_obj, bool value_in);
void random_state_struct_get_h_saved(const void *struct_obj, double *value_out);
void random_state_struct_set_h_saved(void *struct_obj, double value_in);
void random_state_struct_get_engine(const void *struct_obj, int *value_out);
void random_state_struct_set_engine(void *struct_obj, int value_in);
void random_state_struct_get_seed(const void *struct_obj, int *value_out);
void random_state_struct_set_seed(void *struct_obj, int value_in);
void random_state_struct_get_am(const void *struct_obj, double *value_out);
void random_state_struct_set_am(void *struct_obj, double value_in);
void random_state_struct_get_gauss_converter(const void *struct_obj, int *value_out);
void random_state_struct_set_gauss_converter(void *struct_obj, int value_in);
void random_state_struct_get_gauss_sigma_cut(const void *struct_obj, double *value_out);
void random_state_struct_set_gauss_sigma_cut(void *struct_obj, double value_in);
void random_state_struct_get_in_sobseq(const void *struct_obj, int64_t *value_out);
void random_state_struct_set_in_sobseq(void *struct_obj, int64_t value_in);
void random_state_struct_get_ix_sobseq_info(
    const void *s,
    int64_t **d,
    int *bounds,
    bool *is_alloc
);
void random_state_struct_set_ix_sobseq(void *s, const void *d, const int *shape);
void random_state_struct_get_x_sobseq_info(const void *s, double **d, int *bounds, bool *is_alloc);
void random_state_struct_set_x_sobseq(void *s, const void *d, const int *shape);
void bbu_stage_struct_get_ix_ele_lr_wake(const void *struct_obj, int *value_out);
void bbu_stage_struct_set_ix_ele_lr_wake(void *struct_obj, int value_in);
void bbu_stage_struct_get_ix_ele_stage_end(const void *struct_obj, int *value_out);
void bbu_stage_struct_set_ix_ele_stage_end(void *struct_obj, int value_in);
void bbu_stage_struct_get_ix_pass(const void *struct_obj, int *value_out);
void bbu_stage_struct_set_ix_pass(void *struct_obj, int value_in);
void bbu_stage_struct_get_ix_stage_pass1(const void *struct_obj, int *value_out);
void bbu_stage_struct_set_ix_stage_pass1(void *struct_obj, int value_in);
void bbu_stage_struct_get_ix_head_bunch(const void *struct_obj, int *value_out);
void bbu_stage_struct_set_ix_head_bunch(void *struct_obj, int value_in);
void bbu_stage_struct_get_ix_hom_max(const void *struct_obj, int *value_out);
void bbu_stage_struct_set_ix_hom_max(void *struct_obj, int value_in);
void bbu_stage_struct_get_hom_voltage_max(const void *struct_obj, double *value_out);
void bbu_stage_struct_set_hom_voltage_max(void *struct_obj, double value_in);
void bbu_stage_struct_get_time_at_wake_ele(const void *struct_obj, double *value_out);
void bbu_stage_struct_set_time_at_wake_ele(void *struct_obj, double value_in);
void bbu_stage_struct_get_ave_orb_info(const void *s, double **d, int *bounds, bool *is_alloc);
void bbu_stage_struct_set_ave_orb(void *s, const void *d, const int *shape);
void bbu_stage_struct_get_rms_orb_info(const void *s, double **d, int *bounds, bool *is_alloc);
void bbu_stage_struct_set_rms_orb(void *s, const void *d, const int *shape);
void bbu_stage_struct_get_min_orb_info(const void *s, double **d, int *bounds, bool *is_alloc);
void bbu_stage_struct_set_min_orb(void *s, const void *d, const int *shape);
void bbu_stage_struct_get_max_orb_info(const void *s, double **d, int *bounds, bool *is_alloc);
void bbu_stage_struct_set_max_orb(void *s, const void *d, const int *shape);
void bbu_stage_struct_get_n_orb(const void *struct_obj, int *value_out);
void bbu_stage_struct_set_n_orb(void *struct_obj, int value_in);

void bbu_beam_struct_get_bunch_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bbu_beam_struct_get_stage_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void bbu_beam_struct_get_ix_ele_bunch_info(const void *s, int **d, int *bounds, bool *is_alloc);
void bbu_beam_struct_set_ix_ele_bunch(void *s, const void *d, const int *shape);
void bbu_beam_struct_get_ix_bunch_head(const void *struct_obj, int *value_out);
void bbu_beam_struct_set_ix_bunch_head(void *struct_obj, int value_in);
void bbu_beam_struct_get_ix_bunch_end(const void *struct_obj, int *value_out);
void bbu_beam_struct_set_ix_bunch_end(void *struct_obj, int value_in);
void bbu_beam_struct_get_n_bunch_in_lat(const void *struct_obj, int *value_out);
void bbu_beam_struct_set_n_bunch_in_lat(void *struct_obj, int value_in);
void bbu_beam_struct_get_ix_stage_voltage_max(const void *struct_obj, int *value_out);
void bbu_beam_struct_set_ix_stage_voltage_max(void *struct_obj, int value_in);
void bbu_beam_struct_get_hom_voltage_max(const void *struct_obj, double *value_out);
void bbu_beam_struct_set_hom_voltage_max(void *struct_obj, double value_in);
void bbu_beam_struct_get_time_now(const void *struct_obj, double *value_out);
void bbu_beam_struct_set_time_now(void *struct_obj, double value_in);
void bbu_beam_struct_get_one_turn_time(const void *struct_obj, double *value_out);
void bbu_beam_struct_set_one_turn_time(void *struct_obj, double value_in);
void bbu_beam_struct_get_rf_wavelength_max(const void *struct_obj, double *value_out);
void bbu_beam_struct_set_rf_wavelength_max(void *struct_obj, double value_in);
void bbu_param_struct_get_lat_filename_info(const void *s, char **d, int *bounds, bool *a);
void bbu_param_struct_set_lat_filename(void *struct_obj, const char *str_ptr, int str_len);
void bbu_param_struct_get_lat2_filename_info(const void *s, char **d, int *bounds, bool *a);
void bbu_param_struct_set_lat2_filename(void *struct_obj, const char *str_ptr, int str_len);
void bbu_param_struct_get_bunch_by_bunch_info_file_info(
    const void *s,
    char **d,
    int *bounds,
    bool *a
);
void bbu_param_struct_set_bunch_by_bunch_info_file(
    void *struct_obj,
    const char *str_ptr,
    int str_len
);
void bbu_param_struct_get_hybridize(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_hybridize(void *struct_obj, bool value_in);
void bbu_param_struct_get_write_digested_hybrid_lat(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_write_digested_hybrid_lat(void *struct_obj, bool value_in);
void bbu_param_struct_get_write_voltage_vs_time_dat(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_write_voltage_vs_time_dat(void *struct_obj, bool value_in);
void bbu_param_struct_get_keep_overlays_and_groups(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_keep_overlays_and_groups(void *struct_obj, bool value_in);
void bbu_param_struct_get_keep_all_lcavities(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_keep_all_lcavities(void *struct_obj, bool value_in);
void bbu_param_struct_get_use_taylor_for_hybrids(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_use_taylor_for_hybrids(void *struct_obj, bool value_in);
void bbu_param_struct_get_stable_orbit_anal(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_stable_orbit_anal(void *struct_obj, bool value_in);
void bbu_param_struct_get_limit_factor(const void *struct_obj, double *value_out);
void bbu_param_struct_set_limit_factor(void *struct_obj, double value_in);
void bbu_param_struct_get_simulation_turns_max(const void *struct_obj, double *value_out);
void bbu_param_struct_set_simulation_turns_max(void *struct_obj, double value_in);
void bbu_param_struct_get_bunch_freq(const void *struct_obj, double *value_out);
void bbu_param_struct_set_bunch_freq(void *struct_obj, double value_in);
void bbu_param_struct_get_init_particle_offset(const void *struct_obj, double *value_out);
void bbu_param_struct_set_init_particle_offset(void *struct_obj, double value_in);
void bbu_param_struct_get_current(const void *struct_obj, double *value_out);
void bbu_param_struct_set_current(void *struct_obj, double value_in);
void bbu_param_struct_get_rel_tol(const void *struct_obj, double *value_out);
void bbu_param_struct_set_rel_tol(void *struct_obj, double value_in);
void bbu_param_struct_get_drscan(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_drscan(void *struct_obj, bool value_in);
void bbu_param_struct_get_use_interpolated_threshold(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_use_interpolated_threshold(void *struct_obj, bool value_in);
void bbu_param_struct_get_write_hom_info(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_write_hom_info(void *struct_obj, bool value_in);
void bbu_param_struct_get_elindex(const void *struct_obj, int *value_out);
void bbu_param_struct_set_elindex(void *struct_obj, int value_in);
void bbu_param_struct_get_elname_info(const void *s, char **d, int *bounds, bool *a);
void bbu_param_struct_set_elname(void *struct_obj, const char *str_ptr, int str_len);
void bbu_param_struct_get_nstep(const void *struct_obj, int *value_out);
void bbu_param_struct_set_nstep(void *struct_obj, int value_in);
void bbu_param_struct_get_begdr(const void *struct_obj, double *value_out);
void bbu_param_struct_set_begdr(void *struct_obj, double value_in);
void bbu_param_struct_get_enddr(const void *struct_obj, double *value_out);
void bbu_param_struct_set_enddr(void *struct_obj, double value_in);
void bbu_param_struct_get_nrep(const void *struct_obj, int *value_out);
void bbu_param_struct_set_nrep(void *struct_obj, int value_in);
void bbu_param_struct_get_ran_seed(const void *struct_obj, int *value_out);
void bbu_param_struct_set_ran_seed(void *struct_obj, int value_in);
void bbu_param_struct_get_hom_order_cutoff(const void *struct_obj, int *value_out);
void bbu_param_struct_set_hom_order_cutoff(void *struct_obj, int value_in);
void bbu_param_struct_get_ran_gauss_sigma_cut(const void *struct_obj, double *value_out);
void bbu_param_struct_set_ran_gauss_sigma_cut(void *struct_obj, double value_in);
void bbu_param_struct_get_ele_track_end_info(const void *s, char **d, int *bounds, bool *a);
void bbu_param_struct_set_ele_track_end(void *struct_obj, const char *str_ptr, int str_len);
void bbu_param_struct_get_ix_ele_track_end(const void *struct_obj, int *value_out);
void bbu_param_struct_set_ix_ele_track_end(void *struct_obj, int value_in);
void bbu_param_struct_get_regression(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_regression(void *struct_obj, bool value_in);
void bbu_param_struct_get_normalize_z_to_rf(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_normalize_z_to_rf(void *struct_obj, bool value_in);
void bbu_param_struct_get_ramp_on(const void *struct_obj, bool *value_out);
void bbu_param_struct_set_ramp_on(void *struct_obj, bool value_in);
void bbu_param_struct_get_ramp_pattern_info(const void *s, double **d, int *bounds, bool *is_alloc);
void bbu_param_struct_set_ramp_pattern(void *s, const void *d, const int *shape);
void bbu_param_struct_get_ramp_n_start(const void *struct_obj, int *value_out);
void bbu_param_struct_set_ramp_n_start(void *struct_obj, int value_in);
void bbu_param_struct_get_n_ramp_pattern(const void *struct_obj, int *value_out);
void bbu_param_struct_set_n_ramp_pattern(void *struct_obj, int value_in);
void fibre_get_DIR(const void *struct_obj, int **ptr_out);
void fibre_set_DIR(void *struct_obj, int value_in);
void fibre_get_PREVIOUS(const void *struct_obj, void **ptr_out);
void fibre_set_PREVIOUS(void *struct_obj, const void *src_ptr);
void fibre_get_NEXT(const void *struct_obj, void **ptr_out);
void fibre_set_NEXT(void *struct_obj, const void *src_ptr);
void fibre_get_PARENT_LAYOUT(const void *struct_obj, void **ptr_out);
void fibre_set_PARENT_LAYOUT(void *struct_obj, const void *src_ptr);
void fibre_get_pos(const void *struct_obj, int **ptr_out);
void fibre_set_pos(void *struct_obj, int value_in);
void fibre_get_BETA0(const void *struct_obj, double **ptr_out);
void fibre_set_BETA0(void *struct_obj, double value_in);
void fibre_get_GAMMA0I(const void *struct_obj, double **ptr_out);
void fibre_set_GAMMA0I(void *struct_obj, double value_in);
void fibre_get_GAMBET(const void *struct_obj, double **ptr_out);
void fibre_set_GAMBET(void *struct_obj, double value_in);
void fibre_get_MASS(const void *struct_obj, double **ptr_out);
void fibre_set_MASS(void *struct_obj, double value_in);
void fibre_get_CHARGE(const void *struct_obj, double **ptr_out);
void fibre_set_CHARGE(void *struct_obj, double value_in);
void fibre_get_AG(const void *struct_obj, double **ptr_out);
void fibre_set_AG(void *struct_obj, double value_in);
void fibre_get_P(const void *struct_obj, void **ptr_out);
void fibre_set_P(void *struct_obj, const void *src_ptr);
void fibre_get_N(const void *struct_obj, void **ptr_out);
void fibre_set_N(void *struct_obj, const void *src_ptr);
void fibre_get_loc(const void *struct_obj, int **ptr_out);
void fibre_set_loc(void *struct_obj, int value_in);

void layout_get_NAME_info(const void *s, char **d, int *len, bool *is_alloc);

void layout_set_NAME(void *struct_obj, const char *str_ptr, int str_len);
void layout_get_INDEX(const void *struct_obj, int **ptr_out);
void layout_set_INDEX(void *struct_obj, int value_in);
void layout_get_HARMONIC_NUMBER(const void *struct_obj, double **ptr_out);
void layout_set_HARMONIC_NUMBER(void *struct_obj, double value_in);
void layout_get_CLOSED(const void *struct_obj, bool **ptr_out);
void layout_set_CLOSED(void *struct_obj, bool value_in);
void layout_get_N(const void *struct_obj, int **ptr_out);
void layout_set_N(void *struct_obj, int value_in);
void layout_get_NTHIN(const void *struct_obj, int **ptr_out);
void layout_set_NTHIN(void *struct_obj, int value_in);
void layout_get_THIN(const void *struct_obj, double **ptr_out);
void layout_set_THIN(void *struct_obj, double value_in);
void layout_get_LASTPOS(const void *struct_obj, int **ptr_out);
void layout_set_LASTPOS(void *struct_obj, int value_in);
void layout_get_LAST(const void *struct_obj, void **ptr_out);
void layout_set_LAST(void *struct_obj, const void *src_ptr);
void layout_get_END(const void *struct_obj, void **ptr_out);
void layout_set_END(void *struct_obj, const void *src_ptr);
void layout_get_START(const void *struct_obj, void **ptr_out);
void layout_set_START(void *struct_obj, const void *src_ptr);
void layout_get_START_GROUND(const void *struct_obj, void **ptr_out);
void layout_set_START_GROUND(void *struct_obj, const void *src_ptr);
void layout_get_END_GROUND(const void *struct_obj, void **ptr_out);
void layout_set_END_GROUND(void *struct_obj, const void *src_ptr);
void layout_get_NEXT(const void *struct_obj, void **ptr_out);
void layout_set_NEXT(void *struct_obj, const void *src_ptr);
void layout_get_PREVIOUS(const void *struct_obj, void **ptr_out);
void layout_set_PREVIOUS(void *struct_obj, const void *src_ptr);
void all_encompassing_struct_get_real_rp_0d(const void *struct_obj, double *value_out);
void all_encompassing_struct_set_real_rp_0d(void *struct_obj, double value_in);
void all_encompassing_struct_get_real_rp_1d_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_real_rp_1d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_rp_2d_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_rp_2d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_rp_3d_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_rp_3d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_rp_0d_ptr(const void *struct_obj, double **ptr_out);
void all_encompassing_struct_set_real_rp_0d_ptr(void *struct_obj, double value_in);
void all_encompassing_struct_get_real_rp_1d_ptr_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_real_rp_1d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_rp_2d_ptr_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_rp_2d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_rp_3d_ptr_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_rp_3d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_rp_1d_alloc_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_real_rp_1d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_rp_2d_alloc_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_rp_2d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_rp_3d_alloc_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_rp_3d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_dp_0d(const void *struct_obj, double *value_out);
void all_encompassing_struct_set_real_dp_0d(void *struct_obj, double value_in);
void all_encompassing_struct_get_real_dp_1d_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_real_dp_1d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_dp_2d_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_dp_2d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_dp_3d_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_dp_3d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_dp_0d_ptr(const void *struct_obj, double **ptr_out);
void all_encompassing_struct_set_real_dp_0d_ptr(void *struct_obj, double value_in);
void all_encompassing_struct_get_real_dp_1d_ptr_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_real_dp_1d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_dp_2d_ptr_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_dp_2d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_dp_3d_ptr_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_dp_3d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_dp_1d_alloc_info(
    const void *s,
    double **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_real_dp_1d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_dp_2d_alloc_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_dp_2d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_real_dp_3d_alloc_info(
    const void *s,
    double **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_real_dp_3d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_complex_dp_0d(
    const void *struct_obj,
    std::complex<double> *value_out
);
void all_encompassing_struct_set_complex_dp_0d(void *struct_obj, std::complex<double> value_in);
void all_encompassing_struct_get_complex_dp_1d_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_complex_dp_1d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_complex_dp_2d_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_complex_dp_2d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_complex_dp_3d_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_complex_dp_3d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_complex_dp_1d_ptr_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_complex_dp_1d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_complex_dp_2d_ptr_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_complex_dp_2d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_complex_dp_3d_ptr_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_complex_dp_3d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_complex_dp_1d_alloc_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_complex_dp_1d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_complex_dp_2d_alloc_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_complex_dp_2d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_complex_dp_3d_alloc_info(
    const void *s,
    std::complex<double> **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_complex_dp_3d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int_0d(const void *struct_obj, int *value_out);
void all_encompassing_struct_set_int_0d(void *struct_obj, int value_in);
void all_encompassing_struct_get_int_1d_info(const void *s, int **d, int *bounds, bool *is_alloc);
void all_encompassing_struct_set_int_1d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int_2d_info(
    const void *s,
    int **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int_2d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int_3d_info(
    const void *s,
    int **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int_3d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int_0d_ptr(const void *struct_obj, int **ptr_out);
void all_encompassing_struct_set_int_0d_ptr(void *struct_obj, int value_in);
void all_encompassing_struct_get_int_1d_ptr_info(
    const void *s,
    int **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_int_1d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int_2d_ptr_info(
    const void *s,
    int **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int_2d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int_3d_ptr_info(
    const void *s,
    int **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int_3d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int_1d_alloc_info(
    const void *s,
    int **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_int_1d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int_2d_alloc_info(
    const void *s,
    int **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int_2d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int_3d_alloc_info(
    const void *s,
    int **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int_3d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int8_0d(const void *struct_obj, int64_t *value_out);
void all_encompassing_struct_set_int8_0d(void *struct_obj, int64_t value_in);
void all_encompassing_struct_get_int8_1d_info(
    const void *s,
    int64_t **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_int8_1d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int8_2d_info(
    const void *s,
    int64_t **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int8_2d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int8_3d_info(
    const void *s,
    int64_t **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int8_3d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int8_0d_ptr(const void *struct_obj, int64_t **ptr_out);
void all_encompassing_struct_set_int8_0d_ptr(void *struct_obj, int64_t value_in);
void all_encompassing_struct_get_int8_1d_ptr_info(
    const void *s,
    int64_t **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_int8_1d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int8_2d_ptr_info(
    const void *s,
    int64_t **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int8_2d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int8_3d_ptr_info(
    const void *s,
    int64_t **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int8_3d_ptr(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int8_1d_alloc_info(
    const void *s,
    int64_t **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_int8_1d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int8_2d_alloc_info(
    const void *s,
    int64_t **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int8_2d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_int8_3d_alloc_info(
    const void *s,
    int64_t **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_int8_3d_alloc(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_logical_0d(const void *struct_obj, bool *value_out);
void all_encompassing_struct_set_logical_0d(void *struct_obj, bool value_in);
void all_encompassing_struct_get_logical_1d_info(
    const void *s,
    bool **d,
    int *bounds,
    bool *is_alloc
);
void all_encompassing_struct_set_logical_1d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_logical_2d_info(
    const void *s,
    bool **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_logical_2d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_logical_3d_info(
    const void *s,
    bool **d,
    int *bounds,
    int *strides,
    bool *is_alloc
);
void all_encompassing_struct_set_logical_3d(void *s, const void *d, const int *shape);
void all_encompassing_struct_get_logical_0d_ptr(const void *struct_obj, bool **ptr_out);
void all_encompassing_struct_set_logical_0d_ptr(void *struct_obj, bool value_in);
void all_encompassing_struct_get_type_0d(const void *struct_obj, void **ptr_out);
void all_encompassing_struct_set_type_0d(void *struct_obj, const void *src_ptr);

void all_encompassing_struct_get_type_1d_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void all_encompassing_struct_get_type_2d_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void all_encompassing_struct_get_type_3d_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void all_encompassing_struct_get_type_0d_ptr(const void *struct_obj, void **ptr_out);
void all_encompassing_struct_set_type_0d_ptr(void *struct_obj, const void *src_ptr);

void all_encompassing_struct_get_type_1d_ptr_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void all_encompassing_struct_get_type_2d_ptr_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void all_encompassing_struct_get_type_3d_ptr_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void all_encompassing_struct_get_type_1d_alloc_info(
    const void *s,
    void **d,
    int *bounds,
    bool *is_alloc,
    size_t *el_size
);

void all_encompassing_struct_get_type_2d_alloc_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void all_encompassing_struct_get_type_3d_alloc_info(
    const void *s,
    void **d,
    int *bounds,
    int *strides,
    bool *a,
    size_t *es
);

void test_sub_struct_get_sr(const void *struct_obj, void **ptr_out);
void test_sub_struct_set_sr(void *struct_obj, const void *src_ptr);
void test_sub_sub_struct_get_aaa(const void *struct_obj, int64_t *value_out);
void test_sub_sub_struct_set_aaa(void *struct_obj, int64_t value_in);
void test_sub_sub_struct_get_bbb(const void *struct_obj, int *value_out);
void test_sub_sub_struct_set_bbb(void *struct_obj, int value_in);
void test_sub_sub_struct_get_file_info(const void *s, char **d, int *bounds, bool *a);
void test_sub_sub_struct_set_file(void *struct_obj, const char *str_ptr, int str_len);
void test_sub_sub_struct_get_t_ref(const void *struct_obj, double *value_out);
void test_sub_sub_struct_set_t_ref(void *struct_obj, double value_in);
void test_sub_sub_struct_get_freq_spread(const void *struct_obj, double *value_out);
void test_sub_sub_struct_set_freq_spread(void *struct_obj, double value_in);
}

namespace Bmad {

extern "C" {

void *allocate_real_container();
void reallocate_real_container_data(void *, int, size_t) noexcept;
void deallocate_real_container(void *) noexcept;
void access_real_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_real16_container();
void reallocate_real16_container_data(void *, int, size_t) noexcept;
void deallocate_real16_container(void *) noexcept;
void access_real16_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_integer_container();
void reallocate_integer_container_data(void *, int, size_t) noexcept;
void deallocate_integer_container(void *) noexcept;
void access_integer_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_integer8_container();
void reallocate_integer8_container_data(void *, int, size_t) noexcept;
void deallocate_integer8_container(void *) noexcept;
void access_integer8_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_logical_container();
void reallocate_logical_container_data(void *, int, size_t) noexcept;
void deallocate_logical_container(void *) noexcept;
void access_logical_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_complex_container();
void reallocate_complex_container_data(void *, int, size_t) noexcept;
void deallocate_complex_container(void *) noexcept;
void access_complex_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_spline_struct(int n, size_t *element_size);
void deallocate_fortran_spline_struct(void *ptr, int n) noexcept;
void copy_fortran_spline_struct(const void *src, void *dst);

void *allocate_spline_struct_container();
void reallocate_spline_struct_container_data(void *, int, size_t) noexcept;
void deallocate_spline_struct_container(void *) noexcept;
void access_spline_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_spin_polar_struct(int n, size_t *element_size);
void deallocate_fortran_spin_polar_struct(void *ptr, int n) noexcept;
void copy_fortran_spin_polar_struct(const void *src, void *dst);

void *allocate_spin_polar_struct_container();
void reallocate_spin_polar_struct_container_data(void *, int, size_t) noexcept;
void deallocate_spin_polar_struct_container(void *) noexcept;
void access_spin_polar_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_ac_kicker_time_struct(int n, size_t *element_size);
void deallocate_fortran_ac_kicker_time_struct(void *ptr, int n) noexcept;
void copy_fortran_ac_kicker_time_struct(const void *src, void *dst);

void *allocate_ac_kicker_time_struct_container();
void reallocate_ac_kicker_time_struct_container_data(void *, int, size_t) noexcept;
void deallocate_ac_kicker_time_struct_container(void *) noexcept;
void access_ac_kicker_time_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_ac_kicker_freq_struct(int n, size_t *element_size);
void deallocate_fortran_ac_kicker_freq_struct(void *ptr, int n) noexcept;
void copy_fortran_ac_kicker_freq_struct(const void *src, void *dst);

void *allocate_ac_kicker_freq_struct_container();
void reallocate_ac_kicker_freq_struct_container_data(void *, int, size_t) noexcept;
void deallocate_ac_kicker_freq_struct_container(void *) noexcept;
void access_ac_kicker_freq_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_ac_kicker_struct(int n, size_t *element_size);
void deallocate_fortran_ac_kicker_struct(void *ptr, int n) noexcept;
void copy_fortran_ac_kicker_struct(const void *src, void *dst);

void *allocate_ac_kicker_struct_container();
void reallocate_ac_kicker_struct_container_data(void *, int, size_t) noexcept;
void deallocate_ac_kicker_struct_container(void *) noexcept;
void access_ac_kicker_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_interval1_coef_struct(int n, size_t *element_size);
void deallocate_fortran_interval1_coef_struct(void *ptr, int n) noexcept;
void copy_fortran_interval1_coef_struct(const void *src, void *dst);

void *allocate_interval1_coef_struct_container();
void reallocate_interval1_coef_struct_container_data(void *, int, size_t) noexcept;
void deallocate_interval1_coef_struct_container(void *) noexcept;
void access_interval1_coef_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_photon_reflect_table_struct(int n, size_t *element_size);
void deallocate_fortran_photon_reflect_table_struct(void *ptr, int n) noexcept;
void copy_fortran_photon_reflect_table_struct(const void *src, void *dst);

void *allocate_photon_reflect_table_struct_container();
void reallocate_photon_reflect_table_struct_container_data(void *, int, size_t) noexcept;
void deallocate_photon_reflect_table_struct_container(void *) noexcept;
void access_photon_reflect_table_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_photon_reflect_surface_struct(int n, size_t *element_size);
void deallocate_fortran_photon_reflect_surface_struct(void *ptr, int n) noexcept;
void copy_fortran_photon_reflect_surface_struct(const void *src, void *dst);

void *allocate_photon_reflect_surface_struct_container();
void reallocate_photon_reflect_surface_struct_container_data(void *, int, size_t) noexcept;
void deallocate_photon_reflect_surface_struct_container(void *) noexcept;
void access_photon_reflect_surface_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_coord_struct(int n, size_t *element_size);
void deallocate_fortran_coord_struct(void *ptr, int n) noexcept;
void copy_fortran_coord_struct(const void *src, void *dst);

void *allocate_coord_struct_container();
void reallocate_coord_struct_container_data(void *, int, size_t) noexcept;
void deallocate_coord_struct_container(void *) noexcept;
void access_coord_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_coord_array_struct(int n, size_t *element_size);
void deallocate_fortran_coord_array_struct(void *ptr, int n) noexcept;
void copy_fortran_coord_array_struct(const void *src, void *dst);

void *allocate_coord_array_struct_container();
void reallocate_coord_array_struct_container_data(void *, int, size_t) noexcept;
void deallocate_coord_array_struct_container(void *) noexcept;
void access_coord_array_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bpm_phase_coupling_struct(int n, size_t *element_size);
void deallocate_fortran_bpm_phase_coupling_struct(void *ptr, int n) noexcept;
void copy_fortran_bpm_phase_coupling_struct(const void *src, void *dst);

void *allocate_bpm_phase_coupling_struct_container();
void reallocate_bpm_phase_coupling_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bpm_phase_coupling_struct_container(void *) noexcept;
void access_bpm_phase_coupling_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_expression_atom_struct(int n, size_t *element_size);
void deallocate_fortran_expression_atom_struct(void *ptr, int n) noexcept;
void copy_fortran_expression_atom_struct(const void *src, void *dst);

void *allocate_expression_atom_struct_container();
void reallocate_expression_atom_struct_container_data(void *, int, size_t) noexcept;
void deallocate_expression_atom_struct_container(void *) noexcept;
void access_expression_atom_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_wake_sr_z_long_struct(int n, size_t *element_size);
void deallocate_fortran_wake_sr_z_long_struct(void *ptr, int n) noexcept;
void copy_fortran_wake_sr_z_long_struct(const void *src, void *dst);

void *allocate_wake_sr_z_long_struct_container();
void reallocate_wake_sr_z_long_struct_container_data(void *, int, size_t) noexcept;
void deallocate_wake_sr_z_long_struct_container(void *) noexcept;
void access_wake_sr_z_long_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_wake_sr_mode_struct(int n, size_t *element_size);
void deallocate_fortran_wake_sr_mode_struct(void *ptr, int n) noexcept;
void copy_fortran_wake_sr_mode_struct(const void *src, void *dst);

void *allocate_wake_sr_mode_struct_container();
void reallocate_wake_sr_mode_struct_container_data(void *, int, size_t) noexcept;
void deallocate_wake_sr_mode_struct_container(void *) noexcept;
void access_wake_sr_mode_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_wake_sr_struct(int n, size_t *element_size);
void deallocate_fortran_wake_sr_struct(void *ptr, int n) noexcept;
void copy_fortran_wake_sr_struct(const void *src, void *dst);

void *allocate_wake_sr_struct_container();
void reallocate_wake_sr_struct_container_data(void *, int, size_t) noexcept;
void deallocate_wake_sr_struct_container(void *) noexcept;
void access_wake_sr_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_wake_lr_mode_struct(int n, size_t *element_size);
void deallocate_fortran_wake_lr_mode_struct(void *ptr, int n) noexcept;
void copy_fortran_wake_lr_mode_struct(const void *src, void *dst);

void *allocate_wake_lr_mode_struct_container();
void reallocate_wake_lr_mode_struct_container_data(void *, int, size_t) noexcept;
void deallocate_wake_lr_mode_struct_container(void *) noexcept;
void access_wake_lr_mode_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_wake_lr_struct(int n, size_t *element_size);
void deallocate_fortran_wake_lr_struct(void *ptr, int n) noexcept;
void copy_fortran_wake_lr_struct(const void *src, void *dst);

void *allocate_wake_lr_struct_container();
void reallocate_wake_lr_struct_container_data(void *, int, size_t) noexcept;
void deallocate_wake_lr_struct_container(void *) noexcept;
void access_wake_lr_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_lat_ele_loc_struct(int n, size_t *element_size);
void deallocate_fortran_lat_ele_loc_struct(void *ptr, int n) noexcept;
void copy_fortran_lat_ele_loc_struct(const void *src, void *dst);

void *allocate_lat_ele_loc_struct_container();
void reallocate_lat_ele_loc_struct_container_data(void *, int, size_t) noexcept;
void deallocate_lat_ele_loc_struct_container(void *) noexcept;
void access_lat_ele_loc_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_wake_struct(int n, size_t *element_size);
void deallocate_fortran_wake_struct(void *ptr, int n) noexcept;
void copy_fortran_wake_struct(const void *src, void *dst);

void *allocate_wake_struct_container();
void reallocate_wake_struct_container_data(void *, int, size_t) noexcept;
void deallocate_wake_struct_container(void *) noexcept;
void access_wake_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_taylor_term_struct(int n, size_t *element_size);
void deallocate_fortran_taylor_term_struct(void *ptr, int n) noexcept;
void copy_fortran_taylor_term_struct(const void *src, void *dst);

void *allocate_taylor_term_struct_container();
void reallocate_taylor_term_struct_container_data(void *, int, size_t) noexcept;
void deallocate_taylor_term_struct_container(void *) noexcept;
void access_taylor_term_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_taylor_struct(int n, size_t *element_size);
void deallocate_fortran_taylor_struct(void *ptr, int n) noexcept;
void copy_fortran_taylor_struct(const void *src, void *dst);

void *allocate_taylor_struct_container();
void reallocate_taylor_struct_container_data(void *, int, size_t) noexcept;
void deallocate_taylor_struct_container(void *) noexcept;
void access_taylor_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_em_taylor_term_struct(int n, size_t *element_size);
void deallocate_fortran_em_taylor_term_struct(void *ptr, int n) noexcept;
void copy_fortran_em_taylor_term_struct(const void *src, void *dst);

void *allocate_em_taylor_term_struct_container();
void reallocate_em_taylor_term_struct_container_data(void *, int, size_t) noexcept;
void deallocate_em_taylor_term_struct_container(void *) noexcept;
void access_em_taylor_term_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_em_taylor_struct(int n, size_t *element_size);
void deallocate_fortran_em_taylor_struct(void *ptr, int n) noexcept;
void copy_fortran_em_taylor_struct(const void *src, void *dst);

void *allocate_em_taylor_struct_container();
void reallocate_em_taylor_struct_container_data(void *, int, size_t) noexcept;
void deallocate_em_taylor_struct_container(void *) noexcept;
void access_em_taylor_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_cartesian_map_term1_struct(int n, size_t *element_size);
void deallocate_fortran_cartesian_map_term1_struct(void *ptr, int n) noexcept;
void copy_fortran_cartesian_map_term1_struct(const void *src, void *dst);

void *allocate_cartesian_map_term1_struct_container();
void reallocate_cartesian_map_term1_struct_container_data(void *, int, size_t) noexcept;
void deallocate_cartesian_map_term1_struct_container(void *) noexcept;
void access_cartesian_map_term1_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_cartesian_map_term_struct(int n, size_t *element_size);
void deallocate_fortran_cartesian_map_term_struct(void *ptr, int n) noexcept;
void copy_fortran_cartesian_map_term_struct(const void *src, void *dst);

void *allocate_cartesian_map_term_struct_container();
void reallocate_cartesian_map_term_struct_container_data(void *, int, size_t) noexcept;
void deallocate_cartesian_map_term_struct_container(void *) noexcept;
void access_cartesian_map_term_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_cartesian_map_struct(int n, size_t *element_size);
void deallocate_fortran_cartesian_map_struct(void *ptr, int n) noexcept;
void copy_fortran_cartesian_map_struct(const void *src, void *dst);

void *allocate_cartesian_map_struct_container();
void reallocate_cartesian_map_struct_container_data(void *, int, size_t) noexcept;
void deallocate_cartesian_map_struct_container(void *) noexcept;
void access_cartesian_map_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_cylindrical_map_term1_struct(int n, size_t *element_size);
void deallocate_fortran_cylindrical_map_term1_struct(void *ptr, int n) noexcept;
void copy_fortran_cylindrical_map_term1_struct(const void *src, void *dst);

void *allocate_cylindrical_map_term1_struct_container();
void reallocate_cylindrical_map_term1_struct_container_data(void *, int, size_t) noexcept;
void deallocate_cylindrical_map_term1_struct_container(void *) noexcept;
void access_cylindrical_map_term1_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_cylindrical_map_term_struct(int n, size_t *element_size);
void deallocate_fortran_cylindrical_map_term_struct(void *ptr, int n) noexcept;
void copy_fortran_cylindrical_map_term_struct(const void *src, void *dst);

void *allocate_cylindrical_map_term_struct_container();
void reallocate_cylindrical_map_term_struct_container_data(void *, int, size_t) noexcept;
void deallocate_cylindrical_map_term_struct_container(void *) noexcept;
void access_cylindrical_map_term_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_cylindrical_map_struct(int n, size_t *element_size);
void deallocate_fortran_cylindrical_map_struct(void *ptr, int n) noexcept;
void copy_fortran_cylindrical_map_struct(const void *src, void *dst);

void *allocate_cylindrical_map_struct_container();
void reallocate_cylindrical_map_struct_container_data(void *, int, size_t) noexcept;
void deallocate_cylindrical_map_struct_container(void *) noexcept;
void access_cylindrical_map_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bicubic_cmplx_coef_struct(int n, size_t *element_size);
void deallocate_fortran_bicubic_cmplx_coef_struct(void *ptr, int n) noexcept;
void copy_fortran_bicubic_cmplx_coef_struct(const void *src, void *dst);

void *allocate_bicubic_cmplx_coef_struct_container();
void reallocate_bicubic_cmplx_coef_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bicubic_cmplx_coef_struct_container(void *) noexcept;
void access_bicubic_cmplx_coef_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tricubic_cmplx_coef_struct(int n, size_t *element_size);
void deallocate_fortran_tricubic_cmplx_coef_struct(void *ptr, int n) noexcept;
void copy_fortran_tricubic_cmplx_coef_struct(const void *src, void *dst);

void *allocate_tricubic_cmplx_coef_struct_container();
void reallocate_tricubic_cmplx_coef_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tricubic_cmplx_coef_struct_container(void *) noexcept;
void access_tricubic_cmplx_coef_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_grid_field_pt1_struct(int n, size_t *element_size);
void deallocate_fortran_grid_field_pt1_struct(void *ptr, int n) noexcept;
void copy_fortran_grid_field_pt1_struct(const void *src, void *dst);

void *allocate_grid_field_pt1_struct_container();
void reallocate_grid_field_pt1_struct_container_data(void *, int, size_t) noexcept;
void deallocate_grid_field_pt1_struct_container(void *) noexcept;
void access_grid_field_pt1_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_grid_field_pt_struct(int n, size_t *element_size);
void deallocate_fortran_grid_field_pt_struct(void *ptr, int n) noexcept;
void copy_fortran_grid_field_pt_struct(const void *src, void *dst);

void *allocate_grid_field_pt_struct_container();
void reallocate_grid_field_pt_struct_container_data(void *, int, size_t) noexcept;
void deallocate_grid_field_pt_struct_container(void *) noexcept;
void access_grid_field_pt_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_grid_field_struct(int n, size_t *element_size);
void deallocate_fortran_grid_field_struct(void *ptr, int n) noexcept;
void copy_fortran_grid_field_struct(const void *src, void *dst);

void *allocate_grid_field_struct_container();
void reallocate_grid_field_struct_container_data(void *, int, size_t) noexcept;
void deallocate_grid_field_struct_container(void *) noexcept;
void access_grid_field_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_floor_position_struct(int n, size_t *element_size);
void deallocate_fortran_floor_position_struct(void *ptr, int n) noexcept;
void copy_fortran_floor_position_struct(const void *src, void *dst);

void *allocate_floor_position_struct_container();
void reallocate_floor_position_struct_container_data(void *, int, size_t) noexcept;
void deallocate_floor_position_struct_container(void *) noexcept;
void access_floor_position_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_high_energy_space_charge_struct(int n, size_t *element_size);
void deallocate_fortran_high_energy_space_charge_struct(void *ptr, int n) noexcept;
void copy_fortran_high_energy_space_charge_struct(const void *src, void *dst);

void *allocate_high_energy_space_charge_struct_container();
void reallocate_high_energy_space_charge_struct_container_data(void *, int, size_t) noexcept;
void deallocate_high_energy_space_charge_struct_container(void *) noexcept;
void access_high_energy_space_charge_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_xy_disp_struct(int n, size_t *element_size);
void deallocate_fortran_xy_disp_struct(void *ptr, int n) noexcept;
void copy_fortran_xy_disp_struct(const void *src, void *dst);

void *allocate_xy_disp_struct_container();
void reallocate_xy_disp_struct_container_data(void *, int, size_t) noexcept;
void deallocate_xy_disp_struct_container(void *) noexcept;
void access_xy_disp_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_twiss_struct(int n, size_t *element_size);
void deallocate_fortran_twiss_struct(void *ptr, int n) noexcept;
void copy_fortran_twiss_struct(const void *src, void *dst);

void *allocate_twiss_struct_container();
void reallocate_twiss_struct_container_data(void *, int, size_t) noexcept;
void deallocate_twiss_struct_container(void *) noexcept;
void access_twiss_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_mode3_struct(int n, size_t *element_size);
void deallocate_fortran_mode3_struct(void *ptr, int n) noexcept;
void copy_fortran_mode3_struct(const void *src, void *dst);

void *allocate_mode3_struct_container();
void reallocate_mode3_struct_container_data(void *, int, size_t) noexcept;
void deallocate_mode3_struct_container(void *) noexcept;
void access_mode3_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bookkeeping_state_struct(int n, size_t *element_size);
void deallocate_fortran_bookkeeping_state_struct(void *ptr, int n) noexcept;
void copy_fortran_bookkeeping_state_struct(const void *src, void *dst);

void *allocate_bookkeeping_state_struct_container();
void reallocate_bookkeeping_state_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bookkeeping_state_struct_container(void *) noexcept;
void access_bookkeeping_state_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_rad_map_struct(int n, size_t *element_size);
void deallocate_fortran_rad_map_struct(void *ptr, int n) noexcept;
void copy_fortran_rad_map_struct(const void *src, void *dst);

void *allocate_rad_map_struct_container();
void reallocate_rad_map_struct_container_data(void *, int, size_t) noexcept;
void deallocate_rad_map_struct_container(void *) noexcept;
void access_rad_map_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_rad_map_ele_struct(int n, size_t *element_size);
void deallocate_fortran_rad_map_ele_struct(void *ptr, int n) noexcept;
void copy_fortran_rad_map_ele_struct(const void *src, void *dst);

void *allocate_rad_map_ele_struct_container();
void reallocate_rad_map_ele_struct_container_data(void *, int, size_t) noexcept;
void deallocate_rad_map_ele_struct_container(void *) noexcept;
void access_rad_map_ele_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_gen_grad1_struct(int n, size_t *element_size);
void deallocate_fortran_gen_grad1_struct(void *ptr, int n) noexcept;
void copy_fortran_gen_grad1_struct(const void *src, void *dst);

void *allocate_gen_grad1_struct_container();
void reallocate_gen_grad1_struct_container_data(void *, int, size_t) noexcept;
void deallocate_gen_grad1_struct_container(void *) noexcept;
void access_gen_grad1_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_gen_grad_map_struct(int n, size_t *element_size);
void deallocate_fortran_gen_grad_map_struct(void *ptr, int n) noexcept;
void copy_fortran_gen_grad_map_struct(const void *src, void *dst);

void *allocate_gen_grad_map_struct_container();
void reallocate_gen_grad_map_struct_container_data(void *, int, size_t) noexcept;
void deallocate_gen_grad_map_struct_container(void *) noexcept;
void access_gen_grad_map_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_surface_segmented_pt_struct(int n, size_t *element_size);
void deallocate_fortran_surface_segmented_pt_struct(void *ptr, int n) noexcept;
void copy_fortran_surface_segmented_pt_struct(const void *src, void *dst);

void *allocate_surface_segmented_pt_struct_container();
void reallocate_surface_segmented_pt_struct_container_data(void *, int, size_t) noexcept;
void deallocate_surface_segmented_pt_struct_container(void *) noexcept;
void access_surface_segmented_pt_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_surface_segmented_struct(int n, size_t *element_size);
void deallocate_fortran_surface_segmented_struct(void *ptr, int n) noexcept;
void copy_fortran_surface_segmented_struct(const void *src, void *dst);

void *allocate_surface_segmented_struct_container();
void reallocate_surface_segmented_struct_container_data(void *, int, size_t) noexcept;
void deallocate_surface_segmented_struct_container(void *) noexcept;
void access_surface_segmented_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_surface_h_misalign_pt_struct(int n, size_t *element_size);
void deallocate_fortran_surface_h_misalign_pt_struct(void *ptr, int n) noexcept;
void copy_fortran_surface_h_misalign_pt_struct(const void *src, void *dst);

void *allocate_surface_h_misalign_pt_struct_container();
void reallocate_surface_h_misalign_pt_struct_container_data(void *, int, size_t) noexcept;
void deallocate_surface_h_misalign_pt_struct_container(void *) noexcept;
void access_surface_h_misalign_pt_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_surface_h_misalign_struct(int n, size_t *element_size);
void deallocate_fortran_surface_h_misalign_struct(void *ptr, int n) noexcept;
void copy_fortran_surface_h_misalign_struct(const void *src, void *dst);

void *allocate_surface_h_misalign_struct_container();
void reallocate_surface_h_misalign_struct_container_data(void *, int, size_t) noexcept;
void deallocate_surface_h_misalign_struct_container(void *) noexcept;
void access_surface_h_misalign_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_surface_displacement_pt_struct(int n, size_t *element_size);
void deallocate_fortran_surface_displacement_pt_struct(void *ptr, int n) noexcept;
void copy_fortran_surface_displacement_pt_struct(const void *src, void *dst);

void *allocate_surface_displacement_pt_struct_container();
void reallocate_surface_displacement_pt_struct_container_data(void *, int, size_t) noexcept;
void deallocate_surface_displacement_pt_struct_container(void *) noexcept;
void access_surface_displacement_pt_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_surface_displacement_struct(int n, size_t *element_size);
void deallocate_fortran_surface_displacement_struct(void *ptr, int n) noexcept;
void copy_fortran_surface_displacement_struct(const void *src, void *dst);

void *allocate_surface_displacement_struct_container();
void reallocate_surface_displacement_struct_container_data(void *, int, size_t) noexcept;
void deallocate_surface_displacement_struct_container(void *) noexcept;
void access_surface_displacement_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_target_point_struct(int n, size_t *element_size);
void deallocate_fortran_target_point_struct(void *ptr, int n) noexcept;
void copy_fortran_target_point_struct(const void *src, void *dst);

void *allocate_target_point_struct_container();
void reallocate_target_point_struct_container_data(void *, int, size_t) noexcept;
void deallocate_target_point_struct_container(void *) noexcept;
void access_target_point_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_surface_curvature_struct(int n, size_t *element_size);
void deallocate_fortran_surface_curvature_struct(void *ptr, int n) noexcept;
void copy_fortran_surface_curvature_struct(const void *src, void *dst);

void *allocate_surface_curvature_struct_container();
void reallocate_surface_curvature_struct_container_data(void *, int, size_t) noexcept;
void deallocate_surface_curvature_struct_container(void *) noexcept;
void access_surface_curvature_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_photon_target_struct(int n, size_t *element_size);
void deallocate_fortran_photon_target_struct(void *ptr, int n) noexcept;
void copy_fortran_photon_target_struct(const void *src, void *dst);

void *allocate_photon_target_struct_container();
void reallocate_photon_target_struct_container_data(void *, int, size_t) noexcept;
void deallocate_photon_target_struct_container(void *) noexcept;
void access_photon_target_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_photon_material_struct(int n, size_t *element_size);
void deallocate_fortran_photon_material_struct(void *ptr, int n) noexcept;
void copy_fortran_photon_material_struct(const void *src, void *dst);

void *allocate_photon_material_struct_container();
void reallocate_photon_material_struct_container_data(void *, int, size_t) noexcept;
void deallocate_photon_material_struct_container(void *) noexcept;
void access_photon_material_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_pixel_pt_struct(int n, size_t *element_size);
void deallocate_fortran_pixel_pt_struct(void *ptr, int n) noexcept;
void copy_fortran_pixel_pt_struct(const void *src, void *dst);

void *allocate_pixel_pt_struct_container();
void reallocate_pixel_pt_struct_container_data(void *, int, size_t) noexcept;
void deallocate_pixel_pt_struct_container(void *) noexcept;
void access_pixel_pt_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_pixel_detec_struct(int n, size_t *element_size);
void deallocate_fortran_pixel_detec_struct(void *ptr, int n) noexcept;
void copy_fortran_pixel_detec_struct(const void *src, void *dst);

void *allocate_pixel_detec_struct_container();
void reallocate_pixel_detec_struct_container_data(void *, int, size_t) noexcept;
void deallocate_pixel_detec_struct_container(void *) noexcept;
void access_pixel_detec_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_photon_element_struct(int n, size_t *element_size);
void deallocate_fortran_photon_element_struct(void *ptr, int n) noexcept;
void copy_fortran_photon_element_struct(const void *src, void *dst);

void *allocate_photon_element_struct_container();
void reallocate_photon_element_struct_container_data(void *, int, size_t) noexcept;
void deallocate_photon_element_struct_container(void *) noexcept;
void access_photon_element_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_wall3d_vertex_struct(int n, size_t *element_size);
void deallocate_fortran_wall3d_vertex_struct(void *ptr, int n) noexcept;
void copy_fortran_wall3d_vertex_struct(const void *src, void *dst);

void *allocate_wall3d_vertex_struct_container();
void reallocate_wall3d_vertex_struct_container_data(void *, int, size_t) noexcept;
void deallocate_wall3d_vertex_struct_container(void *) noexcept;
void access_wall3d_vertex_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_wall3d_section_struct(int n, size_t *element_size);
void deallocate_fortran_wall3d_section_struct(void *ptr, int n) noexcept;
void copy_fortran_wall3d_section_struct(const void *src, void *dst);

void *allocate_wall3d_section_struct_container();
void reallocate_wall3d_section_struct_container_data(void *, int, size_t) noexcept;
void deallocate_wall3d_section_struct_container(void *) noexcept;
void access_wall3d_section_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_wall3d_struct(int n, size_t *element_size);
void deallocate_fortran_wall3d_struct(void *ptr, int n) noexcept;
void copy_fortran_wall3d_struct(const void *src, void *dst);

void *allocate_wall3d_struct_container();
void reallocate_wall3d_struct_container_data(void *, int, size_t) noexcept;
void deallocate_wall3d_struct_container(void *) noexcept;
void access_wall3d_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_ramper_lord_struct(int n, size_t *element_size);
void deallocate_fortran_ramper_lord_struct(void *ptr, int n) noexcept;
void copy_fortran_ramper_lord_struct(const void *src, void *dst);

void *allocate_ramper_lord_struct_container();
void reallocate_ramper_lord_struct_container_data(void *, int, size_t) noexcept;
void deallocate_ramper_lord_struct_container(void *) noexcept;
void access_ramper_lord_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_control_struct(int n, size_t *element_size);
void deallocate_fortran_control_struct(void *ptr, int n) noexcept;
void copy_fortran_control_struct(const void *src, void *dst);

void *allocate_control_struct_container();
void reallocate_control_struct_container_data(void *, int, size_t) noexcept;
void deallocate_control_struct_container(void *) noexcept;
void access_control_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_control_var1_struct(int n, size_t *element_size);
void deallocate_fortran_control_var1_struct(void *ptr, int n) noexcept;
void copy_fortran_control_var1_struct(const void *src, void *dst);

void *allocate_control_var1_struct_container();
void reallocate_control_var1_struct_container_data(void *, int, size_t) noexcept;
void deallocate_control_var1_struct_container(void *) noexcept;
void access_control_var1_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_control_ramp1_struct(int n, size_t *element_size);
void deallocate_fortran_control_ramp1_struct(void *ptr, int n) noexcept;
void copy_fortran_control_ramp1_struct(const void *src, void *dst);

void *allocate_control_ramp1_struct_container();
void reallocate_control_ramp1_struct_container_data(void *, int, size_t) noexcept;
void deallocate_control_ramp1_struct_container(void *) noexcept;
void access_control_ramp1_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_controller_struct(int n, size_t *element_size);
void deallocate_fortran_controller_struct(void *ptr, int n) noexcept;
void copy_fortran_controller_struct(const void *src, void *dst);

void *allocate_controller_struct_container();
void reallocate_controller_struct_container_data(void *, int, size_t) noexcept;
void deallocate_controller_struct_container(void *) noexcept;
void access_controller_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_ellipse_beam_init_struct(int n, size_t *element_size);
void deallocate_fortran_ellipse_beam_init_struct(void *ptr, int n) noexcept;
void copy_fortran_ellipse_beam_init_struct(const void *src, void *dst);

void *allocate_ellipse_beam_init_struct_container();
void reallocate_ellipse_beam_init_struct_container_data(void *, int, size_t) noexcept;
void deallocate_ellipse_beam_init_struct_container(void *) noexcept;
void access_ellipse_beam_init_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_kv_beam_init_struct(int n, size_t *element_size);
void deallocate_fortran_kv_beam_init_struct(void *ptr, int n) noexcept;
void copy_fortran_kv_beam_init_struct(const void *src, void *dst);

void *allocate_kv_beam_init_struct_container();
void reallocate_kv_beam_init_struct_container_data(void *, int, size_t) noexcept;
void deallocate_kv_beam_init_struct_container(void *) noexcept;
void access_kv_beam_init_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_grid_beam_init_struct(int n, size_t *element_size);
void deallocate_fortran_grid_beam_init_struct(void *ptr, int n) noexcept;
void copy_fortran_grid_beam_init_struct(const void *src, void *dst);

void *allocate_grid_beam_init_struct_container();
void reallocate_grid_beam_init_struct_container_data(void *, int, size_t) noexcept;
void deallocate_grid_beam_init_struct_container(void *) noexcept;
void access_grid_beam_init_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_beam_init_struct(int n, size_t *element_size);
void deallocate_fortran_beam_init_struct(void *ptr, int n) noexcept;
void copy_fortran_beam_init_struct(const void *src, void *dst);

void *allocate_beam_init_struct_container();
void reallocate_beam_init_struct_container_data(void *, int, size_t) noexcept;
void deallocate_beam_init_struct_container(void *) noexcept;
void access_beam_init_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_lat_param_struct(int n, size_t *element_size);
void deallocate_fortran_lat_param_struct(void *ptr, int n) noexcept;
void copy_fortran_lat_param_struct(const void *src, void *dst);

void *allocate_lat_param_struct_container();
void reallocate_lat_param_struct_container_data(void *, int, size_t) noexcept;
void deallocate_lat_param_struct_container(void *) noexcept;
void access_lat_param_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_mode_info_struct(int n, size_t *element_size);
void deallocate_fortran_mode_info_struct(void *ptr, int n) noexcept;
void copy_fortran_mode_info_struct(const void *src, void *dst);

void *allocate_mode_info_struct_container();
void reallocate_mode_info_struct_container_data(void *, int, size_t) noexcept;
void deallocate_mode_info_struct_container(void *) noexcept;
void access_mode_info_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_pre_tracker_struct(int n, size_t *element_size);
void deallocate_fortran_pre_tracker_struct(void *ptr, int n) noexcept;
void copy_fortran_pre_tracker_struct(const void *src, void *dst);

void *allocate_pre_tracker_struct_container();
void reallocate_pre_tracker_struct_container_data(void *, int, size_t) noexcept;
void deallocate_pre_tracker_struct_container(void *) noexcept;
void access_pre_tracker_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_anormal_mode_struct(int n, size_t *element_size);
void deallocate_fortran_anormal_mode_struct(void *ptr, int n) noexcept;
void copy_fortran_anormal_mode_struct(const void *src, void *dst);

void *allocate_anormal_mode_struct_container();
void reallocate_anormal_mode_struct_container_data(void *, int, size_t) noexcept;
void deallocate_anormal_mode_struct_container(void *) noexcept;
void access_anormal_mode_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_linac_normal_mode_struct(int n, size_t *element_size);
void deallocate_fortran_linac_normal_mode_struct(void *ptr, int n) noexcept;
void copy_fortran_linac_normal_mode_struct(const void *src, void *dst);

void *allocate_linac_normal_mode_struct_container();
void reallocate_linac_normal_mode_struct_container_data(void *, int, size_t) noexcept;
void deallocate_linac_normal_mode_struct_container(void *) noexcept;
void access_linac_normal_mode_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_normal_modes_struct(int n, size_t *element_size);
void deallocate_fortran_normal_modes_struct(void *ptr, int n) noexcept;
void copy_fortran_normal_modes_struct(const void *src, void *dst);

void *allocate_normal_modes_struct_container();
void reallocate_normal_modes_struct_container_data(void *, int, size_t) noexcept;
void deallocate_normal_modes_struct_container(void *) noexcept;
void access_normal_modes_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_em_field_struct(int n, size_t *element_size);
void deallocate_fortran_em_field_struct(void *ptr, int n) noexcept;
void copy_fortran_em_field_struct(const void *src, void *dst);

void *allocate_em_field_struct_container();
void reallocate_em_field_struct_container_data(void *, int, size_t) noexcept;
void deallocate_em_field_struct_container(void *) noexcept;
void access_em_field_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_strong_beam_struct(int n, size_t *element_size);
void deallocate_fortran_strong_beam_struct(void *ptr, int n) noexcept;
void copy_fortran_strong_beam_struct(const void *src, void *dst);

void *allocate_strong_beam_struct_container();
void reallocate_strong_beam_struct_container_data(void *, int, size_t) noexcept;
void deallocate_strong_beam_struct_container(void *) noexcept;
void access_strong_beam_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_track_point_struct(int n, size_t *element_size);
void deallocate_fortran_track_point_struct(void *ptr, int n) noexcept;
void copy_fortran_track_point_struct(const void *src, void *dst);

void *allocate_track_point_struct_container();
void reallocate_track_point_struct_container_data(void *, int, size_t) noexcept;
void deallocate_track_point_struct_container(void *) noexcept;
void access_track_point_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_track_struct(int n, size_t *element_size);
void deallocate_fortran_track_struct(void *ptr, int n) noexcept;
void copy_fortran_track_struct(const void *src, void *dst);

void *allocate_track_struct_container();
void reallocate_track_struct_container_data(void *, int, size_t) noexcept;
void deallocate_track_struct_container(void *) noexcept;
void access_track_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_space_charge_common_struct(int n, size_t *element_size);
void deallocate_fortran_space_charge_common_struct(void *ptr, int n) noexcept;
void copy_fortran_space_charge_common_struct(const void *src, void *dst);

void *allocate_space_charge_common_struct_container();
void reallocate_space_charge_common_struct_container_data(void *, int, size_t) noexcept;
void deallocate_space_charge_common_struct_container(void *) noexcept;
void access_space_charge_common_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bmad_common_struct(int n, size_t *element_size);
void deallocate_fortran_bmad_common_struct(void *ptr, int n) noexcept;
void copy_fortran_bmad_common_struct(const void *src, void *dst);

void *allocate_bmad_common_struct_container();
void reallocate_bmad_common_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bmad_common_struct_container(void *) noexcept;
void access_bmad_common_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_rad_int1_struct(int n, size_t *element_size);
void deallocate_fortran_rad_int1_struct(void *ptr, int n) noexcept;
void copy_fortran_rad_int1_struct(const void *src, void *dst);

void *allocate_rad_int1_struct_container();
void reallocate_rad_int1_struct_container_data(void *, int, size_t) noexcept;
void deallocate_rad_int1_struct_container(void *) noexcept;
void access_rad_int1_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_rad_int_branch_struct(int n, size_t *element_size);
void deallocate_fortran_rad_int_branch_struct(void *ptr, int n) noexcept;
void copy_fortran_rad_int_branch_struct(const void *src, void *dst);

void *allocate_rad_int_branch_struct_container();
void reallocate_rad_int_branch_struct_container_data(void *, int, size_t) noexcept;
void deallocate_rad_int_branch_struct_container(void *) noexcept;
void access_rad_int_branch_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_rad_int_all_ele_struct(int n, size_t *element_size);
void deallocate_fortran_rad_int_all_ele_struct(void *ptr, int n) noexcept;
void copy_fortran_rad_int_all_ele_struct(const void *src, void *dst);

void *allocate_rad_int_all_ele_struct_container();
void reallocate_rad_int_all_ele_struct_container_data(void *, int, size_t) noexcept;
void deallocate_rad_int_all_ele_struct_container(void *) noexcept;
void access_rad_int_all_ele_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_rf_stair_step_struct(int n, size_t *element_size);
void deallocate_fortran_rf_stair_step_struct(void *ptr, int n) noexcept;
void copy_fortran_rf_stair_step_struct(const void *src, void *dst);

void *allocate_rf_stair_step_struct_container();
void reallocate_rf_stair_step_struct_container_data(void *, int, size_t) noexcept;
void deallocate_rf_stair_step_struct_container(void *) noexcept;
void access_rf_stair_step_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_rf_ele_struct(int n, size_t *element_size);
void deallocate_fortran_rf_ele_struct(void *ptr, int n) noexcept;
void copy_fortran_rf_ele_struct(const void *src, void *dst);

void *allocate_rf_ele_struct_container();
void reallocate_rf_ele_struct_container_data(void *, int, size_t) noexcept;
void deallocate_rf_ele_struct_container(void *) noexcept;
void access_rf_ele_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_ele_struct(int n, size_t *element_size);
void deallocate_fortran_ele_struct(void *ptr, int n) noexcept;
void copy_fortran_ele_struct(const void *src, void *dst);

void *allocate_ele_struct_container();
void reallocate_ele_struct_container_data(void *, int, size_t) noexcept;
void deallocate_ele_struct_container(void *) noexcept;
void access_ele_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_complex_taylor_term_struct(int n, size_t *element_size);
void deallocate_fortran_complex_taylor_term_struct(void *ptr, int n) noexcept;
void copy_fortran_complex_taylor_term_struct(const void *src, void *dst);

void *allocate_complex_taylor_term_struct_container();
void reallocate_complex_taylor_term_struct_container_data(void *, int, size_t) noexcept;
void deallocate_complex_taylor_term_struct_container(void *) noexcept;
void access_complex_taylor_term_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_complex_taylor_struct(int n, size_t *element_size);
void deallocate_fortran_complex_taylor_struct(void *ptr, int n) noexcept;
void copy_fortran_complex_taylor_struct(const void *src, void *dst);

void *allocate_complex_taylor_struct_container();
void reallocate_complex_taylor_struct_container_data(void *, int, size_t) noexcept;
void deallocate_complex_taylor_struct_container(void *) noexcept;
void access_complex_taylor_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_branch_struct(int n, size_t *element_size);
void deallocate_fortran_branch_struct(void *ptr, int n) noexcept;
void copy_fortran_branch_struct(const void *src, void *dst);

void *allocate_branch_struct_container();
void reallocate_branch_struct_container_data(void *, int, size_t) noexcept;
void deallocate_branch_struct_container(void *) noexcept;
void access_branch_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_lat_struct(int n, size_t *element_size);
void deallocate_fortran_lat_struct(void *ptr, int n) noexcept;
void copy_fortran_lat_struct(const void *src, void *dst);

void *allocate_lat_struct_container();
void reallocate_lat_struct_container_data(void *, int, size_t) noexcept;
void deallocate_lat_struct_container(void *) noexcept;
void access_lat_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bunch_struct(int n, size_t *element_size);
void deallocate_fortran_bunch_struct(void *ptr, int n) noexcept;
void copy_fortran_bunch_struct(const void *src, void *dst);

void *allocate_bunch_struct_container();
void reallocate_bunch_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bunch_struct_container(void *) noexcept;
void access_bunch_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bunch_params_struct(int n, size_t *element_size);
void deallocate_fortran_bunch_params_struct(void *ptr, int n) noexcept;
void copy_fortran_bunch_params_struct(const void *src, void *dst);

void *allocate_bunch_params_struct_container();
void reallocate_bunch_params_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bunch_params_struct_container(void *) noexcept;
void access_bunch_params_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_beam_struct(int n, size_t *element_size);
void deallocate_fortran_beam_struct(void *ptr, int n) noexcept;
void copy_fortran_beam_struct(const void *src, void *dst);

void *allocate_beam_struct_container();
void reallocate_beam_struct_container_data(void *, int, size_t) noexcept;
void deallocate_beam_struct_container(void *) noexcept;
void access_beam_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_aperture_point_struct(int n, size_t *element_size);
void deallocate_fortran_aperture_point_struct(void *ptr, int n) noexcept;
void copy_fortran_aperture_point_struct(const void *src, void *dst);

void *allocate_aperture_point_struct_container();
void reallocate_aperture_point_struct_container_data(void *, int, size_t) noexcept;
void deallocate_aperture_point_struct_container(void *) noexcept;
void access_aperture_point_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_aperture_param_struct(int n, size_t *element_size);
void deallocate_fortran_aperture_param_struct(void *ptr, int n) noexcept;
void copy_fortran_aperture_param_struct(const void *src, void *dst);

void *allocate_aperture_param_struct_container();
void reallocate_aperture_param_struct_container_data(void *, int, size_t) noexcept;
void deallocate_aperture_param_struct_container(void *) noexcept;
void access_aperture_param_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_aperture_scan_struct(int n, size_t *element_size);
void deallocate_fortran_aperture_scan_struct(void *ptr, int n) noexcept;
void copy_fortran_aperture_scan_struct(const void *src, void *dst);

void *allocate_aperture_scan_struct_container();
void reallocate_aperture_scan_struct_container_data(void *, int, size_t) noexcept;
void deallocate_aperture_scan_struct_container(void *) noexcept;
void access_aperture_scan_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_ele_pointer_struct(int n, size_t *element_size);
void deallocate_fortran_ele_pointer_struct(void *ptr, int n) noexcept;
void copy_fortran_ele_pointer_struct(const void *src, void *dst);

void *allocate_ele_pointer_struct_container();
void reallocate_ele_pointer_struct_container_data(void *, int, size_t) noexcept;
void deallocate_ele_pointer_struct_container(void *) noexcept;
void access_ele_pointer_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_expression_tree_struct(int n, size_t *element_size);
void deallocate_fortran_expression_tree_struct(void *ptr, int n) noexcept;
void copy_fortran_expression_tree_struct(const void *src, void *dst);

void *allocate_expression_tree_struct_container();
void reallocate_expression_tree_struct_container_data(void *, int, size_t) noexcept;
void deallocate_expression_tree_struct_container(void *) noexcept;
void access_expression_tree_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_nametable_struct(int n, size_t *element_size);
void deallocate_fortran_nametable_struct(void *ptr, int n) noexcept;
void copy_fortran_nametable_struct(const void *src, void *dst);

void *allocate_nametable_struct_container();
void reallocate_nametable_struct_container_data(void *, int, size_t) noexcept;
void deallocate_nametable_struct_container(void *) noexcept;
void access_nametable_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_spin_dn_dpz_struct(int n, size_t *element_size);
void deallocate_fortran_tao_spin_dn_dpz_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_spin_dn_dpz_struct(const void *src, void *dst);

void *allocate_tao_spin_dn_dpz_struct_container();
void reallocate_tao_spin_dn_dpz_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_spin_dn_dpz_struct_container(void *) noexcept;
void access_tao_spin_dn_dpz_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_resonance_h_struct(int n, size_t *element_size);
void deallocate_fortran_resonance_h_struct(void *ptr, int n) noexcept;
void copy_fortran_resonance_h_struct(const void *src, void *dst);

void *allocate_resonance_h_struct_container();
void reallocate_resonance_h_struct_container_data(void *, int, size_t) noexcept;
void deallocate_resonance_h_struct_container(void *) noexcept;
void access_resonance_h_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_spin_orbit_map1_struct(int n, size_t *element_size);
void deallocate_fortran_spin_orbit_map1_struct(void *ptr, int n) noexcept;
void copy_fortran_spin_orbit_map1_struct(const void *src, void *dst);

void *allocate_spin_orbit_map1_struct_container();
void reallocate_spin_orbit_map1_struct_container_data(void *, int, size_t) noexcept;
void deallocate_spin_orbit_map1_struct_container(void *) noexcept;
void access_spin_orbit_map1_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_spin_axis_struct(int n, size_t *element_size);
void deallocate_fortran_spin_axis_struct(void *ptr, int n) noexcept;
void copy_fortran_spin_axis_struct(const void *src, void *dst);

void *allocate_spin_axis_struct_container();
void reallocate_spin_axis_struct_container_data(void *, int, size_t) noexcept;
void deallocate_spin_axis_struct_container(void *) noexcept;
void access_spin_axis_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_ptc_normal_form_struct(int n, size_t *element_size);
void deallocate_fortran_ptc_normal_form_struct(void *ptr, int n) noexcept;
void copy_fortran_ptc_normal_form_struct(const void *src, void *dst);

void *allocate_ptc_normal_form_struct_container();
void reallocate_ptc_normal_form_struct_container_data(void *, int, size_t) noexcept;
void deallocate_ptc_normal_form_struct_container(void *) noexcept;
void access_ptc_normal_form_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bmad_normal_form_struct(int n, size_t *element_size);
void deallocate_fortran_bmad_normal_form_struct(void *ptr, int n) noexcept;
void copy_fortran_bmad_normal_form_struct(const void *src, void *dst);

void *allocate_bmad_normal_form_struct_container();
void reallocate_bmad_normal_form_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bmad_normal_form_struct_container(void *) noexcept;
void access_bmad_normal_form_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bunch_track_struct(int n, size_t *element_size);
void deallocate_fortran_bunch_track_struct(void *ptr, int n) noexcept;
void copy_fortran_bunch_track_struct(const void *src, void *dst);

void *allocate_bunch_track_struct_container();
void reallocate_bunch_track_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bunch_track_struct_container(void *) noexcept;
void access_bunch_track_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_summation_rdt_struct(int n, size_t *element_size);
void deallocate_fortran_summation_rdt_struct(void *ptr, int n) noexcept;
void copy_fortran_summation_rdt_struct(const void *src, void *dst);

void *allocate_summation_rdt_struct_container();
void reallocate_summation_rdt_struct_container_data(void *, int, size_t) noexcept;
void deallocate_summation_rdt_struct_container(void *) noexcept;
void access_summation_rdt_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_ele_shape_struct(int n, size_t *element_size);
void deallocate_fortran_tao_ele_shape_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_ele_shape_struct(const void *src, void *dst);

void *allocate_tao_ele_shape_struct_container();
void reallocate_tao_ele_shape_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_ele_shape_struct_container(void *) noexcept;
void access_tao_ele_shape_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_ele_pointer_struct(int n, size_t *element_size);
void deallocate_fortran_tao_ele_pointer_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_ele_pointer_struct(const void *src, void *dst);

void *allocate_tao_ele_pointer_struct_container();
void reallocate_tao_ele_pointer_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_ele_pointer_struct_container(void *) noexcept;
void access_tao_ele_pointer_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_curve_struct(int n, size_t *element_size);
void deallocate_fortran_tao_curve_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_curve_struct(const void *src, void *dst);

void *allocate_tao_curve_struct_container();
void reallocate_tao_curve_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_curve_struct_container(void *) noexcept;
void access_tao_curve_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_curve_color_struct(int n, size_t *element_size);
void deallocate_fortran_tao_curve_color_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_curve_color_struct(const void *src, void *dst);

void *allocate_tao_curve_color_struct_container();
void reallocate_tao_curve_color_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_curve_color_struct_container(void *) noexcept;
void access_tao_curve_color_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_curve_orbit_struct(int n, size_t *element_size);
void deallocate_fortran_tao_curve_orbit_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_curve_orbit_struct(const void *src, void *dst);

void *allocate_tao_curve_orbit_struct_container();
void reallocate_tao_curve_orbit_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_curve_orbit_struct_container(void *) noexcept;
void access_tao_curve_orbit_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_histogram_struct(int n, size_t *element_size);
void deallocate_fortran_tao_histogram_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_histogram_struct(const void *src, void *dst);

void *allocate_tao_histogram_struct_container();
void reallocate_tao_histogram_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_histogram_struct_container(void *) noexcept;
void access_tao_histogram_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_lat_ele_order1_struct(int n, size_t *element_size);
void deallocate_fortran_lat_ele_order1_struct(void *ptr, int n) noexcept;
void copy_fortran_lat_ele_order1_struct(const void *src, void *dst);

void *allocate_lat_ele_order1_struct_container();
void reallocate_lat_ele_order1_struct_container_data(void *, int, size_t) noexcept;
void deallocate_lat_ele_order1_struct_container(void *) noexcept;
void access_lat_ele_order1_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_lat_ele_order_array_struct(int n, size_t *element_size);
void deallocate_fortran_lat_ele_order_array_struct(void *ptr, int n) noexcept;
void copy_fortran_lat_ele_order_array_struct(const void *src, void *dst);

void *allocate_lat_ele_order_array_struct_container();
void reallocate_lat_ele_order_array_struct_container_data(void *, int, size_t) noexcept;
void deallocate_lat_ele_order_array_struct_container(void *) noexcept;
void access_lat_ele_order_array_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_lat_sigma_struct(int n, size_t *element_size);
void deallocate_fortran_tao_lat_sigma_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_lat_sigma_struct(const void *src, void *dst);

void *allocate_tao_lat_sigma_struct_container();
void reallocate_tao_lat_sigma_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_lat_sigma_struct_container(void *) noexcept;
void access_tao_lat_sigma_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_spin_ele_struct(int n, size_t *element_size);
void deallocate_fortran_tao_spin_ele_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_spin_ele_struct(const void *src, void *dst);

void *allocate_tao_spin_ele_struct_container();
void reallocate_tao_spin_ele_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_spin_ele_struct_container(void *) noexcept;
void access_tao_spin_ele_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_plot_cache_struct(int n, size_t *element_size);
void deallocate_fortran_tao_plot_cache_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_plot_cache_struct(const void *src, void *dst);

void *allocate_tao_plot_cache_struct_container();
void reallocate_tao_plot_cache_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_plot_cache_struct_container(void *) noexcept;
void access_tao_plot_cache_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_spin_polarization_struct(int n, size_t *element_size);
void deallocate_fortran_tao_spin_polarization_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_spin_polarization_struct(const void *src, void *dst);

void *allocate_tao_spin_polarization_struct_container();
void reallocate_tao_spin_polarization_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_spin_polarization_struct_container(void *) noexcept;
void access_tao_spin_polarization_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_lattice_branch_struct(int n, size_t *element_size);
void deallocate_fortran_tao_lattice_branch_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_lattice_branch_struct(const void *src, void *dst);

void *allocate_tao_lattice_branch_struct_container();
void reallocate_tao_lattice_branch_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_lattice_branch_struct_container(void *) noexcept;
void access_tao_lattice_branch_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_model_element_struct(int n, size_t *element_size);
void deallocate_fortran_tao_model_element_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_model_element_struct(const void *src, void *dst);

void *allocate_tao_model_element_struct_container();
void reallocate_tao_model_element_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_model_element_struct_container(void *) noexcept;
void access_tao_model_element_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_beam_branch_struct(int n, size_t *element_size);
void deallocate_fortran_tao_beam_branch_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_beam_branch_struct(const void *src, void *dst);

void *allocate_tao_beam_branch_struct_container();
void reallocate_tao_beam_branch_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_beam_branch_struct_container(void *) noexcept;
void access_tao_beam_branch_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_d1_data_struct(int n, size_t *element_size);
void deallocate_fortran_tao_d1_data_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_d1_data_struct(const void *src, void *dst);

void *allocate_tao_d1_data_struct_container();
void reallocate_tao_d1_data_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_d1_data_struct_container(void *) noexcept;
void access_tao_d1_data_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_d2_data_struct(int n, size_t *element_size);
void deallocate_fortran_tao_d2_data_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_d2_data_struct(const void *src, void *dst);

void *allocate_tao_d2_data_struct_container();
void reallocate_tao_d2_data_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_d2_data_struct_container(void *) noexcept;
void access_tao_d2_data_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_data_var_component_struct(int n, size_t *element_size);
void deallocate_fortran_tao_data_var_component_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_data_var_component_struct(const void *src, void *dst);

void *allocate_tao_data_var_component_struct_container();
void reallocate_tao_data_var_component_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_data_var_component_struct_container(void *) noexcept;
void access_tao_data_var_component_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_graph_struct(int n, size_t *element_size);
void deallocate_fortran_tao_graph_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_graph_struct(const void *src, void *dst);

void *allocate_tao_graph_struct_container();
void reallocate_tao_graph_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_graph_struct_container(void *) noexcept;
void access_tao_graph_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_plot_struct(int n, size_t *element_size);
void deallocate_fortran_tao_plot_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_plot_struct(const void *src, void *dst);

void *allocate_tao_plot_struct_container();
void reallocate_tao_plot_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_plot_struct_container(void *) noexcept;
void access_tao_plot_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_plot_region_struct(int n, size_t *element_size);
void deallocate_fortran_tao_plot_region_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_plot_region_struct(const void *src, void *dst);

void *allocate_tao_plot_region_struct_container();
void reallocate_tao_plot_region_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_plot_region_struct_container(void *) noexcept;
void access_tao_plot_region_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_universe_pointer_struct(int n, size_t *element_size);
void deallocate_fortran_tao_universe_pointer_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_universe_pointer_struct(const void *src, void *dst);

void *allocate_tao_universe_pointer_struct_container();
void reallocate_tao_universe_pointer_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_universe_pointer_struct_container(void *) noexcept;
void access_tao_universe_pointer_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_super_universe_struct(int n, size_t *element_size);
void deallocate_fortran_tao_super_universe_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_super_universe_struct(const void *src, void *dst);

void *allocate_tao_super_universe_struct_container();
void reallocate_tao_super_universe_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_super_universe_struct_container(void *) noexcept;
void access_tao_super_universe_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_var_struct(int n, size_t *element_size);
void deallocate_fortran_tao_var_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_var_struct(const void *src, void *dst);

void *allocate_tao_var_struct_container();
void reallocate_tao_var_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_var_struct_container(void *) noexcept;
void access_tao_var_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_var_slave_struct(int n, size_t *element_size);
void deallocate_fortran_tao_var_slave_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_var_slave_struct(const void *src, void *dst);

void *allocate_tao_var_slave_struct_container();
void reallocate_tao_var_slave_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_var_slave_struct_container(void *) noexcept;
void access_tao_var_slave_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_lattice_struct(int n, size_t *element_size);
void deallocate_fortran_tao_lattice_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_lattice_struct(const void *src, void *dst);

void *allocate_tao_lattice_struct_container();
void reallocate_tao_lattice_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_lattice_struct_container(void *) noexcept;
void access_tao_lattice_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_beam_uni_struct(int n, size_t *element_size);
void deallocate_fortran_tao_beam_uni_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_beam_uni_struct(const void *src, void *dst);

void *allocate_tao_beam_uni_struct_container();
void reallocate_tao_beam_uni_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_beam_uni_struct_container(void *) noexcept;
void access_tao_beam_uni_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_dynamic_aperture_struct(int n, size_t *element_size);
void deallocate_fortran_tao_dynamic_aperture_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_dynamic_aperture_struct(const void *src, void *dst);

void *allocate_tao_dynamic_aperture_struct_container();
void reallocate_tao_dynamic_aperture_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_dynamic_aperture_struct_container(void *) noexcept;
void access_tao_dynamic_aperture_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_model_branch_struct(int n, size_t *element_size);
void deallocate_fortran_tao_model_branch_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_model_branch_struct(const void *src, void *dst);

void *allocate_tao_model_branch_struct_container();
void reallocate_tao_model_branch_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_model_branch_struct_container(void *) noexcept;
void access_tao_model_branch_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_spin_map_struct(int n, size_t *element_size);
void deallocate_fortran_tao_spin_map_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_spin_map_struct(const void *src, void *dst);

void *allocate_tao_spin_map_struct_container();
void reallocate_tao_spin_map_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_spin_map_struct_container(void *) noexcept;
void access_tao_spin_map_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_data_struct(int n, size_t *element_size);
void deallocate_fortran_tao_data_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_data_struct(const void *src, void *dst);

void *allocate_tao_data_struct_container();
void reallocate_tao_data_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_data_struct_container(void *) noexcept;
void access_tao_data_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_ping_scale_struct(int n, size_t *element_size);
void deallocate_fortran_tao_ping_scale_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_ping_scale_struct(const void *src, void *dst);

void *allocate_tao_ping_scale_struct_container();
void reallocate_tao_ping_scale_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_ping_scale_struct_container(void *) noexcept;
void access_tao_ping_scale_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_universe_calc_struct(int n, size_t *element_size);
void deallocate_fortran_tao_universe_calc_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_universe_calc_struct(const void *src, void *dst);

void *allocate_tao_universe_calc_struct_container();
void reallocate_tao_universe_calc_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_universe_calc_struct_container(void *) noexcept;
void access_tao_universe_calc_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_lat_ele_order_struct(int n, size_t *element_size);
void deallocate_fortran_lat_ele_order_struct(void *ptr, int n) noexcept;
void copy_fortran_lat_ele_order_struct(const void *src, void *dst);

void *allocate_lat_ele_order_struct_container();
void reallocate_lat_ele_order_struct_container_data(void *, int, size_t) noexcept;
void deallocate_lat_ele_order_struct_container(void *) noexcept;
void access_lat_ele_order_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_expression_info_struct(int n, size_t *element_size);
void deallocate_fortran_tao_expression_info_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_expression_info_struct(const void *src, void *dst);

void *allocate_tao_expression_info_struct_container();
void reallocate_tao_expression_info_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_expression_info_struct_container(void *) noexcept;
void access_tao_expression_info_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_eval_node_struct(int n, size_t *element_size);
void deallocate_fortran_tao_eval_node_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_eval_node_struct(const void *src, void *dst);

void *allocate_tao_eval_node_struct_container();
void reallocate_tao_eval_node_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_eval_node_struct_container(void *) noexcept;
void access_tao_eval_node_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_title_struct(int n, size_t *element_size);
void deallocate_fortran_tao_title_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_title_struct(const void *src, void *dst);

void *allocate_tao_title_struct_container();
void reallocate_tao_title_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_title_struct_container(void *) noexcept;
void access_tao_title_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_qp_rect_struct(int n, size_t *element_size);
void deallocate_fortran_qp_rect_struct(void *ptr, int n) noexcept;
void copy_fortran_qp_rect_struct(const void *src, void *dst);

void *allocate_qp_rect_struct_container();
void reallocate_qp_rect_struct_container_data(void *, int, size_t) noexcept;
void deallocate_qp_rect_struct_container(void *) noexcept;
void access_qp_rect_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_drawing_struct(int n, size_t *element_size);
void deallocate_fortran_tao_drawing_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_drawing_struct(const void *src, void *dst);

void *allocate_tao_drawing_struct_container();
void reallocate_tao_drawing_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_drawing_struct_container(void *) noexcept;
void access_tao_drawing_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_shape_pattern_struct(int n, size_t *element_size);
void deallocate_fortran_tao_shape_pattern_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_shape_pattern_struct(const void *src, void *dst);

void *allocate_tao_shape_pattern_struct_container();
void reallocate_tao_shape_pattern_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_shape_pattern_struct_container(void *) noexcept;
void access_tao_shape_pattern_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_shape_pattern_point_struct(int n, size_t *element_size);
void deallocate_fortran_tao_shape_pattern_point_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_shape_pattern_point_struct(const void *src, void *dst);

void *allocate_tao_shape_pattern_point_struct_container();
void reallocate_tao_shape_pattern_point_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_shape_pattern_point_struct_container(void *) noexcept;
void access_tao_shape_pattern_point_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_qp_axis_struct(int n, size_t *element_size);
void deallocate_fortran_qp_axis_struct(void *ptr, int n) noexcept;
void copy_fortran_qp_axis_struct(const void *src, void *dst);

void *allocate_qp_axis_struct_container();
void reallocate_qp_axis_struct_container_data(void *, int, size_t) noexcept;
void deallocate_qp_axis_struct_container(void *) noexcept;
void access_qp_axis_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_qp_legend_struct(int n, size_t *element_size);
void deallocate_fortran_qp_legend_struct(void *ptr, int n) noexcept;
void copy_fortran_qp_legend_struct(const void *src, void *dst);

void *allocate_qp_legend_struct_container();
void reallocate_qp_legend_struct_container_data(void *, int, size_t) noexcept;
void deallocate_qp_legend_struct_container(void *) noexcept;
void access_qp_legend_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_qp_point_struct(int n, size_t *element_size);
void deallocate_fortran_qp_point_struct(void *ptr, int n) noexcept;
void copy_fortran_qp_point_struct(const void *src, void *dst);

void *allocate_qp_point_struct_container();
void reallocate_qp_point_struct_container_data(void *, int, size_t) noexcept;
void deallocate_qp_point_struct_container(void *) noexcept;
void access_qp_point_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_qp_line_struct(int n, size_t *element_size);
void deallocate_fortran_qp_line_struct(void *ptr, int n) noexcept;
void copy_fortran_qp_line_struct(const void *src, void *dst);

void *allocate_qp_line_struct_container();
void reallocate_qp_line_struct_container_data(void *, int, size_t) noexcept;
void deallocate_qp_line_struct_container(void *) noexcept;
void access_qp_line_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_qp_symbol_struct(int n, size_t *element_size);
void deallocate_fortran_qp_symbol_struct(void *ptr, int n) noexcept;
void copy_fortran_qp_symbol_struct(const void *src, void *dst);

void *allocate_qp_symbol_struct_container();
void reallocate_qp_symbol_struct_container_data(void *, int, size_t) noexcept;
void deallocate_qp_symbol_struct_container(void *) noexcept;
void access_qp_symbol_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_floor_plan_struct(int n, size_t *element_size);
void deallocate_fortran_tao_floor_plan_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_floor_plan_struct(const void *src, void *dst);

void *allocate_tao_floor_plan_struct_container();
void reallocate_tao_floor_plan_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_floor_plan_struct_container(void *) noexcept;
void access_tao_floor_plan_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_v1_var_struct(int n, size_t *element_size);
void deallocate_fortran_tao_v1_var_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_v1_var_struct(const void *src, void *dst);

void *allocate_tao_v1_var_struct_container();
void reallocate_tao_v1_var_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_v1_var_struct_container(void *) noexcept;
void access_tao_v1_var_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_global_struct(int n, size_t *element_size);
void deallocate_fortran_tao_global_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_global_struct(const void *src, void *dst);

void *allocate_tao_global_struct_container();
void reallocate_tao_global_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_global_struct_container(void *) noexcept;
void access_tao_global_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_init_struct(int n, size_t *element_size);
void deallocate_fortran_tao_init_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_init_struct(const void *src, void *dst);

void *allocate_tao_init_struct_container();
void reallocate_tao_init_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_init_struct_container(void *) noexcept;
void access_tao_init_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_common_struct(int n, size_t *element_size);
void deallocate_fortran_tao_common_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_common_struct(const void *src, void *dst);

void *allocate_tao_common_struct_container();
void reallocate_tao_common_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_common_struct_container(void *) noexcept;
void access_tao_common_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_plot_page_struct(int n, size_t *element_size);
void deallocate_fortran_tao_plot_page_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_plot_page_struct(const void *src, void *dst);

void *allocate_tao_plot_page_struct_container();
void reallocate_tao_plot_page_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_plot_page_struct_container(void *) noexcept;
void access_tao_plot_page_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_building_wall_struct(int n, size_t *element_size);
void deallocate_fortran_tao_building_wall_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_building_wall_struct(const void *src, void *dst);

void *allocate_tao_building_wall_struct_container();
void reallocate_tao_building_wall_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_building_wall_struct_container(void *) noexcept;
void access_tao_building_wall_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_building_wall_orientation_struct(int n, size_t *element_size);
void deallocate_fortran_tao_building_wall_orientation_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_building_wall_orientation_struct(const void *src, void *dst);

void *allocate_tao_building_wall_orientation_struct_container();
void reallocate_tao_building_wall_orientation_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_building_wall_orientation_struct_container(void *) noexcept;
void access_tao_building_wall_orientation_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_building_wall_section_struct(int n, size_t *element_size);
void deallocate_fortran_tao_building_wall_section_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_building_wall_section_struct(const void *src, void *dst);

void *allocate_tao_building_wall_section_struct_container();
void reallocate_tao_building_wall_section_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_building_wall_section_struct_container(void *) noexcept;
void access_tao_building_wall_section_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_building_wall_point_struct(int n, size_t *element_size);
void deallocate_fortran_tao_building_wall_point_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_building_wall_point_struct(const void *src, void *dst);

void *allocate_tao_building_wall_point_struct_container();
void reallocate_tao_building_wall_point_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_building_wall_point_struct_container(void *) noexcept;
void access_tao_building_wall_point_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_wave_struct(int n, size_t *element_size);
void deallocate_fortran_tao_wave_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_wave_struct(const void *src, void *dst);

void *allocate_tao_wave_struct_container();
void reallocate_tao_wave_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_wave_struct_container(void *) noexcept;
void access_tao_wave_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_wave_kick_pt_struct(int n, size_t *element_size);
void deallocate_fortran_tao_wave_kick_pt_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_wave_kick_pt_struct(const void *src, void *dst);

void *allocate_tao_wave_kick_pt_struct_container();
void reallocate_tao_wave_kick_pt_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_wave_kick_pt_struct_container(void *) noexcept;
void access_tao_wave_kick_pt_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_cmd_history_struct(int n, size_t *element_size);
void deallocate_fortran_tao_cmd_history_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_cmd_history_struct(const void *src, void *dst);

void *allocate_tao_cmd_history_struct_container();
void reallocate_tao_cmd_history_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_cmd_history_struct_container(void *) noexcept;
void access_tao_cmd_history_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_tao_universe_struct(int n, size_t *element_size);
void deallocate_fortran_tao_universe_struct(void *ptr, int n) noexcept;
void copy_fortran_tao_universe_struct(const void *src, void *dst);

void *allocate_tao_universe_struct_container();
void reallocate_tao_universe_struct_container_data(void *, int, size_t) noexcept;
void deallocate_tao_universe_struct_container(void *) noexcept;
void access_tao_universe_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_mad_energy_struct(int n, size_t *element_size);
void deallocate_fortran_mad_energy_struct(void *ptr, int n) noexcept;
void copy_fortran_mad_energy_struct(const void *src, void *dst);

void *allocate_mad_energy_struct_container();
void reallocate_mad_energy_struct_container_data(void *, int, size_t) noexcept;
void deallocate_mad_energy_struct_container(void *) noexcept;
void access_mad_energy_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_mad_map_struct(int n, size_t *element_size);
void deallocate_fortran_mad_map_struct(void *ptr, int n) noexcept;
void copy_fortran_mad_map_struct(const void *src, void *dst);

void *allocate_mad_map_struct_container();
void reallocate_mad_map_struct_container_data(void *, int, size_t) noexcept;
void deallocate_mad_map_struct_container(void *) noexcept;
void access_mad_map_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_random_state_struct(int n, size_t *element_size);
void deallocate_fortran_random_state_struct(void *ptr, int n) noexcept;
void copy_fortran_random_state_struct(const void *src, void *dst);

void *allocate_random_state_struct_container();
void reallocate_random_state_struct_container_data(void *, int, size_t) noexcept;
void deallocate_random_state_struct_container(void *) noexcept;
void access_random_state_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bbu_stage_struct(int n, size_t *element_size);
void deallocate_fortran_bbu_stage_struct(void *ptr, int n) noexcept;
void copy_fortran_bbu_stage_struct(const void *src, void *dst);

void *allocate_bbu_stage_struct_container();
void reallocate_bbu_stage_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bbu_stage_struct_container(void *) noexcept;
void access_bbu_stage_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bbu_beam_struct(int n, size_t *element_size);
void deallocate_fortran_bbu_beam_struct(void *ptr, int n) noexcept;
void copy_fortran_bbu_beam_struct(const void *src, void *dst);

void *allocate_bbu_beam_struct_container();
void reallocate_bbu_beam_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bbu_beam_struct_container(void *) noexcept;
void access_bbu_beam_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_bbu_param_struct(int n, size_t *element_size);
void deallocate_fortran_bbu_param_struct(void *ptr, int n) noexcept;
void copy_fortran_bbu_param_struct(const void *src, void *dst);

void *allocate_bbu_param_struct_container();
void reallocate_bbu_param_struct_container_data(void *, int, size_t) noexcept;
void deallocate_bbu_param_struct_container(void *) noexcept;
void access_bbu_param_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_fibre(int n, size_t *element_size);
void deallocate_fortran_fibre(void *ptr, int n) noexcept;
void copy_fortran_fibre(const void *src, void *dst);

void *allocate_fibre_container();
void reallocate_fibre_container_data(void *, int, size_t) noexcept;
void deallocate_fibre_container(void *) noexcept;
void access_fibre_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_layout(int n, size_t *element_size);
void deallocate_fortran_layout(void *ptr, int n) noexcept;
void copy_fortran_layout(const void *src, void *dst);

void *allocate_layout_container();
void reallocate_layout_container_data(void *, int, size_t) noexcept;
void deallocate_layout_container(void *) noexcept;
void access_layout_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_all_encompassing_struct(int n, size_t *element_size);
void deallocate_fortran_all_encompassing_struct(void *ptr, int n) noexcept;
void copy_fortran_all_encompassing_struct(const void *src, void *dst);

void *allocate_all_encompassing_struct_container();
void reallocate_all_encompassing_struct_container_data(void *, int, size_t) noexcept;
void deallocate_all_encompassing_struct_container(void *) noexcept;
void access_all_encompassing_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_test_sub_struct(int n, size_t *element_size);
void deallocate_fortran_test_sub_struct(void *ptr, int n) noexcept;
void copy_fortran_test_sub_struct(const void *src, void *dst);

void *allocate_test_sub_struct_container();
void reallocate_test_sub_struct_container_data(void *, int, size_t) noexcept;
void deallocate_test_sub_struct_container(void *) noexcept;
void access_test_sub_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);

void *allocate_fortran_test_sub_sub_struct(int n, size_t *element_size);
void deallocate_fortran_test_sub_sub_struct(void *ptr, int n) noexcept;
void copy_fortran_test_sub_sub_struct(const void *src, void *dst);

void *allocate_test_sub_sub_struct_container();
void reallocate_test_sub_sub_struct_container_data(void *, int, size_t) noexcept;
void deallocate_test_sub_sub_struct_container(void *) noexcept;
void access_test_sub_sub_struct_container(
    void *handle,
    void **data,
    int *lbound,
    int *size,
    size_t *elem_size,
    bool *alloc
);
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

class SplineStruct;

using SplineStructArray1D =
    FTypeArray1D<SplineStruct, allocate_fortran_spline_struct, deallocate_fortran_spline_struct>;
using SplineStructArray2D = FTypeArray2D<SplineStruct>;
using SplineStructArray3D = FTypeArray3D<SplineStruct>;

using SplineStructAlloc1D = FTypeAlloc1D<
    SplineStructArray1D,
    allocate_spline_struct_container,
    deallocate_spline_struct_container,
    reallocate_spline_struct_container_data,
    access_spline_struct_container>;

class SpinPolarStruct;

using SpinPolarStructArray1D = FTypeArray1D<
    SpinPolarStruct,
    allocate_fortran_spin_polar_struct,
    deallocate_fortran_spin_polar_struct>;
using SpinPolarStructArray2D = FTypeArray2D<SpinPolarStruct>;
using SpinPolarStructArray3D = FTypeArray3D<SpinPolarStruct>;

using SpinPolarStructAlloc1D = FTypeAlloc1D<
    SpinPolarStructArray1D,
    allocate_spin_polar_struct_container,
    deallocate_spin_polar_struct_container,
    reallocate_spin_polar_struct_container_data,
    access_spin_polar_struct_container>;

class AcKickerTimeStruct;

using AcKickerTimeStructArray1D = FTypeArray1D<
    AcKickerTimeStruct,
    allocate_fortran_ac_kicker_time_struct,
    deallocate_fortran_ac_kicker_time_struct>;
using AcKickerTimeStructArray2D = FTypeArray2D<AcKickerTimeStruct>;
using AcKickerTimeStructArray3D = FTypeArray3D<AcKickerTimeStruct>;

using AcKickerTimeStructAlloc1D = FTypeAlloc1D<
    AcKickerTimeStructArray1D,
    allocate_ac_kicker_time_struct_container,
    deallocate_ac_kicker_time_struct_container,
    reallocate_ac_kicker_time_struct_container_data,
    access_ac_kicker_time_struct_container>;

class AcKickerFreqStruct;

using AcKickerFreqStructArray1D = FTypeArray1D<
    AcKickerFreqStruct,
    allocate_fortran_ac_kicker_freq_struct,
    deallocate_fortran_ac_kicker_freq_struct>;
using AcKickerFreqStructArray2D = FTypeArray2D<AcKickerFreqStruct>;
using AcKickerFreqStructArray3D = FTypeArray3D<AcKickerFreqStruct>;

using AcKickerFreqStructAlloc1D = FTypeAlloc1D<
    AcKickerFreqStructArray1D,
    allocate_ac_kicker_freq_struct_container,
    deallocate_ac_kicker_freq_struct_container,
    reallocate_ac_kicker_freq_struct_container_data,
    access_ac_kicker_freq_struct_container>;

class AcKickerStruct;

using AcKickerStructArray1D = FTypeArray1D<
    AcKickerStruct,
    allocate_fortran_ac_kicker_struct,
    deallocate_fortran_ac_kicker_struct>;
using AcKickerStructArray2D = FTypeArray2D<AcKickerStruct>;
using AcKickerStructArray3D = FTypeArray3D<AcKickerStruct>;

using AcKickerStructAlloc1D = FTypeAlloc1D<
    AcKickerStructArray1D,
    allocate_ac_kicker_struct_container,
    deallocate_ac_kicker_struct_container,
    reallocate_ac_kicker_struct_container_data,
    access_ac_kicker_struct_container>;

class Interval1CoefStruct;

using Interval1CoefStructArray1D = FTypeArray1D<
    Interval1CoefStruct,
    allocate_fortran_interval1_coef_struct,
    deallocate_fortran_interval1_coef_struct>;
using Interval1CoefStructArray2D = FTypeArray2D<Interval1CoefStruct>;
using Interval1CoefStructArray3D = FTypeArray3D<Interval1CoefStruct>;

using Interval1CoefStructAlloc1D = FTypeAlloc1D<
    Interval1CoefStructArray1D,
    allocate_interval1_coef_struct_container,
    deallocate_interval1_coef_struct_container,
    reallocate_interval1_coef_struct_container_data,
    access_interval1_coef_struct_container>;

class PhotonReflectTableStruct;

using PhotonReflectTableStructArray1D = FTypeArray1D<
    PhotonReflectTableStruct,
    allocate_fortran_photon_reflect_table_struct,
    deallocate_fortran_photon_reflect_table_struct>;
using PhotonReflectTableStructArray2D = FTypeArray2D<PhotonReflectTableStruct>;
using PhotonReflectTableStructArray3D = FTypeArray3D<PhotonReflectTableStruct>;

using PhotonReflectTableStructAlloc1D = FTypeAlloc1D<
    PhotonReflectTableStructArray1D,
    allocate_photon_reflect_table_struct_container,
    deallocate_photon_reflect_table_struct_container,
    reallocate_photon_reflect_table_struct_container_data,
    access_photon_reflect_table_struct_container>;

class PhotonReflectSurfaceStruct;

using PhotonReflectSurfaceStructArray1D = FTypeArray1D<
    PhotonReflectSurfaceStruct,
    allocate_fortran_photon_reflect_surface_struct,
    deallocate_fortran_photon_reflect_surface_struct>;
using PhotonReflectSurfaceStructArray2D = FTypeArray2D<PhotonReflectSurfaceStruct>;
using PhotonReflectSurfaceStructArray3D = FTypeArray3D<PhotonReflectSurfaceStruct>;

using PhotonReflectSurfaceStructAlloc1D = FTypeAlloc1D<
    PhotonReflectSurfaceStructArray1D,
    allocate_photon_reflect_surface_struct_container,
    deallocate_photon_reflect_surface_struct_container,
    reallocate_photon_reflect_surface_struct_container_data,
    access_photon_reflect_surface_struct_container>;

class CoordStruct;

using CoordStructArray1D =
    FTypeArray1D<CoordStruct, allocate_fortran_coord_struct, deallocate_fortran_coord_struct>;
using CoordStructArray2D = FTypeArray2D<CoordStruct>;
using CoordStructArray3D = FTypeArray3D<CoordStruct>;

using CoordStructAlloc1D = FTypeAlloc1D<
    CoordStructArray1D,
    allocate_coord_struct_container,
    deallocate_coord_struct_container,
    reallocate_coord_struct_container_data,
    access_coord_struct_container>;

class CoordArrayStruct;

using CoordArrayStructArray1D = FTypeArray1D<
    CoordArrayStruct,
    allocate_fortran_coord_array_struct,
    deallocate_fortran_coord_array_struct>;
using CoordArrayStructArray2D = FTypeArray2D<CoordArrayStruct>;
using CoordArrayStructArray3D = FTypeArray3D<CoordArrayStruct>;

using CoordArrayStructAlloc1D = FTypeAlloc1D<
    CoordArrayStructArray1D,
    allocate_coord_array_struct_container,
    deallocate_coord_array_struct_container,
    reallocate_coord_array_struct_container_data,
    access_coord_array_struct_container>;

class BpmPhaseCouplingStruct;

using BpmPhaseCouplingStructArray1D = FTypeArray1D<
    BpmPhaseCouplingStruct,
    allocate_fortran_bpm_phase_coupling_struct,
    deallocate_fortran_bpm_phase_coupling_struct>;
using BpmPhaseCouplingStructArray2D = FTypeArray2D<BpmPhaseCouplingStruct>;
using BpmPhaseCouplingStructArray3D = FTypeArray3D<BpmPhaseCouplingStruct>;

using BpmPhaseCouplingStructAlloc1D = FTypeAlloc1D<
    BpmPhaseCouplingStructArray1D,
    allocate_bpm_phase_coupling_struct_container,
    deallocate_bpm_phase_coupling_struct_container,
    reallocate_bpm_phase_coupling_struct_container_data,
    access_bpm_phase_coupling_struct_container>;

class ExpressionAtomStruct;

using ExpressionAtomStructArray1D = FTypeArray1D<
    ExpressionAtomStruct,
    allocate_fortran_expression_atom_struct,
    deallocate_fortran_expression_atom_struct>;
using ExpressionAtomStructArray2D = FTypeArray2D<ExpressionAtomStruct>;
using ExpressionAtomStructArray3D = FTypeArray3D<ExpressionAtomStruct>;

using ExpressionAtomStructAlloc1D = FTypeAlloc1D<
    ExpressionAtomStructArray1D,
    allocate_expression_atom_struct_container,
    deallocate_expression_atom_struct_container,
    reallocate_expression_atom_struct_container_data,
    access_expression_atom_struct_container>;

class WakeSrZLongStruct;

using WakeSrZLongStructArray1D = FTypeArray1D<
    WakeSrZLongStruct,
    allocate_fortran_wake_sr_z_long_struct,
    deallocate_fortran_wake_sr_z_long_struct>;
using WakeSrZLongStructArray2D = FTypeArray2D<WakeSrZLongStruct>;
using WakeSrZLongStructArray3D = FTypeArray3D<WakeSrZLongStruct>;

using WakeSrZLongStructAlloc1D = FTypeAlloc1D<
    WakeSrZLongStructArray1D,
    allocate_wake_sr_z_long_struct_container,
    deallocate_wake_sr_z_long_struct_container,
    reallocate_wake_sr_z_long_struct_container_data,
    access_wake_sr_z_long_struct_container>;

class WakeSrModeStruct;

using WakeSrModeStructArray1D = FTypeArray1D<
    WakeSrModeStruct,
    allocate_fortran_wake_sr_mode_struct,
    deallocate_fortran_wake_sr_mode_struct>;
using WakeSrModeStructArray2D = FTypeArray2D<WakeSrModeStruct>;
using WakeSrModeStructArray3D = FTypeArray3D<WakeSrModeStruct>;

using WakeSrModeStructAlloc1D = FTypeAlloc1D<
    WakeSrModeStructArray1D,
    allocate_wake_sr_mode_struct_container,
    deallocate_wake_sr_mode_struct_container,
    reallocate_wake_sr_mode_struct_container_data,
    access_wake_sr_mode_struct_container>;

class WakeSrStruct;

using WakeSrStructArray1D =
    FTypeArray1D<WakeSrStruct, allocate_fortran_wake_sr_struct, deallocate_fortran_wake_sr_struct>;
using WakeSrStructArray2D = FTypeArray2D<WakeSrStruct>;
using WakeSrStructArray3D = FTypeArray3D<WakeSrStruct>;

using WakeSrStructAlloc1D = FTypeAlloc1D<
    WakeSrStructArray1D,
    allocate_wake_sr_struct_container,
    deallocate_wake_sr_struct_container,
    reallocate_wake_sr_struct_container_data,
    access_wake_sr_struct_container>;

class WakeLrModeStruct;

using WakeLrModeStructArray1D = FTypeArray1D<
    WakeLrModeStruct,
    allocate_fortran_wake_lr_mode_struct,
    deallocate_fortran_wake_lr_mode_struct>;
using WakeLrModeStructArray2D = FTypeArray2D<WakeLrModeStruct>;
using WakeLrModeStructArray3D = FTypeArray3D<WakeLrModeStruct>;

using WakeLrModeStructAlloc1D = FTypeAlloc1D<
    WakeLrModeStructArray1D,
    allocate_wake_lr_mode_struct_container,
    deallocate_wake_lr_mode_struct_container,
    reallocate_wake_lr_mode_struct_container_data,
    access_wake_lr_mode_struct_container>;

class WakeLrStruct;

using WakeLrStructArray1D =
    FTypeArray1D<WakeLrStruct, allocate_fortran_wake_lr_struct, deallocate_fortran_wake_lr_struct>;
using WakeLrStructArray2D = FTypeArray2D<WakeLrStruct>;
using WakeLrStructArray3D = FTypeArray3D<WakeLrStruct>;

using WakeLrStructAlloc1D = FTypeAlloc1D<
    WakeLrStructArray1D,
    allocate_wake_lr_struct_container,
    deallocate_wake_lr_struct_container,
    reallocate_wake_lr_struct_container_data,
    access_wake_lr_struct_container>;

class LatEleLocStruct;

using LatEleLocStructArray1D = FTypeArray1D<
    LatEleLocStruct,
    allocate_fortran_lat_ele_loc_struct,
    deallocate_fortran_lat_ele_loc_struct>;
using LatEleLocStructArray2D = FTypeArray2D<LatEleLocStruct>;
using LatEleLocStructArray3D = FTypeArray3D<LatEleLocStruct>;

using LatEleLocStructAlloc1D = FTypeAlloc1D<
    LatEleLocStructArray1D,
    allocate_lat_ele_loc_struct_container,
    deallocate_lat_ele_loc_struct_container,
    reallocate_lat_ele_loc_struct_container_data,
    access_lat_ele_loc_struct_container>;

class WakeStruct;

using WakeStructArray1D =
    FTypeArray1D<WakeStruct, allocate_fortran_wake_struct, deallocate_fortran_wake_struct>;
using WakeStructArray2D = FTypeArray2D<WakeStruct>;
using WakeStructArray3D = FTypeArray3D<WakeStruct>;

using WakeStructAlloc1D = FTypeAlloc1D<
    WakeStructArray1D,
    allocate_wake_struct_container,
    deallocate_wake_struct_container,
    reallocate_wake_struct_container_data,
    access_wake_struct_container>;

class TaylorTermStruct;

using TaylorTermStructArray1D = FTypeArray1D<
    TaylorTermStruct,
    allocate_fortran_taylor_term_struct,
    deallocate_fortran_taylor_term_struct>;
using TaylorTermStructArray2D = FTypeArray2D<TaylorTermStruct>;
using TaylorTermStructArray3D = FTypeArray3D<TaylorTermStruct>;

using TaylorTermStructAlloc1D = FTypeAlloc1D<
    TaylorTermStructArray1D,
    allocate_taylor_term_struct_container,
    deallocate_taylor_term_struct_container,
    reallocate_taylor_term_struct_container_data,
    access_taylor_term_struct_container>;

class TaylorStruct;

using TaylorStructArray1D =
    FTypeArray1D<TaylorStruct, allocate_fortran_taylor_struct, deallocate_fortran_taylor_struct>;
using TaylorStructArray2D = FTypeArray2D<TaylorStruct>;
using TaylorStructArray3D = FTypeArray3D<TaylorStruct>;

using TaylorStructAlloc1D = FTypeAlloc1D<
    TaylorStructArray1D,
    allocate_taylor_struct_container,
    deallocate_taylor_struct_container,
    reallocate_taylor_struct_container_data,
    access_taylor_struct_container>;

class EmTaylorTermStruct;

using EmTaylorTermStructArray1D = FTypeArray1D<
    EmTaylorTermStruct,
    allocate_fortran_em_taylor_term_struct,
    deallocate_fortran_em_taylor_term_struct>;
using EmTaylorTermStructArray2D = FTypeArray2D<EmTaylorTermStruct>;
using EmTaylorTermStructArray3D = FTypeArray3D<EmTaylorTermStruct>;

using EmTaylorTermStructAlloc1D = FTypeAlloc1D<
    EmTaylorTermStructArray1D,
    allocate_em_taylor_term_struct_container,
    deallocate_em_taylor_term_struct_container,
    reallocate_em_taylor_term_struct_container_data,
    access_em_taylor_term_struct_container>;

class EmTaylorStruct;

using EmTaylorStructArray1D = FTypeArray1D<
    EmTaylorStruct,
    allocate_fortran_em_taylor_struct,
    deallocate_fortran_em_taylor_struct>;
using EmTaylorStructArray2D = FTypeArray2D<EmTaylorStruct>;
using EmTaylorStructArray3D = FTypeArray3D<EmTaylorStruct>;

using EmTaylorStructAlloc1D = FTypeAlloc1D<
    EmTaylorStructArray1D,
    allocate_em_taylor_struct_container,
    deallocate_em_taylor_struct_container,
    reallocate_em_taylor_struct_container_data,
    access_em_taylor_struct_container>;

class CartesianMapTerm1Struct;

using CartesianMapTerm1StructArray1D = FTypeArray1D<
    CartesianMapTerm1Struct,
    allocate_fortran_cartesian_map_term1_struct,
    deallocate_fortran_cartesian_map_term1_struct>;
using CartesianMapTerm1StructArray2D = FTypeArray2D<CartesianMapTerm1Struct>;
using CartesianMapTerm1StructArray3D = FTypeArray3D<CartesianMapTerm1Struct>;

using CartesianMapTerm1StructAlloc1D = FTypeAlloc1D<
    CartesianMapTerm1StructArray1D,
    allocate_cartesian_map_term1_struct_container,
    deallocate_cartesian_map_term1_struct_container,
    reallocate_cartesian_map_term1_struct_container_data,
    access_cartesian_map_term1_struct_container>;

class CartesianMapTermStruct;

using CartesianMapTermStructArray1D = FTypeArray1D<
    CartesianMapTermStruct,
    allocate_fortran_cartesian_map_term_struct,
    deallocate_fortran_cartesian_map_term_struct>;
using CartesianMapTermStructArray2D = FTypeArray2D<CartesianMapTermStruct>;
using CartesianMapTermStructArray3D = FTypeArray3D<CartesianMapTermStruct>;

using CartesianMapTermStructAlloc1D = FTypeAlloc1D<
    CartesianMapTermStructArray1D,
    allocate_cartesian_map_term_struct_container,
    deallocate_cartesian_map_term_struct_container,
    reallocate_cartesian_map_term_struct_container_data,
    access_cartesian_map_term_struct_container>;

class CartesianMapStruct;

using CartesianMapStructArray1D = FTypeArray1D<
    CartesianMapStruct,
    allocate_fortran_cartesian_map_struct,
    deallocate_fortran_cartesian_map_struct>;
using CartesianMapStructArray2D = FTypeArray2D<CartesianMapStruct>;
using CartesianMapStructArray3D = FTypeArray3D<CartesianMapStruct>;

using CartesianMapStructAlloc1D = FTypeAlloc1D<
    CartesianMapStructArray1D,
    allocate_cartesian_map_struct_container,
    deallocate_cartesian_map_struct_container,
    reallocate_cartesian_map_struct_container_data,
    access_cartesian_map_struct_container>;

class CylindricalMapTerm1Struct;

using CylindricalMapTerm1StructArray1D = FTypeArray1D<
    CylindricalMapTerm1Struct,
    allocate_fortran_cylindrical_map_term1_struct,
    deallocate_fortran_cylindrical_map_term1_struct>;
using CylindricalMapTerm1StructArray2D = FTypeArray2D<CylindricalMapTerm1Struct>;
using CylindricalMapTerm1StructArray3D = FTypeArray3D<CylindricalMapTerm1Struct>;

using CylindricalMapTerm1StructAlloc1D = FTypeAlloc1D<
    CylindricalMapTerm1StructArray1D,
    allocate_cylindrical_map_term1_struct_container,
    deallocate_cylindrical_map_term1_struct_container,
    reallocate_cylindrical_map_term1_struct_container_data,
    access_cylindrical_map_term1_struct_container>;

class CylindricalMapTermStruct;

using CylindricalMapTermStructArray1D = FTypeArray1D<
    CylindricalMapTermStruct,
    allocate_fortran_cylindrical_map_term_struct,
    deallocate_fortran_cylindrical_map_term_struct>;
using CylindricalMapTermStructArray2D = FTypeArray2D<CylindricalMapTermStruct>;
using CylindricalMapTermStructArray3D = FTypeArray3D<CylindricalMapTermStruct>;

using CylindricalMapTermStructAlloc1D = FTypeAlloc1D<
    CylindricalMapTermStructArray1D,
    allocate_cylindrical_map_term_struct_container,
    deallocate_cylindrical_map_term_struct_container,
    reallocate_cylindrical_map_term_struct_container_data,
    access_cylindrical_map_term_struct_container>;

class CylindricalMapStruct;

using CylindricalMapStructArray1D = FTypeArray1D<
    CylindricalMapStruct,
    allocate_fortran_cylindrical_map_struct,
    deallocate_fortran_cylindrical_map_struct>;
using CylindricalMapStructArray2D = FTypeArray2D<CylindricalMapStruct>;
using CylindricalMapStructArray3D = FTypeArray3D<CylindricalMapStruct>;

using CylindricalMapStructAlloc1D = FTypeAlloc1D<
    CylindricalMapStructArray1D,
    allocate_cylindrical_map_struct_container,
    deallocate_cylindrical_map_struct_container,
    reallocate_cylindrical_map_struct_container_data,
    access_cylindrical_map_struct_container>;

class BicubicCmplxCoefStruct;

using BicubicCmplxCoefStructArray1D = FTypeArray1D<
    BicubicCmplxCoefStruct,
    allocate_fortran_bicubic_cmplx_coef_struct,
    deallocate_fortran_bicubic_cmplx_coef_struct>;
using BicubicCmplxCoefStructArray2D = FTypeArray2D<BicubicCmplxCoefStruct>;
using BicubicCmplxCoefStructArray3D = FTypeArray3D<BicubicCmplxCoefStruct>;

using BicubicCmplxCoefStructAlloc1D = FTypeAlloc1D<
    BicubicCmplxCoefStructArray1D,
    allocate_bicubic_cmplx_coef_struct_container,
    deallocate_bicubic_cmplx_coef_struct_container,
    reallocate_bicubic_cmplx_coef_struct_container_data,
    access_bicubic_cmplx_coef_struct_container>;

class TricubicCmplxCoefStruct;

using TricubicCmplxCoefStructArray1D = FTypeArray1D<
    TricubicCmplxCoefStruct,
    allocate_fortran_tricubic_cmplx_coef_struct,
    deallocate_fortran_tricubic_cmplx_coef_struct>;
using TricubicCmplxCoefStructArray2D = FTypeArray2D<TricubicCmplxCoefStruct>;
using TricubicCmplxCoefStructArray3D = FTypeArray3D<TricubicCmplxCoefStruct>;

using TricubicCmplxCoefStructAlloc1D = FTypeAlloc1D<
    TricubicCmplxCoefStructArray1D,
    allocate_tricubic_cmplx_coef_struct_container,
    deallocate_tricubic_cmplx_coef_struct_container,
    reallocate_tricubic_cmplx_coef_struct_container_data,
    access_tricubic_cmplx_coef_struct_container>;

class GridFieldPt1Struct;

using GridFieldPt1StructArray1D = FTypeArray1D<
    GridFieldPt1Struct,
    allocate_fortran_grid_field_pt1_struct,
    deallocate_fortran_grid_field_pt1_struct>;
using GridFieldPt1StructArray2D = FTypeArray2D<GridFieldPt1Struct>;
using GridFieldPt1StructArray3D = FTypeArray3D<GridFieldPt1Struct>;

using GridFieldPt1StructAlloc1D = FTypeAlloc1D<
    GridFieldPt1StructArray1D,
    allocate_grid_field_pt1_struct_container,
    deallocate_grid_field_pt1_struct_container,
    reallocate_grid_field_pt1_struct_container_data,
    access_grid_field_pt1_struct_container>;

class GridFieldPtStruct;

using GridFieldPtStructArray1D = FTypeArray1D<
    GridFieldPtStruct,
    allocate_fortran_grid_field_pt_struct,
    deallocate_fortran_grid_field_pt_struct>;
using GridFieldPtStructArray2D = FTypeArray2D<GridFieldPtStruct>;
using GridFieldPtStructArray3D = FTypeArray3D<GridFieldPtStruct>;

using GridFieldPtStructAlloc1D = FTypeAlloc1D<
    GridFieldPtStructArray1D,
    allocate_grid_field_pt_struct_container,
    deallocate_grid_field_pt_struct_container,
    reallocate_grid_field_pt_struct_container_data,
    access_grid_field_pt_struct_container>;

class GridFieldStruct;

using GridFieldStructArray1D = FTypeArray1D<
    GridFieldStruct,
    allocate_fortran_grid_field_struct,
    deallocate_fortran_grid_field_struct>;
using GridFieldStructArray2D = FTypeArray2D<GridFieldStruct>;
using GridFieldStructArray3D = FTypeArray3D<GridFieldStruct>;

using GridFieldStructAlloc1D = FTypeAlloc1D<
    GridFieldStructArray1D,
    allocate_grid_field_struct_container,
    deallocate_grid_field_struct_container,
    reallocate_grid_field_struct_container_data,
    access_grid_field_struct_container>;

class FloorPositionStruct;

using FloorPositionStructArray1D = FTypeArray1D<
    FloorPositionStruct,
    allocate_fortran_floor_position_struct,
    deallocate_fortran_floor_position_struct>;
using FloorPositionStructArray2D = FTypeArray2D<FloorPositionStruct>;
using FloorPositionStructArray3D = FTypeArray3D<FloorPositionStruct>;

using FloorPositionStructAlloc1D = FTypeAlloc1D<
    FloorPositionStructArray1D,
    allocate_floor_position_struct_container,
    deallocate_floor_position_struct_container,
    reallocate_floor_position_struct_container_data,
    access_floor_position_struct_container>;

class HighEnergySpaceChargeStruct;

using HighEnergySpaceChargeStructArray1D = FTypeArray1D<
    HighEnergySpaceChargeStruct,
    allocate_fortran_high_energy_space_charge_struct,
    deallocate_fortran_high_energy_space_charge_struct>;
using HighEnergySpaceChargeStructArray2D = FTypeArray2D<HighEnergySpaceChargeStruct>;
using HighEnergySpaceChargeStructArray3D = FTypeArray3D<HighEnergySpaceChargeStruct>;

using HighEnergySpaceChargeStructAlloc1D = FTypeAlloc1D<
    HighEnergySpaceChargeStructArray1D,
    allocate_high_energy_space_charge_struct_container,
    deallocate_high_energy_space_charge_struct_container,
    reallocate_high_energy_space_charge_struct_container_data,
    access_high_energy_space_charge_struct_container>;

class XyDispStruct;

using XyDispStructArray1D =
    FTypeArray1D<XyDispStruct, allocate_fortran_xy_disp_struct, deallocate_fortran_xy_disp_struct>;
using XyDispStructArray2D = FTypeArray2D<XyDispStruct>;
using XyDispStructArray3D = FTypeArray3D<XyDispStruct>;

using XyDispStructAlloc1D = FTypeAlloc1D<
    XyDispStructArray1D,
    allocate_xy_disp_struct_container,
    deallocate_xy_disp_struct_container,
    reallocate_xy_disp_struct_container_data,
    access_xy_disp_struct_container>;

class TwissStruct;

using TwissStructArray1D =
    FTypeArray1D<TwissStruct, allocate_fortran_twiss_struct, deallocate_fortran_twiss_struct>;
using TwissStructArray2D = FTypeArray2D<TwissStruct>;
using TwissStructArray3D = FTypeArray3D<TwissStruct>;

using TwissStructAlloc1D = FTypeAlloc1D<
    TwissStructArray1D,
    allocate_twiss_struct_container,
    deallocate_twiss_struct_container,
    reallocate_twiss_struct_container_data,
    access_twiss_struct_container>;

class Mode3Struct;

using Mode3StructArray1D =
    FTypeArray1D<Mode3Struct, allocate_fortran_mode3_struct, deallocate_fortran_mode3_struct>;
using Mode3StructArray2D = FTypeArray2D<Mode3Struct>;
using Mode3StructArray3D = FTypeArray3D<Mode3Struct>;

using Mode3StructAlloc1D = FTypeAlloc1D<
    Mode3StructArray1D,
    allocate_mode3_struct_container,
    deallocate_mode3_struct_container,
    reallocate_mode3_struct_container_data,
    access_mode3_struct_container>;

class BookkeepingStateStruct;

using BookkeepingStateStructArray1D = FTypeArray1D<
    BookkeepingStateStruct,
    allocate_fortran_bookkeeping_state_struct,
    deallocate_fortran_bookkeeping_state_struct>;
using BookkeepingStateStructArray2D = FTypeArray2D<BookkeepingStateStruct>;
using BookkeepingStateStructArray3D = FTypeArray3D<BookkeepingStateStruct>;

using BookkeepingStateStructAlloc1D = FTypeAlloc1D<
    BookkeepingStateStructArray1D,
    allocate_bookkeeping_state_struct_container,
    deallocate_bookkeeping_state_struct_container,
    reallocate_bookkeeping_state_struct_container_data,
    access_bookkeeping_state_struct_container>;

class RadMapStruct;

using RadMapStructArray1D =
    FTypeArray1D<RadMapStruct, allocate_fortran_rad_map_struct, deallocate_fortran_rad_map_struct>;
using RadMapStructArray2D = FTypeArray2D<RadMapStruct>;
using RadMapStructArray3D = FTypeArray3D<RadMapStruct>;

using RadMapStructAlloc1D = FTypeAlloc1D<
    RadMapStructArray1D,
    allocate_rad_map_struct_container,
    deallocate_rad_map_struct_container,
    reallocate_rad_map_struct_container_data,
    access_rad_map_struct_container>;

class RadMapEleStruct;

using RadMapEleStructArray1D = FTypeArray1D<
    RadMapEleStruct,
    allocate_fortran_rad_map_ele_struct,
    deallocate_fortran_rad_map_ele_struct>;
using RadMapEleStructArray2D = FTypeArray2D<RadMapEleStruct>;
using RadMapEleStructArray3D = FTypeArray3D<RadMapEleStruct>;

using RadMapEleStructAlloc1D = FTypeAlloc1D<
    RadMapEleStructArray1D,
    allocate_rad_map_ele_struct_container,
    deallocate_rad_map_ele_struct_container,
    reallocate_rad_map_ele_struct_container_data,
    access_rad_map_ele_struct_container>;

class GenGrad1Struct;

using GenGrad1StructArray1D = FTypeArray1D<
    GenGrad1Struct,
    allocate_fortran_gen_grad1_struct,
    deallocate_fortran_gen_grad1_struct>;
using GenGrad1StructArray2D = FTypeArray2D<GenGrad1Struct>;
using GenGrad1StructArray3D = FTypeArray3D<GenGrad1Struct>;

using GenGrad1StructAlloc1D = FTypeAlloc1D<
    GenGrad1StructArray1D,
    allocate_gen_grad1_struct_container,
    deallocate_gen_grad1_struct_container,
    reallocate_gen_grad1_struct_container_data,
    access_gen_grad1_struct_container>;

class GenGradMapStruct;

using GenGradMapStructArray1D = FTypeArray1D<
    GenGradMapStruct,
    allocate_fortran_gen_grad_map_struct,
    deallocate_fortran_gen_grad_map_struct>;
using GenGradMapStructArray2D = FTypeArray2D<GenGradMapStruct>;
using GenGradMapStructArray3D = FTypeArray3D<GenGradMapStruct>;

using GenGradMapStructAlloc1D = FTypeAlloc1D<
    GenGradMapStructArray1D,
    allocate_gen_grad_map_struct_container,
    deallocate_gen_grad_map_struct_container,
    reallocate_gen_grad_map_struct_container_data,
    access_gen_grad_map_struct_container>;

class SurfaceSegmentedPtStruct;

using SurfaceSegmentedPtStructArray1D = FTypeArray1D<
    SurfaceSegmentedPtStruct,
    allocate_fortran_surface_segmented_pt_struct,
    deallocate_fortran_surface_segmented_pt_struct>;
using SurfaceSegmentedPtStructArray2D = FTypeArray2D<SurfaceSegmentedPtStruct>;
using SurfaceSegmentedPtStructArray3D = FTypeArray3D<SurfaceSegmentedPtStruct>;

using SurfaceSegmentedPtStructAlloc1D = FTypeAlloc1D<
    SurfaceSegmentedPtStructArray1D,
    allocate_surface_segmented_pt_struct_container,
    deallocate_surface_segmented_pt_struct_container,
    reallocate_surface_segmented_pt_struct_container_data,
    access_surface_segmented_pt_struct_container>;

class SurfaceSegmentedStruct;

using SurfaceSegmentedStructArray1D = FTypeArray1D<
    SurfaceSegmentedStruct,
    allocate_fortran_surface_segmented_struct,
    deallocate_fortran_surface_segmented_struct>;
using SurfaceSegmentedStructArray2D = FTypeArray2D<SurfaceSegmentedStruct>;
using SurfaceSegmentedStructArray3D = FTypeArray3D<SurfaceSegmentedStruct>;

using SurfaceSegmentedStructAlloc1D = FTypeAlloc1D<
    SurfaceSegmentedStructArray1D,
    allocate_surface_segmented_struct_container,
    deallocate_surface_segmented_struct_container,
    reallocate_surface_segmented_struct_container_data,
    access_surface_segmented_struct_container>;

class SurfaceHMisalignPtStruct;

using SurfaceHMisalignPtStructArray1D = FTypeArray1D<
    SurfaceHMisalignPtStruct,
    allocate_fortran_surface_h_misalign_pt_struct,
    deallocate_fortran_surface_h_misalign_pt_struct>;
using SurfaceHMisalignPtStructArray2D = FTypeArray2D<SurfaceHMisalignPtStruct>;
using SurfaceHMisalignPtStructArray3D = FTypeArray3D<SurfaceHMisalignPtStruct>;

using SurfaceHMisalignPtStructAlloc1D = FTypeAlloc1D<
    SurfaceHMisalignPtStructArray1D,
    allocate_surface_h_misalign_pt_struct_container,
    deallocate_surface_h_misalign_pt_struct_container,
    reallocate_surface_h_misalign_pt_struct_container_data,
    access_surface_h_misalign_pt_struct_container>;

class SurfaceHMisalignStruct;

using SurfaceHMisalignStructArray1D = FTypeArray1D<
    SurfaceHMisalignStruct,
    allocate_fortran_surface_h_misalign_struct,
    deallocate_fortran_surface_h_misalign_struct>;
using SurfaceHMisalignStructArray2D = FTypeArray2D<SurfaceHMisalignStruct>;
using SurfaceHMisalignStructArray3D = FTypeArray3D<SurfaceHMisalignStruct>;

using SurfaceHMisalignStructAlloc1D = FTypeAlloc1D<
    SurfaceHMisalignStructArray1D,
    allocate_surface_h_misalign_struct_container,
    deallocate_surface_h_misalign_struct_container,
    reallocate_surface_h_misalign_struct_container_data,
    access_surface_h_misalign_struct_container>;

class SurfaceDisplacementPtStruct;

using SurfaceDisplacementPtStructArray1D = FTypeArray1D<
    SurfaceDisplacementPtStruct,
    allocate_fortran_surface_displacement_pt_struct,
    deallocate_fortran_surface_displacement_pt_struct>;
using SurfaceDisplacementPtStructArray2D = FTypeArray2D<SurfaceDisplacementPtStruct>;
using SurfaceDisplacementPtStructArray3D = FTypeArray3D<SurfaceDisplacementPtStruct>;

using SurfaceDisplacementPtStructAlloc1D = FTypeAlloc1D<
    SurfaceDisplacementPtStructArray1D,
    allocate_surface_displacement_pt_struct_container,
    deallocate_surface_displacement_pt_struct_container,
    reallocate_surface_displacement_pt_struct_container_data,
    access_surface_displacement_pt_struct_container>;

class SurfaceDisplacementStruct;

using SurfaceDisplacementStructArray1D = FTypeArray1D<
    SurfaceDisplacementStruct,
    allocate_fortran_surface_displacement_struct,
    deallocate_fortran_surface_displacement_struct>;
using SurfaceDisplacementStructArray2D = FTypeArray2D<SurfaceDisplacementStruct>;
using SurfaceDisplacementStructArray3D = FTypeArray3D<SurfaceDisplacementStruct>;

using SurfaceDisplacementStructAlloc1D = FTypeAlloc1D<
    SurfaceDisplacementStructArray1D,
    allocate_surface_displacement_struct_container,
    deallocate_surface_displacement_struct_container,
    reallocate_surface_displacement_struct_container_data,
    access_surface_displacement_struct_container>;

class TargetPointStruct;

using TargetPointStructArray1D = FTypeArray1D<
    TargetPointStruct,
    allocate_fortran_target_point_struct,
    deallocate_fortran_target_point_struct>;
using TargetPointStructArray2D = FTypeArray2D<TargetPointStruct>;
using TargetPointStructArray3D = FTypeArray3D<TargetPointStruct>;

using TargetPointStructAlloc1D = FTypeAlloc1D<
    TargetPointStructArray1D,
    allocate_target_point_struct_container,
    deallocate_target_point_struct_container,
    reallocate_target_point_struct_container_data,
    access_target_point_struct_container>;

class SurfaceCurvatureStruct;

using SurfaceCurvatureStructArray1D = FTypeArray1D<
    SurfaceCurvatureStruct,
    allocate_fortran_surface_curvature_struct,
    deallocate_fortran_surface_curvature_struct>;
using SurfaceCurvatureStructArray2D = FTypeArray2D<SurfaceCurvatureStruct>;
using SurfaceCurvatureStructArray3D = FTypeArray3D<SurfaceCurvatureStruct>;

using SurfaceCurvatureStructAlloc1D = FTypeAlloc1D<
    SurfaceCurvatureStructArray1D,
    allocate_surface_curvature_struct_container,
    deallocate_surface_curvature_struct_container,
    reallocate_surface_curvature_struct_container_data,
    access_surface_curvature_struct_container>;

class PhotonTargetStruct;

using PhotonTargetStructArray1D = FTypeArray1D<
    PhotonTargetStruct,
    allocate_fortran_photon_target_struct,
    deallocate_fortran_photon_target_struct>;
using PhotonTargetStructArray2D = FTypeArray2D<PhotonTargetStruct>;
using PhotonTargetStructArray3D = FTypeArray3D<PhotonTargetStruct>;

using PhotonTargetStructAlloc1D = FTypeAlloc1D<
    PhotonTargetStructArray1D,
    allocate_photon_target_struct_container,
    deallocate_photon_target_struct_container,
    reallocate_photon_target_struct_container_data,
    access_photon_target_struct_container>;

class PhotonMaterialStruct;

using PhotonMaterialStructArray1D = FTypeArray1D<
    PhotonMaterialStruct,
    allocate_fortran_photon_material_struct,
    deallocate_fortran_photon_material_struct>;
using PhotonMaterialStructArray2D = FTypeArray2D<PhotonMaterialStruct>;
using PhotonMaterialStructArray3D = FTypeArray3D<PhotonMaterialStruct>;

using PhotonMaterialStructAlloc1D = FTypeAlloc1D<
    PhotonMaterialStructArray1D,
    allocate_photon_material_struct_container,
    deallocate_photon_material_struct_container,
    reallocate_photon_material_struct_container_data,
    access_photon_material_struct_container>;

class PixelPtStruct;

using PixelPtStructArray1D = FTypeArray1D<
    PixelPtStruct,
    allocate_fortran_pixel_pt_struct,
    deallocate_fortran_pixel_pt_struct>;
using PixelPtStructArray2D = FTypeArray2D<PixelPtStruct>;
using PixelPtStructArray3D = FTypeArray3D<PixelPtStruct>;

using PixelPtStructAlloc1D = FTypeAlloc1D<
    PixelPtStructArray1D,
    allocate_pixel_pt_struct_container,
    deallocate_pixel_pt_struct_container,
    reallocate_pixel_pt_struct_container_data,
    access_pixel_pt_struct_container>;

class PixelDetecStruct;

using PixelDetecStructArray1D = FTypeArray1D<
    PixelDetecStruct,
    allocate_fortran_pixel_detec_struct,
    deallocate_fortran_pixel_detec_struct>;
using PixelDetecStructArray2D = FTypeArray2D<PixelDetecStruct>;
using PixelDetecStructArray3D = FTypeArray3D<PixelDetecStruct>;

using PixelDetecStructAlloc1D = FTypeAlloc1D<
    PixelDetecStructArray1D,
    allocate_pixel_detec_struct_container,
    deallocate_pixel_detec_struct_container,
    reallocate_pixel_detec_struct_container_data,
    access_pixel_detec_struct_container>;

class PhotonElementStruct;

using PhotonElementStructArray1D = FTypeArray1D<
    PhotonElementStruct,
    allocate_fortran_photon_element_struct,
    deallocate_fortran_photon_element_struct>;
using PhotonElementStructArray2D = FTypeArray2D<PhotonElementStruct>;
using PhotonElementStructArray3D = FTypeArray3D<PhotonElementStruct>;

using PhotonElementStructAlloc1D = FTypeAlloc1D<
    PhotonElementStructArray1D,
    allocate_photon_element_struct_container,
    deallocate_photon_element_struct_container,
    reallocate_photon_element_struct_container_data,
    access_photon_element_struct_container>;

class Wall3dVertexStruct;

using Wall3dVertexStructArray1D = FTypeArray1D<
    Wall3dVertexStruct,
    allocate_fortran_wall3d_vertex_struct,
    deallocate_fortran_wall3d_vertex_struct>;
using Wall3dVertexStructArray2D = FTypeArray2D<Wall3dVertexStruct>;
using Wall3dVertexStructArray3D = FTypeArray3D<Wall3dVertexStruct>;

using Wall3dVertexStructAlloc1D = FTypeAlloc1D<
    Wall3dVertexStructArray1D,
    allocate_wall3d_vertex_struct_container,
    deallocate_wall3d_vertex_struct_container,
    reallocate_wall3d_vertex_struct_container_data,
    access_wall3d_vertex_struct_container>;

class Wall3dSectionStruct;

using Wall3dSectionStructArray1D = FTypeArray1D<
    Wall3dSectionStruct,
    allocate_fortran_wall3d_section_struct,
    deallocate_fortran_wall3d_section_struct>;
using Wall3dSectionStructArray2D = FTypeArray2D<Wall3dSectionStruct>;
using Wall3dSectionStructArray3D = FTypeArray3D<Wall3dSectionStruct>;

using Wall3dSectionStructAlloc1D = FTypeAlloc1D<
    Wall3dSectionStructArray1D,
    allocate_wall3d_section_struct_container,
    deallocate_wall3d_section_struct_container,
    reallocate_wall3d_section_struct_container_data,
    access_wall3d_section_struct_container>;

class Wall3dStruct;

using Wall3dStructArray1D =
    FTypeArray1D<Wall3dStruct, allocate_fortran_wall3d_struct, deallocate_fortran_wall3d_struct>;
using Wall3dStructArray2D = FTypeArray2D<Wall3dStruct>;
using Wall3dStructArray3D = FTypeArray3D<Wall3dStruct>;

using Wall3dStructAlloc1D = FTypeAlloc1D<
    Wall3dStructArray1D,
    allocate_wall3d_struct_container,
    deallocate_wall3d_struct_container,
    reallocate_wall3d_struct_container_data,
    access_wall3d_struct_container>;

class RamperLordStruct;

using RamperLordStructArray1D = FTypeArray1D<
    RamperLordStruct,
    allocate_fortran_ramper_lord_struct,
    deallocate_fortran_ramper_lord_struct>;
using RamperLordStructArray2D = FTypeArray2D<RamperLordStruct>;
using RamperLordStructArray3D = FTypeArray3D<RamperLordStruct>;

using RamperLordStructAlloc1D = FTypeAlloc1D<
    RamperLordStructArray1D,
    allocate_ramper_lord_struct_container,
    deallocate_ramper_lord_struct_container,
    reallocate_ramper_lord_struct_container_data,
    access_ramper_lord_struct_container>;

class ControlStruct;

using ControlStructArray1D =
    FTypeArray1D<ControlStruct, allocate_fortran_control_struct, deallocate_fortran_control_struct>;
using ControlStructArray2D = FTypeArray2D<ControlStruct>;
using ControlStructArray3D = FTypeArray3D<ControlStruct>;

using ControlStructAlloc1D = FTypeAlloc1D<
    ControlStructArray1D,
    allocate_control_struct_container,
    deallocate_control_struct_container,
    reallocate_control_struct_container_data,
    access_control_struct_container>;

class ControlVar1Struct;

using ControlVar1StructArray1D = FTypeArray1D<
    ControlVar1Struct,
    allocate_fortran_control_var1_struct,
    deallocate_fortran_control_var1_struct>;
using ControlVar1StructArray2D = FTypeArray2D<ControlVar1Struct>;
using ControlVar1StructArray3D = FTypeArray3D<ControlVar1Struct>;

using ControlVar1StructAlloc1D = FTypeAlloc1D<
    ControlVar1StructArray1D,
    allocate_control_var1_struct_container,
    deallocate_control_var1_struct_container,
    reallocate_control_var1_struct_container_data,
    access_control_var1_struct_container>;

class ControlRamp1Struct;

using ControlRamp1StructArray1D = FTypeArray1D<
    ControlRamp1Struct,
    allocate_fortran_control_ramp1_struct,
    deallocate_fortran_control_ramp1_struct>;
using ControlRamp1StructArray2D = FTypeArray2D<ControlRamp1Struct>;
using ControlRamp1StructArray3D = FTypeArray3D<ControlRamp1Struct>;

using ControlRamp1StructAlloc1D = FTypeAlloc1D<
    ControlRamp1StructArray1D,
    allocate_control_ramp1_struct_container,
    deallocate_control_ramp1_struct_container,
    reallocate_control_ramp1_struct_container_data,
    access_control_ramp1_struct_container>;

class ControllerStruct;

using ControllerStructArray1D = FTypeArray1D<
    ControllerStruct,
    allocate_fortran_controller_struct,
    deallocate_fortran_controller_struct>;
using ControllerStructArray2D = FTypeArray2D<ControllerStruct>;
using ControllerStructArray3D = FTypeArray3D<ControllerStruct>;

using ControllerStructAlloc1D = FTypeAlloc1D<
    ControllerStructArray1D,
    allocate_controller_struct_container,
    deallocate_controller_struct_container,
    reallocate_controller_struct_container_data,
    access_controller_struct_container>;

class EllipseBeamInitStruct;

using EllipseBeamInitStructArray1D = FTypeArray1D<
    EllipseBeamInitStruct,
    allocate_fortran_ellipse_beam_init_struct,
    deallocate_fortran_ellipse_beam_init_struct>;
using EllipseBeamInitStructArray2D = FTypeArray2D<EllipseBeamInitStruct>;
using EllipseBeamInitStructArray3D = FTypeArray3D<EllipseBeamInitStruct>;

using EllipseBeamInitStructAlloc1D = FTypeAlloc1D<
    EllipseBeamInitStructArray1D,
    allocate_ellipse_beam_init_struct_container,
    deallocate_ellipse_beam_init_struct_container,
    reallocate_ellipse_beam_init_struct_container_data,
    access_ellipse_beam_init_struct_container>;

class KvBeamInitStruct;

using KvBeamInitStructArray1D = FTypeArray1D<
    KvBeamInitStruct,
    allocate_fortran_kv_beam_init_struct,
    deallocate_fortran_kv_beam_init_struct>;
using KvBeamInitStructArray2D = FTypeArray2D<KvBeamInitStruct>;
using KvBeamInitStructArray3D = FTypeArray3D<KvBeamInitStruct>;

using KvBeamInitStructAlloc1D = FTypeAlloc1D<
    KvBeamInitStructArray1D,
    allocate_kv_beam_init_struct_container,
    deallocate_kv_beam_init_struct_container,
    reallocate_kv_beam_init_struct_container_data,
    access_kv_beam_init_struct_container>;

class GridBeamInitStruct;

using GridBeamInitStructArray1D = FTypeArray1D<
    GridBeamInitStruct,
    allocate_fortran_grid_beam_init_struct,
    deallocate_fortran_grid_beam_init_struct>;
using GridBeamInitStructArray2D = FTypeArray2D<GridBeamInitStruct>;
using GridBeamInitStructArray3D = FTypeArray3D<GridBeamInitStruct>;

using GridBeamInitStructAlloc1D = FTypeAlloc1D<
    GridBeamInitStructArray1D,
    allocate_grid_beam_init_struct_container,
    deallocate_grid_beam_init_struct_container,
    reallocate_grid_beam_init_struct_container_data,
    access_grid_beam_init_struct_container>;

class BeamInitStruct;

using BeamInitStructArray1D = FTypeArray1D<
    BeamInitStruct,
    allocate_fortran_beam_init_struct,
    deallocate_fortran_beam_init_struct>;
using BeamInitStructArray2D = FTypeArray2D<BeamInitStruct>;
using BeamInitStructArray3D = FTypeArray3D<BeamInitStruct>;

using BeamInitStructAlloc1D = FTypeAlloc1D<
    BeamInitStructArray1D,
    allocate_beam_init_struct_container,
    deallocate_beam_init_struct_container,
    reallocate_beam_init_struct_container_data,
    access_beam_init_struct_container>;

class LatParamStruct;

using LatParamStructArray1D = FTypeArray1D<
    LatParamStruct,
    allocate_fortran_lat_param_struct,
    deallocate_fortran_lat_param_struct>;
using LatParamStructArray2D = FTypeArray2D<LatParamStruct>;
using LatParamStructArray3D = FTypeArray3D<LatParamStruct>;

using LatParamStructAlloc1D = FTypeAlloc1D<
    LatParamStructArray1D,
    allocate_lat_param_struct_container,
    deallocate_lat_param_struct_container,
    reallocate_lat_param_struct_container_data,
    access_lat_param_struct_container>;

class ModeInfoStruct;

using ModeInfoStructArray1D = FTypeArray1D<
    ModeInfoStruct,
    allocate_fortran_mode_info_struct,
    deallocate_fortran_mode_info_struct>;
using ModeInfoStructArray2D = FTypeArray2D<ModeInfoStruct>;
using ModeInfoStructArray3D = FTypeArray3D<ModeInfoStruct>;

using ModeInfoStructAlloc1D = FTypeAlloc1D<
    ModeInfoStructArray1D,
    allocate_mode_info_struct_container,
    deallocate_mode_info_struct_container,
    reallocate_mode_info_struct_container_data,
    access_mode_info_struct_container>;

class PreTrackerStruct;

using PreTrackerStructArray1D = FTypeArray1D<
    PreTrackerStruct,
    allocate_fortran_pre_tracker_struct,
    deallocate_fortran_pre_tracker_struct>;
using PreTrackerStructArray2D = FTypeArray2D<PreTrackerStruct>;
using PreTrackerStructArray3D = FTypeArray3D<PreTrackerStruct>;

using PreTrackerStructAlloc1D = FTypeAlloc1D<
    PreTrackerStructArray1D,
    allocate_pre_tracker_struct_container,
    deallocate_pre_tracker_struct_container,
    reallocate_pre_tracker_struct_container_data,
    access_pre_tracker_struct_container>;

class AnormalModeStruct;

using AnormalModeStructArray1D = FTypeArray1D<
    AnormalModeStruct,
    allocate_fortran_anormal_mode_struct,
    deallocate_fortran_anormal_mode_struct>;
using AnormalModeStructArray2D = FTypeArray2D<AnormalModeStruct>;
using AnormalModeStructArray3D = FTypeArray3D<AnormalModeStruct>;

using AnormalModeStructAlloc1D = FTypeAlloc1D<
    AnormalModeStructArray1D,
    allocate_anormal_mode_struct_container,
    deallocate_anormal_mode_struct_container,
    reallocate_anormal_mode_struct_container_data,
    access_anormal_mode_struct_container>;

class LinacNormalModeStruct;

using LinacNormalModeStructArray1D = FTypeArray1D<
    LinacNormalModeStruct,
    allocate_fortran_linac_normal_mode_struct,
    deallocate_fortran_linac_normal_mode_struct>;
using LinacNormalModeStructArray2D = FTypeArray2D<LinacNormalModeStruct>;
using LinacNormalModeStructArray3D = FTypeArray3D<LinacNormalModeStruct>;

using LinacNormalModeStructAlloc1D = FTypeAlloc1D<
    LinacNormalModeStructArray1D,
    allocate_linac_normal_mode_struct_container,
    deallocate_linac_normal_mode_struct_container,
    reallocate_linac_normal_mode_struct_container_data,
    access_linac_normal_mode_struct_container>;

class NormalModesStruct;

using NormalModesStructArray1D = FTypeArray1D<
    NormalModesStruct,
    allocate_fortran_normal_modes_struct,
    deallocate_fortran_normal_modes_struct>;
using NormalModesStructArray2D = FTypeArray2D<NormalModesStruct>;
using NormalModesStructArray3D = FTypeArray3D<NormalModesStruct>;

using NormalModesStructAlloc1D = FTypeAlloc1D<
    NormalModesStructArray1D,
    allocate_normal_modes_struct_container,
    deallocate_normal_modes_struct_container,
    reallocate_normal_modes_struct_container_data,
    access_normal_modes_struct_container>;

class EmFieldStruct;

using EmFieldStructArray1D = FTypeArray1D<
    EmFieldStruct,
    allocate_fortran_em_field_struct,
    deallocate_fortran_em_field_struct>;
using EmFieldStructArray2D = FTypeArray2D<EmFieldStruct>;
using EmFieldStructArray3D = FTypeArray3D<EmFieldStruct>;

using EmFieldStructAlloc1D = FTypeAlloc1D<
    EmFieldStructArray1D,
    allocate_em_field_struct_container,
    deallocate_em_field_struct_container,
    reallocate_em_field_struct_container_data,
    access_em_field_struct_container>;

class StrongBeamStruct;

using StrongBeamStructArray1D = FTypeArray1D<
    StrongBeamStruct,
    allocate_fortran_strong_beam_struct,
    deallocate_fortran_strong_beam_struct>;
using StrongBeamStructArray2D = FTypeArray2D<StrongBeamStruct>;
using StrongBeamStructArray3D = FTypeArray3D<StrongBeamStruct>;

using StrongBeamStructAlloc1D = FTypeAlloc1D<
    StrongBeamStructArray1D,
    allocate_strong_beam_struct_container,
    deallocate_strong_beam_struct_container,
    reallocate_strong_beam_struct_container_data,
    access_strong_beam_struct_container>;

class TrackPointStruct;

using TrackPointStructArray1D = FTypeArray1D<
    TrackPointStruct,
    allocate_fortran_track_point_struct,
    deallocate_fortran_track_point_struct>;
using TrackPointStructArray2D = FTypeArray2D<TrackPointStruct>;
using TrackPointStructArray3D = FTypeArray3D<TrackPointStruct>;

using TrackPointStructAlloc1D = FTypeAlloc1D<
    TrackPointStructArray1D,
    allocate_track_point_struct_container,
    deallocate_track_point_struct_container,
    reallocate_track_point_struct_container_data,
    access_track_point_struct_container>;

class TrackStruct;

using TrackStructArray1D =
    FTypeArray1D<TrackStruct, allocate_fortran_track_struct, deallocate_fortran_track_struct>;
using TrackStructArray2D = FTypeArray2D<TrackStruct>;
using TrackStructArray3D = FTypeArray3D<TrackStruct>;

using TrackStructAlloc1D = FTypeAlloc1D<
    TrackStructArray1D,
    allocate_track_struct_container,
    deallocate_track_struct_container,
    reallocate_track_struct_container_data,
    access_track_struct_container>;

class SpaceChargeCommonStruct;

using SpaceChargeCommonStructArray1D = FTypeArray1D<
    SpaceChargeCommonStruct,
    allocate_fortran_space_charge_common_struct,
    deallocate_fortran_space_charge_common_struct>;
using SpaceChargeCommonStructArray2D = FTypeArray2D<SpaceChargeCommonStruct>;
using SpaceChargeCommonStructArray3D = FTypeArray3D<SpaceChargeCommonStruct>;

using SpaceChargeCommonStructAlloc1D = FTypeAlloc1D<
    SpaceChargeCommonStructArray1D,
    allocate_space_charge_common_struct_container,
    deallocate_space_charge_common_struct_container,
    reallocate_space_charge_common_struct_container_data,
    access_space_charge_common_struct_container>;

class BmadCommonStruct;

using BmadCommonStructArray1D = FTypeArray1D<
    BmadCommonStruct,
    allocate_fortran_bmad_common_struct,
    deallocate_fortran_bmad_common_struct>;
using BmadCommonStructArray2D = FTypeArray2D<BmadCommonStruct>;
using BmadCommonStructArray3D = FTypeArray3D<BmadCommonStruct>;

using BmadCommonStructAlloc1D = FTypeAlloc1D<
    BmadCommonStructArray1D,
    allocate_bmad_common_struct_container,
    deallocate_bmad_common_struct_container,
    reallocate_bmad_common_struct_container_data,
    access_bmad_common_struct_container>;

class RadInt1Struct;

using RadInt1StructArray1D = FTypeArray1D<
    RadInt1Struct,
    allocate_fortran_rad_int1_struct,
    deallocate_fortran_rad_int1_struct>;
using RadInt1StructArray2D = FTypeArray2D<RadInt1Struct>;
using RadInt1StructArray3D = FTypeArray3D<RadInt1Struct>;

using RadInt1StructAlloc1D = FTypeAlloc1D<
    RadInt1StructArray1D,
    allocate_rad_int1_struct_container,
    deallocate_rad_int1_struct_container,
    reallocate_rad_int1_struct_container_data,
    access_rad_int1_struct_container>;

class RadIntBranchStruct;

using RadIntBranchStructArray1D = FTypeArray1D<
    RadIntBranchStruct,
    allocate_fortran_rad_int_branch_struct,
    deallocate_fortran_rad_int_branch_struct>;
using RadIntBranchStructArray2D = FTypeArray2D<RadIntBranchStruct>;
using RadIntBranchStructArray3D = FTypeArray3D<RadIntBranchStruct>;

using RadIntBranchStructAlloc1D = FTypeAlloc1D<
    RadIntBranchStructArray1D,
    allocate_rad_int_branch_struct_container,
    deallocate_rad_int_branch_struct_container,
    reallocate_rad_int_branch_struct_container_data,
    access_rad_int_branch_struct_container>;

class RadIntAllEleStruct;

using RadIntAllEleStructArray1D = FTypeArray1D<
    RadIntAllEleStruct,
    allocate_fortran_rad_int_all_ele_struct,
    deallocate_fortran_rad_int_all_ele_struct>;
using RadIntAllEleStructArray2D = FTypeArray2D<RadIntAllEleStruct>;
using RadIntAllEleStructArray3D = FTypeArray3D<RadIntAllEleStruct>;

using RadIntAllEleStructAlloc1D = FTypeAlloc1D<
    RadIntAllEleStructArray1D,
    allocate_rad_int_all_ele_struct_container,
    deallocate_rad_int_all_ele_struct_container,
    reallocate_rad_int_all_ele_struct_container_data,
    access_rad_int_all_ele_struct_container>;

class RfStairStepStruct;

using RfStairStepStructArray1D = FTypeArray1D<
    RfStairStepStruct,
    allocate_fortran_rf_stair_step_struct,
    deallocate_fortran_rf_stair_step_struct>;
using RfStairStepStructArray2D = FTypeArray2D<RfStairStepStruct>;
using RfStairStepStructArray3D = FTypeArray3D<RfStairStepStruct>;

using RfStairStepStructAlloc1D = FTypeAlloc1D<
    RfStairStepStructArray1D,
    allocate_rf_stair_step_struct_container,
    deallocate_rf_stair_step_struct_container,
    reallocate_rf_stair_step_struct_container_data,
    access_rf_stair_step_struct_container>;

class RfEleStruct;

using RfEleStructArray1D =
    FTypeArray1D<RfEleStruct, allocate_fortran_rf_ele_struct, deallocate_fortran_rf_ele_struct>;
using RfEleStructArray2D = FTypeArray2D<RfEleStruct>;
using RfEleStructArray3D = FTypeArray3D<RfEleStruct>;

using RfEleStructAlloc1D = FTypeAlloc1D<
    RfEleStructArray1D,
    allocate_rf_ele_struct_container,
    deallocate_rf_ele_struct_container,
    reallocate_rf_ele_struct_container_data,
    access_rf_ele_struct_container>;

class EleStruct;

using EleStructArray1D =
    FTypeArray1D<EleStruct, allocate_fortran_ele_struct, deallocate_fortran_ele_struct>;
using EleStructArray2D = FTypeArray2D<EleStruct>;
using EleStructArray3D = FTypeArray3D<EleStruct>;

using EleStructAlloc1D = FTypeAlloc1D<
    EleStructArray1D,
    allocate_ele_struct_container,
    deallocate_ele_struct_container,
    reallocate_ele_struct_container_data,
    access_ele_struct_container>;

class ComplexTaylorTermStruct;

using ComplexTaylorTermStructArray1D = FTypeArray1D<
    ComplexTaylorTermStruct,
    allocate_fortran_complex_taylor_term_struct,
    deallocate_fortran_complex_taylor_term_struct>;
using ComplexTaylorTermStructArray2D = FTypeArray2D<ComplexTaylorTermStruct>;
using ComplexTaylorTermStructArray3D = FTypeArray3D<ComplexTaylorTermStruct>;

using ComplexTaylorTermStructAlloc1D = FTypeAlloc1D<
    ComplexTaylorTermStructArray1D,
    allocate_complex_taylor_term_struct_container,
    deallocate_complex_taylor_term_struct_container,
    reallocate_complex_taylor_term_struct_container_data,
    access_complex_taylor_term_struct_container>;

class ComplexTaylorStruct;

using ComplexTaylorStructArray1D = FTypeArray1D<
    ComplexTaylorStruct,
    allocate_fortran_complex_taylor_struct,
    deallocate_fortran_complex_taylor_struct>;
using ComplexTaylorStructArray2D = FTypeArray2D<ComplexTaylorStruct>;
using ComplexTaylorStructArray3D = FTypeArray3D<ComplexTaylorStruct>;

using ComplexTaylorStructAlloc1D = FTypeAlloc1D<
    ComplexTaylorStructArray1D,
    allocate_complex_taylor_struct_container,
    deallocate_complex_taylor_struct_container,
    reallocate_complex_taylor_struct_container_data,
    access_complex_taylor_struct_container>;

class BranchStruct;

using BranchStructArray1D =
    FTypeArray1D<BranchStruct, allocate_fortran_branch_struct, deallocate_fortran_branch_struct>;
using BranchStructArray2D = FTypeArray2D<BranchStruct>;
using BranchStructArray3D = FTypeArray3D<BranchStruct>;

using BranchStructAlloc1D = FTypeAlloc1D<
    BranchStructArray1D,
    allocate_branch_struct_container,
    deallocate_branch_struct_container,
    reallocate_branch_struct_container_data,
    access_branch_struct_container>;

class LatStruct;

using LatStructArray1D =
    FTypeArray1D<LatStruct, allocate_fortran_lat_struct, deallocate_fortran_lat_struct>;
using LatStructArray2D = FTypeArray2D<LatStruct>;
using LatStructArray3D = FTypeArray3D<LatStruct>;

using LatStructAlloc1D = FTypeAlloc1D<
    LatStructArray1D,
    allocate_lat_struct_container,
    deallocate_lat_struct_container,
    reallocate_lat_struct_container_data,
    access_lat_struct_container>;

class BunchStruct;

using BunchStructArray1D =
    FTypeArray1D<BunchStruct, allocate_fortran_bunch_struct, deallocate_fortran_bunch_struct>;
using BunchStructArray2D = FTypeArray2D<BunchStruct>;
using BunchStructArray3D = FTypeArray3D<BunchStruct>;

using BunchStructAlloc1D = FTypeAlloc1D<
    BunchStructArray1D,
    allocate_bunch_struct_container,
    deallocate_bunch_struct_container,
    reallocate_bunch_struct_container_data,
    access_bunch_struct_container>;

class BunchParamsStruct;

using BunchParamsStructArray1D = FTypeArray1D<
    BunchParamsStruct,
    allocate_fortran_bunch_params_struct,
    deallocate_fortran_bunch_params_struct>;
using BunchParamsStructArray2D = FTypeArray2D<BunchParamsStruct>;
using BunchParamsStructArray3D = FTypeArray3D<BunchParamsStruct>;

using BunchParamsStructAlloc1D = FTypeAlloc1D<
    BunchParamsStructArray1D,
    allocate_bunch_params_struct_container,
    deallocate_bunch_params_struct_container,
    reallocate_bunch_params_struct_container_data,
    access_bunch_params_struct_container>;

class BeamStruct;

using BeamStructArray1D =
    FTypeArray1D<BeamStruct, allocate_fortran_beam_struct, deallocate_fortran_beam_struct>;
using BeamStructArray2D = FTypeArray2D<BeamStruct>;
using BeamStructArray3D = FTypeArray3D<BeamStruct>;

using BeamStructAlloc1D = FTypeAlloc1D<
    BeamStructArray1D,
    allocate_beam_struct_container,
    deallocate_beam_struct_container,
    reallocate_beam_struct_container_data,
    access_beam_struct_container>;

class AperturePointStruct;

using AperturePointStructArray1D = FTypeArray1D<
    AperturePointStruct,
    allocate_fortran_aperture_point_struct,
    deallocate_fortran_aperture_point_struct>;
using AperturePointStructArray2D = FTypeArray2D<AperturePointStruct>;
using AperturePointStructArray3D = FTypeArray3D<AperturePointStruct>;

using AperturePointStructAlloc1D = FTypeAlloc1D<
    AperturePointStructArray1D,
    allocate_aperture_point_struct_container,
    deallocate_aperture_point_struct_container,
    reallocate_aperture_point_struct_container_data,
    access_aperture_point_struct_container>;

class ApertureParamStruct;

using ApertureParamStructArray1D = FTypeArray1D<
    ApertureParamStruct,
    allocate_fortran_aperture_param_struct,
    deallocate_fortran_aperture_param_struct>;
using ApertureParamStructArray2D = FTypeArray2D<ApertureParamStruct>;
using ApertureParamStructArray3D = FTypeArray3D<ApertureParamStruct>;

using ApertureParamStructAlloc1D = FTypeAlloc1D<
    ApertureParamStructArray1D,
    allocate_aperture_param_struct_container,
    deallocate_aperture_param_struct_container,
    reallocate_aperture_param_struct_container_data,
    access_aperture_param_struct_container>;

class ApertureScanStruct;

using ApertureScanStructArray1D = FTypeArray1D<
    ApertureScanStruct,
    allocate_fortran_aperture_scan_struct,
    deallocate_fortran_aperture_scan_struct>;
using ApertureScanStructArray2D = FTypeArray2D<ApertureScanStruct>;
using ApertureScanStructArray3D = FTypeArray3D<ApertureScanStruct>;

using ApertureScanStructAlloc1D = FTypeAlloc1D<
    ApertureScanStructArray1D,
    allocate_aperture_scan_struct_container,
    deallocate_aperture_scan_struct_container,
    reallocate_aperture_scan_struct_container_data,
    access_aperture_scan_struct_container>;

class ElePointerStruct;

using ElePointerStructArray1D = FTypeArray1D<
    ElePointerStruct,
    allocate_fortran_ele_pointer_struct,
    deallocate_fortran_ele_pointer_struct>;
using ElePointerStructArray2D = FTypeArray2D<ElePointerStruct>;
using ElePointerStructArray3D = FTypeArray3D<ElePointerStruct>;

using ElePointerStructAlloc1D = FTypeAlloc1D<
    ElePointerStructArray1D,
    allocate_ele_pointer_struct_container,
    deallocate_ele_pointer_struct_container,
    reallocate_ele_pointer_struct_container_data,
    access_ele_pointer_struct_container>;

class ExpressionTreeStruct;

using ExpressionTreeStructArray1D = FTypeArray1D<
    ExpressionTreeStruct,
    allocate_fortran_expression_tree_struct,
    deallocate_fortran_expression_tree_struct>;
using ExpressionTreeStructArray2D = FTypeArray2D<ExpressionTreeStruct>;
using ExpressionTreeStructArray3D = FTypeArray3D<ExpressionTreeStruct>;

using ExpressionTreeStructAlloc1D = FTypeAlloc1D<
    ExpressionTreeStructArray1D,
    allocate_expression_tree_struct_container,
    deallocate_expression_tree_struct_container,
    reallocate_expression_tree_struct_container_data,
    access_expression_tree_struct_container>;

class NametableStruct;

using NametableStructArray1D = FTypeArray1D<
    NametableStruct,
    allocate_fortran_nametable_struct,
    deallocate_fortran_nametable_struct>;
using NametableStructArray2D = FTypeArray2D<NametableStruct>;
using NametableStructArray3D = FTypeArray3D<NametableStruct>;

using NametableStructAlloc1D = FTypeAlloc1D<
    NametableStructArray1D,
    allocate_nametable_struct_container,
    deallocate_nametable_struct_container,
    reallocate_nametable_struct_container_data,
    access_nametable_struct_container>;

class TaoSpinDnDpzStruct;

using TaoSpinDnDpzStructArray1D = FTypeArray1D<
    TaoSpinDnDpzStruct,
    allocate_fortran_tao_spin_dn_dpz_struct,
    deallocate_fortran_tao_spin_dn_dpz_struct>;
using TaoSpinDnDpzStructArray2D = FTypeArray2D<TaoSpinDnDpzStruct>;
using TaoSpinDnDpzStructArray3D = FTypeArray3D<TaoSpinDnDpzStruct>;

using TaoSpinDnDpzStructAlloc1D = FTypeAlloc1D<
    TaoSpinDnDpzStructArray1D,
    allocate_tao_spin_dn_dpz_struct_container,
    deallocate_tao_spin_dn_dpz_struct_container,
    reallocate_tao_spin_dn_dpz_struct_container_data,
    access_tao_spin_dn_dpz_struct_container>;

class ResonanceHStruct;

using ResonanceHStructArray1D = FTypeArray1D<
    ResonanceHStruct,
    allocate_fortran_resonance_h_struct,
    deallocate_fortran_resonance_h_struct>;
using ResonanceHStructArray2D = FTypeArray2D<ResonanceHStruct>;
using ResonanceHStructArray3D = FTypeArray3D<ResonanceHStruct>;

using ResonanceHStructAlloc1D = FTypeAlloc1D<
    ResonanceHStructArray1D,
    allocate_resonance_h_struct_container,
    deallocate_resonance_h_struct_container,
    reallocate_resonance_h_struct_container_data,
    access_resonance_h_struct_container>;

class SpinOrbitMap1Struct;

using SpinOrbitMap1StructArray1D = FTypeArray1D<
    SpinOrbitMap1Struct,
    allocate_fortran_spin_orbit_map1_struct,
    deallocate_fortran_spin_orbit_map1_struct>;
using SpinOrbitMap1StructArray2D = FTypeArray2D<SpinOrbitMap1Struct>;
using SpinOrbitMap1StructArray3D = FTypeArray3D<SpinOrbitMap1Struct>;

using SpinOrbitMap1StructAlloc1D = FTypeAlloc1D<
    SpinOrbitMap1StructArray1D,
    allocate_spin_orbit_map1_struct_container,
    deallocate_spin_orbit_map1_struct_container,
    reallocate_spin_orbit_map1_struct_container_data,
    access_spin_orbit_map1_struct_container>;

class SpinAxisStruct;

using SpinAxisStructArray1D = FTypeArray1D<
    SpinAxisStruct,
    allocate_fortran_spin_axis_struct,
    deallocate_fortran_spin_axis_struct>;
using SpinAxisStructArray2D = FTypeArray2D<SpinAxisStruct>;
using SpinAxisStructArray3D = FTypeArray3D<SpinAxisStruct>;

using SpinAxisStructAlloc1D = FTypeAlloc1D<
    SpinAxisStructArray1D,
    allocate_spin_axis_struct_container,
    deallocate_spin_axis_struct_container,
    reallocate_spin_axis_struct_container_data,
    access_spin_axis_struct_container>;

class PtcNormalFormStruct;

using PtcNormalFormStructArray1D = FTypeArray1D<
    PtcNormalFormStruct,
    allocate_fortran_ptc_normal_form_struct,
    deallocate_fortran_ptc_normal_form_struct>;
using PtcNormalFormStructArray2D = FTypeArray2D<PtcNormalFormStruct>;
using PtcNormalFormStructArray3D = FTypeArray3D<PtcNormalFormStruct>;

using PtcNormalFormStructAlloc1D = FTypeAlloc1D<
    PtcNormalFormStructArray1D,
    allocate_ptc_normal_form_struct_container,
    deallocate_ptc_normal_form_struct_container,
    reallocate_ptc_normal_form_struct_container_data,
    access_ptc_normal_form_struct_container>;

class BmadNormalFormStruct;

using BmadNormalFormStructArray1D = FTypeArray1D<
    BmadNormalFormStruct,
    allocate_fortran_bmad_normal_form_struct,
    deallocate_fortran_bmad_normal_form_struct>;
using BmadNormalFormStructArray2D = FTypeArray2D<BmadNormalFormStruct>;
using BmadNormalFormStructArray3D = FTypeArray3D<BmadNormalFormStruct>;

using BmadNormalFormStructAlloc1D = FTypeAlloc1D<
    BmadNormalFormStructArray1D,
    allocate_bmad_normal_form_struct_container,
    deallocate_bmad_normal_form_struct_container,
    reallocate_bmad_normal_form_struct_container_data,
    access_bmad_normal_form_struct_container>;

class BunchTrackStruct;

using BunchTrackStructArray1D = FTypeArray1D<
    BunchTrackStruct,
    allocate_fortran_bunch_track_struct,
    deallocate_fortran_bunch_track_struct>;
using BunchTrackStructArray2D = FTypeArray2D<BunchTrackStruct>;
using BunchTrackStructArray3D = FTypeArray3D<BunchTrackStruct>;

using BunchTrackStructAlloc1D = FTypeAlloc1D<
    BunchTrackStructArray1D,
    allocate_bunch_track_struct_container,
    deallocate_bunch_track_struct_container,
    reallocate_bunch_track_struct_container_data,
    access_bunch_track_struct_container>;

class SummationRdtStruct;

using SummationRdtStructArray1D = FTypeArray1D<
    SummationRdtStruct,
    allocate_fortran_summation_rdt_struct,
    deallocate_fortran_summation_rdt_struct>;
using SummationRdtStructArray2D = FTypeArray2D<SummationRdtStruct>;
using SummationRdtStructArray3D = FTypeArray3D<SummationRdtStruct>;

using SummationRdtStructAlloc1D = FTypeAlloc1D<
    SummationRdtStructArray1D,
    allocate_summation_rdt_struct_container,
    deallocate_summation_rdt_struct_container,
    reallocate_summation_rdt_struct_container_data,
    access_summation_rdt_struct_container>;

class TaoEleShapeStruct;

using TaoEleShapeStructArray1D = FTypeArray1D<
    TaoEleShapeStruct,
    allocate_fortran_tao_ele_shape_struct,
    deallocate_fortran_tao_ele_shape_struct>;
using TaoEleShapeStructArray2D = FTypeArray2D<TaoEleShapeStruct>;
using TaoEleShapeStructArray3D = FTypeArray3D<TaoEleShapeStruct>;

using TaoEleShapeStructAlloc1D = FTypeAlloc1D<
    TaoEleShapeStructArray1D,
    allocate_tao_ele_shape_struct_container,
    deallocate_tao_ele_shape_struct_container,
    reallocate_tao_ele_shape_struct_container_data,
    access_tao_ele_shape_struct_container>;

class TaoElePointerStruct;

using TaoElePointerStructArray1D = FTypeArray1D<
    TaoElePointerStruct,
    allocate_fortran_tao_ele_pointer_struct,
    deallocate_fortran_tao_ele_pointer_struct>;
using TaoElePointerStructArray2D = FTypeArray2D<TaoElePointerStruct>;
using TaoElePointerStructArray3D = FTypeArray3D<TaoElePointerStruct>;

using TaoElePointerStructAlloc1D = FTypeAlloc1D<
    TaoElePointerStructArray1D,
    allocate_tao_ele_pointer_struct_container,
    deallocate_tao_ele_pointer_struct_container,
    reallocate_tao_ele_pointer_struct_container_data,
    access_tao_ele_pointer_struct_container>;

class TaoCurveStruct;

using TaoCurveStructArray1D = FTypeArray1D<
    TaoCurveStruct,
    allocate_fortran_tao_curve_struct,
    deallocate_fortran_tao_curve_struct>;
using TaoCurveStructArray2D = FTypeArray2D<TaoCurveStruct>;
using TaoCurveStructArray3D = FTypeArray3D<TaoCurveStruct>;

using TaoCurveStructAlloc1D = FTypeAlloc1D<
    TaoCurveStructArray1D,
    allocate_tao_curve_struct_container,
    deallocate_tao_curve_struct_container,
    reallocate_tao_curve_struct_container_data,
    access_tao_curve_struct_container>;

class TaoCurveColorStruct;

using TaoCurveColorStructArray1D = FTypeArray1D<
    TaoCurveColorStruct,
    allocate_fortran_tao_curve_color_struct,
    deallocate_fortran_tao_curve_color_struct>;
using TaoCurveColorStructArray2D = FTypeArray2D<TaoCurveColorStruct>;
using TaoCurveColorStructArray3D = FTypeArray3D<TaoCurveColorStruct>;

using TaoCurveColorStructAlloc1D = FTypeAlloc1D<
    TaoCurveColorStructArray1D,
    allocate_tao_curve_color_struct_container,
    deallocate_tao_curve_color_struct_container,
    reallocate_tao_curve_color_struct_container_data,
    access_tao_curve_color_struct_container>;

class TaoCurveOrbitStruct;

using TaoCurveOrbitStructArray1D = FTypeArray1D<
    TaoCurveOrbitStruct,
    allocate_fortran_tao_curve_orbit_struct,
    deallocate_fortran_tao_curve_orbit_struct>;
using TaoCurveOrbitStructArray2D = FTypeArray2D<TaoCurveOrbitStruct>;
using TaoCurveOrbitStructArray3D = FTypeArray3D<TaoCurveOrbitStruct>;

using TaoCurveOrbitStructAlloc1D = FTypeAlloc1D<
    TaoCurveOrbitStructArray1D,
    allocate_tao_curve_orbit_struct_container,
    deallocate_tao_curve_orbit_struct_container,
    reallocate_tao_curve_orbit_struct_container_data,
    access_tao_curve_orbit_struct_container>;

class TaoHistogramStruct;

using TaoHistogramStructArray1D = FTypeArray1D<
    TaoHistogramStruct,
    allocate_fortran_tao_histogram_struct,
    deallocate_fortran_tao_histogram_struct>;
using TaoHistogramStructArray2D = FTypeArray2D<TaoHistogramStruct>;
using TaoHistogramStructArray3D = FTypeArray3D<TaoHistogramStruct>;

using TaoHistogramStructAlloc1D = FTypeAlloc1D<
    TaoHistogramStructArray1D,
    allocate_tao_histogram_struct_container,
    deallocate_tao_histogram_struct_container,
    reallocate_tao_histogram_struct_container_data,
    access_tao_histogram_struct_container>;

class LatEleOrder1Struct;

using LatEleOrder1StructArray1D = FTypeArray1D<
    LatEleOrder1Struct,
    allocate_fortran_lat_ele_order1_struct,
    deallocate_fortran_lat_ele_order1_struct>;
using LatEleOrder1StructArray2D = FTypeArray2D<LatEleOrder1Struct>;
using LatEleOrder1StructArray3D = FTypeArray3D<LatEleOrder1Struct>;

using LatEleOrder1StructAlloc1D = FTypeAlloc1D<
    LatEleOrder1StructArray1D,
    allocate_lat_ele_order1_struct_container,
    deallocate_lat_ele_order1_struct_container,
    reallocate_lat_ele_order1_struct_container_data,
    access_lat_ele_order1_struct_container>;

class LatEleOrderArrayStruct;

using LatEleOrderArrayStructArray1D = FTypeArray1D<
    LatEleOrderArrayStruct,
    allocate_fortran_lat_ele_order_array_struct,
    deallocate_fortran_lat_ele_order_array_struct>;
using LatEleOrderArrayStructArray2D = FTypeArray2D<LatEleOrderArrayStruct>;
using LatEleOrderArrayStructArray3D = FTypeArray3D<LatEleOrderArrayStruct>;

using LatEleOrderArrayStructAlloc1D = FTypeAlloc1D<
    LatEleOrderArrayStructArray1D,
    allocate_lat_ele_order_array_struct_container,
    deallocate_lat_ele_order_array_struct_container,
    reallocate_lat_ele_order_array_struct_container_data,
    access_lat_ele_order_array_struct_container>;

class TaoLatSigmaStruct;

using TaoLatSigmaStructArray1D = FTypeArray1D<
    TaoLatSigmaStruct,
    allocate_fortran_tao_lat_sigma_struct,
    deallocate_fortran_tao_lat_sigma_struct>;
using TaoLatSigmaStructArray2D = FTypeArray2D<TaoLatSigmaStruct>;
using TaoLatSigmaStructArray3D = FTypeArray3D<TaoLatSigmaStruct>;

using TaoLatSigmaStructAlloc1D = FTypeAlloc1D<
    TaoLatSigmaStructArray1D,
    allocate_tao_lat_sigma_struct_container,
    deallocate_tao_lat_sigma_struct_container,
    reallocate_tao_lat_sigma_struct_container_data,
    access_tao_lat_sigma_struct_container>;

class TaoSpinEleStruct;

using TaoSpinEleStructArray1D = FTypeArray1D<
    TaoSpinEleStruct,
    allocate_fortran_tao_spin_ele_struct,
    deallocate_fortran_tao_spin_ele_struct>;
using TaoSpinEleStructArray2D = FTypeArray2D<TaoSpinEleStruct>;
using TaoSpinEleStructArray3D = FTypeArray3D<TaoSpinEleStruct>;

using TaoSpinEleStructAlloc1D = FTypeAlloc1D<
    TaoSpinEleStructArray1D,
    allocate_tao_spin_ele_struct_container,
    deallocate_tao_spin_ele_struct_container,
    reallocate_tao_spin_ele_struct_container_data,
    access_tao_spin_ele_struct_container>;

class TaoPlotCacheStruct;

using TaoPlotCacheStructArray1D = FTypeArray1D<
    TaoPlotCacheStruct,
    allocate_fortran_tao_plot_cache_struct,
    deallocate_fortran_tao_plot_cache_struct>;
using TaoPlotCacheStructArray2D = FTypeArray2D<TaoPlotCacheStruct>;
using TaoPlotCacheStructArray3D = FTypeArray3D<TaoPlotCacheStruct>;

using TaoPlotCacheStructAlloc1D = FTypeAlloc1D<
    TaoPlotCacheStructArray1D,
    allocate_tao_plot_cache_struct_container,
    deallocate_tao_plot_cache_struct_container,
    reallocate_tao_plot_cache_struct_container_data,
    access_tao_plot_cache_struct_container>;

class TaoSpinPolarizationStruct;

using TaoSpinPolarizationStructArray1D = FTypeArray1D<
    TaoSpinPolarizationStruct,
    allocate_fortran_tao_spin_polarization_struct,
    deallocate_fortran_tao_spin_polarization_struct>;
using TaoSpinPolarizationStructArray2D = FTypeArray2D<TaoSpinPolarizationStruct>;
using TaoSpinPolarizationStructArray3D = FTypeArray3D<TaoSpinPolarizationStruct>;

using TaoSpinPolarizationStructAlloc1D = FTypeAlloc1D<
    TaoSpinPolarizationStructArray1D,
    allocate_tao_spin_polarization_struct_container,
    deallocate_tao_spin_polarization_struct_container,
    reallocate_tao_spin_polarization_struct_container_data,
    access_tao_spin_polarization_struct_container>;

class TaoLatticeBranchStruct;

using TaoLatticeBranchStructArray1D = FTypeArray1D<
    TaoLatticeBranchStruct,
    allocate_fortran_tao_lattice_branch_struct,
    deallocate_fortran_tao_lattice_branch_struct>;
using TaoLatticeBranchStructArray2D = FTypeArray2D<TaoLatticeBranchStruct>;
using TaoLatticeBranchStructArray3D = FTypeArray3D<TaoLatticeBranchStruct>;

using TaoLatticeBranchStructAlloc1D = FTypeAlloc1D<
    TaoLatticeBranchStructArray1D,
    allocate_tao_lattice_branch_struct_container,
    deallocate_tao_lattice_branch_struct_container,
    reallocate_tao_lattice_branch_struct_container_data,
    access_tao_lattice_branch_struct_container>;

class TaoModelElementStruct;

using TaoModelElementStructArray1D = FTypeArray1D<
    TaoModelElementStruct,
    allocate_fortran_tao_model_element_struct,
    deallocate_fortran_tao_model_element_struct>;
using TaoModelElementStructArray2D = FTypeArray2D<TaoModelElementStruct>;
using TaoModelElementStructArray3D = FTypeArray3D<TaoModelElementStruct>;

using TaoModelElementStructAlloc1D = FTypeAlloc1D<
    TaoModelElementStructArray1D,
    allocate_tao_model_element_struct_container,
    deallocate_tao_model_element_struct_container,
    reallocate_tao_model_element_struct_container_data,
    access_tao_model_element_struct_container>;

class TaoBeamBranchStruct;

using TaoBeamBranchStructArray1D = FTypeArray1D<
    TaoBeamBranchStruct,
    allocate_fortran_tao_beam_branch_struct,
    deallocate_fortran_tao_beam_branch_struct>;
using TaoBeamBranchStructArray2D = FTypeArray2D<TaoBeamBranchStruct>;
using TaoBeamBranchStructArray3D = FTypeArray3D<TaoBeamBranchStruct>;

using TaoBeamBranchStructAlloc1D = FTypeAlloc1D<
    TaoBeamBranchStructArray1D,
    allocate_tao_beam_branch_struct_container,
    deallocate_tao_beam_branch_struct_container,
    reallocate_tao_beam_branch_struct_container_data,
    access_tao_beam_branch_struct_container>;

class TaoD1DataStruct;

using TaoD1DataStructArray1D = FTypeArray1D<
    TaoD1DataStruct,
    allocate_fortran_tao_d1_data_struct,
    deallocate_fortran_tao_d1_data_struct>;
using TaoD1DataStructArray2D = FTypeArray2D<TaoD1DataStruct>;
using TaoD1DataStructArray3D = FTypeArray3D<TaoD1DataStruct>;

using TaoD1DataStructAlloc1D = FTypeAlloc1D<
    TaoD1DataStructArray1D,
    allocate_tao_d1_data_struct_container,
    deallocate_tao_d1_data_struct_container,
    reallocate_tao_d1_data_struct_container_data,
    access_tao_d1_data_struct_container>;

class TaoD2DataStruct;

using TaoD2DataStructArray1D = FTypeArray1D<
    TaoD2DataStruct,
    allocate_fortran_tao_d2_data_struct,
    deallocate_fortran_tao_d2_data_struct>;
using TaoD2DataStructArray2D = FTypeArray2D<TaoD2DataStruct>;
using TaoD2DataStructArray3D = FTypeArray3D<TaoD2DataStruct>;

using TaoD2DataStructAlloc1D = FTypeAlloc1D<
    TaoD2DataStructArray1D,
    allocate_tao_d2_data_struct_container,
    deallocate_tao_d2_data_struct_container,
    reallocate_tao_d2_data_struct_container_data,
    access_tao_d2_data_struct_container>;

class TaoDataVarComponentStruct;

using TaoDataVarComponentStructArray1D = FTypeArray1D<
    TaoDataVarComponentStruct,
    allocate_fortran_tao_data_var_component_struct,
    deallocate_fortran_tao_data_var_component_struct>;
using TaoDataVarComponentStructArray2D = FTypeArray2D<TaoDataVarComponentStruct>;
using TaoDataVarComponentStructArray3D = FTypeArray3D<TaoDataVarComponentStruct>;

using TaoDataVarComponentStructAlloc1D = FTypeAlloc1D<
    TaoDataVarComponentStructArray1D,
    allocate_tao_data_var_component_struct_container,
    deallocate_tao_data_var_component_struct_container,
    reallocate_tao_data_var_component_struct_container_data,
    access_tao_data_var_component_struct_container>;

class TaoGraphStruct;

using TaoGraphStructArray1D = FTypeArray1D<
    TaoGraphStruct,
    allocate_fortran_tao_graph_struct,
    deallocate_fortran_tao_graph_struct>;
using TaoGraphStructArray2D = FTypeArray2D<TaoGraphStruct>;
using TaoGraphStructArray3D = FTypeArray3D<TaoGraphStruct>;

using TaoGraphStructAlloc1D = FTypeAlloc1D<
    TaoGraphStructArray1D,
    allocate_tao_graph_struct_container,
    deallocate_tao_graph_struct_container,
    reallocate_tao_graph_struct_container_data,
    access_tao_graph_struct_container>;

class TaoPlotStruct;

using TaoPlotStructArray1D = FTypeArray1D<
    TaoPlotStruct,
    allocate_fortran_tao_plot_struct,
    deallocate_fortran_tao_plot_struct>;
using TaoPlotStructArray2D = FTypeArray2D<TaoPlotStruct>;
using TaoPlotStructArray3D = FTypeArray3D<TaoPlotStruct>;

using TaoPlotStructAlloc1D = FTypeAlloc1D<
    TaoPlotStructArray1D,
    allocate_tao_plot_struct_container,
    deallocate_tao_plot_struct_container,
    reallocate_tao_plot_struct_container_data,
    access_tao_plot_struct_container>;

class TaoPlotRegionStruct;

using TaoPlotRegionStructArray1D = FTypeArray1D<
    TaoPlotRegionStruct,
    allocate_fortran_tao_plot_region_struct,
    deallocate_fortran_tao_plot_region_struct>;
using TaoPlotRegionStructArray2D = FTypeArray2D<TaoPlotRegionStruct>;
using TaoPlotRegionStructArray3D = FTypeArray3D<TaoPlotRegionStruct>;

using TaoPlotRegionStructAlloc1D = FTypeAlloc1D<
    TaoPlotRegionStructArray1D,
    allocate_tao_plot_region_struct_container,
    deallocate_tao_plot_region_struct_container,
    reallocate_tao_plot_region_struct_container_data,
    access_tao_plot_region_struct_container>;

class TaoUniversePointerStruct;

using TaoUniversePointerStructArray1D = FTypeArray1D<
    TaoUniversePointerStruct,
    allocate_fortran_tao_universe_pointer_struct,
    deallocate_fortran_tao_universe_pointer_struct>;
using TaoUniversePointerStructArray2D = FTypeArray2D<TaoUniversePointerStruct>;
using TaoUniversePointerStructArray3D = FTypeArray3D<TaoUniversePointerStruct>;

using TaoUniversePointerStructAlloc1D = FTypeAlloc1D<
    TaoUniversePointerStructArray1D,
    allocate_tao_universe_pointer_struct_container,
    deallocate_tao_universe_pointer_struct_container,
    reallocate_tao_universe_pointer_struct_container_data,
    access_tao_universe_pointer_struct_container>;

class TaoSuperUniverseStruct;

using TaoSuperUniverseStructArray1D = FTypeArray1D<
    TaoSuperUniverseStruct,
    allocate_fortran_tao_super_universe_struct,
    deallocate_fortran_tao_super_universe_struct>;
using TaoSuperUniverseStructArray2D = FTypeArray2D<TaoSuperUniverseStruct>;
using TaoSuperUniverseStructArray3D = FTypeArray3D<TaoSuperUniverseStruct>;

using TaoSuperUniverseStructAlloc1D = FTypeAlloc1D<
    TaoSuperUniverseStructArray1D,
    allocate_tao_super_universe_struct_container,
    deallocate_tao_super_universe_struct_container,
    reallocate_tao_super_universe_struct_container_data,
    access_tao_super_universe_struct_container>;

class TaoVarStruct;

using TaoVarStructArray1D =
    FTypeArray1D<TaoVarStruct, allocate_fortran_tao_var_struct, deallocate_fortran_tao_var_struct>;
using TaoVarStructArray2D = FTypeArray2D<TaoVarStruct>;
using TaoVarStructArray3D = FTypeArray3D<TaoVarStruct>;

using TaoVarStructAlloc1D = FTypeAlloc1D<
    TaoVarStructArray1D,
    allocate_tao_var_struct_container,
    deallocate_tao_var_struct_container,
    reallocate_tao_var_struct_container_data,
    access_tao_var_struct_container>;

class TaoVarSlaveStruct;

using TaoVarSlaveStructArray1D = FTypeArray1D<
    TaoVarSlaveStruct,
    allocate_fortran_tao_var_slave_struct,
    deallocate_fortran_tao_var_slave_struct>;
using TaoVarSlaveStructArray2D = FTypeArray2D<TaoVarSlaveStruct>;
using TaoVarSlaveStructArray3D = FTypeArray3D<TaoVarSlaveStruct>;

using TaoVarSlaveStructAlloc1D = FTypeAlloc1D<
    TaoVarSlaveStructArray1D,
    allocate_tao_var_slave_struct_container,
    deallocate_tao_var_slave_struct_container,
    reallocate_tao_var_slave_struct_container_data,
    access_tao_var_slave_struct_container>;

class TaoLatticeStruct;

using TaoLatticeStructArray1D = FTypeArray1D<
    TaoLatticeStruct,
    allocate_fortran_tao_lattice_struct,
    deallocate_fortran_tao_lattice_struct>;
using TaoLatticeStructArray2D = FTypeArray2D<TaoLatticeStruct>;
using TaoLatticeStructArray3D = FTypeArray3D<TaoLatticeStruct>;

using TaoLatticeStructAlloc1D = FTypeAlloc1D<
    TaoLatticeStructArray1D,
    allocate_tao_lattice_struct_container,
    deallocate_tao_lattice_struct_container,
    reallocate_tao_lattice_struct_container_data,
    access_tao_lattice_struct_container>;

class TaoBeamUniStruct;

using TaoBeamUniStructArray1D = FTypeArray1D<
    TaoBeamUniStruct,
    allocate_fortran_tao_beam_uni_struct,
    deallocate_fortran_tao_beam_uni_struct>;
using TaoBeamUniStructArray2D = FTypeArray2D<TaoBeamUniStruct>;
using TaoBeamUniStructArray3D = FTypeArray3D<TaoBeamUniStruct>;

using TaoBeamUniStructAlloc1D = FTypeAlloc1D<
    TaoBeamUniStructArray1D,
    allocate_tao_beam_uni_struct_container,
    deallocate_tao_beam_uni_struct_container,
    reallocate_tao_beam_uni_struct_container_data,
    access_tao_beam_uni_struct_container>;

class TaoDynamicApertureStruct;

using TaoDynamicApertureStructArray1D = FTypeArray1D<
    TaoDynamicApertureStruct,
    allocate_fortran_tao_dynamic_aperture_struct,
    deallocate_fortran_tao_dynamic_aperture_struct>;
using TaoDynamicApertureStructArray2D = FTypeArray2D<TaoDynamicApertureStruct>;
using TaoDynamicApertureStructArray3D = FTypeArray3D<TaoDynamicApertureStruct>;

using TaoDynamicApertureStructAlloc1D = FTypeAlloc1D<
    TaoDynamicApertureStructArray1D,
    allocate_tao_dynamic_aperture_struct_container,
    deallocate_tao_dynamic_aperture_struct_container,
    reallocate_tao_dynamic_aperture_struct_container_data,
    access_tao_dynamic_aperture_struct_container>;

class TaoModelBranchStruct;

using TaoModelBranchStructArray1D = FTypeArray1D<
    TaoModelBranchStruct,
    allocate_fortran_tao_model_branch_struct,
    deallocate_fortran_tao_model_branch_struct>;
using TaoModelBranchStructArray2D = FTypeArray2D<TaoModelBranchStruct>;
using TaoModelBranchStructArray3D = FTypeArray3D<TaoModelBranchStruct>;

using TaoModelBranchStructAlloc1D = FTypeAlloc1D<
    TaoModelBranchStructArray1D,
    allocate_tao_model_branch_struct_container,
    deallocate_tao_model_branch_struct_container,
    reallocate_tao_model_branch_struct_container_data,
    access_tao_model_branch_struct_container>;

class TaoSpinMapStruct;

using TaoSpinMapStructArray1D = FTypeArray1D<
    TaoSpinMapStruct,
    allocate_fortran_tao_spin_map_struct,
    deallocate_fortran_tao_spin_map_struct>;
using TaoSpinMapStructArray2D = FTypeArray2D<TaoSpinMapStruct>;
using TaoSpinMapStructArray3D = FTypeArray3D<TaoSpinMapStruct>;

using TaoSpinMapStructAlloc1D = FTypeAlloc1D<
    TaoSpinMapStructArray1D,
    allocate_tao_spin_map_struct_container,
    deallocate_tao_spin_map_struct_container,
    reallocate_tao_spin_map_struct_container_data,
    access_tao_spin_map_struct_container>;

class TaoDataStruct;

using TaoDataStructArray1D = FTypeArray1D<
    TaoDataStruct,
    allocate_fortran_tao_data_struct,
    deallocate_fortran_tao_data_struct>;
using TaoDataStructArray2D = FTypeArray2D<TaoDataStruct>;
using TaoDataStructArray3D = FTypeArray3D<TaoDataStruct>;

using TaoDataStructAlloc1D = FTypeAlloc1D<
    TaoDataStructArray1D,
    allocate_tao_data_struct_container,
    deallocate_tao_data_struct_container,
    reallocate_tao_data_struct_container_data,
    access_tao_data_struct_container>;

class TaoPingScaleStruct;

using TaoPingScaleStructArray1D = FTypeArray1D<
    TaoPingScaleStruct,
    allocate_fortran_tao_ping_scale_struct,
    deallocate_fortran_tao_ping_scale_struct>;
using TaoPingScaleStructArray2D = FTypeArray2D<TaoPingScaleStruct>;
using TaoPingScaleStructArray3D = FTypeArray3D<TaoPingScaleStruct>;

using TaoPingScaleStructAlloc1D = FTypeAlloc1D<
    TaoPingScaleStructArray1D,
    allocate_tao_ping_scale_struct_container,
    deallocate_tao_ping_scale_struct_container,
    reallocate_tao_ping_scale_struct_container_data,
    access_tao_ping_scale_struct_container>;

class TaoUniverseCalcStruct;

using TaoUniverseCalcStructArray1D = FTypeArray1D<
    TaoUniverseCalcStruct,
    allocate_fortran_tao_universe_calc_struct,
    deallocate_fortran_tao_universe_calc_struct>;
using TaoUniverseCalcStructArray2D = FTypeArray2D<TaoUniverseCalcStruct>;
using TaoUniverseCalcStructArray3D = FTypeArray3D<TaoUniverseCalcStruct>;

using TaoUniverseCalcStructAlloc1D = FTypeAlloc1D<
    TaoUniverseCalcStructArray1D,
    allocate_tao_universe_calc_struct_container,
    deallocate_tao_universe_calc_struct_container,
    reallocate_tao_universe_calc_struct_container_data,
    access_tao_universe_calc_struct_container>;

class LatEleOrderStruct;

using LatEleOrderStructArray1D = FTypeArray1D<
    LatEleOrderStruct,
    allocate_fortran_lat_ele_order_struct,
    deallocate_fortran_lat_ele_order_struct>;
using LatEleOrderStructArray2D = FTypeArray2D<LatEleOrderStruct>;
using LatEleOrderStructArray3D = FTypeArray3D<LatEleOrderStruct>;

using LatEleOrderStructAlloc1D = FTypeAlloc1D<
    LatEleOrderStructArray1D,
    allocate_lat_ele_order_struct_container,
    deallocate_lat_ele_order_struct_container,
    reallocate_lat_ele_order_struct_container_data,
    access_lat_ele_order_struct_container>;

class TaoExpressionInfoStruct;

using TaoExpressionInfoStructArray1D = FTypeArray1D<
    TaoExpressionInfoStruct,
    allocate_fortran_tao_expression_info_struct,
    deallocate_fortran_tao_expression_info_struct>;
using TaoExpressionInfoStructArray2D = FTypeArray2D<TaoExpressionInfoStruct>;
using TaoExpressionInfoStructArray3D = FTypeArray3D<TaoExpressionInfoStruct>;

using TaoExpressionInfoStructAlloc1D = FTypeAlloc1D<
    TaoExpressionInfoStructArray1D,
    allocate_tao_expression_info_struct_container,
    deallocate_tao_expression_info_struct_container,
    reallocate_tao_expression_info_struct_container_data,
    access_tao_expression_info_struct_container>;

class TaoEvalNodeStruct;

using TaoEvalNodeStructArray1D = FTypeArray1D<
    TaoEvalNodeStruct,
    allocate_fortran_tao_eval_node_struct,
    deallocate_fortran_tao_eval_node_struct>;
using TaoEvalNodeStructArray2D = FTypeArray2D<TaoEvalNodeStruct>;
using TaoEvalNodeStructArray3D = FTypeArray3D<TaoEvalNodeStruct>;

using TaoEvalNodeStructAlloc1D = FTypeAlloc1D<
    TaoEvalNodeStructArray1D,
    allocate_tao_eval_node_struct_container,
    deallocate_tao_eval_node_struct_container,
    reallocate_tao_eval_node_struct_container_data,
    access_tao_eval_node_struct_container>;

class TaoTitleStruct;

using TaoTitleStructArray1D = FTypeArray1D<
    TaoTitleStruct,
    allocate_fortran_tao_title_struct,
    deallocate_fortran_tao_title_struct>;
using TaoTitleStructArray2D = FTypeArray2D<TaoTitleStruct>;
using TaoTitleStructArray3D = FTypeArray3D<TaoTitleStruct>;

using TaoTitleStructAlloc1D = FTypeAlloc1D<
    TaoTitleStructArray1D,
    allocate_tao_title_struct_container,
    deallocate_tao_title_struct_container,
    reallocate_tao_title_struct_container_data,
    access_tao_title_struct_container>;

class QpRectStruct;

using QpRectStructArray1D =
    FTypeArray1D<QpRectStruct, allocate_fortran_qp_rect_struct, deallocate_fortran_qp_rect_struct>;
using QpRectStructArray2D = FTypeArray2D<QpRectStruct>;
using QpRectStructArray3D = FTypeArray3D<QpRectStruct>;

using QpRectStructAlloc1D = FTypeAlloc1D<
    QpRectStructArray1D,
    allocate_qp_rect_struct_container,
    deallocate_qp_rect_struct_container,
    reallocate_qp_rect_struct_container_data,
    access_qp_rect_struct_container>;

class TaoDrawingStruct;

using TaoDrawingStructArray1D = FTypeArray1D<
    TaoDrawingStruct,
    allocate_fortran_tao_drawing_struct,
    deallocate_fortran_tao_drawing_struct>;
using TaoDrawingStructArray2D = FTypeArray2D<TaoDrawingStruct>;
using TaoDrawingStructArray3D = FTypeArray3D<TaoDrawingStruct>;

using TaoDrawingStructAlloc1D = FTypeAlloc1D<
    TaoDrawingStructArray1D,
    allocate_tao_drawing_struct_container,
    deallocate_tao_drawing_struct_container,
    reallocate_tao_drawing_struct_container_data,
    access_tao_drawing_struct_container>;

class TaoShapePatternStruct;

using TaoShapePatternStructArray1D = FTypeArray1D<
    TaoShapePatternStruct,
    allocate_fortran_tao_shape_pattern_struct,
    deallocate_fortran_tao_shape_pattern_struct>;
using TaoShapePatternStructArray2D = FTypeArray2D<TaoShapePatternStruct>;
using TaoShapePatternStructArray3D = FTypeArray3D<TaoShapePatternStruct>;

using TaoShapePatternStructAlloc1D = FTypeAlloc1D<
    TaoShapePatternStructArray1D,
    allocate_tao_shape_pattern_struct_container,
    deallocate_tao_shape_pattern_struct_container,
    reallocate_tao_shape_pattern_struct_container_data,
    access_tao_shape_pattern_struct_container>;

class TaoShapePatternPointStruct;

using TaoShapePatternPointStructArray1D = FTypeArray1D<
    TaoShapePatternPointStruct,
    allocate_fortran_tao_shape_pattern_point_struct,
    deallocate_fortran_tao_shape_pattern_point_struct>;
using TaoShapePatternPointStructArray2D = FTypeArray2D<TaoShapePatternPointStruct>;
using TaoShapePatternPointStructArray3D = FTypeArray3D<TaoShapePatternPointStruct>;

using TaoShapePatternPointStructAlloc1D = FTypeAlloc1D<
    TaoShapePatternPointStructArray1D,
    allocate_tao_shape_pattern_point_struct_container,
    deallocate_tao_shape_pattern_point_struct_container,
    reallocate_tao_shape_pattern_point_struct_container_data,
    access_tao_shape_pattern_point_struct_container>;

class QpAxisStruct;

using QpAxisStructArray1D =
    FTypeArray1D<QpAxisStruct, allocate_fortran_qp_axis_struct, deallocate_fortran_qp_axis_struct>;
using QpAxisStructArray2D = FTypeArray2D<QpAxisStruct>;
using QpAxisStructArray3D = FTypeArray3D<QpAxisStruct>;

using QpAxisStructAlloc1D = FTypeAlloc1D<
    QpAxisStructArray1D,
    allocate_qp_axis_struct_container,
    deallocate_qp_axis_struct_container,
    reallocate_qp_axis_struct_container_data,
    access_qp_axis_struct_container>;

class QpLegendStruct;

using QpLegendStructArray1D = FTypeArray1D<
    QpLegendStruct,
    allocate_fortran_qp_legend_struct,
    deallocate_fortran_qp_legend_struct>;
using QpLegendStructArray2D = FTypeArray2D<QpLegendStruct>;
using QpLegendStructArray3D = FTypeArray3D<QpLegendStruct>;

using QpLegendStructAlloc1D = FTypeAlloc1D<
    QpLegendStructArray1D,
    allocate_qp_legend_struct_container,
    deallocate_qp_legend_struct_container,
    reallocate_qp_legend_struct_container_data,
    access_qp_legend_struct_container>;

class QpPointStruct;

using QpPointStructArray1D = FTypeArray1D<
    QpPointStruct,
    allocate_fortran_qp_point_struct,
    deallocate_fortran_qp_point_struct>;
using QpPointStructArray2D = FTypeArray2D<QpPointStruct>;
using QpPointStructArray3D = FTypeArray3D<QpPointStruct>;

using QpPointStructAlloc1D = FTypeAlloc1D<
    QpPointStructArray1D,
    allocate_qp_point_struct_container,
    deallocate_qp_point_struct_container,
    reallocate_qp_point_struct_container_data,
    access_qp_point_struct_container>;

class QpLineStruct;

using QpLineStructArray1D =
    FTypeArray1D<QpLineStruct, allocate_fortran_qp_line_struct, deallocate_fortran_qp_line_struct>;
using QpLineStructArray2D = FTypeArray2D<QpLineStruct>;
using QpLineStructArray3D = FTypeArray3D<QpLineStruct>;

using QpLineStructAlloc1D = FTypeAlloc1D<
    QpLineStructArray1D,
    allocate_qp_line_struct_container,
    deallocate_qp_line_struct_container,
    reallocate_qp_line_struct_container_data,
    access_qp_line_struct_container>;

class QpSymbolStruct;

using QpSymbolStructArray1D = FTypeArray1D<
    QpSymbolStruct,
    allocate_fortran_qp_symbol_struct,
    deallocate_fortran_qp_symbol_struct>;
using QpSymbolStructArray2D = FTypeArray2D<QpSymbolStruct>;
using QpSymbolStructArray3D = FTypeArray3D<QpSymbolStruct>;

using QpSymbolStructAlloc1D = FTypeAlloc1D<
    QpSymbolStructArray1D,
    allocate_qp_symbol_struct_container,
    deallocate_qp_symbol_struct_container,
    reallocate_qp_symbol_struct_container_data,
    access_qp_symbol_struct_container>;

class TaoFloorPlanStruct;

using TaoFloorPlanStructArray1D = FTypeArray1D<
    TaoFloorPlanStruct,
    allocate_fortran_tao_floor_plan_struct,
    deallocate_fortran_tao_floor_plan_struct>;
using TaoFloorPlanStructArray2D = FTypeArray2D<TaoFloorPlanStruct>;
using TaoFloorPlanStructArray3D = FTypeArray3D<TaoFloorPlanStruct>;

using TaoFloorPlanStructAlloc1D = FTypeAlloc1D<
    TaoFloorPlanStructArray1D,
    allocate_tao_floor_plan_struct_container,
    deallocate_tao_floor_plan_struct_container,
    reallocate_tao_floor_plan_struct_container_data,
    access_tao_floor_plan_struct_container>;

class TaoV1VarStruct;

using TaoV1VarStructArray1D = FTypeArray1D<
    TaoV1VarStruct,
    allocate_fortran_tao_v1_var_struct,
    deallocate_fortran_tao_v1_var_struct>;
using TaoV1VarStructArray2D = FTypeArray2D<TaoV1VarStruct>;
using TaoV1VarStructArray3D = FTypeArray3D<TaoV1VarStruct>;

using TaoV1VarStructAlloc1D = FTypeAlloc1D<
    TaoV1VarStructArray1D,
    allocate_tao_v1_var_struct_container,
    deallocate_tao_v1_var_struct_container,
    reallocate_tao_v1_var_struct_container_data,
    access_tao_v1_var_struct_container>;

class TaoGlobalStruct;

using TaoGlobalStructArray1D = FTypeArray1D<
    TaoGlobalStruct,
    allocate_fortran_tao_global_struct,
    deallocate_fortran_tao_global_struct>;
using TaoGlobalStructArray2D = FTypeArray2D<TaoGlobalStruct>;
using TaoGlobalStructArray3D = FTypeArray3D<TaoGlobalStruct>;

using TaoGlobalStructAlloc1D = FTypeAlloc1D<
    TaoGlobalStructArray1D,
    allocate_tao_global_struct_container,
    deallocate_tao_global_struct_container,
    reallocate_tao_global_struct_container_data,
    access_tao_global_struct_container>;

class TaoInitStruct;

using TaoInitStructArray1D = FTypeArray1D<
    TaoInitStruct,
    allocate_fortran_tao_init_struct,
    deallocate_fortran_tao_init_struct>;
using TaoInitStructArray2D = FTypeArray2D<TaoInitStruct>;
using TaoInitStructArray3D = FTypeArray3D<TaoInitStruct>;

using TaoInitStructAlloc1D = FTypeAlloc1D<
    TaoInitStructArray1D,
    allocate_tao_init_struct_container,
    deallocate_tao_init_struct_container,
    reallocate_tao_init_struct_container_data,
    access_tao_init_struct_container>;

class TaoCommonStruct;

using TaoCommonStructArray1D = FTypeArray1D<
    TaoCommonStruct,
    allocate_fortran_tao_common_struct,
    deallocate_fortran_tao_common_struct>;
using TaoCommonStructArray2D = FTypeArray2D<TaoCommonStruct>;
using TaoCommonStructArray3D = FTypeArray3D<TaoCommonStruct>;

using TaoCommonStructAlloc1D = FTypeAlloc1D<
    TaoCommonStructArray1D,
    allocate_tao_common_struct_container,
    deallocate_tao_common_struct_container,
    reallocate_tao_common_struct_container_data,
    access_tao_common_struct_container>;

class TaoPlotPageStruct;

using TaoPlotPageStructArray1D = FTypeArray1D<
    TaoPlotPageStruct,
    allocate_fortran_tao_plot_page_struct,
    deallocate_fortran_tao_plot_page_struct>;
using TaoPlotPageStructArray2D = FTypeArray2D<TaoPlotPageStruct>;
using TaoPlotPageStructArray3D = FTypeArray3D<TaoPlotPageStruct>;

using TaoPlotPageStructAlloc1D = FTypeAlloc1D<
    TaoPlotPageStructArray1D,
    allocate_tao_plot_page_struct_container,
    deallocate_tao_plot_page_struct_container,
    reallocate_tao_plot_page_struct_container_data,
    access_tao_plot_page_struct_container>;

class TaoBuildingWallStruct;

using TaoBuildingWallStructArray1D = FTypeArray1D<
    TaoBuildingWallStruct,
    allocate_fortran_tao_building_wall_struct,
    deallocate_fortran_tao_building_wall_struct>;
using TaoBuildingWallStructArray2D = FTypeArray2D<TaoBuildingWallStruct>;
using TaoBuildingWallStructArray3D = FTypeArray3D<TaoBuildingWallStruct>;

using TaoBuildingWallStructAlloc1D = FTypeAlloc1D<
    TaoBuildingWallStructArray1D,
    allocate_tao_building_wall_struct_container,
    deallocate_tao_building_wall_struct_container,
    reallocate_tao_building_wall_struct_container_data,
    access_tao_building_wall_struct_container>;

class TaoBuildingWallOrientationStruct;

using TaoBuildingWallOrientationStructArray1D = FTypeArray1D<
    TaoBuildingWallOrientationStruct,
    allocate_fortran_tao_building_wall_orientation_struct,
    deallocate_fortran_tao_building_wall_orientation_struct>;
using TaoBuildingWallOrientationStructArray2D = FTypeArray2D<TaoBuildingWallOrientationStruct>;
using TaoBuildingWallOrientationStructArray3D = FTypeArray3D<TaoBuildingWallOrientationStruct>;

using TaoBuildingWallOrientationStructAlloc1D = FTypeAlloc1D<
    TaoBuildingWallOrientationStructArray1D,
    allocate_tao_building_wall_orientation_struct_container,
    deallocate_tao_building_wall_orientation_struct_container,
    reallocate_tao_building_wall_orientation_struct_container_data,
    access_tao_building_wall_orientation_struct_container>;

class TaoBuildingWallSectionStruct;

using TaoBuildingWallSectionStructArray1D = FTypeArray1D<
    TaoBuildingWallSectionStruct,
    allocate_fortran_tao_building_wall_section_struct,
    deallocate_fortran_tao_building_wall_section_struct>;
using TaoBuildingWallSectionStructArray2D = FTypeArray2D<TaoBuildingWallSectionStruct>;
using TaoBuildingWallSectionStructArray3D = FTypeArray3D<TaoBuildingWallSectionStruct>;

using TaoBuildingWallSectionStructAlloc1D = FTypeAlloc1D<
    TaoBuildingWallSectionStructArray1D,
    allocate_tao_building_wall_section_struct_container,
    deallocate_tao_building_wall_section_struct_container,
    reallocate_tao_building_wall_section_struct_container_data,
    access_tao_building_wall_section_struct_container>;

class TaoBuildingWallPointStruct;

using TaoBuildingWallPointStructArray1D = FTypeArray1D<
    TaoBuildingWallPointStruct,
    allocate_fortran_tao_building_wall_point_struct,
    deallocate_fortran_tao_building_wall_point_struct>;
using TaoBuildingWallPointStructArray2D = FTypeArray2D<TaoBuildingWallPointStruct>;
using TaoBuildingWallPointStructArray3D = FTypeArray3D<TaoBuildingWallPointStruct>;

using TaoBuildingWallPointStructAlloc1D = FTypeAlloc1D<
    TaoBuildingWallPointStructArray1D,
    allocate_tao_building_wall_point_struct_container,
    deallocate_tao_building_wall_point_struct_container,
    reallocate_tao_building_wall_point_struct_container_data,
    access_tao_building_wall_point_struct_container>;

class TaoWaveStruct;

using TaoWaveStructArray1D = FTypeArray1D<
    TaoWaveStruct,
    allocate_fortran_tao_wave_struct,
    deallocate_fortran_tao_wave_struct>;
using TaoWaveStructArray2D = FTypeArray2D<TaoWaveStruct>;
using TaoWaveStructArray3D = FTypeArray3D<TaoWaveStruct>;

using TaoWaveStructAlloc1D = FTypeAlloc1D<
    TaoWaveStructArray1D,
    allocate_tao_wave_struct_container,
    deallocate_tao_wave_struct_container,
    reallocate_tao_wave_struct_container_data,
    access_tao_wave_struct_container>;

class TaoWaveKickPtStruct;

using TaoWaveKickPtStructArray1D = FTypeArray1D<
    TaoWaveKickPtStruct,
    allocate_fortran_tao_wave_kick_pt_struct,
    deallocate_fortran_tao_wave_kick_pt_struct>;
using TaoWaveKickPtStructArray2D = FTypeArray2D<TaoWaveKickPtStruct>;
using TaoWaveKickPtStructArray3D = FTypeArray3D<TaoWaveKickPtStruct>;

using TaoWaveKickPtStructAlloc1D = FTypeAlloc1D<
    TaoWaveKickPtStructArray1D,
    allocate_tao_wave_kick_pt_struct_container,
    deallocate_tao_wave_kick_pt_struct_container,
    reallocate_tao_wave_kick_pt_struct_container_data,
    access_tao_wave_kick_pt_struct_container>;

class TaoCmdHistoryStruct;

using TaoCmdHistoryStructArray1D = FTypeArray1D<
    TaoCmdHistoryStruct,
    allocate_fortran_tao_cmd_history_struct,
    deallocate_fortran_tao_cmd_history_struct>;
using TaoCmdHistoryStructArray2D = FTypeArray2D<TaoCmdHistoryStruct>;
using TaoCmdHistoryStructArray3D = FTypeArray3D<TaoCmdHistoryStruct>;

using TaoCmdHistoryStructAlloc1D = FTypeAlloc1D<
    TaoCmdHistoryStructArray1D,
    allocate_tao_cmd_history_struct_container,
    deallocate_tao_cmd_history_struct_container,
    reallocate_tao_cmd_history_struct_container_data,
    access_tao_cmd_history_struct_container>;

class TaoUniverseStruct;

using TaoUniverseStructArray1D = FTypeArray1D<
    TaoUniverseStruct,
    allocate_fortran_tao_universe_struct,
    deallocate_fortran_tao_universe_struct>;
using TaoUniverseStructArray2D = FTypeArray2D<TaoUniverseStruct>;
using TaoUniverseStructArray3D = FTypeArray3D<TaoUniverseStruct>;

using TaoUniverseStructAlloc1D = FTypeAlloc1D<
    TaoUniverseStructArray1D,
    allocate_tao_universe_struct_container,
    deallocate_tao_universe_struct_container,
    reallocate_tao_universe_struct_container_data,
    access_tao_universe_struct_container>;

class MadEnergyStruct;

using MadEnergyStructArray1D = FTypeArray1D<
    MadEnergyStruct,
    allocate_fortran_mad_energy_struct,
    deallocate_fortran_mad_energy_struct>;
using MadEnergyStructArray2D = FTypeArray2D<MadEnergyStruct>;
using MadEnergyStructArray3D = FTypeArray3D<MadEnergyStruct>;

using MadEnergyStructAlloc1D = FTypeAlloc1D<
    MadEnergyStructArray1D,
    allocate_mad_energy_struct_container,
    deallocate_mad_energy_struct_container,
    reallocate_mad_energy_struct_container_data,
    access_mad_energy_struct_container>;

class MadMapStruct;

using MadMapStructArray1D =
    FTypeArray1D<MadMapStruct, allocate_fortran_mad_map_struct, deallocate_fortran_mad_map_struct>;
using MadMapStructArray2D = FTypeArray2D<MadMapStruct>;
using MadMapStructArray3D = FTypeArray3D<MadMapStruct>;

using MadMapStructAlloc1D = FTypeAlloc1D<
    MadMapStructArray1D,
    allocate_mad_map_struct_container,
    deallocate_mad_map_struct_container,
    reallocate_mad_map_struct_container_data,
    access_mad_map_struct_container>;

class RandomStateStruct;

using RandomStateStructArray1D = FTypeArray1D<
    RandomStateStruct,
    allocate_fortran_random_state_struct,
    deallocate_fortran_random_state_struct>;
using RandomStateStructArray2D = FTypeArray2D<RandomStateStruct>;
using RandomStateStructArray3D = FTypeArray3D<RandomStateStruct>;

using RandomStateStructAlloc1D = FTypeAlloc1D<
    RandomStateStructArray1D,
    allocate_random_state_struct_container,
    deallocate_random_state_struct_container,
    reallocate_random_state_struct_container_data,
    access_random_state_struct_container>;

class BbuStageStruct;

using BbuStageStructArray1D = FTypeArray1D<
    BbuStageStruct,
    allocate_fortran_bbu_stage_struct,
    deallocate_fortran_bbu_stage_struct>;
using BbuStageStructArray2D = FTypeArray2D<BbuStageStruct>;
using BbuStageStructArray3D = FTypeArray3D<BbuStageStruct>;

using BbuStageStructAlloc1D = FTypeAlloc1D<
    BbuStageStructArray1D,
    allocate_bbu_stage_struct_container,
    deallocate_bbu_stage_struct_container,
    reallocate_bbu_stage_struct_container_data,
    access_bbu_stage_struct_container>;

class BbuBeamStruct;

using BbuBeamStructArray1D = FTypeArray1D<
    BbuBeamStruct,
    allocate_fortran_bbu_beam_struct,
    deallocate_fortran_bbu_beam_struct>;
using BbuBeamStructArray2D = FTypeArray2D<BbuBeamStruct>;
using BbuBeamStructArray3D = FTypeArray3D<BbuBeamStruct>;

using BbuBeamStructAlloc1D = FTypeAlloc1D<
    BbuBeamStructArray1D,
    allocate_bbu_beam_struct_container,
    deallocate_bbu_beam_struct_container,
    reallocate_bbu_beam_struct_container_data,
    access_bbu_beam_struct_container>;

class BbuParamStruct;

using BbuParamStructArray1D = FTypeArray1D<
    BbuParamStruct,
    allocate_fortran_bbu_param_struct,
    deallocate_fortran_bbu_param_struct>;
using BbuParamStructArray2D = FTypeArray2D<BbuParamStruct>;
using BbuParamStructArray3D = FTypeArray3D<BbuParamStruct>;

using BbuParamStructAlloc1D = FTypeAlloc1D<
    BbuParamStructArray1D,
    allocate_bbu_param_struct_container,
    deallocate_bbu_param_struct_container,
    reallocate_bbu_param_struct_container_data,
    access_bbu_param_struct_container>;

class Fibre;

using FibreArray1D = FTypeArray1D<Fibre, allocate_fortran_fibre, deallocate_fortran_fibre>;
using FibreArray2D = FTypeArray2D<Fibre>;
using FibreArray3D = FTypeArray3D<Fibre>;

using FibreAlloc1D = FTypeAlloc1D<
    FibreArray1D,
    allocate_fibre_container,
    deallocate_fibre_container,
    reallocate_fibre_container_data,
    access_fibre_container>;

class Layout;

using LayoutArray1D = FTypeArray1D<Layout, allocate_fortran_layout, deallocate_fortran_layout>;
using LayoutArray2D = FTypeArray2D<Layout>;
using LayoutArray3D = FTypeArray3D<Layout>;

using LayoutAlloc1D = FTypeAlloc1D<
    LayoutArray1D,
    allocate_layout_container,
    deallocate_layout_container,
    reallocate_layout_container_data,
    access_layout_container>;

class AllEncompassingStruct;

using AllEncompassingStructArray1D = FTypeArray1D<
    AllEncompassingStruct,
    allocate_fortran_all_encompassing_struct,
    deallocate_fortran_all_encompassing_struct>;
using AllEncompassingStructArray2D = FTypeArray2D<AllEncompassingStruct>;
using AllEncompassingStructArray3D = FTypeArray3D<AllEncompassingStruct>;

using AllEncompassingStructAlloc1D = FTypeAlloc1D<
    AllEncompassingStructArray1D,
    allocate_all_encompassing_struct_container,
    deallocate_all_encompassing_struct_container,
    reallocate_all_encompassing_struct_container_data,
    access_all_encompassing_struct_container>;

class TestSubStruct;

using TestSubStructArray1D = FTypeArray1D<
    TestSubStruct,
    allocate_fortran_test_sub_struct,
    deallocate_fortran_test_sub_struct>;
using TestSubStructArray2D = FTypeArray2D<TestSubStruct>;
using TestSubStructArray3D = FTypeArray3D<TestSubStruct>;

using TestSubStructAlloc1D = FTypeAlloc1D<
    TestSubStructArray1D,
    allocate_test_sub_struct_container,
    deallocate_test_sub_struct_container,
    reallocate_test_sub_struct_container_data,
    access_test_sub_struct_container>;

class TestSubSubStruct;

using TestSubSubStructArray1D = FTypeArray1D<
    TestSubSubStruct,
    allocate_fortran_test_sub_sub_struct,
    deallocate_fortran_test_sub_sub_struct>;
using TestSubSubStructArray2D = FTypeArray2D<TestSubSubStruct>;
using TestSubSubStructArray3D = FTypeArray3D<TestSubSubStruct>;

using TestSubSubStructAlloc1D = FTypeAlloc1D<
    TestSubSubStructArray1D,
    allocate_test_sub_sub_struct_container,
    deallocate_test_sub_sub_struct_container,
    reallocate_test_sub_sub_struct_container_data,
    access_test_sub_sub_struct_container>;

template <>
struct FortranTraits<SplineStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_spline_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_spline_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_spline_struct(src, dst); }
  static constexpr std::string_view type_name() { return "spline_struct"; }
};

class SplineStruct : public FortranProxy<SplineStruct> {
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
  void set_coef(const std::vector<double> &v);
};

template <>
struct FortranTraits<SpinPolarStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_spin_polar_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_spin_polar_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_spin_polar_struct(src, dst); }
  static constexpr std::string_view type_name() { return "spin_polar_struct"; }
};

class SpinPolarStruct : public FortranProxy<SpinPolarStruct> {
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
struct FortranTraits<AcKickerTimeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_ac_kicker_time_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_ac_kicker_time_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_ac_kicker_time_struct(src, dst); }
  static constexpr std::string_view type_name() { return "ac_kicker_time_struct"; }
};

class AcKickerTimeStruct : public FortranProxy<AcKickerTimeStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double amp() const; // 0D_NOT_real
  void set_amp(double value);
  double time() const; // 0D_NOT_real
  void set_time(double value);
  SplineStruct spline() const; // 0D_NOT_type
  void set_spline(const SplineStruct &src);
};

template <>
struct FortranTraits<AcKickerFreqStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_ac_kicker_freq_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_ac_kicker_freq_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_ac_kicker_freq_struct(src, dst); }
  static constexpr std::string_view type_name() { return "ac_kicker_freq_struct"; }
};

class AcKickerFreqStruct : public FortranProxy<AcKickerFreqStruct> {
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
struct FortranTraits<AcKickerStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_ac_kicker_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_ac_kicker_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_ac_kicker_struct(src, dst); }
  static constexpr std::string_view type_name() { return "ac_kicker_struct"; }
};

class AcKickerStruct : public FortranProxy<AcKickerStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  AcKickerTimeStructArray1D amp_vs_time() const; // 1D_ALLOC_type
  AcKickerFreqStructArray1D frequency() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<Interval1CoefStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_interval1_coef_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_interval1_coef_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_interval1_coef_struct(src, dst); }
  static constexpr std::string_view type_name() { return "interval1_coef_struct"; }
};

class Interval1CoefStruct : public FortranProxy<Interval1CoefStruct> {
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
struct FortranTraits<PhotonReflectTableStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_photon_reflect_table_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_photon_reflect_table_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_photon_reflect_table_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "photon_reflect_table_struct"; }
};

class PhotonReflectTableStruct : public FortranProxy<PhotonReflectTableStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> angle() const; // 1D_ALLOC_real
  void set_angle(const std::vector<double> &v);
  FArray1D<double> energy() const; // 1D_ALLOC_real
  void set_energy(const std::vector<double> &v);
  Interval1CoefStructArray1D int1() const; // 1D_ALLOC_type
  FArray2D<double> p_reflect() const; // 2D_ALLOC_real
  void set_p_reflect(const std::vector<std::vector<double>> &v);
  double max_energy() const; // 0D_NOT_real
  void set_max_energy(double value);
  FArray1D<double> p_reflect_scratch() const; // 1D_ALLOC_real
  void set_p_reflect_scratch(const std::vector<double> &v);
  FArray1D<double> bragg_angle() const; // 1D_ALLOC_real
  void set_bragg_angle(const std::vector<double> &v);
};

template <>
struct FortranTraits<PhotonReflectSurfaceStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_photon_reflect_surface_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_photon_reflect_surface_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_photon_reflect_surface_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "photon_reflect_surface_struct"; }
};

class PhotonReflectSurfaceStruct : public FortranProxy<PhotonReflectSurfaceStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  std::string description() const; // 0D_NOT_character
  void set_description(const std::string &value);
  std::string reflectivity_file() const; // 0D_NOT_character
  void set_reflectivity_file(const std::string &value);
  PhotonReflectTableStructArray1D table() const; // 1D_ALLOC_type
  double surface_roughness_rms() const; // 0D_NOT_real
  void set_surface_roughness_rms(double value);
  double roughness_correlation_len() const; // 0D_NOT_real
  void set_roughness_correlation_len(double value);
  int ix_surface() const; // 0D_NOT_integer
  void set_ix_surface(int value);
};

template <>
struct FortranTraits<CoordStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_coord_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_coord_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_coord_struct(src, dst); }
  static constexpr std::string_view type_name() { return "coord_struct"; }
};

class CoordStruct : public FortranProxy<CoordStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> vec() const; // 1D_NOT_real
  void set_vec(const std::vector<double> &v);
  double s() const; // 0D_NOT_real
  void set_s(double value);
  long double t() const; // 0D_NOT_real16
  void set_t(long double value);
  FArray1D<double> spin() const; // 1D_NOT_real
  void set_spin(const std::vector<double> &v);
  FArray1D<double> field() const; // 1D_NOT_real
  void set_field(const std::vector<double> &v);
  FArray1D<double> phase() const; // 1D_NOT_real
  void set_phase(const std::vector<double> &v);
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
struct FortranTraits<CoordArrayStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_coord_array_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_coord_array_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_coord_array_struct(src, dst); }
  static constexpr std::string_view type_name() { return "coord_array_struct"; }
};

class CoordArrayStruct : public FortranProxy<CoordArrayStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  CoordStructArray1D orbit() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<BpmPhaseCouplingStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bpm_phase_coupling_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_bpm_phase_coupling_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_bpm_phase_coupling_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bpm_phase_coupling_struct"; }
};

class BpmPhaseCouplingStruct : public FortranProxy<BpmPhaseCouplingStruct> {
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
struct FortranTraits<ExpressionAtomStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_expression_atom_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_expression_atom_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_expression_atom_struct(src, dst); }
  static constexpr std::string_view type_name() { return "expression_atom_struct"; }
};

class ExpressionAtomStruct : public FortranProxy<ExpressionAtomStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  int type() const; // 0D_NOT_integer
  void set_type(int value);
  double value() const; // 0D_NOT_real
  void set_value(double value);
};

template <>
struct FortranTraits<WakeSrZLongStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_wake_sr_z_long_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_wake_sr_z_long_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_wake_sr_z_long_struct(src, dst); }
  static constexpr std::string_view type_name() { return "wake_sr_z_long_struct"; }
};

class WakeSrZLongStruct : public FortranProxy<WakeSrZLongStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> w() const; // 1D_ALLOC_real
  void set_w(const std::vector<double> &v);
  FArray1D<std::complex<double>> fw() const; // 1D_ALLOC_complex
  void set_fw(const std::vector<std::complex<double>> &v);
  FArray1D<std::complex<double>> fbunch() const; // 1D_ALLOC_complex
  void set_fbunch(const std::vector<std::complex<double>> &v);
  FArray1D<std::complex<double>> w_out() const; // 1D_ALLOC_complex
  void set_w_out(const std::vector<std::complex<double>> &v);
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
struct FortranTraits<WakeSrModeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_wake_sr_mode_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_wake_sr_mode_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_wake_sr_mode_struct(src, dst); }
  static constexpr std::string_view type_name() { return "wake_sr_mode_struct"; }
};

class WakeSrModeStruct : public FortranProxy<WakeSrModeStruct> {
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
struct FortranTraits<WakeSrStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_wake_sr_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_wake_sr_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_wake_sr_struct(src, dst); }
  static constexpr std::string_view type_name() { return "wake_sr_struct"; }
};

class WakeSrStruct : public FortranProxy<WakeSrStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string &value);
  WakeSrZLongStruct z_long() const; // 0D_NOT_type
  void set_z_long(const WakeSrZLongStruct &src);
  WakeSrModeStructArray1D long_wake() const; // 1D_ALLOC_type
  WakeSrModeStructArray1D trans_wake() const; // 1D_ALLOC_type
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
struct FortranTraits<WakeLrModeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_wake_lr_mode_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_wake_lr_mode_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_wake_lr_mode_struct(src, dst); }
  static constexpr std::string_view type_name() { return "wake_lr_mode_struct"; }
};

class WakeLrModeStruct : public FortranProxy<WakeLrModeStruct> {
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
struct FortranTraits<WakeLrStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_wake_lr_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_wake_lr_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_wake_lr_struct(src, dst); }
  static constexpr std::string_view type_name() { return "wake_lr_struct"; }
};

class WakeLrStruct : public FortranProxy<WakeLrStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string &value);
  WakeLrModeStructArray1D mode() const; // 1D_ALLOC_type
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
struct FortranTraits<LatEleLocStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_lat_ele_loc_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_lat_ele_loc_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_lat_ele_loc_struct(src, dst); }
  static constexpr std::string_view type_name() { return "lat_ele_loc_struct"; }
};

class LatEleLocStruct : public FortranProxy<LatEleLocStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
};

template <>
struct FortranTraits<WakeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_wake_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_wake_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_wake_struct(src, dst); }
  static constexpr std::string_view type_name() { return "wake_struct"; }
};

class WakeStruct : public FortranProxy<WakeStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  WakeSrStruct sr() const; // 0D_NOT_type
  void set_sr(const WakeSrStruct &src);
  WakeLrStruct lr() const; // 0D_NOT_type
  void set_lr(const WakeLrStruct &src);
};

template <>
struct FortranTraits<TaylorTermStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_taylor_term_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_taylor_term_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_taylor_term_struct(src, dst); }
  static constexpr std::string_view type_name() { return "taylor_term_struct"; }
};

class TaylorTermStruct : public FortranProxy<TaylorTermStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double coef() const; // 0D_NOT_real
  void set_coef(double value);
  FArray1D<int> expn() const; // 1D_NOT_integer
  void set_expn(const std::vector<int> &v);
};

template <>
struct FortranTraits<TaylorStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_taylor_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_taylor_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_taylor_struct(src, dst); }
  static constexpr std::string_view type_name() { return "taylor_struct"; }
};

class TaylorStruct : public FortranProxy<TaylorStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double ref() const; // 0D_NOT_real
  void set_ref(double value);
  TaylorTermStructArray1D term() const; // 1D_PTR_type
};

template <>
struct FortranTraits<EmTaylorTermStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_em_taylor_term_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_em_taylor_term_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_em_taylor_term_struct(src, dst); }
  static constexpr std::string_view type_name() { return "em_taylor_term_struct"; }
};

class EmTaylorTermStruct : public FortranProxy<EmTaylorTermStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double coef() const; // 0D_NOT_real
  void set_coef(double value);
  FArray1D<int> expn() const; // 1D_NOT_integer
  void set_expn(const std::vector<int> &v);
};

template <>
struct FortranTraits<EmTaylorStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_em_taylor_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_em_taylor_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_em_taylor_struct(src, dst); }
  static constexpr std::string_view type_name() { return "em_taylor_struct"; }
};

class EmTaylorStruct : public FortranProxy<EmTaylorStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double ref() const; // 0D_NOT_real
  void set_ref(double value);
  EmTaylorTermStructArray1D term() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<CartesianMapTerm1Struct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_cartesian_map_term1_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_cartesian_map_term1_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_cartesian_map_term1_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "cartesian_map_term1_struct"; }
};

class CartesianMapTerm1Struct : public FortranProxy<CartesianMapTerm1Struct> {
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
struct FortranTraits<CartesianMapTermStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_cartesian_map_term_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_cartesian_map_term_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_cartesian_map_term_struct(src, dst); }
  static constexpr std::string_view type_name() { return "cartesian_map_term_struct"; }
};

class CartesianMapTermStruct : public FortranProxy<CartesianMapTermStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string &value);
  int n_link() const; // 0D_NOT_integer
  void set_n_link(int value);
  CartesianMapTerm1StructArray1D term() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<CartesianMapStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_cartesian_map_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_cartesian_map_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_cartesian_map_struct(src, dst); }
  static constexpr std::string_view type_name() { return "cartesian_map_struct"; }
};

class CartesianMapStruct : public FortranProxy<CartesianMapStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double field_scale() const; // 0D_NOT_real
  void set_field_scale(double value);
  FArray1D<double> r0() const; // 1D_NOT_real
  void set_r0(const std::vector<double> &v);
  int master_parameter() const; // 0D_NOT_integer
  void set_master_parameter(int value);
  int ele_anchor_pt() const; // 0D_NOT_integer
  void set_ele_anchor_pt(int value);
  int field_type() const; // 0D_NOT_integer
  void set_field_type(int value);
  std::optional<CartesianMapTermStruct> ptr() const; // 0D_PTR_type
  void set_ptr(const CartesianMapTermStruct &src);
};

template <>
struct FortranTraits<CylindricalMapTerm1Struct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_cylindrical_map_term1_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_cylindrical_map_term1_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_cylindrical_map_term1_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "cylindrical_map_term1_struct"; }
};

class CylindricalMapTerm1Struct : public FortranProxy<CylindricalMapTerm1Struct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::complex<double> e_coef() const; // 0D_NOT_complex
  void set_e_coef(std::complex<double> value);
  std::complex<double> b_coef() const; // 0D_NOT_complex
  void set_b_coef(std::complex<double> value);
};

template <>
struct FortranTraits<CylindricalMapTermStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_cylindrical_map_term_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_cylindrical_map_term_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_cylindrical_map_term_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "cylindrical_map_term_struct"; }
};

class CylindricalMapTermStruct : public FortranProxy<CylindricalMapTermStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string &value);
  int n_link() const; // 0D_NOT_integer
  void set_n_link(int value);
  CylindricalMapTerm1StructArray1D term() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<CylindricalMapStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_cylindrical_map_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_cylindrical_map_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_cylindrical_map_struct(src, dst); }
  static constexpr std::string_view type_name() { return "cylindrical_map_struct"; }
};

class CylindricalMapStruct : public FortranProxy<CylindricalMapStruct> {
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
  void set_r0(const std::vector<double> &v);
  std::optional<CylindricalMapTermStruct> ptr() const; // 0D_PTR_type
  void set_ptr(const CylindricalMapTermStruct &src);
};

template <>
struct FortranTraits<BicubicCmplxCoefStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bicubic_cmplx_coef_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_bicubic_cmplx_coef_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_bicubic_cmplx_coef_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bicubic_cmplx_coef_struct"; }
};

class BicubicCmplxCoefStruct : public FortranProxy<BicubicCmplxCoefStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray2D<std::complex<double>> coef() const; // 2D_NOT_complex
  void set_coef(const std::vector<std::vector<std::complex<double>>> &v);
  FArray1D<int> i_box() const; // 1D_NOT_integer
  void set_i_box(const std::vector<int> &v);
};

template <>
struct FortranTraits<TricubicCmplxCoefStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tricubic_cmplx_coef_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tricubic_cmplx_coef_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_tricubic_cmplx_coef_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "tricubic_cmplx_coef_struct"; }
};

class TricubicCmplxCoefStruct : public FortranProxy<TricubicCmplxCoefStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray3D<std::complex<double>> coef() const; // 3D_NOT_complex
  void set_coef(const std::vector<std::vector<std::vector<std::complex<double>>>> &v);
  FArray1D<int> i_box() const; // 1D_NOT_integer
  void set_i_box(const std::vector<int> &v);
};

template <>
struct FortranTraits<GridFieldPt1Struct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_grid_field_pt1_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_grid_field_pt1_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_grid_field_pt1_struct(src, dst); }
  static constexpr std::string_view type_name() { return "grid_field_pt1_struct"; }
};

class GridFieldPt1Struct : public FortranProxy<GridFieldPt1Struct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<std::complex<double>> E() const; // 1D_NOT_complex
  void set_E(const std::vector<std::complex<double>> &v);
  FArray1D<std::complex<double>> B() const; // 1D_NOT_complex
  void set_B(const std::vector<std::complex<double>> &v);
};

template <>
struct FortranTraits<GridFieldPtStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_grid_field_pt_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_grid_field_pt_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_grid_field_pt_struct(src, dst); }
  static constexpr std::string_view type_name() { return "grid_field_pt_struct"; }
};

class GridFieldPtStruct : public FortranProxy<GridFieldPtStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string &value);
  int n_link() const; // 0D_NOT_integer
  void set_n_link(int value);
  GridFieldPt1StructArray3D pt() const; // 3D_ALLOC_type
};

template <>
struct FortranTraits<GridFieldStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_grid_field_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_grid_field_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_grid_field_struct(src, dst); }
  static constexpr std::string_view type_name() { return "grid_field_struct"; }
};

class GridFieldStruct : public FortranProxy<GridFieldStruct> {
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
  void set_dr(const std::vector<double> &v);
  FArray1D<double> r0() const; // 1D_NOT_real
  void set_r0(const std::vector<double> &v);
  bool curved_ref_frame() const; // 0D_NOT_logical
  void set_curved_ref_frame(bool value);
  std::optional<GridFieldPtStruct> ptr() const; // 0D_PTR_type
  void set_ptr(const GridFieldPtStruct &src);
  BicubicCmplxCoefStructArray3D bi_coef() const; // 3D_NOT_type
  TricubicCmplxCoefStructArray3D tri_coef() const; // 3D_NOT_type
};

template <>
struct FortranTraits<FloorPositionStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_floor_position_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_floor_position_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_floor_position_struct(src, dst); }
  static constexpr std::string_view type_name() { return "floor_position_struct"; }
};

class FloorPositionStruct : public FortranProxy<FloorPositionStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> r() const; // 1D_NOT_real
  void set_r(const std::vector<double> &v);
  FArray2D<double> w() const; // 2D_NOT_real
  void set_w(const std::vector<std::vector<double>> &v);
  double theta() const; // 0D_NOT_real
  void set_theta(double value);
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
  double psi() const; // 0D_NOT_real
  void set_psi(double value);
};

template <>
struct FortranTraits<HighEnergySpaceChargeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_high_energy_space_charge_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_high_energy_space_charge_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_high_energy_space_charge_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "high_energy_space_charge_struct"; }
};

class HighEnergySpaceChargeStruct : public FortranProxy<HighEnergySpaceChargeStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  CoordStruct closed_orb() const; // 0D_NOT_type
  void set_closed_orb(const CoordStruct &src);
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
struct FortranTraits<XyDispStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_xy_disp_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_xy_disp_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_xy_disp_struct(src, dst); }
  static constexpr std::string_view type_name() { return "xy_disp_struct"; }
};

class XyDispStruct : public FortranProxy<XyDispStruct> {
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
struct FortranTraits<TwissStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_twiss_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_twiss_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_twiss_struct(src, dst); }
  static constexpr std::string_view type_name() { return "twiss_struct"; }
};

class TwissStruct : public FortranProxy<TwissStruct> {
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
struct FortranTraits<Mode3Struct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_mode3_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_mode3_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_mode3_struct(src, dst); }
  static constexpr std::string_view type_name() { return "mode3_struct"; }
};

class Mode3Struct : public FortranProxy<Mode3Struct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray2D<double> v() const; // 2D_NOT_real
  void set_v(const std::vector<std::vector<double>> &v);
  TwissStruct a() const; // 0D_NOT_type
  void set_a(const TwissStruct &src);
  TwissStruct b() const; // 0D_NOT_type
  void set_b(const TwissStruct &src);
  TwissStruct c() const; // 0D_NOT_type
  void set_c(const TwissStruct &src);
  TwissStruct x() const; // 0D_NOT_type
  void set_x(const TwissStruct &src);
  TwissStruct y() const; // 0D_NOT_type
  void set_y(const TwissStruct &src);
};

template <>
struct FortranTraits<BookkeepingStateStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bookkeeping_state_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_bookkeeping_state_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_bookkeeping_state_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bookkeeping_state_struct"; }
};

class BookkeepingStateStruct : public FortranProxy<BookkeepingStateStruct> {
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
struct FortranTraits<RadMapStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_rad_map_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_rad_map_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_rad_map_struct(src, dst); }
  static constexpr std::string_view type_name() { return "rad_map_struct"; }
};

class RadMapStruct : public FortranProxy<RadMapStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> ref_orb() const; // 1D_NOT_real
  void set_ref_orb(const std::vector<double> &v);
  FArray2D<double> damp_dmat() const; // 2D_NOT_real
  void set_damp_dmat(const std::vector<std::vector<double>> &v);
  FArray1D<double> xfer_damp_vec() const; // 1D_NOT_real
  void set_xfer_damp_vec(const std::vector<double> &v);
  FArray2D<double> xfer_damp_mat() const; // 2D_NOT_real
  void set_xfer_damp_mat(const std::vector<std::vector<double>> &v);
  FArray2D<double> stoc_mat() const; // 2D_NOT_real
  void set_stoc_mat(const std::vector<std::vector<double>> &v);
};

template <>
struct FortranTraits<RadMapEleStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_rad_map_ele_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_rad_map_ele_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_rad_map_ele_struct(src, dst); }
  static constexpr std::string_view type_name() { return "rad_map_ele_struct"; }
};

class RadMapEleStruct : public FortranProxy<RadMapEleStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  RadMapStruct rm0() const; // 0D_NOT_type
  void set_rm0(const RadMapStruct &src);
  RadMapStruct rm1() const; // 0D_NOT_type
  void set_rm1(const RadMapStruct &src);
  bool stale() const; // 0D_NOT_logical
  void set_stale(bool value);
};

template <>
struct FortranTraits<GenGrad1Struct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_gen_grad1_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_gen_grad1_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_gen_grad1_struct(src, dst); }
  static constexpr std::string_view type_name() { return "gen_grad1_struct"; }
};

class GenGrad1Struct : public FortranProxy<GenGrad1Struct> {
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
  void set_deriv(const std::vector<std::vector<double>> &v);
};

template <>
struct FortranTraits<GenGradMapStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_gen_grad_map_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_gen_grad_map_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_gen_grad_map_struct(src, dst); }
  static constexpr std::string_view type_name() { return "gen_grad_map_struct"; }
};

class GenGradMapStruct : public FortranProxy<GenGradMapStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string file() const; // 0D_NOT_character
  void set_file(const std::string &value);
  GenGrad1StructArray1D gg() const; // 1D_ALLOC_type
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
  void set_r0(const std::vector<double> &v);
  double field_scale() const; // 0D_NOT_real
  void set_field_scale(double value);
  int master_parameter() const; // 0D_NOT_integer
  void set_master_parameter(int value);
  bool curved_ref_frame() const; // 0D_NOT_logical
  void set_curved_ref_frame(bool value);
};

template <>
struct FortranTraits<SurfaceSegmentedPtStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_surface_segmented_pt_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_surface_segmented_pt_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_surface_segmented_pt_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "surface_segmented_pt_struct"; }
};

class SurfaceSegmentedPtStruct : public FortranProxy<SurfaceSegmentedPtStruct> {
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
struct FortranTraits<SurfaceSegmentedStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_surface_segmented_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_surface_segmented_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_surface_segmented_struct(src, dst); }
  static constexpr std::string_view type_name() { return "surface_segmented_struct"; }
};

class SurfaceSegmentedStruct : public FortranProxy<SurfaceSegmentedStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool active() const; // 0D_NOT_logical
  void set_active(bool value);
  FArray1D<double> dr() const; // 1D_NOT_real
  void set_dr(const std::vector<double> &v);
  FArray1D<double> r0() const; // 1D_NOT_real
  void set_r0(const std::vector<double> &v);
  SurfaceSegmentedPtStructArray2D pt() const; // 2D_ALLOC_type
};

template <>
struct FortranTraits<SurfaceHMisalignPtStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_surface_h_misalign_pt_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_surface_h_misalign_pt_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_surface_h_misalign_pt_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "surface_h_misalign_pt_struct"; }
};

class SurfaceHMisalignPtStruct : public FortranProxy<SurfaceHMisalignPtStruct> {
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
struct FortranTraits<SurfaceHMisalignStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_surface_h_misalign_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_surface_h_misalign_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_surface_h_misalign_struct(src, dst); }
  static constexpr std::string_view type_name() { return "surface_h_misalign_struct"; }
};

class SurfaceHMisalignStruct : public FortranProxy<SurfaceHMisalignStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool active() const; // 0D_NOT_logical
  void set_active(bool value);
  FArray1D<double> dr() const; // 1D_NOT_real
  void set_dr(const std::vector<double> &v);
  FArray1D<double> r0() const; // 1D_NOT_real
  void set_r0(const std::vector<double> &v);
  SurfaceHMisalignPtStructArray2D pt() const; // 2D_ALLOC_type
};

template <>
struct FortranTraits<SurfaceDisplacementPtStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_surface_displacement_pt_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_surface_displacement_pt_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_surface_displacement_pt_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "surface_displacement_pt_struct"; }
};

class SurfaceDisplacementPtStruct : public FortranProxy<SurfaceDisplacementPtStruct> {
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
struct FortranTraits<SurfaceDisplacementStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_surface_displacement_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_surface_displacement_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_surface_displacement_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "surface_displacement_struct"; }
};

class SurfaceDisplacementStruct : public FortranProxy<SurfaceDisplacementStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool active() const; // 0D_NOT_logical
  void set_active(bool value);
  FArray1D<double> dr() const; // 1D_NOT_real
  void set_dr(const std::vector<double> &v);
  FArray1D<double> r0() const; // 1D_NOT_real
  void set_r0(const std::vector<double> &v);
  SurfaceDisplacementPtStructArray2D pt() const; // 2D_ALLOC_type
};

template <>
struct FortranTraits<TargetPointStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_target_point_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_target_point_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_target_point_struct(src, dst); }
  static constexpr std::string_view type_name() { return "target_point_struct"; }
};

class TargetPointStruct : public FortranProxy<TargetPointStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> r() const; // 1D_NOT_real
  void set_r(const std::vector<double> &v);
};

template <>
struct FortranTraits<SurfaceCurvatureStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_surface_curvature_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_surface_curvature_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_surface_curvature_struct(src, dst); }
  static constexpr std::string_view type_name() { return "surface_curvature_struct"; }
};

class SurfaceCurvatureStruct : public FortranProxy<SurfaceCurvatureStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray2D<double> xy() const; // 2D_NOT_real
  void set_xy(const std::vector<std::vector<double>> &v);
  double spherical() const; // 0D_NOT_real
  void set_spherical(double value);
  FArray1D<double> elliptical() const; // 1D_NOT_real
  void set_elliptical(const std::vector<double> &v);
  bool has_curvature() const; // 0D_NOT_logical
  void set_has_curvature(bool value);
};

template <>
struct FortranTraits<PhotonTargetStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_photon_target_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_photon_target_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_photon_target_struct(src, dst); }
  static constexpr std::string_view type_name() { return "photon_target_struct"; }
};

class PhotonTargetStruct : public FortranProxy<PhotonTargetStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int type() const; // 0D_NOT_integer
  void set_type(int value);
  int n_corner() const; // 0D_NOT_integer
  void set_n_corner(int value);
  LatEleLocStruct ele_loc() const; // 0D_NOT_type
  void set_ele_loc(const LatEleLocStruct &src);
  TargetPointStructArray1D corner() const; // 1D_NOT_type
  TargetPointStruct center() const; // 0D_NOT_type
  void set_center(const TargetPointStruct &src);
};

template <>
struct FortranTraits<PhotonMaterialStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_photon_material_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_photon_material_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_photon_material_struct(src, dst); }
  static constexpr std::string_view type_name() { return "photon_material_struct"; }
};

class PhotonMaterialStruct : public FortranProxy<PhotonMaterialStruct> {
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
  void set_h_norm(const std::vector<double> &v);
  FArray1D<double> l_ref() const; // 1D_NOT_real
  void set_l_ref(const std::vector<double> &v);
};

template <>
struct FortranTraits<PixelPtStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_pixel_pt_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_pixel_pt_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_pixel_pt_struct(src, dst); }
  static constexpr std::string_view type_name() { return "pixel_pt_struct"; }
};

class PixelPtStruct : public FortranProxy<PixelPtStruct> {
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
  void set_orbit(const std::vector<double> &v);
  FArray1D<double> orbit_rms() const; // 1D_NOT_real
  void set_orbit_rms(const std::vector<double> &v);
  FArray1D<double> init_orbit() const; // 1D_NOT_real
  void set_init_orbit(const std::vector<double> &v);
  FArray1D<double> init_orbit_rms() const; // 1D_NOT_real
  void set_init_orbit_rms(const std::vector<double> &v);
};

template <>
struct FortranTraits<PixelDetecStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_pixel_detec_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_pixel_detec_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_pixel_detec_struct(src, dst); }
  static constexpr std::string_view type_name() { return "pixel_detec_struct"; }
};

class PixelDetecStruct : public FortranProxy<PixelDetecStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> dr() const; // 1D_NOT_real
  void set_dr(const std::vector<double> &v);
  FArray1D<double> r0() const; // 1D_NOT_real
  void set_r0(const std::vector<double> &v);
  int64_t n_track_tot() const; // 0D_NOT_integer8
  void set_n_track_tot(int64_t value);
  int64_t n_hit_detec() const; // 0D_NOT_integer8
  void set_n_hit_detec(int64_t value);
  int64_t n_hit_pixel() const; // 0D_NOT_integer8
  void set_n_hit_pixel(int64_t value);
  PixelPtStructArray2D pt() const; // 2D_ALLOC_type
};

template <>
struct FortranTraits<PhotonElementStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_photon_element_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_photon_element_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_photon_element_struct(src, dst); }
  static constexpr std::string_view type_name() { return "photon_element_struct"; }
};

class PhotonElementStruct : public FortranProxy<PhotonElementStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  SurfaceCurvatureStruct curvature() const; // 0D_NOT_type
  void set_curvature(const SurfaceCurvatureStruct &src);
  PhotonTargetStruct target() const; // 0D_NOT_type
  void set_target(const PhotonTargetStruct &src);
  PhotonMaterialStruct material() const; // 0D_NOT_type
  void set_material(const PhotonMaterialStruct &src);
  SurfaceSegmentedStruct segmented() const; // 0D_NOT_type
  void set_segmented(const SurfaceSegmentedStruct &src);
  SurfaceHMisalignStruct h_misalign() const; // 0D_NOT_type
  void set_h_misalign(const SurfaceHMisalignStruct &src);
  SurfaceDisplacementStruct displacement() const; // 0D_NOT_type
  void set_displacement(const SurfaceDisplacementStruct &src);
  PixelDetecStruct pixel() const; // 0D_NOT_type
  void set_pixel(const PixelDetecStruct &src);
  int reflectivity_table_type() const; // 0D_NOT_integer
  void set_reflectivity_table_type(int value);
  PhotonReflectTableStruct reflectivity_table_sigma() const; // 0D_NOT_type
  void set_reflectivity_table_sigma(const PhotonReflectTableStruct &src);
  PhotonReflectTableStruct reflectivity_table_pi() const; // 0D_NOT_type
  void set_reflectivity_table_pi(const PhotonReflectTableStruct &src);
  SplineStructArray1D init_energy_prob() const; // 1D_ALLOC_type
  FArray1D<double> integrated_init_energy_prob() const; // 1D_ALLOC_real
  void set_integrated_init_energy_prob(const std::vector<double> &v);
};

template <>
struct FortranTraits<Wall3dVertexStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_wall3d_vertex_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_wall3d_vertex_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_wall3d_vertex_struct(src, dst); }
  static constexpr std::string_view type_name() { return "wall3d_vertex_struct"; }
};

class Wall3dVertexStruct : public FortranProxy<Wall3dVertexStruct> {
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
struct FortranTraits<Wall3dSectionStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_wall3d_section_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_wall3d_section_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_wall3d_section_struct(src, dst); }
  static constexpr std::string_view type_name() { return "wall3d_section_struct"; }
};

class Wall3dSectionStruct : public FortranProxy<Wall3dSectionStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  std::string material() const; // 0D_NOT_character
  void set_material(const std::string &value);
  Wall3dVertexStructArray1D v() const; // 1D_ALLOC_type
  std::optional<PhotonReflectSurfaceStruct> surface() const; // 0D_PTR_type
  void set_surface(const PhotonReflectSurfaceStruct &src);
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
  void set_r0(const std::vector<double> &v);
  double dx0_ds() const; // 0D_NOT_real
  void set_dx0_ds(double value);
  double dy0_ds() const; // 0D_NOT_real
  void set_dy0_ds(double value);
  FArray1D<double> x0_coef() const; // 1D_NOT_real
  void set_x0_coef(const std::vector<double> &v);
  FArray1D<double> y0_coef() const; // 1D_NOT_real
  void set_y0_coef(const std::vector<double> &v);
  double dr_ds() const; // 0D_NOT_real
  void set_dr_ds(double value);
  FArray1D<double> p1_coef() const; // 1D_NOT_real
  void set_p1_coef(const std::vector<double> &v);
  FArray1D<double> p2_coef() const; // 1D_NOT_real
  void set_p2_coef(const std::vector<double> &v);
};

template <>
struct FortranTraits<Wall3dStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_wall3d_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_wall3d_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_wall3d_struct(src, dst); }
  static constexpr std::string_view type_name() { return "wall3d_struct"; }
};

class Wall3dStruct : public FortranProxy<Wall3dStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  int type() const; // 0D_NOT_integer
  void set_type(int value);
  int ix_wall3d() const; // 0D_NOT_integer
  void set_ix_wall3d(int value);
  int n_link() const; // 0D_NOT_integer
  void set_n_link(int value);
  double thickness() const; // 0D_NOT_real
  void set_thickness(double value);
  std::string clear_material() const; // 0D_NOT_character
  void set_clear_material(const std::string &value);
  std::string opaque_material() const; // 0D_NOT_character
  void set_opaque_material(const std::string &value);
  bool superimpose() const; // 0D_NOT_logical
  void set_superimpose(bool value);
  int ele_anchor_pt() const; // 0D_NOT_integer
  void set_ele_anchor_pt(int value);
  Wall3dSectionStructArray1D section() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<RamperLordStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_ramper_lord_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_ramper_lord_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_ramper_lord_struct(src, dst); }
  static constexpr std::string_view type_name() { return "ramper_lord_struct"; }
};

class RamperLordStruct : public FortranProxy<RamperLordStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_con() const; // 0D_NOT_integer
  void set_ix_con(int value);
  double *attrib_ptr() const; // 0D_PTR_real
  void set_attrib_ptr(double value);
};

template <>
struct FortranTraits<ControlStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_control_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_control_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_control_struct(src, dst); }
  static constexpr std::string_view type_name() { return "control_struct"; }
};

class ControlStruct : public FortranProxy<ControlStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double value() const; // 0D_NOT_real
  void set_value(double value);
  FArray1D<double> y_knot() const; // 1D_ALLOC_real
  void set_y_knot(const std::vector<double> &v);
  ExpressionAtomStructArray1D stack() const; // 1D_ALLOC_type
  LatEleLocStruct slave() const; // 0D_NOT_type
  void set_slave(const LatEleLocStruct &src);
  LatEleLocStruct lord() const; // 0D_NOT_type
  void set_lord(const LatEleLocStruct &src);
  std::string slave_name() const; // 0D_NOT_character
  void set_slave_name(const std::string &value);
  std::string attribute() const; // 0D_NOT_character
  void set_attribute(const std::string &value);
  int ix_attrib() const; // 0D_NOT_integer
  void set_ix_attrib(int value);
};

template <>
struct FortranTraits<ControlVar1Struct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_control_var1_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_control_var1_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_control_var1_struct(src, dst); }
  static constexpr std::string_view type_name() { return "control_var1_struct"; }
};

class ControlVar1Struct : public FortranProxy<ControlVar1Struct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  double value() const; // 0D_NOT_real
  void set_value(double value);
  double old_value() const; // 0D_NOT_real
  void set_old_value(double value);
};

template <>
struct FortranTraits<ControlRamp1Struct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_control_ramp1_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_control_ramp1_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_control_ramp1_struct(src, dst); }
  static constexpr std::string_view type_name() { return "control_ramp1_struct"; }
};

class ControlRamp1Struct : public FortranProxy<ControlRamp1Struct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> y_knot() const; // 1D_ALLOC_real
  void set_y_knot(const std::vector<double> &v);
  ExpressionAtomStructArray1D stack() const; // 1D_ALLOC_type
  std::string attribute() const; // 0D_NOT_character
  void set_attribute(const std::string &value);
  std::string slave_name() const; // 0D_NOT_character
  void set_slave_name(const std::string &value);
  bool is_controller() const; // 0D_NOT_logical
  void set_is_controller(bool value);
};

template <>
struct FortranTraits<ControllerStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_controller_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_controller_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_controller_struct(src, dst); }
  static constexpr std::string_view type_name() { return "controller_struct"; }
};

class ControllerStruct : public FortranProxy<ControllerStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  ControlVar1StructArray1D var() const; // 1D_ALLOC_type
  ControlRamp1StructArray1D ramp() const; // 1D_ALLOC_type
  RamperLordStructArray1D ramper_lord() const; // 1D_ALLOC_type
  FArray1D<double> x_knot() const; // 1D_ALLOC_real
  void set_x_knot(const std::vector<double> &v);
};

template <>
struct FortranTraits<EllipseBeamInitStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_ellipse_beam_init_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_ellipse_beam_init_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_ellipse_beam_init_struct(src, dst); }
  static constexpr std::string_view type_name() { return "ellipse_beam_init_struct"; }
};

class EllipseBeamInitStruct : public FortranProxy<EllipseBeamInitStruct> {
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
struct FortranTraits<KvBeamInitStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_kv_beam_init_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_kv_beam_init_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_kv_beam_init_struct(src, dst); }
  static constexpr std::string_view type_name() { return "kv_beam_init_struct"; }
};

class KvBeamInitStruct : public FortranProxy<KvBeamInitStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<int> part_per_phi() const; // 1D_NOT_integer
  void set_part_per_phi(const std::vector<int> &v);
  int n_I2() const; // 0D_NOT_integer
  void set_n_I2(int value);
  double A() const; // 0D_NOT_real
  void set_A(double value);
};

template <>
struct FortranTraits<GridBeamInitStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_grid_beam_init_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_grid_beam_init_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_grid_beam_init_struct(src, dst); }
  static constexpr std::string_view type_name() { return "grid_beam_init_struct"; }
};

class GridBeamInitStruct : public FortranProxy<GridBeamInitStruct> {
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
struct FortranTraits<BeamInitStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_beam_init_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_beam_init_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_beam_init_struct(src, dst); }
  static constexpr std::string_view type_name() { return "beam_init_struct"; }
};

class BeamInitStruct : public FortranProxy<BeamInitStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string position_file() const; // 0D_NOT_character
  void set_position_file(const std::string &value);
  FCharArray1D distribution_type() const; // 1D_NOT_character
  FArray1D<double> spin() const; // 1D_NOT_real
  void set_spin(const std::vector<double> &v);
  EllipseBeamInitStructArray1D ellipse() const; // 1D_NOT_type
  KvBeamInitStruct KV() const; // 0D_NOT_type
  void set_KV(const KvBeamInitStruct &src);
  GridBeamInitStructArray1D grid() const; // 1D_NOT_type
  FArray1D<double> center_jitter() const; // 1D_NOT_real
  void set_center_jitter(const std::vector<double> &v);
  FArray1D<double> emit_jitter() const; // 1D_NOT_real
  void set_emit_jitter(const std::vector<double> &v);
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
  void set_random_engine(const std::string &value);
  std::string random_gauss_converter() const; // 0D_NOT_character
  void set_random_gauss_converter(const std::string &value);
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
  void set_center(const std::vector<double> &v);
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
  void set_species(const std::string &value);
  bool full_6D_coupling_calc() const; // 0D_NOT_logical
  void set_full_6D_coupling_calc(bool value);
  bool use_particle_start() const; // 0D_NOT_logical
  void set_use_particle_start(bool value);
  bool use_t_coords() const; // 0D_NOT_logical
  void set_use_t_coords(bool value);
  bool use_z_as_t() const; // 0D_NOT_logical
  void set_use_z_as_t(bool value);
  std::string file_name() const; // 0D_NOT_character
  void set_file_name(const std::string &value);
};

template <>
struct FortranTraits<LatParamStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_lat_param_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_lat_param_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_lat_param_struct(src, dst); }
  static constexpr std::string_view type_name() { return "lat_param_struct"; }
};

class LatParamStruct : public FortranProxy<LatParamStruct> {
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
  void set_t1_with_RF(const std::vector<std::vector<double>> &v);
  FArray2D<double> t1_no_RF() const; // 2D_NOT_real
  void set_t1_no_RF(const std::vector<std::vector<double>> &v);
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
  BookkeepingStateStruct bookkeeping_state() const; // 0D_NOT_type
  void set_bookkeeping_state(const BookkeepingStateStruct &src);
  BeamInitStruct beam_init() const; // 0D_NOT_type
  void set_beam_init(const BeamInitStruct &src);
};

template <>
struct FortranTraits<ModeInfoStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_mode_info_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_mode_info_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_mode_info_struct(src, dst); }
  static constexpr std::string_view type_name() { return "mode_info_struct"; }
};

class ModeInfoStruct : public FortranProxy<ModeInfoStruct> {
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
struct FortranTraits<PreTrackerStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_pre_tracker_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_pre_tracker_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_pre_tracker_struct(src, dst); }
  static constexpr std::string_view type_name() { return "pre_tracker_struct"; }
};

class PreTrackerStruct : public FortranProxy<PreTrackerStruct> {
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
  void set_input_file(const std::string &value);
};

template <>
struct FortranTraits<AnormalModeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_anormal_mode_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_anormal_mode_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_anormal_mode_struct(src, dst); }
  static constexpr std::string_view type_name() { return "anormal_mode_struct"; }
};

class AnormalModeStruct : public FortranProxy<AnormalModeStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double emittance() const; // 0D_NOT_real
  void set_emittance(double value);
  double emittance_no_vert() const; // 0D_NOT_real
  void set_emittance_no_vert(double value);
  FArray1D<double> synch_int() const; // 1D_NOT_real
  void set_synch_int(const std::vector<double> &v);
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
struct FortranTraits<LinacNormalModeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_linac_normal_mode_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_linac_normal_mode_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_linac_normal_mode_struct(src, dst); }
  static constexpr std::string_view type_name() { return "linac_normal_mode_struct"; }
};

class LinacNormalModeStruct : public FortranProxy<LinacNormalModeStruct> {
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
struct FortranTraits<NormalModesStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_normal_modes_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_normal_modes_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_normal_modes_struct(src, dst); }
  static constexpr std::string_view type_name() { return "normal_modes_struct"; }
};

class NormalModesStruct : public FortranProxy<NormalModesStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> synch_int() const; // 1D_NOT_real
  void set_synch_int(const std::vector<double> &v);
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
  AnormalModeStruct a() const; // 0D_NOT_type
  void set_a(const AnormalModeStruct &src);
  AnormalModeStruct b() const; // 0D_NOT_type
  void set_b(const AnormalModeStruct &src);
  AnormalModeStruct z() const; // 0D_NOT_type
  void set_z(const AnormalModeStruct &src);
  LinacNormalModeStruct lin() const; // 0D_NOT_type
  void set_lin(const LinacNormalModeStruct &src);
};

template <>
struct FortranTraits<EmFieldStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_em_field_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_em_field_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_em_field_struct(src, dst); }
  static constexpr std::string_view type_name() { return "em_field_struct"; }
};

class EmFieldStruct : public FortranProxy<EmFieldStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> E() const; // 1D_NOT_real
  void set_E(const std::vector<double> &v);
  FArray1D<double> B() const; // 1D_NOT_real
  void set_B(const std::vector<double> &v);
  FArray2D<double> dE() const; // 2D_NOT_real
  void set_dE(const std::vector<std::vector<double>> &v);
  FArray2D<double> dB() const; // 2D_NOT_real
  void set_dB(const std::vector<std::vector<double>> &v);
  double phi() const; // 0D_NOT_real
  void set_phi(double value);
  double phi_B() const; // 0D_NOT_real
  void set_phi_B(double value);
  FArray1D<double> A() const; // 1D_NOT_real
  void set_A(const std::vector<double> &v);
};

template <>
struct FortranTraits<StrongBeamStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_strong_beam_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_strong_beam_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_strong_beam_struct(src, dst); }
  static constexpr std::string_view type_name() { return "strong_beam_struct"; }
};

class StrongBeamStruct : public FortranProxy<StrongBeamStruct> {
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
struct FortranTraits<TrackPointStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_track_point_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_track_point_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_track_point_struct(src, dst); }
  static constexpr std::string_view type_name() { return "track_point_struct"; }
};

class TrackPointStruct : public FortranProxy<TrackPointStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double s_lab() const; // 0D_NOT_real
  void set_s_lab(double value);
  double s_body() const; // 0D_NOT_real
  void set_s_body(double value);
  CoordStruct orb() const; // 0D_NOT_type
  void set_orb(const CoordStruct &src);
  EmFieldStruct field() const; // 0D_NOT_type
  void set_field(const EmFieldStruct &src);
  StrongBeamStruct strong_beam() const; // 0D_NOT_type
  void set_strong_beam(const StrongBeamStruct &src);
  FArray1D<double> vec0() const; // 1D_NOT_real
  void set_vec0(const std::vector<double> &v);
  FArray2D<double> mat6() const; // 2D_NOT_real
  void set_mat6(const std::vector<std::vector<double>> &v);
};

template <>
struct FortranTraits<TrackStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_track_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_track_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_track_struct(src, dst); }
  static constexpr std::string_view type_name() { return "track_struct"; }
};

class TrackStruct : public FortranProxy<TrackStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TrackPointStructArray1D pt() const; // 1D_ALLOC_type
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
struct FortranTraits<SpaceChargeCommonStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_space_charge_common_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_space_charge_common_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_space_charge_common_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "space_charge_common_struct"; }
};

class SpaceChargeCommonStruct : public FortranProxy<SpaceChargeCommonStruct> {
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
  void set_space_charge_mesh_size(const std::vector<int> &v);
  FArray1D<int> csr3d_mesh_size() const; // 1D_NOT_integer
  void set_csr3d_mesh_size(const std::vector<int> &v);
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
  void set_diagnostic_output_file(const std::string &value);
};

template <>
struct FortranTraits<BmadCommonStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bmad_common_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_bmad_common_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_bmad_common_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bmad_common_struct"; }
};

class BmadCommonStruct : public FortranProxy<BmadCommonStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double max_aperture_limit() const; // 0D_NOT_real
  void set_max_aperture_limit(double value);
  FArray1D<double> d_orb() const; // 1D_NOT_real
  void set_d_orb(const std::vector<double> &v);
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
struct FortranTraits<RadInt1Struct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_rad_int1_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_rad_int1_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_rad_int1_struct(src, dst); }
  static constexpr std::string_view type_name() { return "rad_int1_struct"; }
};

class RadInt1Struct : public FortranProxy<RadInt1Struct> {
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
struct FortranTraits<RadIntBranchStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_rad_int_branch_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_rad_int_branch_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_rad_int_branch_struct(src, dst); }
  static constexpr std::string_view type_name() { return "rad_int_branch_struct"; }
};

class RadIntBranchStruct : public FortranProxy<RadIntBranchStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  RadInt1StructArray1D ele() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<RadIntAllEleStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_rad_int_all_ele_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_rad_int_all_ele_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_rad_int_all_ele_struct(src, dst); }
  static constexpr std::string_view type_name() { return "rad_int_all_ele_struct"; }
};

class RadIntAllEleStruct : public FortranProxy<RadIntAllEleStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  RadIntBranchStructArray1D branch() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<RfStairStepStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_rf_stair_step_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_rf_stair_step_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_rf_stair_step_struct(src, dst); }
  static constexpr std::string_view type_name() { return "rf_stair_step_struct"; }
};

class RfStairStepStruct : public FortranProxy<RfStairStepStruct> {
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
struct FortranTraits<RfEleStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_rf_ele_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_rf_ele_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_rf_ele_struct(src, dst); }
  static constexpr std::string_view type_name() { return "rf_ele_struct"; }
};

class RfEleStruct : public FortranProxy<RfEleStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  RfStairStepStructArray1D steps() const; // 1D_ALLOC_type
  double ds_step() const; // 0D_NOT_real
  void set_ds_step(double value);
};

template <>
struct FortranTraits<EleStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_ele_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_ele_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_ele_struct(src, dst); }
  static constexpr std::string_view type_name() { return "ele_struct"; }
};

class EleStruct : public FortranProxy<EleStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  std::string type() const; // 0D_NOT_character
  void set_type(const std::string &value);
  std::string alias() const; // 0D_NOT_character
  void set_alias(const std::string &value);
  std::string component_name() const; // 0D_NOT_character
  void set_component_name(const std::string &value);
  std::string descrip() const; // 0D_PTR_character
  void set_descrip(const std::string &value);
  TwissStruct a() const; // 0D_NOT_type
  void set_a(const TwissStruct &src);
  TwissStruct b() const; // 0D_NOT_type
  void set_b(const TwissStruct &src);
  TwissStruct z() const; // 0D_NOT_type
  void set_z(const TwissStruct &src);
  XyDispStruct x() const; // 0D_NOT_type
  void set_x(const XyDispStruct &src);
  XyDispStruct y() const; // 0D_NOT_type
  void set_y(const XyDispStruct &src);
  std::optional<AcKickerStruct> ac_kick() const; // 0D_PTR_type
  void set_ac_kick(const AcKickerStruct &src);
  BookkeepingStateStruct bookkeeping_state() const; // 0D_NOT_type
  void set_bookkeeping_state(const BookkeepingStateStruct &src);
  std::optional<BranchStruct> branch() const; // 0D_PTR_type
  void set_branch(const BranchStruct &src);
  std::optional<ControllerStruct> control() const; // 0D_PTR_type
  void set_control(const ControllerStruct &src);
  std::optional<RfEleStruct> rf() const; // 0D_PTR_type
  void set_rf(const RfEleStruct &src);
  std::optional<EleStruct> lord() const; // 0D_PTR_type
  void set_lord(const EleStruct &src);
  std::optional<Fibre> ptc_fibre() const; // 0D_PTR_type
  void set_ptc_fibre(const Fibre &src);
  FloorPositionStruct floor() const; // 0D_NOT_type
  void set_floor(const FloorPositionStruct &src);
  std::optional<HighEnergySpaceChargeStruct> high_energy_space_charge() const; // 0D_PTR_type
  void set_high_energy_space_charge(const HighEnergySpaceChargeStruct &src);
  std::optional<Mode3Struct> mode3() const; // 0D_PTR_type
  void set_mode3(const Mode3Struct &src);
  std::optional<PhotonElementStruct> photon() const; // 0D_PTR_type
  void set_photon(const PhotonElementStruct &src);
  std::optional<RadMapEleStruct> rad_map() const; // 0D_PTR_type
  void set_rad_map(const RadMapEleStruct &src);
  TaylorStructArray1D taylor() const; // 1D_NOT_type
  FArray1D<double> spin_taylor_ref_orb_in() const; // 1D_NOT_real
  void set_spin_taylor_ref_orb_in(const std::vector<double> &v);
  TaylorStructArray1D spin_taylor() const; // 1D_NOT_type
  std::optional<WakeStruct> wake() const; // 0D_PTR_type
  void set_wake(const WakeStruct &src);
  Wall3dStructArray1D wall3d() const; // 1D_PTR_type
  CartesianMapStructArray1D cartesian_map() const; // 1D_PTR_type
  CylindricalMapStructArray1D cylindrical_map() const; // 1D_PTR_type
  GenGradMapStructArray1D gen_grad_map() const; // 1D_PTR_type
  GridFieldStructArray1D grid_field() const; // 1D_PTR_type
  CoordStruct map_ref_orb_in() const; // 0D_NOT_type
  void set_map_ref_orb_in(const CoordStruct &src);
  CoordStruct map_ref_orb_out() const; // 0D_NOT_type
  void set_map_ref_orb_out(const CoordStruct &src);
  CoordStruct time_ref_orb_in() const; // 0D_NOT_type
  void set_time_ref_orb_in(const CoordStruct &src);
  CoordStruct time_ref_orb_out() const; // 0D_NOT_type
  void set_time_ref_orb_out(const CoordStruct &src);
  FArray1D<double> value() const; // 1D_NOT_real
  void set_value(const std::vector<double> &v);
  FArray1D<double> old_value() const; // 1D_NOT_real
  void set_old_value(const std::vector<double> &v);
  FArray2D<double> spin_q() const; // 2D_NOT_real
  void set_spin_q(const std::vector<std::vector<double>> &v);
  FArray1D<double> vec0() const; // 1D_NOT_real
  void set_vec0(const std::vector<double> &v);
  FArray2D<double> mat6() const; // 2D_NOT_real
  void set_mat6(const std::vector<std::vector<double>> &v);
  FArray2D<double> c_mat() const; // 2D_NOT_real
  void set_c_mat(const std::vector<std::vector<double>> &v);
  FArray2D<double> dc_mat_dpz() const; // 2D_NOT_real
  void set_dc_mat_dpz(const std::vector<std::vector<double>> &v);
  double gamma_c() const; // 0D_NOT_real
  void set_gamma_c(double value);
  double s_start() const; // 0D_NOT_real
  void set_s_start(double value);
  double s() const; // 0D_NOT_real
  void set_s(double value);
  double ref_time() const; // 0D_NOT_real
  void set_ref_time(double value);
  FArray1D<double> a_pole() const; // 1D_PTR_real
  void set_a_pole(const std::vector<double> &v);
  FArray1D<double> b_pole() const; // 1D_PTR_real
  void set_b_pole(const std::vector<double> &v);
  FArray1D<double> a_pole_elec() const; // 1D_PTR_real
  void set_a_pole_elec(const std::vector<double> &v);
  FArray1D<double> b_pole_elec() const; // 1D_PTR_real
  void set_b_pole_elec(const std::vector<double> &v);
  FArray1D<double> custom() const; // 1D_PTR_real
  void set_custom(const std::vector<double> &v);
  FArray3D<double> r() const; // 3D_PTR_real
  void set_r(const std::vector<std::vector<std::vector<double>>> &v);
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
struct FortranTraits<ComplexTaylorTermStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_complex_taylor_term_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_complex_taylor_term_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_complex_taylor_term_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "complex_taylor_term_struct"; }
};

class ComplexTaylorTermStruct : public FortranProxy<ComplexTaylorTermStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::complex<double> coef() const; // 0D_NOT_complex
  void set_coef(std::complex<double> value);
  FArray1D<int> expn() const; // 1D_NOT_integer
  void set_expn(const std::vector<int> &v);
};

template <>
struct FortranTraits<ComplexTaylorStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_complex_taylor_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_complex_taylor_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_complex_taylor_struct(src, dst); }
  static constexpr std::string_view type_name() { return "complex_taylor_struct"; }
};

class ComplexTaylorStruct : public FortranProxy<ComplexTaylorStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::complex<double> ref() const; // 0D_NOT_complex
  void set_ref(std::complex<double> value);
  ComplexTaylorTermStructArray1D term() const; // 1D_PTR_type
};

template <>
struct FortranTraits<BranchStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_branch_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_branch_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_branch_struct(src, dst); }
  static constexpr std::string_view type_name() { return "branch_struct"; }
};

class BranchStruct : public FortranProxy<BranchStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
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
  std::optional<LatStruct> lat() const; // 0D_PTR_type
  void set_lat(const LatStruct &src);
  ModeInfoStruct a() const; // 0D_NOT_type
  void set_a(const ModeInfoStruct &src);
  ModeInfoStruct b() const; // 0D_NOT_type
  void set_b(const ModeInfoStruct &src);
  ModeInfoStruct z() const; // 0D_NOT_type
  void set_z(const ModeInfoStruct &src);
  EleStructArray1D ele() const; // 1D_PTR_type
  LatParamStruct param() const; // 0D_NOT_type
  void set_param(const LatParamStruct &src);
  CoordStruct particle_start() const; // 0D_NOT_type
  void set_particle_start(const CoordStruct &src);
  Wall3dStructArray1D wall3d() const; // 1D_PTR_type
};

template <>
struct FortranTraits<LatStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_lat_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_lat_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_lat_struct(src, dst); }
  static constexpr std::string_view type_name() { return "lat_struct"; }
};

class LatStruct : public FortranProxy<LatStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string use_name() const; // 0D_NOT_character
  void set_use_name(const std::string &value);
  std::string lattice() const; // 0D_NOT_character
  void set_lattice(const std::string &value);
  std::string machine() const; // 0D_NOT_character
  void set_machine(const std::string &value);
  std::string input_file_name() const; // 0D_NOT_character
  void set_input_file_name(const std::string &value);
  std::string title() const; // 0D_NOT_character
  void set_title(const std::string &value);
  FCharArray1D print_str() const; // 1D_ALLOC_character
  ExpressionAtomStructArray1D constant() const; // 1D_ALLOC_type
  std::optional<ModeInfoStruct> a() const; // 0D_PTR_type
  void set_a(const ModeInfoStruct &src);
  std::optional<ModeInfoStruct> b() const; // 0D_PTR_type
  void set_b(const ModeInfoStruct &src);
  std::optional<ModeInfoStruct> z() const; // 0D_PTR_type
  void set_z(const ModeInfoStruct &src);
  std::optional<LatParamStruct> param() const; // 0D_PTR_type
  void set_param(const LatParamStruct &src);
  BookkeepingStateStruct lord_state() const; // 0D_NOT_type
  void set_lord_state(const BookkeepingStateStruct &src);
  EleStruct ele_init() const; // 0D_NOT_type
  void set_ele_init(const EleStruct &src);
  EleStructArray1D ele() const; // 1D_PTR_type
  BranchStructArray1D branch() const; // 1D_ALLOC_type
  ControlStructArray1D control() const; // 1D_ALLOC_type
  std::optional<CoordStruct> particle_start() const; // 0D_PTR_type
  void set_particle_start(const CoordStruct &src);
  BeamInitStruct beam_init() const; // 0D_NOT_type
  void set_beam_init(const BeamInitStruct &src);
  PreTrackerStruct pre_tracker() const; // 0D_NOT_type
  void set_pre_tracker(const PreTrackerStruct &src);
  FArray1D<double> custom() const; // 1D_ALLOC_real
  void set_custom(const std::vector<double> &v);
  int version() const; // 0D_NOT_integer
  void set_version(int value);
  int *n_ele_track() const; // 0D_PTR_integer
  void set_n_ele_track(int value);
  int *n_ele_max() const; // 0D_PTR_integer
  void set_n_ele_max(int value);
  int n_control_max() const; // 0D_NOT_integer
  void set_n_control_max(int value);
  int n_ic_max() const; // 0D_NOT_integer
  void set_n_ic_max(int value);
  int input_taylor_order() const; // 0D_NOT_integer
  void set_input_taylor_order(int value);
  FArray1D<int> ic() const; // 1D_ALLOC_integer
  void set_ic(const std::vector<int> &v);
  int photon_type() const; // 0D_NOT_integer
  void set_photon_type(int value);
  int creation_hash() const; // 0D_NOT_integer
  void set_creation_hash(int value);
  int ramper_slave_bookkeeping() const; // 0D_NOT_integer
  void set_ramper_slave_bookkeeping(int value);
};

template <>
struct FortranTraits<BunchStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bunch_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_bunch_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_bunch_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bunch_struct"; }
};

class BunchStruct : public FortranProxy<BunchStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  CoordStructArray1D particle() const; // 1D_ALLOC_type
  FArray1D<int> ix_z() const; // 1D_ALLOC_integer
  void set_ix_z(const std::vector<int> &v);
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
struct FortranTraits<BunchParamsStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bunch_params_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_bunch_params_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_bunch_params_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bunch_params_struct"; }
};

class BunchParamsStruct : public FortranProxy<BunchParamsStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  CoordStruct centroid() const; // 0D_NOT_type
  void set_centroid(const CoordStruct &src);
  TwissStruct x() const; // 0D_NOT_type
  void set_x(const TwissStruct &src);
  TwissStruct y() const; // 0D_NOT_type
  void set_y(const TwissStruct &src);
  TwissStruct z() const; // 0D_NOT_type
  void set_z(const TwissStruct &src);
  TwissStruct a() const; // 0D_NOT_type
  void set_a(const TwissStruct &src);
  TwissStruct b() const; // 0D_NOT_type
  void set_b(const TwissStruct &src);
  TwissStruct c() const; // 0D_NOT_type
  void set_c(const TwissStruct &src);
  FArray2D<double> sigma() const; // 2D_NOT_real
  void set_sigma(const std::vector<std::vector<double>> &v);
  FArray1D<double> rel_max() const; // 1D_NOT_real
  void set_rel_max(const std::vector<double> &v);
  FArray1D<double> rel_min() const; // 1D_NOT_real
  void set_rel_min(const std::vector<double> &v);
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
struct FortranTraits<BeamStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_beam_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_beam_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_beam_struct(src, dst); }
  static constexpr std::string_view type_name() { return "beam_struct"; }
};

class BeamStruct : public FortranProxy<BeamStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  BunchStructArray1D bunch() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<AperturePointStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_aperture_point_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_aperture_point_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_aperture_point_struct(src, dst); }
  static constexpr std::string_view type_name() { return "aperture_point_struct"; }
};

class AperturePointStruct : public FortranProxy<AperturePointStruct> {
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
struct FortranTraits<ApertureParamStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_aperture_param_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_aperture_param_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_aperture_param_struct(src, dst); }
  static constexpr std::string_view type_name() { return "aperture_param_struct"; }
};

class ApertureParamStruct : public FortranProxy<ApertureParamStruct> {
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
  void set_start_ele(const std::string &value);
};

template <>
struct FortranTraits<ApertureScanStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_aperture_scan_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_aperture_scan_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_aperture_scan_struct(src, dst); }
  static constexpr std::string_view type_name() { return "aperture_scan_struct"; }
};

class ApertureScanStruct : public FortranProxy<ApertureScanStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  AperturePointStructArray1D point() const; // 1D_ALLOC_type
  CoordStruct ref_orb() const; // 0D_NOT_type
  void set_ref_orb(const CoordStruct &src);
  double pz_start() const; // 0D_NOT_real
  void set_pz_start(double value);
};

template <>
struct FortranTraits<ElePointerStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_ele_pointer_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_ele_pointer_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_ele_pointer_struct(src, dst); }
  static constexpr std::string_view type_name() { return "ele_pointer_struct"; }
};

class ElePointerStruct : public FortranProxy<ElePointerStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<EleStruct> ele() const; // 0D_PTR_type
  void set_ele(const EleStruct &src);
  LatEleLocStruct loc() const; // 0D_NOT_type
  void set_loc(const LatEleLocStruct &src);
  int id() const; // 0D_NOT_integer
  void set_id(int value);
};

template <>
struct FortranTraits<ExpressionTreeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_expression_tree_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_expression_tree_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_expression_tree_struct(src, dst); }
  static constexpr std::string_view type_name() { return "expression_tree_struct"; }
};

class ExpressionTreeStruct : public FortranProxy<ExpressionTreeStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  int type() const; // 0D_NOT_integer
  void set_type(int value);
  double value() const; // 0D_NOT_real
  void set_value(double value);
  ExpressionTreeStructArray1D node() const; // 1D_PTR_type
};

template <>
struct FortranTraits<NametableStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_nametable_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_nametable_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_nametable_struct(src, dst); }
  static constexpr std::string_view type_name() { return "nametable_struct"; }
};

class NametableStruct : public FortranProxy<NametableStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FCharArray1D name() const; // 1D_ALLOC_character
  FArray1D<int> index() const; // 1D_ALLOC_integer
  void set_index(const std::vector<int> &v);
  int n_min() const; // 0D_NOT_integer
  void set_n_min(int value);
  int n_max() const; // 0D_NOT_integer
  void set_n_max(int value);
};

template <>
struct FortranTraits<TaoSpinDnDpzStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_spin_dn_dpz_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_spin_dn_dpz_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_spin_dn_dpz_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_spin_dn_dpz_struct"; }
};

class TaoSpinDnDpzStruct : public FortranProxy<TaoSpinDnDpzStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> vec() const; // 1D_NOT_real
  void set_vec(const std::vector<double> &v);
  FArray2D<double> partial() const; // 2D_NOT_real
  void set_partial(const std::vector<std::vector<double>> &v);
  FArray2D<double> partial2() const; // 2D_NOT_real
  void set_partial2(const std::vector<std::vector<double>> &v);
};

template <>
struct FortranTraits<ResonanceHStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_resonance_h_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_resonance_h_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_resonance_h_struct(src, dst); }
  static constexpr std::string_view type_name() { return "resonance_h_struct"; }
};

class ResonanceHStruct : public FortranProxy<ResonanceHStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string id() const; // 0D_NOT_character
  void set_id(const std::string &value);
  std::complex<double> c_val() const; // 0D_NOT_complex
  void set_c_val(std::complex<double> value);
};

template <>
struct FortranTraits<SpinOrbitMap1Struct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_spin_orbit_map1_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_spin_orbit_map1_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_spin_orbit_map1_struct(src, dst); }
  static constexpr std::string_view type_name() { return "spin_orbit_map1_struct"; }
};

class SpinOrbitMap1Struct : public FortranProxy<SpinOrbitMap1Struct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray2D<double> orb_mat() const; // 2D_NOT_real
  void set_orb_mat(const std::vector<std::vector<double>> &v);
  FArray1D<double> vec0() const; // 1D_NOT_real
  void set_vec0(const std::vector<double> &v);
  FArray2D<double> spin_q() const; // 2D_NOT_real
  void set_spin_q(const std::vector<std::vector<double>> &v);
};

template <>
struct FortranTraits<SpinAxisStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_spin_axis_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_spin_axis_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_spin_axis_struct(src, dst); }
  static constexpr std::string_view type_name() { return "spin_axis_struct"; }
};

class SpinAxisStruct : public FortranProxy<SpinAxisStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> l() const; // 1D_NOT_real
  void set_l(const std::vector<double> &v);
  FArray1D<double> n0() const; // 1D_NOT_real
  void set_n0(const std::vector<double> &v);
  FArray1D<double> m() const; // 1D_NOT_real
  void set_m(const std::vector<double> &v);
};

template <>
struct FortranTraits<PtcNormalFormStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_ptc_normal_form_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_ptc_normal_form_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_ptc_normal_form_struct(src, dst); }
  static constexpr std::string_view type_name() { return "ptc_normal_form_struct"; }
};

class PtcNormalFormStruct : public FortranProxy<PtcNormalFormStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<EleStruct> ele_origin() const; // 0D_PTR_type
  void set_ele_origin(const EleStruct &src);
  FArray1D<double> orb0() const; // 1D_NOT_real
  void set_orb0(const std::vector<double> &v);
  bool valid_map() const; // 0D_NOT_logical
  void set_valid_map(bool value);
};

template <>
struct FortranTraits<BmadNormalFormStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bmad_normal_form_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_bmad_normal_form_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_bmad_normal_form_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bmad_normal_form_struct"; }
};

class BmadNormalFormStruct : public FortranProxy<BmadNormalFormStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<EleStruct> ele_origin() const; // 0D_PTR_type
  void set_ele_origin(const EleStruct &src);
  TaylorStructArray1D M() const; // 1D_NOT_type
  TaylorStructArray1D A() const; // 1D_NOT_type
  TaylorStructArray1D A_inv() const; // 1D_NOT_type
  TaylorStructArray1D dhdj() const; // 1D_NOT_type
  ComplexTaylorStructArray1D F() const; // 1D_NOT_type
  ComplexTaylorStructArray1D L() const; // 1D_NOT_type
  ResonanceHStructArray1D h() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<BunchTrackStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bunch_track_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_bunch_track_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_bunch_track_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bunch_track_struct"; }
};

class BunchTrackStruct : public FortranProxy<BunchTrackStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  BunchParamsStructArray1D pt() const; // 1D_ALLOC_type
  double ds_save() const; // 0D_NOT_real
  void set_ds_save(double value);
  int n_pt() const; // 0D_NOT_integer
  void set_n_pt(int value);
};

template <>
struct FortranTraits<SummationRdtStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_summation_rdt_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_summation_rdt_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_summation_rdt_struct(src, dst); }
  static constexpr std::string_view type_name() { return "summation_rdt_struct"; }
};

class SummationRdtStruct : public FortranProxy<SummationRdtStruct> {
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
struct FortranTraits<TaoEleShapeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_ele_shape_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_ele_shape_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_ele_shape_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_ele_shape_struct"; }
};

class TaoEleShapeStruct : public FortranProxy<TaoEleShapeStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string ele_id() const; // 0D_NOT_character
  void set_ele_id(const std::string &value);
  std::string shape() const; // 0D_NOT_character
  void set_shape(const std::string &value);
  std::string color() const; // 0D_NOT_character
  void set_color(const std::string &value);
  double size() const; // 0D_NOT_real
  void set_size(double value);
  std::string label() const; // 0D_NOT_character
  void set_label(const std::string &value);
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
  void set_name_ele(const std::string &value);
  TaoElePointerStructArray1D uni() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoElePointerStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_ele_pointer_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_ele_pointer_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_ele_pointer_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_ele_pointer_struct"; }
};

class TaoElePointerStruct : public FortranProxy<TaoElePointerStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  ElePointerStructArray1D eles() const; // 1D_ALLOC_type
  int n_loc() const; // 0D_NOT_integer
  void set_n_loc(int value);
};

template <>
struct FortranTraits<TaoCurveStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_curve_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_curve_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_curve_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_curve_struct"; }
};

class TaoCurveStruct : public FortranProxy<TaoCurveStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  std::string data_source() const; // 0D_NOT_character
  void set_data_source(const std::string &value);
  std::string data_index() const; // 0D_NOT_character
  void set_data_index(const std::string &value);
  std::string data_type_x() const; // 0D_NOT_character
  void set_data_type_x(const std::string &value);
  std::string data_type() const; // 0D_ALLOC_character
  void set_data_type(const std::string &value);
  std::string ele_ref_name() const; // 0D_NOT_character
  void set_ele_ref_name(const std::string &value);
  std::string legend_text() const; // 0D_NOT_character
  void set_legend_text(const std::string &value);
  std::string message_text() const; // 0D_NOT_character
  void set_message_text(const std::string &value);
  std::string component() const; // 0D_NOT_character
  void set_component(const std::string &value);
  std::string why_invalid() const; // 0D_NOT_character
  void set_why_invalid(const std::string &value);
  std::optional<TaoGraphStruct> g() const; // 0D_PTR_type
  void set_g(const TaoGraphStruct &src);
  TaoHistogramStruct hist() const; // 0D_NOT_type
  void set_hist(const TaoHistogramStruct &src);
  TaoCurveColorStruct z_color() const; // 0D_NOT_type
  void set_z_color(const TaoCurveColorStruct &src);
  FArray1D<double> x_line() const; // 1D_ALLOC_real
  void set_x_line(const std::vector<double> &v);
  FArray1D<double> y_line() const; // 1D_ALLOC_real
  void set_y_line(const std::vector<double> &v);
  FArray1D<double> y2_line() const; // 1D_ALLOC_real
  void set_y2_line(const std::vector<double> &v);
  FArray1D<int> ix_line() const; // 1D_ALLOC_integer
  void set_ix_line(const std::vector<int> &v);
  FArray1D<double> x_symb() const; // 1D_ALLOC_real
  void set_x_symb(const std::vector<double> &v);
  FArray1D<double> y_symb() const; // 1D_ALLOC_real
  void set_y_symb(const std::vector<double> &v);
  FArray1D<double> z_symb() const; // 1D_ALLOC_real
  void set_z_symb(const std::vector<double> &v);
  FArray1D<double> err_symb() const; // 1D_ALLOC_real
  void set_err_symb(const std::vector<double> &v);
  FArray1D<double> symb_size() const; // 1D_ALLOC_real
  void set_symb_size(const std::vector<double> &v);
  FArray1D<int> ix_symb() const; // 1D_ALLOC_integer
  void set_ix_symb(const std::vector<int> &v);
  double y_axis_scale_factor() const; // 0D_NOT_real
  void set_y_axis_scale_factor(double value);
  QpLineStruct line() const; // 0D_NOT_type
  void set_line(const QpLineStruct &src);
  QpSymbolStruct symbol() const; // 0D_NOT_type
  void set_symbol(const QpSymbolStruct &src);
  TaoCurveOrbitStruct orbit() const; // 0D_NOT_type
  void set_orbit(const TaoCurveOrbitStruct &src);
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
struct FortranTraits<TaoCurveColorStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_curve_color_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_curve_color_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_curve_color_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_curve_color_struct"; }
};

class TaoCurveColorStruct : public FortranProxy<TaoCurveColorStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string data_type() const; // 0D_NOT_character
  void set_data_type(const std::string &value);
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
struct FortranTraits<TaoCurveOrbitStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_curve_orbit_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_curve_orbit_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_curve_orbit_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_curve_orbit_struct"; }
};

class TaoCurveOrbitStruct : public FortranProxy<TaoCurveOrbitStruct> {
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
struct FortranTraits<TaoHistogramStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_histogram_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_histogram_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_histogram_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_histogram_struct"; }
};

class TaoHistogramStruct : public FortranProxy<TaoHistogramStruct> {
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
struct FortranTraits<LatEleOrder1Struct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_lat_ele_order1_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_lat_ele_order1_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_lat_ele_order1_struct(src, dst); }
  static constexpr std::string_view type_name() { return "lat_ele_order1_struct"; }
};

class LatEleOrder1Struct : public FortranProxy<LatEleOrder1Struct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_order() const; // 0D_NOT_integer
  void set_ix_order(int value);
};

template <>
struct FortranTraits<LatEleOrderArrayStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_lat_ele_order_array_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_lat_ele_order_array_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_lat_ele_order_array_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "lat_ele_order_array_struct"; }
};

class LatEleOrderArrayStruct : public FortranProxy<LatEleOrderArrayStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  LatEleOrder1StructArray1D ele() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoLatSigmaStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_lat_sigma_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_lat_sigma_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_lat_sigma_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_lat_sigma_struct"; }
};

class TaoLatSigmaStruct : public FortranProxy<TaoLatSigmaStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray2D<double> mat() const; // 2D_NOT_real
  void set_mat(const std::vector<std::vector<double>> &v);
};

template <>
struct FortranTraits<TaoSpinEleStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_spin_ele_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_spin_ele_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_spin_ele_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_spin_ele_struct"; }
};

class TaoSpinEleStruct : public FortranProxy<TaoSpinEleStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoSpinDnDpzStruct dn_dpz() const; // 0D_NOT_type
  void set_dn_dpz(const TaoSpinDnDpzStruct &src);
  FArray1D<double> orb_eigen_val() const; // 1D_NOT_real
  void set_orb_eigen_val(const std::vector<double> &v);
  FArray2D<double> orb_eigen_vec() const; // 2D_NOT_real
  void set_orb_eigen_vec(const std::vector<std::vector<double>> &v);
  FArray2D<double> spin_eigen_vec() const; // 2D_NOT_real
  void set_spin_eigen_vec(const std::vector<std::vector<double>> &v);
  bool valid() const; // 0D_NOT_logical
  void set_valid(bool value);
};

template <>
struct FortranTraits<TaoPlotCacheStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_plot_cache_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_plot_cache_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_plot_cache_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_plot_cache_struct"; }
};

class TaoPlotCacheStruct : public FortranProxy<TaoPlotCacheStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  EleStruct ele_to_s() const; // 0D_NOT_type
  void set_ele_to_s(const EleStruct &src);
  CoordStruct orbit() const; // 0D_NOT_type
  void set_orbit(const CoordStruct &src);
  bool err() const; // 0D_NOT_logical
  void set_err(bool value);
};

template <>
struct FortranTraits<TaoSpinPolarizationStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_spin_polarization_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_spin_polarization_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_tao_spin_polarization_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "tao_spin_polarization_struct"; }
};

class TaoSpinPolarizationStruct : public FortranProxy<TaoSpinPolarizationStruct> {
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
  void set_pol_limit_dk_partial(const std::vector<double> &v);
  FArray1D<double> pol_limit_dk_partial2() const; // 1D_NOT_real
  void set_pol_limit_dk_partial2(const std::vector<double> &v);
  double pol_rate_bks() const; // 0D_NOT_real
  void set_pol_rate_bks(double value);
  double depol_rate() const; // 0D_NOT_real
  void set_depol_rate(double value);
  FArray1D<double> depol_rate_partial() const; // 1D_NOT_real
  void set_depol_rate_partial(const std::vector<double> &v);
  FArray1D<double> depol_rate_partial2() const; // 1D_NOT_real
  void set_depol_rate_partial2(const std::vector<double> &v);
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
  SpinOrbitMap1Struct q_1turn() const; // 0D_NOT_type
  void set_q_1turn(const SpinOrbitMap1Struct &src);
  SpinOrbitMap1StructArray1D q_ele() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoLatticeBranchStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_lattice_branch_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_lattice_branch_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_tao_lattice_branch_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_lattice_branch_struct"; }
};

class TaoLatticeBranchStruct : public FortranProxy<TaoLatticeBranchStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<TaoLatticeStruct> tao_lat() const; // 0D_PTR_type
  void set_tao_lat(const TaoLatticeStruct &src);
  TaoLatSigmaStructArray1D lat_sigma() const; // 1D_ALLOC_type
  TaoSpinEleStructArray1D spin_ele() const; // 1D_ALLOC_type
  BunchParamsStructArray1D bunch_params() const; // 1D_ALLOC_type
  BunchTrackStructArray1D bunch_params_comb() const; // 1D_ALLOC_type
  CoordStructArray1D orbit() const; // 1D_ALLOC_type
  TaoPlotCacheStructArray1D plot_cache() const; // 1D_ALLOC_type
  TaoSpinPolarizationStruct spin() const; // 0D_NOT_type
  void set_spin(const TaoSpinPolarizationStruct &src);
  SummationRdtStruct srdt() const; // 0D_NOT_type
  void set_srdt(const SummationRdtStruct &src);
  CoordStruct orb0() const; // 0D_NOT_type
  void set_orb0(const CoordStruct &src);
  NormalModesStruct modes_ri() const; // 0D_NOT_type
  void set_modes_ri(const NormalModesStruct &src);
  NormalModesStruct modes_6d() const; // 0D_NOT_type
  void set_modes_6d(const NormalModesStruct &src);
  PtcNormalFormStruct ptc_normal_form() const; // 0D_NOT_type
  void set_ptc_normal_form(const PtcNormalFormStruct &src);
  BmadNormalFormStruct bmad_normal_form() const; // 0D_NOT_type
  void set_bmad_normal_form(const BmadNormalFormStruct &src);
  CoordStructArray1D high_E_orb() const; // 1D_ALLOC_type
  CoordStructArray1D low_E_orb() const; // 1D_ALLOC_type
  TaylorStructArray1D taylor_save() const; // 1D_NOT_type
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
struct FortranTraits<TaoModelElementStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_model_element_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_model_element_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_tao_model_element_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_model_element_struct"; }
};

class TaoModelElementStruct : public FortranProxy<TaoModelElementStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  BeamStruct beam() const; // 0D_NOT_type
  void set_beam(const BeamStruct &src);
  bool save_beam_internally() const; // 0D_NOT_logical
  void set_save_beam_internally(bool value);
  bool save_beam_to_file() const; // 0D_NOT_logical
  void set_save_beam_to_file(bool value);
};

template <>
struct FortranTraits<TaoBeamBranchStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_beam_branch_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_beam_branch_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_beam_branch_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_beam_branch_struct"; }
};

class TaoBeamBranchStruct : public FortranProxy<TaoBeamBranchStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  BeamStruct beam_at_start() const; // 0D_NOT_type
  void set_beam_at_start(const BeamStruct &src);
  BeamInitStruct beam_init() const; // 0D_NOT_type
  void set_beam_init(const BeamInitStruct &src);
  BeamInitStruct beam_init_used() const; // 0D_NOT_type
  void set_beam_init_used(const BeamInitStruct &src);
  bool init_starting_distribution() const; // 0D_NOT_logical
  void set_init_starting_distribution(bool value);
  std::string track_start() const; // 0D_NOT_character
  void set_track_start(const std::string &value);
  std::string track_end() const; // 0D_NOT_character
  void set_track_end(const std::string &value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_track_start() const; // 0D_NOT_integer
  void set_ix_track_start(int value);
  int ix_track_end() const; // 0D_NOT_integer
  void set_ix_track_end(int value);
};

template <>
struct FortranTraits<TaoD1DataStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_d1_data_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_d1_data_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_d1_data_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_d1_data_struct"; }
};

class TaoD1DataStruct : public FortranProxy<TaoD1DataStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  std::optional<TaoD2DataStruct> d2() const; // 0D_PTR_type
  void set_d2(const TaoD2DataStruct &src);
  TaoDataStructArray1D d() const; // 1D_PTR_type
};

template <>
struct FortranTraits<TaoD2DataStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_d2_data_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_d2_data_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_d2_data_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_d2_data_struct"; }
};

class TaoD2DataStruct : public FortranProxy<TaoD2DataStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  std::string data_file_name() const; // 0D_NOT_character
  void set_data_file_name(const std::string &value);
  std::string ref_file_name() const; // 0D_NOT_character
  void set_ref_file_name(const std::string &value);
  std::string data_date() const; // 0D_NOT_character
  void set_data_date(const std::string &value);
  std::string ref_date() const; // 0D_NOT_character
  void set_ref_date(const std::string &value);
  FCharArray1D descrip() const; // 1D_NOT_character
  TaoD1DataStructArray1D d1() const; // 1D_ALLOC_type
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
struct FortranTraits<TaoDataVarComponentStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_data_var_component_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_data_var_component_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_tao_data_var_component_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "tao_data_var_component_struct"; }
};

class TaoDataVarComponentStruct : public FortranProxy<TaoDataVarComponentStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  double sign() const; // 0D_NOT_real
  void set_sign(double value);
};

template <>
struct FortranTraits<TaoGraphStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_graph_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_graph_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_graph_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_graph_struct"; }
};

class TaoGraphStruct : public FortranProxy<TaoGraphStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  std::string type() const; // 0D_NOT_character
  void set_type(const std::string &value);
  std::string title() const; // 0D_NOT_character
  void set_title(const std::string &value);
  std::string title_suffix() const; // 0D_NOT_character
  void set_title_suffix(const std::string &value);
  FCharArray1D text_legend() const; // 1D_NOT_character
  FCharArray1D text_legend_out() const; // 1D_NOT_character
  std::string why_invalid() const; // 0D_NOT_character
  void set_why_invalid(const std::string &value);
  TaoCurveStructArray1D curve() const; // 1D_ALLOC_type
  std::optional<TaoPlotStruct> p() const; // 0D_PTR_type
  void set_p(const TaoPlotStruct &src);
  TaoFloorPlanStruct floor_plan() const; // 0D_NOT_type
  void set_floor_plan(const TaoFloorPlanStruct &src);
  QpPointStruct text_legend_origin() const; // 0D_NOT_type
  void set_text_legend_origin(const QpPointStruct &src);
  QpPointStruct curve_legend_origin() const; // 0D_NOT_type
  void set_curve_legend_origin(const QpPointStruct &src);
  QpLegendStruct curve_legend() const; // 0D_NOT_type
  void set_curve_legend(const QpLegendStruct &src);
  QpAxisStruct x() const; // 0D_NOT_type
  void set_x(const QpAxisStruct &src);
  QpAxisStruct y() const; // 0D_NOT_type
  void set_y(const QpAxisStruct &src);
  QpAxisStruct x2() const; // 0D_NOT_type
  void set_x2(const QpAxisStruct &src);
  QpAxisStruct y2() const; // 0D_NOT_type
  void set_y2(const QpAxisStruct &src);
  QpRectStruct margin() const; // 0D_NOT_type
  void set_margin(const QpRectStruct &src);
  QpRectStruct scale_margin() const; // 0D_NOT_type
  void set_scale_margin(const QpRectStruct &src);
  double x_axis_scale_factor() const; // 0D_NOT_real
  void set_x_axis_scale_factor(double value);
  double symbol_size_scale() const; // 0D_NOT_real
  void set_symbol_size_scale(double value);
  FArray1D<int> box() const; // 1D_NOT_integer
  void set_box(const std::vector<int> &v);
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
struct FortranTraits<TaoPlotStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_plot_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_plot_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_plot_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_plot_struct"; }
};

class TaoPlotStruct : public FortranProxy<TaoPlotStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  std::string description() const; // 0D_NOT_character
  void set_description(const std::string &value);
  TaoGraphStructArray1D graph() const; // 1D_ALLOC_type
  std::optional<TaoPlotRegionStruct> r() const; // 0D_PTR_type
  void set_r(const TaoPlotRegionStruct &src);
  int ix_plot() const; // 0D_NOT_integer
  void set_ix_plot(int value);
  int n_curve_pts() const; // 0D_NOT_integer
  void set_n_curve_pts(int value);
  std::string type() const; // 0D_NOT_character
  void set_type(const std::string &value);
  std::string x_axis_type() const; // 0D_NOT_character
  void set_x_axis_type(const std::string &value);
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
struct FortranTraits<TaoPlotRegionStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_plot_region_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_plot_region_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_plot_region_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_plot_region_struct"; }
};

class TaoPlotRegionStruct : public FortranProxy<TaoPlotRegionStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  TaoPlotStruct plot() const; // 0D_NOT_type
  void set_plot(const TaoPlotStruct &src);
  FArray1D<double> location() const; // 1D_NOT_real
  void set_location(const std::vector<double> &v);
  bool visible() const; // 0D_NOT_logical
  void set_visible(bool value);
  bool list_with_show_plot_command() const; // 0D_NOT_logical
  void set_list_with_show_plot_command(bool value);
  bool setup_done() const; // 0D_NOT_logical
  void set_setup_done(bool value);
};

template <>
struct FortranTraits<TaoUniversePointerStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_universe_pointer_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_universe_pointer_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_tao_universe_pointer_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "tao_universe_pointer_struct"; }
};

class TaoUniversePointerStruct : public FortranProxy<TaoUniversePointerStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<TaoUniverseStruct> u() const; // 0D_PTR_type
  void set_u(const TaoUniverseStruct &src);
};

template <>
struct FortranTraits<TaoSuperUniverseStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_super_universe_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_super_universe_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_tao_super_universe_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_super_universe_struct"; }
};

class TaoSuperUniverseStruct : public FortranProxy<TaoSuperUniverseStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoGlobalStruct global() const; // 0D_NOT_type
  void set_global(const TaoGlobalStruct &src);
  TaoInitStruct init() const; // 0D_NOT_type
  void set_init(const TaoInitStruct &src);
  TaoCommonStruct com() const; // 0D_NOT_type
  void set_com(const TaoCommonStruct &src);
  TaoPlotPageStruct plot_page() const; // 0D_NOT_type
  void set_plot_page(const TaoPlotPageStruct &src);
  TaoV1VarStructArray1D v1_var() const; // 1D_ALLOC_type
  TaoVarStructArray1D var() const; // 1D_ALLOC_type
  TaoUniverseStructArray1D u() const; // 1D_ALLOC_type
  FArray1D<int> key() const; // 1D_ALLOC_integer
  void set_key(const std::vector<int> &v);
  TaoBuildingWallStruct building_wall() const; // 0D_NOT_type
  void set_building_wall(const TaoBuildingWallStruct &src);
  TaoWaveStruct wave() const; // 0D_NOT_type
  void set_wave(const TaoWaveStruct &src);
  int n_var_used() const; // 0D_NOT_integer
  void set_n_var_used(int value);
  int n_v1_var_used() const; // 0D_NOT_integer
  void set_n_v1_var_used(int value);
  TaoCmdHistoryStructArray1D history() const; // 1D_NOT_type
  bool initialized() const; // 0D_NOT_logical
  void set_initialized(bool value);
};

template <>
struct FortranTraits<TaoVarStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_var_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_var_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_var_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_var_struct"; }
};

class TaoVarStruct : public FortranProxy<TaoVarStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string ele_name() const; // 0D_NOT_character
  void set_ele_name(const std::string &value);
  std::string attrib_name() const; // 0D_NOT_character
  void set_attrib_name(const std::string &value);
  std::string id() const; // 0D_NOT_character
  void set_id(const std::string &value);
  TaoVarSlaveStructArray1D slave() const; // 1D_ALLOC_type
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
  double *model_value() const; // 0D_PTR_real
  void set_model_value(double value);
  double *base_value() const; // 0D_PTR_real
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
  void set_merit_type(const std::string &value);
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
  std::optional<TaoV1VarStruct> v1() const; // 0D_PTR_type
  void set_v1(const TaoV1VarStruct &src);
};

template <>
struct FortranTraits<TaoVarSlaveStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_var_slave_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_var_slave_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_var_slave_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_var_slave_struct"; }
};

class TaoVarSlaveStruct : public FortranProxy<TaoVarSlaveStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int ix_uni() const; // 0D_NOT_integer
  void set_ix_uni(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  double *model_value() const; // 0D_PTR_real
  void set_model_value(double value);
  double *base_value() const; // 0D_PTR_real
  void set_base_value(double value);
};

template <>
struct FortranTraits<TaoLatticeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_lattice_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_lattice_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_lattice_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_lattice_struct"; }
};

class TaoLatticeStruct : public FortranProxy<TaoLatticeStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  LatStruct lat() const; // 0D_NOT_type
  void set_lat(const LatStruct &src);
  LatStruct high_E_lat() const; // 0D_NOT_type
  void set_high_E_lat(const LatStruct &src);
  LatStruct low_E_lat() const; // 0D_NOT_type
  void set_low_E_lat(const LatStruct &src);
  RadIntAllEleStruct rad_int_by_ele_ri() const; // 0D_NOT_type
  void set_rad_int_by_ele_ri(const RadIntAllEleStruct &src);
  RadIntAllEleStruct rad_int_by_ele_6d() const; // 0D_NOT_type
  void set_rad_int_by_ele_6d(const RadIntAllEleStruct &src);
  TaoLatticeBranchStructArray1D tao_branch() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoBeamUniStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_beam_uni_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_beam_uni_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_beam_uni_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_beam_uni_struct"; }
};

class TaoBeamUniStruct : public FortranProxy<TaoBeamUniStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string saved_at() const; // 0D_NOT_character
  void set_saved_at(const std::string &value);
  std::string dump_file() const; // 0D_NOT_character
  void set_dump_file(const std::string &value);
  std::string dump_at() const; // 0D_NOT_character
  void set_dump_at(const std::string &value);
  bool track_beam_in_universe() const; // 0D_NOT_logical
  void set_track_beam_in_universe(bool value);
  bool always_reinit() const; // 0D_NOT_logical
  void set_always_reinit(bool value);
};

template <>
struct FortranTraits<TaoDynamicApertureStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_dynamic_aperture_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_dynamic_aperture_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_tao_dynamic_aperture_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "tao_dynamic_aperture_struct"; }
};

class TaoDynamicApertureStruct : public FortranProxy<TaoDynamicApertureStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  ApertureParamStruct param() const; // 0D_NOT_type
  void set_param(const ApertureParamStruct &src);
  ApertureScanStructArray1D scan() const; // 1D_ALLOC_type
  FArray1D<double> pz() const; // 1D_ALLOC_real
  void set_pz(const std::vector<double> &v);
  double ellipse_scale() const; // 0D_NOT_real
  void set_ellipse_scale(double value);
  double a_emit() const; // 0D_NOT_real
  void set_a_emit(double value);
  double b_emit() const; // 0D_NOT_real
  void set_b_emit(double value);
};

template <>
struct FortranTraits<TaoModelBranchStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_model_branch_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_model_branch_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_model_branch_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_model_branch_struct"; }
};

class TaoModelBranchStruct : public FortranProxy<TaoModelBranchStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoModelElementStructArray1D ele() const; // 1D_ALLOC_type
  TaoBeamBranchStruct beam() const; // 0D_NOT_type
  void set_beam(const TaoBeamBranchStruct &src);
};

template <>
struct FortranTraits<TaoSpinMapStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_spin_map_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_spin_map_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_spin_map_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_spin_map_struct"; }
};

class TaoSpinMapStruct : public FortranProxy<TaoSpinMapStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool valid() const; // 0D_NOT_logical
  void set_valid(bool value);
  SpinOrbitMap1Struct map1() const; // 0D_NOT_type
  void set_map1(const SpinOrbitMap1Struct &src);
  SpinAxisStruct axis_input() const; // 0D_NOT_type
  void set_axis_input(const SpinAxisStruct &src);
  SpinAxisStruct axis0() const; // 0D_NOT_type
  void set_axis0(const SpinAxisStruct &src);
  SpinAxisStruct axis1() const; // 0D_NOT_type
  void set_axis1(const SpinAxisStruct &src);
  int ix_ele() const; // 0D_NOT_integer
  void set_ix_ele(int value);
  int ix_ref() const; // 0D_NOT_integer
  void set_ix_ref(int value);
  int ix_uni() const; // 0D_NOT_integer
  void set_ix_uni(int value);
  int ix_branch() const; // 0D_NOT_integer
  void set_ix_branch(int value);
  FArray2D<double> mat8() const; // 2D_NOT_real
  void set_mat8(const std::vector<std::vector<double>> &v);
};

template <>
struct FortranTraits<TaoDataStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_data_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_data_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_data_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_data_struct"; }
};

class TaoDataStruct : public FortranProxy<TaoDataStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string ele_name() const; // 0D_NOT_character
  void set_ele_name(const std::string &value);
  std::string ele_start_name() const; // 0D_NOT_character
  void set_ele_start_name(const std::string &value);
  std::string ele_ref_name() const; // 0D_NOT_character
  void set_ele_ref_name(const std::string &value);
  std::string data_type() const; // 0D_ALLOC_character
  void set_data_type(const std::string &value);
  std::string merit_type() const; // 0D_NOT_character
  void set_merit_type(const std::string &value);
  std::string id() const; // 0D_NOT_character
  void set_id(const std::string &value);
  std::string data_source() const; // 0D_NOT_character
  void set_data_source(const std::string &value);
  std::string why_invalid() const; // 0D_NOT_character
  void set_why_invalid(const std::string &value);
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
  TaoSpinMapStruct spin_map() const; // 0D_NOT_type
  void set_spin_map(const TaoSpinMapStruct &src);
  std::optional<TaoD1DataStruct> d1() const; // 0D_PTR_type
  void set_d1(const TaoD1DataStruct &src);
};

template <>
struct FortranTraits<TaoPingScaleStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_ping_scale_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_ping_scale_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_ping_scale_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_ping_scale_struct"; }
};

class TaoPingScaleStruct : public FortranProxy<TaoPingScaleStruct> {
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
struct FortranTraits<TaoUniverseCalcStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_universe_calc_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_universe_calc_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_tao_universe_calc_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_universe_calc_struct"; }
};

class TaoUniverseCalcStruct : public FortranProxy<TaoUniverseCalcStruct> {
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
struct FortranTraits<LatEleOrderStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_lat_ele_order_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_lat_ele_order_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_lat_ele_order_struct(src, dst); }
  static constexpr std::string_view type_name() { return "lat_ele_order_struct"; }
};

class LatEleOrderStruct : public FortranProxy<LatEleOrderStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  LatEleOrderArrayStructArray1D branch() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoExpressionInfoStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_expression_info_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_expression_info_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_tao_expression_info_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "tao_expression_info_struct"; }
};

class TaoExpressionInfoStruct : public FortranProxy<TaoExpressionInfoStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  bool good() const; // 0D_NOT_logical
  void set_good(bool value);
  std::optional<EleStruct> ele() const; // 0D_PTR_type
  void set_ele(const EleStruct &src);
  double s() const; // 0D_NOT_real
  void set_s(double value);
};

template <>
struct FortranTraits<TaoEvalNodeStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_eval_node_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_eval_node_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_eval_node_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_eval_node_struct"; }
};

class TaoEvalNodeStruct : public FortranProxy<TaoEvalNodeStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int type() const; // 0D_NOT_integer
  void set_type(int value);
  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  double scale() const; // 0D_NOT_real
  void set_scale(double value);
  FArray1D<double> value() const; // 1D_ALLOC_real
  void set_value(const std::vector<double> &v);
  TaoExpressionInfoStructArray1D info() const; // 1D_ALLOC_type
  TaoEvalNodeStructArray1D node() const; // 1D_PTR_type
};

template <>
struct FortranTraits<TaoTitleStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_title_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_title_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_title_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_title_struct"; }
};

class TaoTitleStruct : public FortranProxy<TaoTitleStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string string() const; // 0D_NOT_character
  void set_string(const std::string &value);
  double x() const; // 0D_NOT_real
  void set_x(double value);
  double y() const; // 0D_NOT_real
  void set_y(double value);
  std::string units() const; // 0D_NOT_character
  void set_units(const std::string &value);
  std::string justify() const; // 0D_NOT_character
  void set_justify(const std::string &value);
  bool draw_it() const; // 0D_NOT_logical
  void set_draw_it(bool value);
};

template <>
struct FortranTraits<QpRectStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_qp_rect_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_qp_rect_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_qp_rect_struct(src, dst); }
  static constexpr std::string_view type_name() { return "qp_rect_struct"; }
};

class QpRectStruct : public FortranProxy<QpRectStruct> {
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
  void set_units(const std::string &value);
};

template <>
struct FortranTraits<TaoDrawingStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_drawing_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_drawing_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_drawing_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_drawing_struct"; }
};

class TaoDrawingStruct : public FortranProxy<TaoDrawingStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoEleShapeStructArray1D ele_shape() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoShapePatternStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_shape_pattern_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_shape_pattern_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_tao_shape_pattern_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_shape_pattern_struct"; }
};

class TaoShapePatternStruct : public FortranProxy<TaoShapePatternStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  QpLineStruct line() const; // 0D_NOT_type
  void set_line(const QpLineStruct &src);
  TaoShapePatternPointStructArray1D pt() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoShapePatternPointStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_shape_pattern_point_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_shape_pattern_point_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_tao_shape_pattern_point_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "tao_shape_pattern_point_struct"; }
};

class TaoShapePatternPointStruct : public FortranProxy<TaoShapePatternPointStruct> {
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
struct FortranTraits<QpAxisStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_qp_axis_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_qp_axis_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_qp_axis_struct(src, dst); }
  static constexpr std::string_view type_name() { return "qp_axis_struct"; }
};

class QpAxisStruct : public FortranProxy<QpAxisStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string label() const; // 0D_NOT_character
  void set_label(const std::string &value);
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
  void set_label_color(const std::string &value);
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
  void set_type(const std::string &value);
  std::string bounds() const; // 0D_NOT_character
  void set_bounds(const std::string &value);
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
struct FortranTraits<QpLegendStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_qp_legend_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_qp_legend_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_qp_legend_struct(src, dst); }
  static constexpr std::string_view type_name() { return "qp_legend_struct"; }
};

class QpLegendStruct : public FortranProxy<QpLegendStruct> {
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
struct FortranTraits<QpPointStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_qp_point_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_qp_point_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_qp_point_struct(src, dst); }
  static constexpr std::string_view type_name() { return "qp_point_struct"; }
};

class QpPointStruct : public FortranProxy<QpPointStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double x() const; // 0D_NOT_real
  void set_x(double value);
  double y() const; // 0D_NOT_real
  void set_y(double value);
  std::string units() const; // 0D_NOT_character
  void set_units(const std::string &value);
};

template <>
struct FortranTraits<QpLineStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_qp_line_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_qp_line_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_qp_line_struct(src, dst); }
  static constexpr std::string_view type_name() { return "qp_line_struct"; }
};

class QpLineStruct : public FortranProxy<QpLineStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int width() const; // 0D_NOT_integer
  void set_width(int value);
  std::string color() const; // 0D_NOT_character
  void set_color(const std::string &value);
  std::string pattern() const; // 0D_NOT_character
  void set_pattern(const std::string &value);
};

template <>
struct FortranTraits<QpSymbolStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_qp_symbol_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_qp_symbol_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_qp_symbol_struct(src, dst); }
  static constexpr std::string_view type_name() { return "qp_symbol_struct"; }
};

class QpSymbolStruct : public FortranProxy<QpSymbolStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string type() const; // 0D_NOT_character
  void set_type(const std::string &value);
  double height() const; // 0D_NOT_real
  void set_height(double value);
  std::string color() const; // 0D_NOT_character
  void set_color(const std::string &value);
  std::string fill_pattern() const; // 0D_NOT_character
  void set_fill_pattern(const std::string &value);
  int line_width() const; // 0D_NOT_integer
  void set_line_width(int value);
};

template <>
struct FortranTraits<TaoFloorPlanStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_floor_plan_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_floor_plan_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_floor_plan_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_floor_plan_struct"; }
};

class TaoFloorPlanStruct : public FortranProxy<TaoFloorPlanStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string view() const; // 0D_NOT_character
  void set_view(const std::string &value);
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
  void set_orbit_color(const std::string &value);
  std::string orbit_pattern() const; // 0D_NOT_character
  void set_orbit_pattern(const std::string &value);
  std::string orbit_lattice() const; // 0D_NOT_character
  void set_orbit_lattice(const std::string &value);
  int orbit_width() const; // 0D_NOT_integer
  void set_orbit_width(int value);
};

template <>
struct FortranTraits<TaoV1VarStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_v1_var_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_v1_var_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_v1_var_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_v1_var_struct"; }
};

class TaoV1VarStruct : public FortranProxy<TaoV1VarStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  int ix_v1_var() const; // 0D_NOT_integer
  void set_ix_v1_var(int value);
  TaoVarStructArray1D v() const; // 1D_PTR_type
};

template <>
struct FortranTraits<TaoGlobalStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_global_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_global_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_global_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_global_struct"; }
};

class TaoGlobalStruct : public FortranProxy<TaoGlobalStruct> {
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
  void set_quiet(const std::string &value);
  std::string random_engine() const; // 0D_NOT_character
  void set_random_engine(const std::string &value);
  std::string random_gauss_converter() const; // 0D_NOT_character
  void set_random_gauss_converter(const std::string &value);
  std::string track_type() const; // 0D_NOT_character
  void set_track_type(const std::string &value);
  std::string lat_sigma_calc_uses_emit_from() const; // 0D_NOT_character
  void set_lat_sigma_calc_uses_emit_from(const std::string &value);
  std::string prompt_string() const; // 0D_NOT_character
  void set_prompt_string(const std::string &value);
  std::string prompt_color() const; // 0D_NOT_character
  void set_prompt_color(const std::string &value);
  std::string optimizer() const; // 0D_NOT_character
  void set_optimizer(const std::string &value);
  std::string print_command() const; // 0D_NOT_character
  void set_print_command(const std::string &value);
  std::string var_out_file() const; // 0D_NOT_character
  void set_var_out_file(const std::string &value);
  std::string history_file() const; // 0D_NOT_character
  void set_history_file(const std::string &value);
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
struct FortranTraits<TaoInitStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_init_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_init_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_init_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_init_struct"; }
};

class TaoInitStruct : public FortranProxy<TaoInitStruct> {
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
  void set_init_name(const std::string &value);
  std::string hook_init_file() const; // 0D_NOT_character
  void set_hook_init_file(const std::string &value);
  std::string hook_lat_file() const; // 0D_NOT_character
  void set_hook_lat_file(const std::string &value);
  std::string hook_beam_file() const; // 0D_NOT_character
  void set_hook_beam_file(const std::string &value);
  std::string hook_data_file() const; // 0D_NOT_character
  void set_hook_data_file(const std::string &value);
  std::string hook_plot_file() const; // 0D_NOT_character
  void set_hook_plot_file(const std::string &value);
  std::string hook_startup_file() const; // 0D_NOT_character
  void set_hook_startup_file(const std::string &value);
  std::string hook_var_file() const; // 0D_NOT_character
  void set_hook_var_file(const std::string &value);
  std::string hook_building_wall_file() const; // 0D_NOT_character
  void set_hook_building_wall_file(const std::string &value);
  std::string init_file_arg_path() const; // 0D_NOT_character
  void set_init_file_arg_path(const std::string &value);
  std::string lattice_file_arg() const; // 0D_NOT_character
  void set_lattice_file_arg(const std::string &value);
  std::string hook_init_file_arg() const; // 0D_NOT_character
  void set_hook_init_file_arg(const std::string &value);
  std::string init_file_arg() const; // 0D_NOT_character
  void set_init_file_arg(const std::string &value);
  std::string beam_file_arg() const; // 0D_NOT_character
  void set_beam_file_arg(const std::string &value);
  std::string beam_init_position_file_arg() const; // 0D_NOT_character
  void set_beam_init_position_file_arg(const std::string &value);
  std::string command_arg() const; // 0D_NOT_character
  void set_command_arg(const std::string &value);
  std::string data_file_arg() const; // 0D_NOT_character
  void set_data_file_arg(const std::string &value);
  std::string plot_file_arg() const; // 0D_NOT_character
  void set_plot_file_arg(const std::string &value);
  std::string startup_file_arg() const; // 0D_NOT_character
  void set_startup_file_arg(const std::string &value);
  std::string var_file_arg() const; // 0D_NOT_character
  void set_var_file_arg(const std::string &value);
  std::string building_wall_file_arg() const; // 0D_NOT_character
  void set_building_wall_file_arg(const std::string &value);
  std::string geometry_arg() const; // 0D_NOT_character
  void set_geometry_arg(const std::string &value);
  std::string slice_lattice_arg() const; // 0D_NOT_character
  void set_slice_lattice_arg(const std::string &value);
  std::string start_branch_at_arg() const; // 0D_NOT_character
  void set_start_branch_at_arg(const std::string &value);
  std::string log_startup_arg() const; // 0D_NOT_character
  void set_log_startup_arg(const std::string &value);
  std::string no_stopping_arg() const; // 0D_NOT_character
  void set_no_stopping_arg(const std::string &value);
  std::string noplot_arg() const; // 0D_NOT_character
  void set_noplot_arg(const std::string &value);
  std::string no_rad_int_arg() const; // 0D_NOT_character
  void set_no_rad_int_arg(const std::string &value);
  std::string reverse_arg() const; // 0D_NOT_character
  void set_reverse_arg(const std::string &value);
  std::string debug_arg() const; // 0D_NOT_character
  void set_debug_arg(const std::string &value);
  std::string disable_smooth_line_calc_arg() const; // 0D_NOT_character
  void set_disable_smooth_line_calc_arg(const std::string &value);
  std::string rf_on_arg() const; // 0D_NOT_character
  void set_rf_on_arg(const std::string &value);
  std::string prompt_color_arg() const; // 0D_NOT_character
  void set_prompt_color_arg(const std::string &value);
  std::string quiet_arg() const; // 0D_NOT_character
  void set_quiet_arg(const std::string &value);
  std::string noinit_arg() const; // 0D_NOT_character
  void set_noinit_arg(const std::string &value);
  std::string nostartup_arg() const; // 0D_NOT_character
  void set_nostartup_arg(const std::string &value);
  std::string symbol_import_arg() const; // 0D_NOT_character
  void set_symbol_import_arg(const std::string &value);
  std::string unique_name_suffix() const; // 0D_NOT_character
  void set_unique_name_suffix(const std::string &value);
};

template <>
struct FortranTraits<TaoCommonStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_common_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_common_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_common_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_common_struct"; }
};

class TaoCommonStruct : public FortranProxy<TaoCommonStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoPlotRegionStructArray1D plot_place_buffer() const; // 1D_ALLOC_type
  FArray2D<double> covar() const; // 2D_ALLOC_real
  void set_covar(const std::vector<std::vector<double>> &v);
  FArray2D<double> alpha() const; // 2D_ALLOC_real
  void set_alpha(const std::vector<std::vector<double>> &v);
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
  FArray1D<bool> is_err_message_printed() const; // 1D_NOT_logical
  void set_is_err_message_printed(const std::vector<bool> &v);
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
  void set_single_mode_buffer(const std::string &value);
  std::string cmd() const; // 0D_NOT_character
  void set_cmd(const std::string &value);
  std::string saved_cmd_line() const; // 0D_NOT_character
  void set_saved_cmd_line(const std::string &value);
};

template <>
struct FortranTraits<TaoPlotPageStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_plot_page_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_plot_page_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_plot_page_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_plot_page_struct"; }
};

class TaoPlotPageStruct : public FortranProxy<TaoPlotPageStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoTitleStruct title() const; // 0D_NOT_type
  void set_title(const TaoTitleStruct &src);
  TaoTitleStruct subtitle() const; // 0D_NOT_type
  void set_subtitle(const TaoTitleStruct &src);
  QpRectStruct border() const; // 0D_NOT_type
  void set_border(const QpRectStruct &src);
  TaoDrawingStruct floor_plan() const; // 0D_NOT_type
  void set_floor_plan(const TaoDrawingStruct &src);
  TaoDrawingStruct lat_layout() const; // 0D_NOT_type
  void set_lat_layout(const TaoDrawingStruct &src);
  TaoShapePatternStructArray1D pattern() const; // 1D_ALLOC_type
  TaoPlotStructArray1D template_() const; // 1D_ALLOC_type
  TaoPlotRegionStructArray1D region() const; // 1D_ALLOC_type
  std::string plot_display_type() const; // 0D_NOT_character
  void set_plot_display_type(const std::string &value);
  FArray1D<double> size() const; // 1D_NOT_real
  void set_size(const std::vector<double> &v);
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
struct FortranTraits<TaoBuildingWallStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_building_wall_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_building_wall_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) { copy_fortran_tao_building_wall_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_building_wall_struct"; }
};

class TaoBuildingWallStruct : public FortranProxy<TaoBuildingWallStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TaoBuildingWallOrientationStruct orientation() const; // 0D_NOT_type
  void set_orientation(const TaoBuildingWallOrientationStruct &src);
  TaoBuildingWallSectionStructArray1D section() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoBuildingWallOrientationStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_building_wall_orientation_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_building_wall_orientation_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_tao_building_wall_orientation_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "tao_building_wall_orientation_struct"; }
};

class TaoBuildingWallOrientationStruct : public FortranProxy<TaoBuildingWallOrientationStruct> {
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
struct FortranTraits<TaoBuildingWallSectionStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_building_wall_section_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_building_wall_section_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_tao_building_wall_section_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "tao_building_wall_section_struct"; }
};

class TaoBuildingWallSectionStruct : public FortranProxy<TaoBuildingWallSectionStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string name() const; // 0D_NOT_character
  void set_name(const std::string &value);
  std::string constraint() const; // 0D_NOT_character
  void set_constraint(const std::string &value);
  TaoBuildingWallPointStructArray1D point() const; // 1D_ALLOC_type
};

template <>
struct FortranTraits<TaoBuildingWallPointStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_building_wall_point_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept {
    deallocate_fortran_tao_building_wall_point_struct(ptr, 0);
  }
  static void copy(const void *src, void *dst) {
    copy_fortran_tao_building_wall_point_struct(src, dst);
  }
  static constexpr std::string_view type_name() { return "tao_building_wall_point_struct"; }
};

class TaoBuildingWallPointStruct : public FortranProxy<TaoBuildingWallPointStruct> {
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
struct FortranTraits<TaoWaveStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_wave_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_wave_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_wave_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_wave_struct"; }
};

class TaoWaveStruct : public FortranProxy<TaoWaveStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string data_type() const; // 0D_NOT_character
  void set_data_type(const std::string &value);
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
  void set_amp_a(const std::vector<double> &v);
  FArray1D<double> amp_b() const; // 1D_NOT_real
  void set_amp_b(const std::vector<double> &v);
  FArray1D<double> amp_ba() const; // 1D_NOT_real
  void set_amp_ba(const std::vector<double> &v);
  FArray1D<double> coef_a() const; // 1D_NOT_real
  void set_coef_a(const std::vector<double> &v);
  FArray1D<double> coef_b() const; // 1D_NOT_real
  void set_coef_b(const std::vector<double> &v);
  FArray1D<double> coef_ba() const; // 1D_NOT_real
  void set_coef_ba(const std::vector<double> &v);
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
  void set_ix_data(const std::vector<int> &v);
  int n_kick() const; // 0D_NOT_integer
  void set_n_kick(int value);
  TaoWaveKickPtStructArray1D kick() const; // 1D_ALLOC_type
  TaoGraphStruct base_graph() const; // 0D_NOT_type
  void set_base_graph(const TaoGraphStruct &src);
  std::optional<TaoPlotRegionStruct> region() const; // 0D_PTR_type
  void set_region(const TaoPlotRegionStruct &src);
  std::optional<TaoD1DataStruct> d1_dat() const; // 0D_PTR_type
  void set_d1_dat(const TaoD1DataStruct &src);
};

template <>
struct FortranTraits<TaoWaveKickPtStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_wave_kick_pt_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_wave_kick_pt_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_wave_kick_pt_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_wave_kick_pt_struct"; }
};

class TaoWaveKickPtStruct : public FortranProxy<TaoWaveKickPtStruct> {
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
  std::optional<EleStruct> ele() const; // 0D_PTR_type
  void set_ele(const EleStruct &src);
};

template <>
struct FortranTraits<TaoCmdHistoryStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_cmd_history_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_cmd_history_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_cmd_history_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_cmd_history_struct"; }
};

class TaoCmdHistoryStruct : public FortranProxy<TaoCmdHistoryStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string cmd() const; // 0D_ALLOC_character
  void set_cmd(const std::string &value);
  int ix() const; // 0D_NOT_integer
  void set_ix(int value);
};

template <>
struct FortranTraits<TaoUniverseStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_tao_universe_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_tao_universe_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_tao_universe_struct(src, dst); }
  static constexpr std::string_view type_name() { return "tao_universe_struct"; }
};

class TaoUniverseStruct : public FortranProxy<TaoUniverseStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::optional<TaoLatticeStruct> model() const; // 0D_PTR_type
  void set_model(const TaoLatticeStruct &src);
  std::optional<TaoLatticeStruct> design() const; // 0D_PTR_type
  void set_design(const TaoLatticeStruct &src);
  std::optional<TaoLatticeStruct> base() const; // 0D_PTR_type
  void set_base(const TaoLatticeStruct &src);
  TaoBeamUniStruct beam() const; // 0D_NOT_type
  void set_beam(const TaoBeamUniStruct &src);
  TaoDynamicApertureStruct dynamic_aperture() const; // 0D_NOT_type
  void set_dynamic_aperture(const TaoDynamicApertureStruct &src);
  TaoModelBranchStructArray1D model_branch() const; // 1D_PTR_type
  TaoD2DataStructArray1D d2_data() const; // 1D_ALLOC_type
  TaoDataStructArray1D data() const; // 1D_ALLOC_type
  TaoPingScaleStruct ping_scale() const; // 0D_NOT_type
  void set_ping_scale(const TaoPingScaleStruct &src);
  LatStruct scratch_lat() const; // 0D_NOT_type
  void set_scratch_lat(const LatStruct &src);
  TaoUniverseCalcStruct calc() const; // 0D_NOT_type
  void set_calc(const TaoUniverseCalcStruct &src);
  LatEleOrderStruct ele_order() const; // 0D_NOT_type
  void set_ele_order(const LatEleOrderStruct &src);
  TaoSpinMapStruct spin_map() const; // 0D_NOT_type
  void set_spin_map(const TaoSpinMapStruct &src);
  FArray2D<double> dModel_dVar() const; // 2D_ALLOC_real
  void set_dModel_dVar(const std::vector<std::vector<double>> &v);
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
struct FortranTraits<MadEnergyStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_mad_energy_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_mad_energy_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_mad_energy_struct(src, dst); }
  static constexpr std::string_view type_name() { return "mad_energy_struct"; }
};

class MadEnergyStruct : public FortranProxy<MadEnergyStruct> {
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
struct FortranTraits<MadMapStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_mad_map_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_mad_map_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_mad_map_struct(src, dst); }
  static constexpr std::string_view type_name() { return "mad_map_struct"; }
};

class MadMapStruct : public FortranProxy<MadMapStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  FArray1D<double> k() const; // 1D_NOT_real
  void set_k(const std::vector<double> &v);
  FArray2D<double> r() const; // 2D_NOT_real
  void set_r(const std::vector<std::vector<double>> &v);
  FArray3D<double> t() const; // 3D_NOT_real
  void set_t(const std::vector<std::vector<std::vector<double>>> &v);
};

template <>
struct FortranTraits<RandomStateStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_random_state_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_random_state_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_random_state_struct(src, dst); }
  static constexpr std::string_view type_name() { return "random_state_struct"; }
};

class RandomStateStruct : public FortranProxy<RandomStateStruct> {
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
  FArray1D<int64_t> ix_sobseq() const; // 1D_NOT_integer8
  void set_ix_sobseq(const std::vector<int64_t> &v);
  FArray1D<double> x_sobseq() const; // 1D_NOT_real
  void set_x_sobseq(const std::vector<double> &v);
};

template <>
struct FortranTraits<BbuStageStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bbu_stage_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_bbu_stage_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_bbu_stage_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bbu_stage_struct"; }
};

class BbuStageStruct : public FortranProxy<BbuStageStruct> {
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
  void set_ave_orb(const std::vector<double> &v);
  FArray1D<double> rms_orb() const; // 1D_NOT_real
  void set_rms_orb(const std::vector<double> &v);
  FArray1D<double> min_orb() const; // 1D_NOT_real
  void set_min_orb(const std::vector<double> &v);
  FArray1D<double> max_orb() const; // 1D_NOT_real
  void set_max_orb(const std::vector<double> &v);
  int n_orb() const; // 0D_NOT_integer
  void set_n_orb(int value);
};

template <>
struct FortranTraits<BbuBeamStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bbu_beam_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_bbu_beam_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_bbu_beam_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bbu_beam_struct"; }
};

class BbuBeamStruct : public FortranProxy<BbuBeamStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  BunchStructArray1D bunch() const; // 1D_ALLOC_type
  BbuStageStructArray1D stage() const; // 1D_ALLOC_type
  FArray1D<int> ix_ele_bunch() const; // 1D_ALLOC_integer
  void set_ix_ele_bunch(const std::vector<int> &v);
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
struct FortranTraits<BbuParamStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_bbu_param_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_bbu_param_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_bbu_param_struct(src, dst); }
  static constexpr std::string_view type_name() { return "bbu_param_struct"; }
};

class BbuParamStruct : public FortranProxy<BbuParamStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string lat_filename() const; // 0D_NOT_character
  void set_lat_filename(const std::string &value);
  std::string lat2_filename() const; // 0D_NOT_character
  void set_lat2_filename(const std::string &value);
  std::string bunch_by_bunch_info_file() const; // 0D_NOT_character
  void set_bunch_by_bunch_info_file(const std::string &value);
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
  void set_elname(const std::string &value);
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
  void set_ele_track_end(const std::string &value);
  int ix_ele_track_end() const; // 0D_NOT_integer
  void set_ix_ele_track_end(int value);
  bool regression() const; // 0D_NOT_logical
  void set_regression(bool value);
  bool normalize_z_to_rf() const; // 0D_NOT_logical
  void set_normalize_z_to_rf(bool value);
  bool ramp_on() const; // 0D_NOT_logical
  void set_ramp_on(bool value);
  FArray1D<double> ramp_pattern() const; // 1D_NOT_real
  void set_ramp_pattern(const std::vector<double> &v);
  int ramp_n_start() const; // 0D_NOT_integer
  void set_ramp_n_start(int value);
  int n_ramp_pattern() const; // 0D_NOT_integer
  void set_n_ramp_pattern(int value);
};

template <>
struct FortranTraits<Fibre> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_fibre(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_fibre(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_fibre(src, dst); }
  static constexpr std::string_view type_name() { return "fibre"; }
};

class Fibre : public FortranProxy<Fibre> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int *DIR() const; // 0D_PTR_integer
  void set_DIR(int value);
  std::optional<Fibre> PREVIOUS() const; // 0D_PTR_type
  void set_PREVIOUS(const Fibre &src);
  std::optional<Fibre> NEXT() const; // 0D_PTR_type
  void set_NEXT(const Fibre &src);
  std::optional<Layout> PARENT_LAYOUT() const; // 0D_PTR_type
  void set_PARENT_LAYOUT(const Layout &src);
  int *pos() const; // 0D_PTR_integer
  void set_pos(int value);
  double *BETA0() const; // 0D_PTR_real
  void set_BETA0(double value);
  double *GAMMA0I() const; // 0D_PTR_real
  void set_GAMMA0I(double value);
  double *GAMBET() const; // 0D_PTR_real
  void set_GAMBET(double value);
  double *MASS() const; // 0D_PTR_real
  void set_MASS(double value);
  double *CHARGE() const; // 0D_PTR_real
  void set_CHARGE(double value);
  double *AG() const; // 0D_PTR_real
  void set_AG(double value);
  std::optional<Fibre> P() const; // 0D_PTR_type
  void set_P(const Fibre &src);
  std::optional<Fibre> N() const; // 0D_PTR_type
  void set_N(const Fibre &src);
  int *loc() const; // 0D_PTR_integer
  void set_loc(int value);
};

template <>
struct FortranTraits<Layout> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_layout(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_layout(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_layout(src, dst); }
  static constexpr std::string_view type_name() { return "layout"; }
};

class Layout : public FortranProxy<Layout> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  std::string NAME() const; // 0D_PTR_character
  void set_NAME(const std::string &value);
  int *INDEX() const; // 0D_PTR_integer
  void set_INDEX(int value);
  double *HARMONIC_NUMBER() const; // 0D_PTR_real
  void set_HARMONIC_NUMBER(double value);
  bool *CLOSED() const; // 0D_PTR_logical
  void set_CLOSED(bool value);
  int *N() const; // 0D_PTR_integer
  void set_N(int value);
  int *NTHIN() const; // 0D_PTR_integer
  void set_NTHIN(int value);
  double *THIN() const; // 0D_PTR_real
  void set_THIN(double value);
  int *LASTPOS() const; // 0D_PTR_integer
  void set_LASTPOS(int value);
  std::optional<Fibre> LAST() const; // 0D_PTR_type
  void set_LAST(const Fibre &src);
  std::optional<Fibre> END() const; // 0D_PTR_type
  void set_END(const Fibre &src);
  std::optional<Fibre> START() const; // 0D_PTR_type
  void set_START(const Fibre &src);
  std::optional<Fibre> START_GROUND() const; // 0D_PTR_type
  void set_START_GROUND(const Fibre &src);
  std::optional<Fibre> END_GROUND() const; // 0D_PTR_type
  void set_END_GROUND(const Fibre &src);
  std::optional<Layout> NEXT() const; // 0D_PTR_type
  void set_NEXT(const Layout &src);
  std::optional<Layout> PREVIOUS() const; // 0D_PTR_type
  void set_PREVIOUS(const Layout &src);
};

template <>
struct FortranTraits<AllEncompassingStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_all_encompassing_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_all_encompassing_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_all_encompassing_struct(src, dst); }
  static constexpr std::string_view type_name() { return "all_encompassing_struct"; }
};

class AllEncompassingStruct : public FortranProxy<AllEncompassingStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  double real_rp_0d() const; // 0D_NOT_real
  void set_real_rp_0d(double value);
  FArray1D<double> real_rp_1d() const; // 1D_NOT_real
  void set_real_rp_1d(const std::vector<double> &v);
  FArray2D<double> real_rp_2d() const; // 2D_NOT_real
  void set_real_rp_2d(const std::vector<std::vector<double>> &v);
  FArray3D<double> real_rp_3d() const; // 3D_NOT_real
  void set_real_rp_3d(const std::vector<std::vector<std::vector<double>>> &v);
  double *real_rp_0d_ptr() const; // 0D_PTR_real
  void set_real_rp_0d_ptr(double value);
  FArray1D<double> real_rp_1d_ptr() const; // 1D_PTR_real
  void set_real_rp_1d_ptr(const std::vector<double> &v);
  FArray2D<double> real_rp_2d_ptr() const; // 2D_PTR_real
  void set_real_rp_2d_ptr(const std::vector<std::vector<double>> &v);
  FArray3D<double> real_rp_3d_ptr() const; // 3D_PTR_real
  void set_real_rp_3d_ptr(const std::vector<std::vector<std::vector<double>>> &v);
  FArray1D<double> real_rp_1d_alloc() const; // 1D_ALLOC_real
  void set_real_rp_1d_alloc(const std::vector<double> &v);
  FArray2D<double> real_rp_2d_alloc() const; // 2D_ALLOC_real
  void set_real_rp_2d_alloc(const std::vector<std::vector<double>> &v);
  FArray3D<double> real_rp_3d_alloc() const; // 3D_ALLOC_real
  void set_real_rp_3d_alloc(const std::vector<std::vector<std::vector<double>>> &v);
  double real_dp_0d() const; // 0D_NOT_real
  void set_real_dp_0d(double value);
  FArray1D<double> real_dp_1d() const; // 1D_NOT_real
  void set_real_dp_1d(const std::vector<double> &v);
  FArray2D<double> real_dp_2d() const; // 2D_NOT_real
  void set_real_dp_2d(const std::vector<std::vector<double>> &v);
  FArray3D<double> real_dp_3d() const; // 3D_NOT_real
  void set_real_dp_3d(const std::vector<std::vector<std::vector<double>>> &v);
  double *real_dp_0d_ptr() const; // 0D_PTR_real
  void set_real_dp_0d_ptr(double value);
  FArray1D<double> real_dp_1d_ptr() const; // 1D_PTR_real
  void set_real_dp_1d_ptr(const std::vector<double> &v);
  FArray2D<double> real_dp_2d_ptr() const; // 2D_PTR_real
  void set_real_dp_2d_ptr(const std::vector<std::vector<double>> &v);
  FArray3D<double> real_dp_3d_ptr() const; // 3D_PTR_real
  void set_real_dp_3d_ptr(const std::vector<std::vector<std::vector<double>>> &v);
  FArray1D<double> real_dp_1d_alloc() const; // 1D_ALLOC_real
  void set_real_dp_1d_alloc(const std::vector<double> &v);
  FArray2D<double> real_dp_2d_alloc() const; // 2D_ALLOC_real
  void set_real_dp_2d_alloc(const std::vector<std::vector<double>> &v);
  FArray3D<double> real_dp_3d_alloc() const; // 3D_ALLOC_real
  void set_real_dp_3d_alloc(const std::vector<std::vector<std::vector<double>>> &v);
  std::complex<double> complex_dp_0d() const; // 0D_NOT_complex
  void set_complex_dp_0d(std::complex<double> value);
  FArray1D<std::complex<double>> complex_dp_1d() const; // 1D_NOT_complex
  void set_complex_dp_1d(const std::vector<std::complex<double>> &v);
  FArray2D<std::complex<double>> complex_dp_2d() const; // 2D_NOT_complex
  void set_complex_dp_2d(const std::vector<std::vector<std::complex<double>>> &v);
  FArray3D<std::complex<double>> complex_dp_3d() const; // 3D_NOT_complex
  void set_complex_dp_3d(const std::vector<std::vector<std::vector<std::complex<double>>>> &v);
  FArray1D<std::complex<double>> complex_dp_1d_ptr() const; // 1D_PTR_complex
  void set_complex_dp_1d_ptr(const std::vector<std::complex<double>> &v);
  FArray2D<std::complex<double>> complex_dp_2d_ptr() const; // 2D_PTR_complex
  void set_complex_dp_2d_ptr(const std::vector<std::vector<std::complex<double>>> &v);
  FArray3D<std::complex<double>> complex_dp_3d_ptr() const; // 3D_PTR_complex
  void set_complex_dp_3d_ptr(const std::vector<std::vector<std::vector<std::complex<double>>>> &v);
  FArray1D<std::complex<double>> complex_dp_1d_alloc() const; // 1D_ALLOC_complex
  void set_complex_dp_1d_alloc(const std::vector<std::complex<double>> &v);
  FArray2D<std::complex<double>> complex_dp_2d_alloc() const; // 2D_ALLOC_complex
  void set_complex_dp_2d_alloc(const std::vector<std::vector<std::complex<double>>> &v);
  FArray3D<std::complex<double>> complex_dp_3d_alloc() const; // 3D_ALLOC_complex
  void set_complex_dp_3d_alloc(const std::vector<std::vector<std::vector<std::complex<double>>>> &v
  );
  int int_0d() const; // 0D_NOT_integer
  void set_int_0d(int value);
  FArray1D<int> int_1d() const; // 1D_NOT_integer
  void set_int_1d(const std::vector<int> &v);
  FArray2D<int> int_2d() const; // 2D_NOT_integer
  void set_int_2d(const std::vector<std::vector<int>> &v);
  FArray3D<int> int_3d() const; // 3D_NOT_integer
  void set_int_3d(const std::vector<std::vector<std::vector<int>>> &v);
  int *int_0d_ptr() const; // 0D_PTR_integer
  void set_int_0d_ptr(int value);
  FArray1D<int> int_1d_ptr() const; // 1D_PTR_integer
  void set_int_1d_ptr(const std::vector<int> &v);
  FArray2D<int> int_2d_ptr() const; // 2D_PTR_integer
  void set_int_2d_ptr(const std::vector<std::vector<int>> &v);
  FArray3D<int> int_3d_ptr() const; // 3D_PTR_integer
  void set_int_3d_ptr(const std::vector<std::vector<std::vector<int>>> &v);
  FArray1D<int> int_1d_alloc() const; // 1D_ALLOC_integer
  void set_int_1d_alloc(const std::vector<int> &v);
  FArray2D<int> int_2d_alloc() const; // 2D_ALLOC_integer
  void set_int_2d_alloc(const std::vector<std::vector<int>> &v);
  FArray3D<int> int_3d_alloc() const; // 3D_ALLOC_integer
  void set_int_3d_alloc(const std::vector<std::vector<std::vector<int>>> &v);
  int64_t int8_0d() const; // 0D_NOT_integer8
  void set_int8_0d(int64_t value);
  FArray1D<int64_t> int8_1d() const; // 1D_NOT_integer8
  void set_int8_1d(const std::vector<int64_t> &v);
  FArray2D<int64_t> int8_2d() const; // 2D_NOT_integer8
  void set_int8_2d(const std::vector<std::vector<int64_t>> &v);
  FArray3D<int64_t> int8_3d() const; // 3D_NOT_integer8
  void set_int8_3d(const std::vector<std::vector<std::vector<int64_t>>> &v);
  int64_t *int8_0d_ptr() const; // 0D_PTR_integer8
  void set_int8_0d_ptr(int64_t value);
  FArray1D<int64_t> int8_1d_ptr() const; // 1D_PTR_integer8
  void set_int8_1d_ptr(const std::vector<int64_t> &v);
  FArray2D<int64_t> int8_2d_ptr() const; // 2D_PTR_integer8
  void set_int8_2d_ptr(const std::vector<std::vector<int64_t>> &v);
  FArray3D<int64_t> int8_3d_ptr() const; // 3D_PTR_integer8
  void set_int8_3d_ptr(const std::vector<std::vector<std::vector<int64_t>>> &v);
  FArray1D<int64_t> int8_1d_alloc() const; // 1D_ALLOC_integer8
  void set_int8_1d_alloc(const std::vector<int64_t> &v);
  FArray2D<int64_t> int8_2d_alloc() const; // 2D_ALLOC_integer8
  void set_int8_2d_alloc(const std::vector<std::vector<int64_t>> &v);
  FArray3D<int64_t> int8_3d_alloc() const; // 3D_ALLOC_integer8
  void set_int8_3d_alloc(const std::vector<std::vector<std::vector<int64_t>>> &v);
  bool logical_0d() const; // 0D_NOT_logical
  void set_logical_0d(bool value);
  FArray1D<bool> logical_1d() const; // 1D_NOT_logical
  void set_logical_1d(const std::vector<bool> &v);
  FArray2D<bool> logical_2d() const; // 2D_NOT_logical
  void set_logical_2d(const std::vector<std::vector<bool>> &v);
  FArray3D<bool> logical_3d() const; // 3D_NOT_logical
  void set_logical_3d(const std::vector<std::vector<std::vector<bool>>> &v);
  bool *logical_0d_ptr() const; // 0D_PTR_logical
  void set_logical_0d_ptr(bool value);
  TestSubStruct type_0d() const; // 0D_NOT_type
  void set_type_0d(const TestSubStruct &src);
  TestSubStructArray1D type_1d() const; // 1D_NOT_type
  TestSubStructArray2D type_2d() const; // 2D_NOT_type
  TestSubStructArray3D type_3d() const; // 3D_NOT_type
  std::optional<TestSubStruct> type_0d_ptr() const; // 0D_PTR_type
  void set_type_0d_ptr(const TestSubStruct &src);
  TestSubStructArray1D type_1d_ptr() const; // 1D_PTR_type
  TestSubStructArray2D type_2d_ptr() const; // 2D_PTR_type
  TestSubStructArray3D type_3d_ptr() const; // 3D_PTR_type
  TestSubStructArray1D type_1d_alloc() const; // 1D_ALLOC_type
  TestSubStructArray2D type_2d_alloc() const; // 2D_ALLOC_type
  TestSubStructArray3D type_3d_alloc() const; // 3D_ALLOC_type
};

template <>
struct FortranTraits<TestSubStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_test_sub_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_test_sub_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_test_sub_struct(src, dst); }
  static constexpr std::string_view type_name() { return "test_sub_struct"; }
};

class TestSubStruct : public FortranProxy<TestSubStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  TestSubSubStruct sr() const; // 0D_NOT_type
  void set_sr(const TestSubSubStruct &src);
};

template <>
struct FortranTraits<TestSubSubStruct> {
  static void *allocate() {
    size_t sz;
    return allocate_fortran_test_sub_sub_struct(0, &sz);
  }
  static void deallocate(void *ptr) noexcept { deallocate_fortran_test_sub_sub_struct(ptr, 0); }
  static void copy(const void *src, void *dst) { copy_fortran_test_sub_sub_struct(src, dst); }
  static constexpr std::string_view type_name() { return "test_sub_sub_struct"; }
};

class TestSubSubStruct : public FortranProxy<TestSubSubStruct> {
public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  int64_t aaa() const; // 0D_NOT_integer8
  void set_aaa(int64_t value);
  int bbb() const; // 0D_NOT_integer
  void set_bbb(int value);
  std::string file() const; // 0D_NOT_character
  void set_file(const std::string &value);
  double t_ref() const; // 0D_NOT_real
  void set_t_ref(double value);
  double freq_spread() const; // 0D_NOT_real
  void set_freq_spread(double value);
};

} // namespace Bmad
