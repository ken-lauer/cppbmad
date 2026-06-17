# Sim Utils

Simulation utility helpers from the sim_utils library.

## Classes (Fortran Structures)

::: pybmad.BicubicCmplxCoefStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BicubicCmplxCoefStruct

Fortran struct: `bicubic_cmplx_coef_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/cubic_interpolation_mod.f90#L70))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | 2D array of complex (shape: 0:3,0:3) | Coefs |
| `i_box` | 1D array of int (shape: 2) | index at lower box corner. |

::: pybmad.NametableStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### NametableStruct

Fortran struct: `nametable_struct` ([`sim_utils/interfaces/sim_utils_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_struct.f90#L34))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | 1D array of str | Array of names. |
| `index` | 1D array of int | Sorted index for names(:) array. names(an_index(i)) is in alphabetical order. |
| `n_min` | int | Set to 0 for use in a lattice. |
| `n_max` | int | Use only names(n_min:n_max) part of array. |

::: pybmad.QpAxisStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### QpAxisStruct

Fortran struct: `qp_axis_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/plot/quick_plot_struct.f90#L58))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `label` | str |  |
| `min` | float | Axis min/max in data units. |
| `max` | float | Axis min/max in data units. |
| `tick_min` | float | Min tick location along axis in data units. |
| `tick_max` | float | Max tick location along axis in data units. |
| `eval_min` | float | For general use. Not set by quick_plot. |
| `eval_max` | float | For general use. Not set by quick_plot. |
| `dtick` | float | Distance between ticks. In data units. Ticks will be drawn between %min and %max. |
| `number_offset` | float | Offset from axis line in inches. |
| `label_offset` | float | Offset from numbers in inches. |
| `major_tick_len` | float | In inches. |
| `minor_tick_len` | float | In inches. |
| `label_color` | str | Color of the label. |
| `major_div` | int | Actual numbrer of major divisions |
| `major_div_nominal` | int | Nominal value. |
| `minor_div` | int | 0 = auto choose. |
| `minor_div_max` | int | Max number for auto choose. |
| `places` | int | Number of places after the decimal point to print. |
| `type` | str | Or 'LOG', or 'CUSTOM' |
| `bounds` | str | Or 'ZERO_AT_END' or 'ZERO_SYMMETRIC' |
| `tick_side` | int | +1 = Draw on the side inside the graph, 0 = both (longer tick), -1 = outside. |
| `number_side` | int | +1 = Draw to the side inside the graph, -1 = outside. |
| `draw_label` | bool |  |
| `draw_numbers` | bool |  |

::: pybmad.QpLegendStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### QpLegendStruct

Fortran struct: `qp_legend_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/plot/quick_plot_struct.f90#L138))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `row_spacing` | float | Spacing between rows. |
| `line_length` | float | Length of the line in points. |
| `text_offset` | float | Horizontal offset in points between the line and the text. |
| `draw_line` | bool | Draw lines? |
| `draw_symbol` | bool | Draw symbols? |
| `draw_text` | bool | Draw text? |

::: pybmad.QpLineStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### QpLineStruct

Fortran struct: `qp_line_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/plot/quick_plot_struct.f90#L116))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `width` | int |  |
| `color` | str |  |
| `pattern` | str |  |

::: pybmad.QpPointStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### QpPointStruct

Fortran struct: `qp_point_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/plot/quick_plot_struct.f90#L100))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x` | float |  |
| `y` | float |  |
| `units` | str |  |

::: pybmad.QpRectStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### QpRectStruct

Fortran struct: `qp_rect_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/plot/quick_plot_struct.f90#L105))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x1` | float |  |
| `x2` | float |  |
| `y1` | float |  |
| `y2` | float |  |
| `units` | str |  |

::: pybmad.QpSymbolStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### QpSymbolStruct

Fortran struct: `qp_symbol_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/plot/quick_plot_struct.f90#L122))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `type` | str |  |
| `height` | float | in points (same as text height) |
| `color` | str |  |
| `fill_pattern` | str |  |
| `line_width` | int |  |

::: pybmad.RandomStateStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RandomStateStruct

Fortran struct: `random_state_struct` ([`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L28))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ix` | int |  |
| `iy` | int |  |
| `number_stored` | bool |  |
| `h_saved` | float |  |
| `engine` | int | Params |
| `seed` | int |  |
| `am` | float |  |
| `gauss_converter` | int |  |
| `gauss_sigma_cut` | float | Only used if positive. |
| `in_sobseq` | int |  |
| `ix_sobseq` | 1D array of int (shape: sobseq_maxdim) |  |
| `x_sobseq` | 1D array of float (shape: sobseq_maxdim) |  |

::: pybmad.SplineStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SplineStruct

Fortran struct: `spline_struct` ([`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/spline_mod.f90#L18))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x0` | float | Point at start of spline |
| `y0` | float | Point at start of spline |
| `x1` | float | Point at end of spline |
| `coef` | 1D array of float (shape: 0:3) | coefficients for cubic spline |

::: pybmad.TricubicCmplxCoefStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TricubicCmplxCoefStruct

Fortran struct: `tricubic_cmplx_coef_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/cubic_interpolation_mod.f90#L126))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | 3D array of complex (shape: 0:3,0:3,0:3) | Coefs |
| `i_box` | 1D array of int (shape: 3) | index at lower box corner. |

## Procedures

### allocate_thread_states

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L1151)

::: pybmad.allocate_thread_states
    options:
      show_root_heading: false
      show_root_toc_entry: false

### anomalous_moment_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1430)

::: pybmad.anomalous_moment_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### antiparticle

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L882)

::: pybmad.antiparticle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apfft

Fortran source: [`sim_utils/math/all_phase_fft.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/all_phase_fft.f90#L135)

::: pybmad.apfft
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apfft_corr

Fortran source: [`sim_utils/math/all_phase_fft.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/all_phase_fft.f90#L39)

::: pybmad.apfft_corr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apfft_ext

Fortran source: [`sim_utils/math/all_phase_fft.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/all_phase_fft.f90#L174)

::: pybmad.apfft_ext
    options:
      show_root_heading: false
      show_root_toc_entry: false

### asinc

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L35)

::: pybmad.asinc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### assert_equal

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L43)

::: pybmad.assert_equal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### atomic_number

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1773)

::: pybmad.atomic_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### atomic_species_id

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1213)

::: pybmad.atomic_species_id
    options:
      show_root_heading: false
      show_root_toc_entry: false

### axis_angle_to_quat

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L385)

::: pybmad.axis_angle_to_quat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### axis_angle_to_w_mat

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L234)

::: pybmad.axis_angle_to_w_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bicubic_cmplx_eval

Fortran source: [`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/cubic_interpolation_mod.f90#L1103)

::: pybmad.bicubic_cmplx_eval
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bin_index

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/bin_mod.f90#L329)

::: pybmad.bin_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bin_x_center

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/bin_mod.f90#L353)

::: pybmad.bin_x_center
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bit_set

Fortran source: [`sim_utils/math/bit_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/bit_set_mod.f90#L22)

::: pybmad.bit_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bracket_index_for_spline

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/spline_mod.f90#L327)

::: pybmad.bracket_index_for_spline
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_file_number

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L80)

::: pybmad.calc_file_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### celbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L420)

::: pybmad.celbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cesr_getarg

Fortran source: [`sim_utils/io/command_line_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/command_line_mod.f90#L57)

::: pybmad.cesr_getarg
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cesr_iargc

Fortran source: [`sim_utils/io/command_line_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/command_line_mod.f90#L23)

::: pybmad.cesr_iargc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### change_file_number

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L88)

::: pybmad.change_file_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### charge_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1505)

::: pybmad.charge_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### charge_to_mass_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1669)

::: pybmad.charge_to_mass_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coarse_frequency_estimate

Fortran source: [`sim_utils/math/fourier_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/fourier_mod.f90#L29)

::: pybmad.coarse_frequency_estimate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_error_function

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L114)

::: pybmad.complex_error_function
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cos_one

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L94)

::: pybmad.cos_one
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cosc

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L27)

::: pybmad.cosc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_a_spline

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/spline_mod.f90#L105)

::: pybmad.create_a_spline
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cross_product

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L123)

::: pybmad.cross_product
    options:
      show_root_heading: false
      show_root_toc_entry: false

### date_and_time_stamp

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L130)

::: pybmad.date_and_time_stamp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### destfixedwindowls

Fortran source: [`sim_utils/math/windowLS.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/windowLS.f90#L145)

::: pybmad.destfixedwindowls
    options:
      show_root_heading: false
      show_root_toc_entry: false

### detab

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L136)

::: pybmad.detab
    options:
      show_root_heading: false
      show_root_toc_entry: false

### display_size_and_resolution

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L148)

::: pybmad.display_size_and_resolution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### dj_bessel

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L155)

::: pybmad.dj_bessel
    options:
      show_root_heading: false
      show_root_toc_entry: false

### djb_hash

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L162)

::: pybmad.djb_hash
    options:
      show_root_heading: false
      show_root_toc_entry: false

### djb_str_hash

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L169)

::: pybmad.djb_str_hash
    options:
      show_root_heading: false
      show_root_toc_entry: false

### downcase_string

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L181)

::: pybmad.downcase_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### elbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L271)

::: pybmad.elbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### elcbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L371)

::: pybmad.elcbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ellipinc

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L53)

::: pybmad.ellipinc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### elsbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L325)

::: pybmad.elsbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### end_akima_spline_calc

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/spline_mod.f90#L573)

::: pybmad.end_akima_spline_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### err_exit

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L192)

::: pybmad.err_exit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### factorial

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L197)

::: pybmad.factorial
    options:
      show_root_heading: false
      show_root_toc_entry: false

### faddeeva_function

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L204)

::: pybmad.faddeeva_function
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fft_1d

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L217)

::: pybmad.fft_1d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### file_directorizer

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L223)

::: pybmad.file_directorizer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### file_get

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L231)

::: pybmad.file_get
    options:
      show_root_heading: false
      show_root_toc_entry: false

### file_get_open

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L238)

::: pybmad.file_get_open
    options:
      show_root_heading: false
      show_root_toc_entry: false

### file_suffixer

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L247)

::: pybmad.file_suffixer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### find_location

Fortran sources (overloaded):

- `find_location_real`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L1033)
- `find_location_int`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L1039)
- `find_location_logic`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L1044)
- `find_location_str`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L1049)

::: pybmad.find_location
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fine_frequency_estimate

Fortran source: [`sim_utils/math/fourier_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/fourier_mod.f90#L113)

::: pybmad.fine_frequency_estimate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fixedwindowls

Fortran source: [`sim_utils/math/windowLS.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/windowLS.f90#L171)

::: pybmad.fixedwindowls
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fourier_amplitude

Fortran source: [`sim_utils/math/fourier_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/fourier_mod.f90#L207)

::: pybmad.fourier_amplitude
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gelbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L140)

::: pybmad.gelbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gen_complete_elliptic

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L255)

::: pybmad.gen_complete_elliptic
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_a_char

Fortran source: [`sim_utils/io/input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/input_mod.f90#L85)

::: pybmad.get_a_char
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_file_number

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L262)

::: pybmad.get_file_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_file_time_stamp

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L278)

::: pybmad.get_file_time_stamp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_tty_char

Fortran source: [`sim_utils/io/input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/input_mod.f90#L36)

::: pybmad.get_tty_char
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hanhan

Fortran source: [`sim_utils/math/all_phase_fft.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/all_phase_fft.f90#L277)

::: pybmad.hanhan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### i_bessel

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L283)

::: pybmad.i_bessel
    options:
      show_root_heading: false
      show_root_toc_entry: false

### i_bessel_extended

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L290)

::: pybmad.i_bessel_extended
    options:
      show_root_heading: false
      show_root_toc_entry: false

### increment_file_number

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L298)

::: pybmad.increment_file_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### index_nocase

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L306)

::: pybmad.index_nocase
    options:
      show_root_heading: false
      show_root_toc_entry: false

### initfixedwindowls

Fortran source: [`sim_utils/math/windowLS.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/windowLS.f90#L76)

::: pybmad.initfixedwindowls
    options:
      show_root_heading: false
      show_root_toc_entry: false

### initial_lmdif

Fortran source: [`sim_utils/optimizers/lmdif_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/optimizers/lmdif_mod.f90#L81)

::: pybmad.initial_lmdif
    options:
      show_root_heading: false
      show_root_toc_entry: false

### int_str

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L313)

::: pybmad.int_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### interpolated_fft

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/naff.f90#L328)

::: pybmad.interpolated_fft
    options:
      show_root_heading: false
      show_root_toc_entry: false

### interpolated_fft_gsl

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/naff.f90#L233)

::: pybmad.interpolated_fft_gsl
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_alphabetic

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L341)

::: pybmad.is_alphabetic
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_decreasing_sequence

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L348)

::: pybmad.is_decreasing_sequence
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_false

Fortran source: [`sim_utils/interfaces/sim_utils_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_struct.f90#L211)

::: pybmad.is_false
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_increasing_sequence

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L356)

::: pybmad.is_increasing_sequence
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_integer

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L364)

::: pybmad.is_integer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_logical

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L372)

::: pybmad.is_logical
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_real

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L379)

::: pybmad.is_real
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_subatomic_species

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1806)

::: pybmad.is_subatomic_species
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_true

Fortran source: [`sim_utils/interfaces/sim_utils_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_struct.f90#L178)

::: pybmad.is_true
    options:
      show_root_heading: false
      show_root_toc_entry: false

### j_bessel

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L396)

::: pybmad.j_bessel
    options:
      show_root_heading: false
      show_root_toc_entry: false

### linear_fit

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L403)

::: pybmad.linear_fit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### linear_fit_2d

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L415)

::: pybmad.linear_fit_2d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### logic_str

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L432)

::: pybmad.logic_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lunget

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L438)

::: pybmad.lunget
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_legal_comment

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L454)

::: pybmad.make_legal_comment
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mass_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1557)

::: pybmad.mass_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### match_reg

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L443)

::: pybmad.match_reg
    options:
      show_root_heading: false
      show_root_toc_entry: false

### match_wild

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L460)

::: pybmad.match_wild
    options:
      show_root_heading: false
      show_root_toc_entry: false

### match_word

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L551)

::: pybmad.match_word
    options:
      show_root_heading: false
      show_root_toc_entry: false

### maximize_projection

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/naff.f90#L164)

::: pybmad.maximize_projection
    options:
      show_root_heading: false
      show_root_toc_entry: false

### milli_sleep

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L449)

::: pybmad.milli_sleep
    options:
      show_root_heading: false
      show_root_toc_entry: false

### modulo2_dp

Fortran source: [`sim_utils/special_functions/modulo2_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/modulo2_mod.f90#L85)

::: pybmad.modulo2_dp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### modulo2_int

Fortran source: [`sim_utils/special_functions/modulo2_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/modulo2_mod.f90#L157)

::: pybmad.modulo2_int
    options:
      show_root_heading: false
      show_root_toc_entry: false

### modulo2_qp

Fortran source: [`sim_utils/special_functions/modulo2_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/modulo2_mod.f90#L121)

::: pybmad.modulo2_qp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### modulo2_sp

Fortran source: [`sim_utils/special_functions/modulo2_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/modulo2_mod.f90#L49)

::: pybmad.modulo2_sp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### n_bins_automatic

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/bin_mod.f90#L386)

::: pybmad.n_bins_automatic
    options:
      show_root_heading: false
      show_root_toc_entry: false

### n_choose_k

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L568)

::: pybmad.n_choose_k
    options:
      show_root_heading: false
      show_root_toc_entry: false

### n_spline_create

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L575)

::: pybmad.n_spline_create
    options:
      show_root_heading: false
      show_root_toc_entry: false

### naff

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/naff.f90#L72)

::: pybmad.naff
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nametable_add

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L583)

::: pybmad.nametable_add
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nametable_bracket_indexx

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L591)

::: pybmad.nametable_bracket_indexx
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nametable_change1

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L600)

::: pybmad.nametable_change1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nametable_init

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L608)

::: pybmad.nametable_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nametable_remove

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L615)

::: pybmad.nametable_remove
    options:
      show_root_heading: false
      show_root_toc_entry: false

### negative_ampsquared

Fortran source: [`sim_utils/math/fourier_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/fourier_mod.f90#L153)

::: pybmad.negative_ampsquared
    options:
      show_root_heading: false
      show_root_toc_entry: false

### negative_dampsquared

Fortran source: [`sim_utils/math/fourier_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/fourier_mod.f90#L173)

::: pybmad.negative_dampsquared
    options:
      show_root_heading: false
      show_root_toc_entry: false

### omega_to_quat

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L314)

::: pybmad.omega_to_quat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### openpmd_species_name

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1394)

::: pybmad.openpmd_species_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ordinal_str

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L633)

::: pybmad.ordinal_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### out_io

Fortran sources (overloaded):

- `out_io_real`: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/output_mod.f90#L249)
- `out_io_int`: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/output_mod.f90#L293)
- `out_io_logical`: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/output_mod.f90#L337)
- `out_io_line12`: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/output_mod.f90#L382)
- `out_io_lines`: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/output_mod.f90#L434)

::: pybmad.out_io
    options:
      show_root_heading: false
      show_root_toc_entry: false

### out_io_buffer_get_line

Fortran source: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/output_mod.f90#L842)

::: pybmad.out_io_buffer_get_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### out_io_buffer_num_lines

Fortran source: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/output_mod.f90#L822)

::: pybmad.out_io_buffer_num_lines
    options:
      show_root_heading: false
      show_root_toc_entry: false

### out_io_buffer_reset

Fortran source: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/output_mod.f90#L799)

::: pybmad.out_io_buffer_reset
    options:
      show_root_heading: false
      show_root_toc_entry: false

### out_io_print_and_capture_setup

Fortran source: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/output_mod.f90#L768)

::: pybmad.out_io_print_and_capture_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_fortran_format

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L646)

::: pybmad.parse_fortran_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_locations

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L652)

::: pybmad.pointer_to_locations
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_ran_state

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L1114)

::: pybmad.pointer_to_ran_state
    options:
      show_root_heading: false
      show_root_toc_entry: false

### poly_eval

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L661)

::: pybmad.poly_eval
    options:
      show_root_heading: false
      show_root_toc_entry: false

### probability_funct

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L668)

::: pybmad.probability_funct
    options:
      show_root_heading: false
      show_root_toc_entry: false

### projdd

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/naff.f90#L132)

::: pybmad.projdd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quadratic_roots

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L675)

::: pybmad.quadratic_roots
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_conj

Fortran sources (overloaded):

- `quat_conj_real`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L411)
- `quat_conj_complex`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L435)

::: pybmad.quat_conj
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_inverse

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L458)

::: pybmad.quat_inverse
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_mul

Fortran sources (overloaded):

- `quat_mul_real`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L484)
- `quat_mul_complex`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L551)

::: pybmad.quat_mul
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_rotate

Fortran sources (overloaded):

- `quat_rotate_real`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L616)
- `quat_rotate_complex`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L652)

::: pybmad.quat_rotate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_to_axis_angle

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L347)

::: pybmad.quat_to_axis_angle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_to_omega

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L281)

::: pybmad.quat_to_omega
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_to_w_mat

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L184)

::: pybmad.quat_to_w_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### query_string

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L684)

::: pybmad.query_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quote

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L693)

::: pybmad.quote
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quoten

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L698)

::: pybmad.quoten
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_default_state

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L891)

::: pybmad.ran_default_state
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_engine

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L622)

::: pybmad.ran_engine
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_gauss_converter

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L700)

::: pybmad.ran_gauss_converter
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_gauss_scalar

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L171)

::: pybmad.ran_gauss_scalar
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_gauss_vector

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L405)

::: pybmad.ran_gauss_vector
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_seed_get

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L860)

::: pybmad.ran_seed_get
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_seed_put

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L780)

::: pybmad.ran_seed_put
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_uniform

Fortran sources (overloaded):

- `ran_uniform_scalar`: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L920)
- `ran_uniform_vector`: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L984)

::: pybmad.ran_uniform
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rcelbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L1185)

::: pybmad.rcelbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_a_line

Fortran source: [`sim_utils/io/input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/input_mod.f90#L144)

::: pybmad.read_a_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### readline_read_history

Fortran source: [`sim_utils/io/input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/input_mod.f90#L267)

::: pybmad.readline_read_history
    options:
      show_root_heading: false
      show_root_toc_entry: false

### readline_write_history

Fortran source: [`sim_utils/io/input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/io/input_mod.f90#L291)

::: pybmad.readline_write_history
    options:
      show_root_heading: false
      show_root_toc_entry: false

### real_num_fortran_format

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L745)

::: pybmad.real_num_fortran_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### real_path

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L779)

::: pybmad.real_path
    options:
      show_root_heading: false
      show_root_toc_entry: false

### real_str

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L785)

::: pybmad.real_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### real_to_string

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L704)

::: pybmad.real_to_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_spline

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/spline_mod.f90#L50)

::: pybmad.reallocate_spline
    options:
      show_root_heading: false
      show_root_toc_entry: false

### relbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L1042)

::: pybmad.relbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### relcbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L1135)

::: pybmad.relcbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### relsbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L1090)

::: pybmad.relsbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rgelbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L206)

::: pybmad.rgelbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rms_value

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L793)

::: pybmad.rms_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rot_2d

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L801)

::: pybmad.rot_2d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_vec

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L737)

::: pybmad.rotate_vec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_vec_given_axis_angle

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L689)

::: pybmad.rotate_vec_given_axis_angle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rp8

Fortran source: [`sim_utils/interfaces/precision_def.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/precision_def.f90#L42)

::: pybmad.rp8
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rserbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L1377)

::: pybmad.rserbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### run_timer

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L807)

::: pybmad.run_timer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### serbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L923)

::: pybmad.serbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_env

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L754)

::: pybmad.set_env
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_parameter

Fortran sources (overloaded):

- `set_parameter_real`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L6)
- `set_parameter_int`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L11)
- `set_parameter_logic`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L16)

::: pybmad.set_parameter
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_species_charge

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1700)

::: pybmad.set_species_charge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sign_of

Fortran sources (overloaded):

- `sign_of_real`: [`sim_utils/special_functions/sign_of_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/sign_of_mod.f90#L41)
- `sign_of_int`: [`sim_utils/special_functions/sign_of_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/sign_of_mod.f90#L72)

::: pybmad.sign_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sinc

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L814)

::: pybmad.sinc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sincc

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L822)

::: pybmad.sincc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sinhx_x

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L830)

::: pybmad.sinhx_x
    options:
      show_root_heading: false
      show_root_toc_entry: false

### skip_header

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L838)

::: pybmad.skip_header
    options:
      show_root_heading: false
      show_root_toc_entry: false

### special_projection

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/naff.f90#L215)

::: pybmad.special_projection
    options:
      show_root_heading: false
      show_root_toc_entry: false

### species_id

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L966)

::: pybmad.species_id
    options:
      show_root_heading: false
      show_root_toc_entry: false

### species_id_from_openpmd

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1356)

::: pybmad.species_id_from_openpmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### species_name

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1246)

::: pybmad.species_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### species_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L929)

::: pybmad.species_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1473)

::: pybmad.spin_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spline1

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/spline_mod.f90#L405)

::: pybmad.spline1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spline_akima

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/spline_mod.f90#L475)

::: pybmad.spline_akima
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spline_akima_interpolate

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/spline_mod.f90#L163)

::: pybmad.spline_akima_interpolate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spline_evaluate

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/spline_mod.f90#L281)

::: pybmad.spline_evaluate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sqrt_alpha

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L852)

::: pybmad.sqrt_alpha
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sqrt_one

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L844)

::: pybmad.sqrt_one
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_count

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L760)

::: pybmad.str_count
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_downcase

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L998)

::: pybmad.str_downcase
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_first_in_set

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L859)

::: pybmad.str_first_in_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_first_not_in_set

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L867)

::: pybmad.str_first_not_in_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_last_in_set

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L874)

::: pybmad.str_last_in_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_last_not_in_set

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L881)

::: pybmad.str_last_not_in_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_match_wild

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L986)

::: pybmad.str_match_wild
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_substitute

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L979)

::: pybmad.str_substitute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_upcase

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L992)

::: pybmad.str_upcase
    options:
      show_root_heading: false
      show_root_toc_entry: false

### string_to_int

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L888)

::: pybmad.string_to_int
    options:
      show_root_heading: false
      show_root_toc_entry: false

### string_to_real

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L897)

::: pybmad.string_to_real
    options:
      show_root_heading: false
      show_root_toc_entry: false

### string_trim

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L1010)

::: pybmad.string_trim
    options:
      show_root_heading: false
      show_root_toc_entry: false

### string_trim2

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L907)

::: pybmad.string_trim2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### suggest_lmdif

Fortran source: [`sim_utils/optimizers/lmdif_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/optimizers/lmdif_mod.f90#L121)

::: pybmad.suggest_lmdif
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_bicubic_coef

Fortran source: [`sim_utils/math/super_recipes_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/super_recipes_mod.f90#L112)

::: pybmad.super_bicubic_coef
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_bicubic_interpolation

Fortran source: [`sim_utils/math/super_recipes_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/super_recipes_mod.f90#L58)

::: pybmad.super_bicubic_interpolation
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_polint

Fortran source: [`sim_utils/math/super_recipes_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/super_recipes_mod.f90#L1506)

::: pybmad.super_polint
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_poly

Fortran source: [`sim_utils/math/super_recipes_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/super_recipes_mod.f90#L1739)

::: pybmad.super_poly
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_sobseq

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L1019)

::: pybmad.super_sobseq
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_sort

Fortran source: [`sim_utils/math/super_recipes_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/super_recipes_mod.f90#L159)

::: pybmad.super_sort
    options:
      show_root_heading: false
      show_root_toc_entry: false

### system_command

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L1004)

::: pybmad.system_command
    options:
      show_root_heading: false
      show_root_toc_entry: false

### test_xgelbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/special_functions/elliptic_integral_mod.f90#L102)

::: pybmad.test_xgelbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_str

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L930)

::: pybmad.to_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tricubic_cmplx_eval

Fortran source: [`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/cubic_interpolation_mod.f90#L1433)

::: pybmad.tricubic_cmplx_eval
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_this_file

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L938)

::: pybmad.type_this_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### upcase_string

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L949)

::: pybmad.upcase_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### virtual_memory_usage

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L1023)

::: pybmad.virtual_memory_usage
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_to_axis_angle

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L102)

::: pybmad.w_mat_to_axis_angle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_to_quat

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/rotation_3d_mod.f90#L129)

::: pybmad.w_mat_to_quat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### word_len

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L954)

::: pybmad.word_len
    options:
      show_root_heading: false
      show_root_toc_entry: false

### word_read

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/sim_utils_interface.f90#L960)

::: pybmad.word_read
    options:
      show_root_heading: false
      show_root_toc_entry: false

### x0_radiation_length

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/interfaces/particle_species_mod.f90#L1739)

::: pybmad.x0_radiation_length
    options:
      show_root_heading: false
      show_root_toc_entry: false

### zig_table_init

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4bb23ee832b1d289c45d4d295e48cb96a52a6a6d/sim_utils/math/random_mod.f90#L136)

::: pybmad.zig_table_init
    options:
      show_root_heading: false
      show_root_toc_entry: false
