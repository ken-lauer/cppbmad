# Sim Utils

Simulation utility helpers from the sim_utils library.

## Classes (Fortran Structures)

::: pybmad.AllPointerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### AllPointerStruct

Fortran struct: `all_pointer_struct` ([`sim_utils/interfaces/sim_utils_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_struct.f90#L45))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `r` | float |  |
| `q` | float |  |
| `i` | int |  |
| `l` | bool |  |
| `r1` | 1D array of float |  |
| `i1` | 1D array of int |  |

::: pybmad.BicubicCmplxCoefStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BicubicCmplxCoefStruct

Fortran struct: `bicubic_cmplx_coef_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L70))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | 2D array of complex (shape: 0:3,0:3) | Coefs |
| `i_box` | 1D array of int (shape: 2) | index at lower box corner. |

::: pybmad.BicubicCoefStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BicubicCoefStruct

Fortran struct: `bicubic_coef_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L43))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | 2D array of float (shape: 0:3,0:3) | Coefs |
| `i_box` | 1D array of int (shape: 2) | index at lower box corner. |

::: pybmad.BinStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BinStruct

Fortran struct: `bin_struct` ([`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L9))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `count` | 1D array of float | Counts (or weight) in each bin |
| `min` | float | Bounds for the bins |
| `max` | float |  |
| `delta` | float | Size of a bin |
| `n` | int | Number of bins |

::: pybmad.CmplxField1At2DPtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CmplxField1At2DPtStruct

Fortran struct: `cmplx_field1_at_2D_pt_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L55))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `f` | complex | Field |
| `df_dx` | complex | Normalized field 1st derivatives |
| `df_dy` | complex | Normalized field 1st derivatives |
| `d2f_dxdy` | complex | Normalized field 2nd derivative |

::: pybmad.CmplxField1At3DPtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CmplxField1At3DPtStruct

Fortran struct: `cmplx_field1_at_3D_pt_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L110))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `f` | complex | Field |
| `df_dx` | complex | Normalized field 1st derivatives |
| `df_dy` | complex | Normalized field 1st derivatives |
| `df_dz` | complex | Normalized field 1st derivatives |
| `d2f_dxdy` | complex | Normalized field 2nd derivatives |
| `d2f_dxdz` | complex | Normalized field 2nd derivatives |
| `d2f_dydz` | complex | Normalized field 2nd derivatives |
| `d3f_dxdydz` | complex | Normalized field 3rd derivative |

::: pybmad.CmplxFieldAt2DBoxStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CmplxFieldAt2DBoxStruct

Fortran struct: `cmplx_field_at_2D_box_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L63))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `pt` | 2D array of CmplxField1At2dPtStruct (shape: 0:1, 0:1) |  |
| `i_box` | 1D array of int (shape: 2) | index at lower box corner. |

::: pybmad.CmplxFieldAt3DBoxStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CmplxFieldAt3DBoxStruct

Fortran struct: `cmplx_field_at_3D_box_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L119))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `pt` | 3D array of CmplxField1At3dPtStruct (shape: 0:1, 0:1, 0:1) |  |
| `i_box` | 1D array of int (shape: 3) | index at lower box corner. |

::: pybmad.Field1At2DPtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Field1At2DPtStruct

Fortran struct: `field1_at_2D_pt_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L28))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `f` | float | Field |
| `df_dx` | float | Normalized field 1st derivatives |
| `df_dy` | float | Normalized field 1st derivatives |
| `d2f_dxdy` | float | Normalized field 2nd derivative |

::: pybmad.Field1At3DPtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Field1At3DPtStruct

Fortran struct: `field1_at_3D_pt_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L82))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `f` | float | Field |
| `df_dx` | float | Normalized field 1st derivatives |
| `df_dy` | float | Normalized field 1st derivatives |
| `df_dz` | float | Normalized field 1st derivatives |
| `d2f_dxdy` | float | Normalized field 2nd derivatives |
| `d2f_dxdz` | float | Normalized field 2nd derivatives |
| `d2f_dydz` | float | Normalized field 2nd derivatives |
| `d3f_dxdydz` | float | Normalized field 3rd derivative |

::: pybmad.FieldAt2DBoxStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### FieldAt2DBoxStruct

Fortran struct: `field_at_2D_box_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L36))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `pt` | 2D array of Field1At2dPtStruct (shape: 0:1, 0:1) |  |
| `i_box` | 1D array of int (shape: 2) | index at lower box corner. |

::: pybmad.FieldAt3DBoxStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### FieldAt3DBoxStruct

Fortran struct: `field_at_3D_box_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L91))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `pt` | 3D array of Field1At3dPtStruct (shape: 0:1, 0:1, 0:1) |  |
| `i_box` | 1D array of int (shape: 3) | index at lower box corner. |

::: pybmad.GeneralBinStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GeneralBinStruct

Fortran struct: `general_bin_struct` ([`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L17))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `count` | 1D array of float | Counts (or weight) in each bin |
| `min` | 1D array of float (shape: 3) | Bounds for the bins |
| `max` | 1D array of float (shape: 3) |  |
| `delta` | 1D array of float (shape: 3) | Size of a bin |
| `dim` | int | Number of dimensions |
| `n` | 1D array of int (shape: 3) | number of bins in each dimension |

::: pybmad.MolecularComponentStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MolecularComponentStruct

Fortran struct: `molecular_component_struct` ([`sim_utils/interfaces/sim_utils_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_struct.f90#L78))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `atom` | str |  |
| `number` | int |  |

::: pybmad.NamedNumberStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### NamedNumberStruct

Fortran struct: `named_number_struct` ([`sim_utils/interfaces/precision_def.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/precision_def.f90#L18))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `value` | float |  |

::: pybmad.NametableStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### NametableStruct

Fortran struct: `nametable_struct` ([`sim_utils/interfaces/sim_utils_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_struct.f90#L34))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | 1D array of str | Array of names. |
| `index` | 1D array of int | Sorted index for names(:) array. names(an_index(i)) is in alphabetical order. |
| `n_min` | int | Set to 0 for use in a lattice. |
| `n_max` | int | Use only names(n_min:n_max) part of array. |

::: pybmad.OutIoOutputDirectStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### OutIoOutputDirectStruct

Fortran struct: `out_io_output_direct_struct` ([`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L28))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `print_and_capture` | 1D array of bool (shape: -1:10) |  |
| `file_unit` | 1D array of int (shape: -1:10) |  |

::: pybmad.QpAxisStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### QpAxisStruct

Fortran struct: `qp_axis_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/plot/quick_plot_struct.f90#L58))

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

Fortran struct: `qp_legend_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/plot/quick_plot_struct.f90#L138))

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

Fortran struct: `qp_line_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/plot/quick_plot_struct.f90#L116))

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

Fortran struct: `qp_point_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/plot/quick_plot_struct.f90#L100))

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

Fortran struct: `qp_rect_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/plot/quick_plot_struct.f90#L105))

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

Fortran struct: `qp_symbol_struct` ([`sim_utils/plot/quick_plot_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/plot/quick_plot_struct.f90#L122))

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

Fortran struct: `random_state_struct` ([`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L28))

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

Fortran struct: `spline_struct` ([`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/spline_mod.f90#L18))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x0` | float | Point at start of spline |
| `y0` | float | Point at start of spline |
| `x1` | float | Point at end of spline |
| `coef` | 1D array of float (shape: 0:3) | coefficients for cubic spline |

::: pybmad.StrIndexStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### StrIndexStruct

Fortran struct: `str_index_struct` ([`sim_utils/interfaces/sim_utils_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_struct.f90#L18))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | [1D array of VarLengthStringStruct](sim_utils.md#varlengthstringstruct) | Array of names. |
| `index` | 1D array of int | Sorted index for names(:) array. names(an_index(i)) is in alphabetical order. |
| `n_min` | int |  |
| `n_max` | int | Use only names(n_min:n_max) part of array. |

::: pybmad.TricubicCmplxCoefStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TricubicCmplxCoefStruct

Fortran struct: `tricubic_cmplx_coef_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L126))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | 3D array of complex (shape: 0:3,0:3,0:3) | Coefs |
| `i_box` | 1D array of int (shape: 3) | index at lower box corner. |

::: pybmad.TricubicCoefStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TricubicCoefStruct

Fortran struct: `tricubic_coef_struct` ([`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L98))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | 3D array of float (shape: 0:3,0:3,0:3) | Coefs |
| `i_box` | 1D array of int (shape: 3) | index at lower box corner. |

::: pybmad.VarLengthStringStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### VarLengthStringStruct

Fortran struct: `var_length_string_struct` ([`sim_utils/interfaces/sim_utils_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_struct.f90#L14))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `str` | str |  |

## Procedures

### all_pointer_to_string

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L731)

::: pybmad.simutils.all_pointer_to_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### allocate_thread_states

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L1152)

::: pybmad.simutils.allocate_thread_states
    options:
      show_root_heading: false
      show_root_toc_entry: false

### anomalous_moment_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1423)

::: pybmad.simutils.anomalous_moment_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### antiparticle

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L882)

::: pybmad.simutils.antiparticle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apfft

Fortran source: [`sim_utils/math/all_phase_fft.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/all_phase_fft.f90#L135)

::: pybmad.simutils.apfft
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apfft_corr

Fortran source: [`sim_utils/math/all_phase_fft.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/all_phase_fft.f90#L39)

::: pybmad.simutils.apfft_corr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apfft_ext

Fortran source: [`sim_utils/math/all_phase_fft.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/all_phase_fft.f90#L174)

::: pybmad.simutils.apfft_ext
    options:
      show_root_heading: false
      show_root_toc_entry: false

### asinc

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L35)

::: pybmad.simutils.asinc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### assert_equal

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L43)

::: pybmad.simutils.assert_equal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### atomic_number

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1766)

::: pybmad.simutils.atomic_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### atomic_species_id

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1213)

::: pybmad.simutils.atomic_species_id
    options:
      show_root_heading: false
      show_root_toc_entry: false

### axis_angle_to_quat

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L385)

::: pybmad.simutils.axis_angle_to_quat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### axis_angle_to_w_mat

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L234)

::: pybmad.simutils.axis_angle_to_w_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bicubic_cmplx_eval

Fortran source: [`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L1103)

::: pybmad.simutils.bicubic_cmplx_eval
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bicubic_eval

Fortran source: [`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L472)

::: pybmad.simutils.bicubic_eval
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bicubic_interpolation_cmplx_coefs

Fortran source: [`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L1049)

::: pybmad.simutils.bicubic_interpolation_cmplx_coefs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bicubic_interpolation_coefs

Fortran source: [`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L418)

::: pybmad.simutils.bicubic_interpolation_coefs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bin_2d

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L254)

::: pybmad.simutils.bin_2d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bin_data

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L107)

::: pybmad.simutils.bin_data
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bin_data_density

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L46)

::: pybmad.simutils.bin_data_density
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bin_data_density_2d

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L162)

::: pybmad.simutils.bin_data_density_2d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bin_index

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L329)

::: pybmad.simutils.bin_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bin_x_center

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L353)

::: pybmad.simutils.bin_x_center
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bit_set

Fortran source: [`sim_utils/math/bit_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bit_set_mod.f90#L22)

::: pybmad.simutils.bit_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bracket_index_for_spline

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/spline_mod.f90#L327)

::: pybmad.simutils.bracket_index_for_spline
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_file_number

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L80)

::: pybmad.simutils.calc_file_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### celbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L420)

::: pybmad.simutils.celbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cesr_getarg

Fortran source: [`sim_utils/io/command_line_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/command_line_mod.f90#L57)

::: pybmad.simutils.cesr_getarg
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cesr_iargc

Fortran source: [`sim_utils/io/command_line_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/command_line_mod.f90#L23)

::: pybmad.simutils.cesr_iargc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### change_file_number

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L88)

::: pybmad.simutils.change_file_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### charge_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1498)

::: pybmad.simutils.charge_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### charge_to_mass_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1662)

::: pybmad.simutils.charge_to_mass_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coarse_frequency_estimate

Fortran source: [`sim_utils/math/fourier_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/fourier_mod.f90#L29)

::: pybmad.simutils.coarse_frequency_estimate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_error_function

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L114)

::: pybmad.simutils.complex_error_function
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cos_one

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L94)

::: pybmad.simutils.cos_one
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cosc

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L27)

::: pybmad.simutils.cosc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### count_at_index

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L369)

::: pybmad.simutils.count_at_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_a_spline

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/spline_mod.f90#L105)

::: pybmad.simutils.create_a_spline
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cross_product

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L123)

::: pybmad.simutils.cross_product
    options:
      show_root_heading: false
      show_root_toc_entry: false

### date_and_time_stamp

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L130)

::: pybmad.simutils.date_and_time_stamp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### destfixedwindowls

Fortran source: [`sim_utils/math/windowLS.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/windowLS.f90#L145)

::: pybmad.simutils.destfixedwindowls
    options:
      show_root_heading: false
      show_root_toc_entry: false

### detab

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L136)

::: pybmad.simutils.detab
    options:
      show_root_heading: false
      show_root_toc_entry: false

### display_size_and_resolution

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L148)

::: pybmad.simutils.display_size_and_resolution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### dj_bessel

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L155)

::: pybmad.simutils.dj_bessel
    options:
      show_root_heading: false
      show_root_toc_entry: false

### djb_hash

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L162)

::: pybmad.simutils.djb_hash
    options:
      show_root_heading: false
      show_root_toc_entry: false

### djb_str_hash

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L169)

::: pybmad.simutils.djb_str_hash
    options:
      show_root_heading: false
      show_root_toc_entry: false

### downcase_string

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L181)

::: pybmad.simutils.downcase_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### elbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L271)

::: pybmad.simutils.elbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### elcbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L371)

::: pybmad.simutils.elcbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ellipinc

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L53)

::: pybmad.simutils.ellipinc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### elsbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L325)

::: pybmad.simutils.elsbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### end_akima_spline_calc

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/spline_mod.f90#L573)

::: pybmad.simutils.end_akima_spline_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### err_exit

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L192)

::: pybmad.simutils.err_exit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### factorial

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L197)

::: pybmad.simutils.factorial
    options:
      show_root_heading: false
      show_root_toc_entry: false

### faddeeva_function

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L204)

::: pybmad.simutils.faddeeva_function
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fft_1d

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L217)

::: pybmad.simutils.fft_1d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### file_directorizer

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L223)

::: pybmad.simutils.file_directorizer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### file_get

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L231)

::: pybmad.simutils.file_get
    options:
      show_root_heading: false
      show_root_toc_entry: false

### file_get_open

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L238)

::: pybmad.simutils.file_get_open
    options:
      show_root_heading: false
      show_root_toc_entry: false

### file_suffixer

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L247)

::: pybmad.simutils.file_suffixer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### find_location

Fortran sources (overloaded):

- `find_location_real`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L1033)
- `find_location_int`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L1039)
- `find_location_logic`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L1044)
- `find_location_str`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L1049)

::: pybmad.simutils.find_location
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fine_frequency_estimate

Fortran source: [`sim_utils/math/fourier_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/fourier_mod.f90#L113)

::: pybmad.simutils.fine_frequency_estimate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fixedwindowls

Fortran source: [`sim_utils/math/windowLS.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/windowLS.f90#L171)

::: pybmad.simutils.fixedwindowls
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fourier_amplitude

Fortran source: [`sim_utils/math/fourier_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/fourier_mod.f90#L207)

::: pybmad.simutils.fourier_amplitude
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gelbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L140)

::: pybmad.simutils.gelbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gen_complete_elliptic

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L255)

::: pybmad.simutils.gen_complete_elliptic
    options:
      show_root_heading: false
      show_root_toc_entry: false

### general_bin_count

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L414)

::: pybmad.simutils.general_bin_count
    options:
      show_root_heading: false
      show_root_toc_entry: false

### general_bin_index

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L399)

::: pybmad.simutils.general_bin_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### general_bin_index_in_bounds

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L438)

::: pybmad.simutils.general_bin_index_in_bounds
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_a_char

Fortran source: [`sim_utils/io/input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/input_mod.f90#L85)

::: pybmad.simutils.get_a_char
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_file_number

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L262)

::: pybmad.simutils.get_file_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_file_time_stamp

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L278)

::: pybmad.simutils.get_file_time_stamp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_tty_char

Fortran source: [`sim_utils/io/input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/input_mod.f90#L36)

::: pybmad.simutils.get_tty_char
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hanhan

Fortran source: [`sim_utils/math/all_phase_fft.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/all_phase_fft.f90#L277)

::: pybmad.simutils.hanhan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### i_bessel

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L283)

::: pybmad.simutils.i_bessel
    options:
      show_root_heading: false
      show_root_toc_entry: false

### i_bessel_extended

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L290)

::: pybmad.simutils.i_bessel_extended
    options:
      show_root_heading: false
      show_root_toc_entry: false

### increment_file_number

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L298)

::: pybmad.simutils.increment_file_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### index_nocase

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L306)

::: pybmad.simutils.index_nocase
    options:
      show_root_heading: false
      show_root_toc_entry: false

### initfixedwindowls

Fortran source: [`sim_utils/math/windowLS.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/windowLS.f90#L76)

::: pybmad.simutils.initfixedwindowls
    options:
      show_root_heading: false
      show_root_toc_entry: false

### initial_lmdif

Fortran source: [`sim_utils/optimizers/lmdif_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/optimizers/lmdif_mod.f90#L81)

::: pybmad.simutils.initial_lmdif
    options:
      show_root_heading: false
      show_root_toc_entry: false

### int_str

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L313)

::: pybmad.simutils.int_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### interpolated_fft

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/naff.f90#L328)

::: pybmad.simutils.interpolated_fft
    options:
      show_root_heading: false
      show_root_toc_entry: false

### interpolated_fft_gsl

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/naff.f90#L233)

::: pybmad.simutils.interpolated_fft_gsl
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_alphabetic

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L341)

::: pybmad.simutils.is_alphabetic
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_decreasing_sequence

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L348)

::: pybmad.simutils.is_decreasing_sequence
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_false

Fortran source: [`sim_utils/interfaces/sim_utils_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_struct.f90#L211)

::: pybmad.simutils.is_false
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_increasing_sequence

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L356)

::: pybmad.simutils.is_increasing_sequence
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_integer

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L364)

::: pybmad.simutils.is_integer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_logical

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L372)

::: pybmad.simutils.is_logical
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_real

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L379)

::: pybmad.simutils.is_real
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_subatomic_species

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1799)

::: pybmad.simutils.is_subatomic_species
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_true

Fortran source: [`sim_utils/interfaces/sim_utils_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_struct.f90#L178)

::: pybmad.simutils.is_true
    options:
      show_root_heading: false
      show_root_toc_entry: false

### j_bessel

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L396)

::: pybmad.simutils.j_bessel
    options:
      show_root_heading: false
      show_root_toc_entry: false

### linear_fit

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L403)

::: pybmad.simutils.linear_fit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### linear_fit_2d

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L415)

::: pybmad.simutils.linear_fit_2d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### logic_str

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L432)

::: pybmad.simutils.logic_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lunget

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L438)

::: pybmad.simutils.lunget
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_legal_comment

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L454)

::: pybmad.simutils.make_legal_comment
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mass_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1550)

::: pybmad.simutils.mass_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### match_reg

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L443)

::: pybmad.simutils.match_reg
    options:
      show_root_heading: false
      show_root_toc_entry: false

### match_wild

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L460)

::: pybmad.simutils.match_wild
    options:
      show_root_heading: false
      show_root_toc_entry: false

### match_word

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L551)

::: pybmad.simutils.match_word
    options:
      show_root_heading: false
      show_root_toc_entry: false

### maximize_projection

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/naff.f90#L164)

::: pybmad.simutils.maximize_projection
    options:
      show_root_heading: false
      show_root_toc_entry: false

### milli_sleep

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L449)

::: pybmad.simutils.milli_sleep
    options:
      show_root_heading: false
      show_root_toc_entry: false

### modulo2_dp

Fortran source: [`sim_utils/special_functions/modulo2_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/modulo2_mod.f90#L85)

::: pybmad.simutils.modulo2_dp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### modulo2_int

Fortran source: [`sim_utils/special_functions/modulo2_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/modulo2_mod.f90#L157)

::: pybmad.simutils.modulo2_int
    options:
      show_root_heading: false
      show_root_toc_entry: false

### modulo2_qp

Fortran source: [`sim_utils/special_functions/modulo2_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/modulo2_mod.f90#L121)

::: pybmad.simutils.modulo2_qp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### modulo2_sp

Fortran source: [`sim_utils/special_functions/modulo2_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/modulo2_mod.f90#L49)

::: pybmad.simutils.modulo2_sp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### molecular_components

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L812)

::: pybmad.simutils.molecular_components
    options:
      show_root_heading: false
      show_root_toc_entry: false

### n_bins_automatic

Fortran source: [`sim_utils/math/bin_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/bin_mod.f90#L386)

::: pybmad.simutils.n_bins_automatic
    options:
      show_root_heading: false
      show_root_toc_entry: false

### n_choose_k

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L568)

::: pybmad.simutils.n_choose_k
    options:
      show_root_heading: false
      show_root_toc_entry: false

### n_spline_create

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L575)

::: pybmad.simutils.n_spline_create
    options:
      show_root_heading: false
      show_root_toc_entry: false

### naff

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/naff.f90#L72)

::: pybmad.simutils.naff
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nametable_add

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L583)

::: pybmad.simutils.nametable_add
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nametable_bracket_indexx

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L591)

::: pybmad.simutils.nametable_bracket_indexx
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nametable_change1

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L600)

::: pybmad.simutils.nametable_change1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nametable_init

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L608)

::: pybmad.simutils.nametable_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nametable_remove

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L615)

::: pybmad.simutils.nametable_remove
    options:
      show_root_heading: false
      show_root_toc_entry: false

### negative_ampsquared

Fortran source: [`sim_utils/math/fourier_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/fourier_mod.f90#L153)

::: pybmad.simutils.negative_ampsquared
    options:
      show_root_heading: false
      show_root_toc_entry: false

### negative_dampsquared

Fortran source: [`sim_utils/math/fourier_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/fourier_mod.f90#L173)

::: pybmad.simutils.negative_dampsquared
    options:
      show_root_heading: false
      show_root_toc_entry: false

### omega_to_quat

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L314)

::: pybmad.simutils.omega_to_quat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### openpmd_species_name

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1392)

::: pybmad.simutils.openpmd_species_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ordinal_str

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L633)

::: pybmad.simutils.ordinal_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### out_io

Fortran sources (overloaded):

- `out_io_real`: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L249)
- `out_io_int`: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L293)
- `out_io_logical`: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L337)
- `out_io_line12`: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L382)
- `out_io_lines`: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L434)

::: pybmad.simutils.out_io
    options:
      show_root_heading: false
      show_root_toc_entry: false

### out_io_buffer_get_line

Fortran source: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L842)

::: pybmad.simutils.out_io_buffer_get_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### out_io_buffer_num_lines

Fortran source: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L822)

::: pybmad.simutils.out_io_buffer_num_lines
    options:
      show_root_heading: false
      show_root_toc_entry: false

### out_io_buffer_reset

Fortran source: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L799)

::: pybmad.simutils.out_io_buffer_reset
    options:
      show_root_heading: false
      show_root_toc_entry: false

### out_io_print_and_capture_setup

Fortran source: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L768)

::: pybmad.simutils.out_io_print_and_capture_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### output_direct

Fortran source: [`sim_utils/io/output_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/output_mod.f90#L161)

::: pybmad.simutils.output_direct
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_fortran_format

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L646)

::: pybmad.simutils.parse_fortran_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_locations

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L652)

::: pybmad.simutils.pointer_to_locations
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_ran_state

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L1115)

::: pybmad.simutils.pointer_to_ran_state
    options:
      show_root_heading: false
      show_root_toc_entry: false

### poly_eval

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L661)

::: pybmad.simutils.poly_eval
    options:
      show_root_heading: false
      show_root_toc_entry: false

### probability_funct

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L668)

::: pybmad.simutils.probability_funct
    options:
      show_root_heading: false
      show_root_toc_entry: false

### projdd

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/naff.f90#L132)

::: pybmad.simutils.projdd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quadratic_roots

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L675)

::: pybmad.simutils.quadratic_roots
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_conj

Fortran sources (overloaded):

- `quat_conj_real`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L411)
- `quat_conj_complex`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L435)

::: pybmad.simutils.quat_conj
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_inverse

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L458)

::: pybmad.simutils.quat_inverse
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_mul

Fortran sources (overloaded):

- `quat_mul_real`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L484)
- `quat_mul_complex`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L551)

::: pybmad.simutils.quat_mul
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_rotate

Fortran sources (overloaded):

- `quat_rotate_real`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L616)
- `quat_rotate_complex`: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L652)

::: pybmad.simutils.quat_rotate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_to_axis_angle

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L347)

::: pybmad.simutils.quat_to_axis_angle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_to_omega

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L281)

::: pybmad.simutils.quat_to_omega
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quat_to_w_mat

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L184)

::: pybmad.simutils.quat_to_w_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### query_string

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L684)

::: pybmad.simutils.query_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quote

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L693)

::: pybmad.simutils.quote
    options:
      show_root_heading: false
      show_root_toc_entry: false

### quoten

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L698)

::: pybmad.simutils.quoten
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_default_state

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L892)

::: pybmad.simutils.ran_default_state
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_engine

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L623)

::: pybmad.simutils.ran_engine
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_gauss_converter

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L701)

::: pybmad.simutils.ran_gauss_converter
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_gauss_scalar

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L172)

::: pybmad.simutils.ran_gauss_scalar
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_gauss_vector

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L406)

::: pybmad.simutils.ran_gauss_vector
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_seed_get

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L861)

::: pybmad.simutils.ran_seed_get
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_seed_put

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L781)

::: pybmad.simutils.ran_seed_put
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ran_uniform

Fortran sources (overloaded):

- `ran_uniform_scalar`: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L921)
- `ran_uniform_vector`: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L985)

::: pybmad.simutils.ran_uniform
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rcelbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L1185)

::: pybmad.simutils.rcelbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_a_line

Fortran source: [`sim_utils/io/input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/input_mod.f90#L144)

::: pybmad.simutils.read_a_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### readline_read_history

Fortran source: [`sim_utils/io/input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/input_mod.f90#L267)

::: pybmad.simutils.readline_read_history
    options:
      show_root_heading: false
      show_root_toc_entry: false

### readline_write_history

Fortran source: [`sim_utils/io/input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/io/input_mod.f90#L291)

::: pybmad.simutils.readline_write_history
    options:
      show_root_heading: false
      show_root_toc_entry: false

### real_num_fortran_format

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L745)

::: pybmad.simutils.real_num_fortran_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### real_path

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L779)

::: pybmad.simutils.real_path
    options:
      show_root_heading: false
      show_root_toc_entry: false

### real_str

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L785)

::: pybmad.simutils.real_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### real_to_string

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L704)

::: pybmad.simutils.real_to_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_spline

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/spline_mod.f90#L50)

::: pybmad.simutils.reallocate_spline
    options:
      show_root_heading: false
      show_root_toc_entry: false

### relbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L1042)

::: pybmad.simutils.relbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### relcbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L1135)

::: pybmad.simutils.relcbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### relsbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L1090)

::: pybmad.simutils.relsbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rgelbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L206)

::: pybmad.simutils.rgelbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rms_value

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L793)

::: pybmad.simutils.rms_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rot_2d

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L801)

::: pybmad.simutils.rot_2d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_vec

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L737)

::: pybmad.simutils.rotate_vec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_vec_given_axis_angle

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L689)

::: pybmad.simutils.rotate_vec_given_axis_angle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rp8

Fortran source: [`sim_utils/interfaces/precision_def.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/precision_def.f90#L42)

::: pybmad.simutils.rp8
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rserbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L1377)

::: pybmad.simutils.rserbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### run_timer

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L807)

::: pybmad.simutils.run_timer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### serbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L923)

::: pybmad.simutils.serbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_all_ptr

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L970)

::: pybmad.simutils.set_all_ptr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_env

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L754)

::: pybmad.simutils.set_env
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_parameter

Fortran sources (overloaded):

- `set_parameter_real`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L6)
- `set_parameter_int`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L11)
- `set_parameter_logic`: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L16)

::: pybmad.simutils.set_parameter
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_species_charge

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1693)

::: pybmad.simutils.set_species_charge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sign_of

Fortran sources (overloaded):

- `sign_of_real`: [`sim_utils/special_functions/sign_of_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/sign_of_mod.f90#L41)
- `sign_of_int`: [`sim_utils/special_functions/sign_of_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/sign_of_mod.f90#L72)

::: pybmad.simutils.sign_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sinc

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L814)

::: pybmad.simutils.sinc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sincc

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L822)

::: pybmad.simutils.sincc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sinhx_x

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L830)

::: pybmad.simutils.sinhx_x
    options:
      show_root_heading: false
      show_root_toc_entry: false

### skip_header

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L838)

::: pybmad.simutils.skip_header
    options:
      show_root_heading: false
      show_root_toc_entry: false

### special_projection

Fortran source: [`sim_utils/math/naff.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/naff.f90#L215)

::: pybmad.simutils.special_projection
    options:
      show_root_heading: false
      show_root_toc_entry: false

### species_id

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L966)

::: pybmad.simutils.species_id
    options:
      show_root_heading: false
      show_root_toc_entry: false

### species_id_from_openpmd

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1356)

::: pybmad.simutils.species_id_from_openpmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### species_name

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1246)

::: pybmad.simutils.species_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### species_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L929)

::: pybmad.simutils.species_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_of

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1466)

::: pybmad.simutils.spin_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spline1

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/spline_mod.f90#L405)

::: pybmad.simutils.spline1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spline_akima

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/spline_mod.f90#L475)

::: pybmad.simutils.spline_akima
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spline_akima_interpolate

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/spline_mod.f90#L163)

::: pybmad.simutils.spline_akima_interpolate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spline_evaluate

Fortran source: [`sim_utils/math/spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/spline_mod.f90#L281)

::: pybmad.simutils.spline_evaluate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sqrt_alpha

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L852)

::: pybmad.simutils.sqrt_alpha
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sqrt_one

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L844)

::: pybmad.simutils.sqrt_one
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_count

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L760)

::: pybmad.simutils.str_count
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_downcase

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L998)

::: pybmad.simutils.str_downcase
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_first_in_set

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L859)

::: pybmad.simutils.str_first_in_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_first_not_in_set

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L867)

::: pybmad.simutils.str_first_not_in_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_last_in_set

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L874)

::: pybmad.simutils.str_last_in_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_last_not_in_set

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L881)

::: pybmad.simutils.str_last_not_in_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_match_wild

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L986)

::: pybmad.simutils.str_match_wild
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_substitute

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L979)

::: pybmad.simutils.str_substitute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### str_upcase

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L992)

::: pybmad.simutils.str_upcase
    options:
      show_root_heading: false
      show_root_toc_entry: false

### string_to_int

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L888)

::: pybmad.simutils.string_to_int
    options:
      show_root_heading: false
      show_root_toc_entry: false

### string_to_real

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L897)

::: pybmad.simutils.string_to_real
    options:
      show_root_heading: false
      show_root_toc_entry: false

### string_trim

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L1010)

::: pybmad.simutils.string_trim
    options:
      show_root_heading: false
      show_root_toc_entry: false

### string_trim2

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L907)

::: pybmad.simutils.string_trim2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### suggest_lmdif

Fortran source: [`sim_utils/optimizers/lmdif_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/optimizers/lmdif_mod.f90#L121)

::: pybmad.simutils.suggest_lmdif
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_bicubic_coef

Fortran source: [`sim_utils/math/super_recipes_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/super_recipes_mod.f90#L112)

::: pybmad.simutils.super_bicubic_coef
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_bicubic_interpolation

Fortran source: [`sim_utils/math/super_recipes_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/super_recipes_mod.f90#L58)

::: pybmad.simutils.super_bicubic_interpolation
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_polint

Fortran source: [`sim_utils/math/super_recipes_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/super_recipes_mod.f90#L1506)

::: pybmad.simutils.super_polint
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_poly

Fortran source: [`sim_utils/math/super_recipes_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/super_recipes_mod.f90#L1739)

::: pybmad.simutils.super_poly
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_sobseq

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L1020)

::: pybmad.simutils.super_sobseq
    options:
      show_root_heading: false
      show_root_toc_entry: false

### super_sort

Fortran source: [`sim_utils/math/super_recipes_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/super_recipes_mod.f90#L159)

::: pybmad.simutils.super_sort
    options:
      show_root_heading: false
      show_root_toc_entry: false

### system_command

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L1004)

::: pybmad.simutils.system_command
    options:
      show_root_heading: false
      show_root_toc_entry: false

### test_xgelbd

Fortran source: [`sim_utils/special_functions/elliptic_integral_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/special_functions/elliptic_integral_mod.f90#L102)

::: pybmad.simutils.test_xgelbd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_str

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L930)

::: pybmad.simutils.to_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tricubic_cmplx_eval

Fortran source: [`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L1433)

::: pybmad.simutils.tricubic_cmplx_eval
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tricubic_eval

Fortran source: [`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L800)

::: pybmad.simutils.tricubic_eval
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tricubic_interpolation_cmplx_coefs

Fortran source: [`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L1373)

::: pybmad.simutils.tricubic_interpolation_cmplx_coefs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tricubic_interpolation_coefs

Fortran source: [`sim_utils/math/cubic_interpolation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/cubic_interpolation_mod.f90#L740)

::: pybmad.simutils.tricubic_interpolation_coefs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_this_file

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L938)

::: pybmad.simutils.type_this_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### upcase_string

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L949)

::: pybmad.simutils.upcase_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### value_of_all_ptr

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L1016)

::: pybmad.simutils.value_of_all_ptr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### virtual_memory_usage

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L1023)

::: pybmad.simutils.virtual_memory_usage
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_to_axis_angle

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L102)

::: pybmad.simutils.w_mat_to_axis_angle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_to_quat

Fortran source: [`sim_utils/math/rotation_3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/rotation_3d_mod.f90#L129)

::: pybmad.simutils.w_mat_to_quat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### word_len

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L954)

::: pybmad.simutils.word_len
    options:
      show_root_heading: false
      show_root_toc_entry: false

### word_read

Fortran source: [`sim_utils/interfaces/sim_utils_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/sim_utils_interface.f90#L960)

::: pybmad.simutils.word_read
    options:
      show_root_heading: false
      show_root_toc_entry: false

### x0_radiation_length

Fortran source: [`sim_utils/interfaces/particle_species_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/interfaces/particle_species_mod.f90#L1732)

::: pybmad.simutils.x0_radiation_length
    options:
      show_root_heading: false
      show_root_toc_entry: false

### zig_table_init

Fortran source: [`sim_utils/math/random_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b41ae6ef92b61c52c3bd27486c07fec6df7dd1e9/sim_utils/math/random_mod.f90#L136)

::: pybmad.simutils.zig_table_init
    options:
      show_root_heading: false
      show_root_toc_entry: false
