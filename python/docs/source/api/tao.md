# Tao

Tao (The Tool for Accelerator Optics)

## Classes (Fortran Structures)

::: pybmad.TaoBeamBranchStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoBeamBranchStruct

Fortran struct: `tao_beam_branch_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L994))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `beam_at_start` | [BeamStruct](bmad.md#beamstruct) | Initial beam |
| `beam_init` | [BeamInitStruct](bmad.md#beaminitstruct) | User set beam distrubution at track start. |
| `beam_init_used` | [BeamInitStruct](bmad.md#beaminitstruct) | beam distribution with emit values set. |
| `init_starting_distribution` | bool | Init beam |
| `track_start` | str | Tracking start element. |
| `track_end` | str |  |
| `ix_branch` | int | Branch tracked. If track_start or track_end is a lord, ix_track_start/end index will be a index of slave. |
| `ix_track_start` | int | Element track start index. |
| `ix_track_end` | int | Element track end index |

::: pybmad.TaoBeamUniStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoBeamUniStruct

Fortran struct: `tao_beam_uni_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1016))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `saved_at` | str |  |
| `dump_file` | str |  |
| `dump_at` | str |  |
| `track_beam_in_universe` | bool | Beam tracking enabled in this universe? |
| `always_reinit` | bool |  |

::: pybmad.TaoBuildingWallOrientationStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoBuildingWallOrientationStruct

Fortran struct: `tao_building_wall_orientation_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L614))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `theta` | float |  |
| `x_offset` | float |  |
| `z_offset` | float |  |

::: pybmad.TaoBuildingWallPointStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoBuildingWallPointStruct

Fortran struct: `tao_building_wall_point_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L620))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `z` | float | Global floor position |
| `x` | float | Global floor position |
| `radius` | float | Arc radius. +r -> CW rotation, same as bends. |
| `z_center` | float | Arc center. |
| `x_center` | float | Arc center. |

::: pybmad.TaoBuildingWallSectionStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoBuildingWallSectionStruct

Fortran struct: `tao_building_wall_section_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L626))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `constraint` | str | "left_side" or "right_side" constraint. |
| `point` | [1D array of TaoBuildingWallPointStruct](tao.md#taobuildingwallpointstruct) |  |

::: pybmad.TaoBuildingWallStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoBuildingWallStruct

Fortran struct: `tao_building_wall_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L632))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `orientation` | [TaoBuildingWallOrientationStruct](tao.md#taobuildingwallorientationstruct) |  |
| `section` | [1D array of TaoBuildingWallSectionStruct](tao.md#taobuildingwallsectionstruct) |  |

::: pybmad.TaoCmdHistoryStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoCmdHistoryStruct

Fortran struct: `tao_cmd_history_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L62))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `cmd` | str | The command |
| `ix` | int | Command index (1st command has ix = 1, etc.) Note: Commands from command files will be assigned an index. |

::: pybmad.TaoCommonStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoCommonStruct

Fortran struct: `tao_common_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L753))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `plot_place_buffer` | [1D array of TaoPlotRegionStruct](tao.md#taoplotregionstruct) | Used when %external_plotting is on. |
| `covar` | 2D array of float |  |
| `alpha` | 2D array of float |  |
| `dummy_target` | float | Dummy varaible |
| `n_alias` | int |  |
| `cmd_file_level` | int | For nested command files. 0 -> no command file. |
| `ix_key_bank` | int | For single mode. |
| `ix_history` | int | Index to latest command in the history circular buffer. |
| `n_history` | int | Number of commands issued from beginning of starting Tao. |
| `lev_loop` | int | in do loop nest level |
| `n_err_messages_printed` | int | Used by tao_set_invalid to limit number of messages. |
| `n_universes` | int |  |
| `ix_beam_track_active_element` | int | Element being tracked through `tao_beam_track`. |
| `cmd_file_paused` | bool |  |
| `use_cmd_here` | bool | Used for commands recalled from the cmd history stack |
| `cmd_from_cmd_file` | bool | was command from a command file? |
| `use_saved_beam_in_tracking` | bool |  |
| `single_mode` | bool |  |
| `combine_consecutive_elements_of_like_name` | bool |  |
| `have_tracked_beam` | bool | Used to catch error when beam plotting without having tracked a beam. |
| `init_plot_needed` | bool | reinitialize plotting? |
| `init_beam` | bool | Used by custom programs to control Tao init |
| `init_var` | bool | Used by custom programs to control Tao init |
| `init_read_lat_info` | bool | Used by custom programs to control Tao init |
| `optimizer_running` | bool |  |
| `have_datums_using_expressions` | bool |  |
| `print_to_terminal` | bool | Print command prompt to the terminal? For use with GUIs. |
| `lattice_calc_done` | bool | Used by GUI for deciding when to refresh. |
| `add_measurement_noise` | bool | Turn off to take data derivatives. |
| `is_err_message_printed` | 1D array of bool (shape: 2) | Used by tao_set_invalid |
| `command_arg_has_been_executed` | bool | Has the -command command line argument been executed? |
| `all_merit_weights_positive` | bool |  |
| `multi_turn_orbit_is_plotted` | bool | Is a multi_turn_orbit being plotted? |
| `force_chrom_calc` | bool | Used by a routine to force a single chromaticity calculation. |
| `force_rad_int_calc` | bool | Used by a routine to force a single radiation integrals calculation |
| `rad_int_ri_calc_on` | bool | "Classical" radiation integrals calculation on/off. |
| `rad_int_6d_calc_on` | bool | 6D Radiation integrals calculation on/off. |
| `valid_plot_who` | 1D array of str (shape: 10) | model, base, ref etc... |
| `single_mode_buffer` | str |  |
| `cmd` | str | Used for the cmd history |

::: pybmad.TaoCurveColorStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoCurveColorStruct

Fortran struct: `tao_curve_color_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L182))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `data_type` | str | Datum type to use for z-axis. |
| `is_on` | bool | On/Off |
| `min` | float | Min and max values for mapping z-axis to color. |
| `max` | float | Min and max values for mapping z-axis to color. |
| `autoscale` | bool | Set %min, %max automatically to the limits of %data_type |

::: pybmad.TaoCurveOrbitStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoCurveOrbitStruct

Fortran struct: `tao_curve_orbit_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L174))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x` | float | Transverse offset |
| `y` | float | Transverse offset |
| `t` | float | Time |

::: pybmad.TaoCurveStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoCurveStruct

Fortran struct: `tao_curve_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L192))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Name identifying the curve. |
| `data_source` | str | 'lat', 'beam', 'data' (deprecated: 'dat'), 'var', 'multi_turn_orbit' |
| `data_index` | str | Used for calculating %ix_symb(:). |
| `data_type_x` | str | Used for data slices and phase space plots. |
| `data_type` | str | 'orbit.x', etc. |
| `ele_ref_name` | str | Reference element. |
| `legend_text` | str | String to draw in a curve legend. |
| `message_text` | str | Informational message to draw with graph. |
| `component` | str | Who to plot. Eg: 'meas - design' |
| `why_invalid` | str | Informative string to print. |
| `g` | [TaoGraphStruct](tao.md#taographstruct) | pointer to parent graph |
| `hist` | [TaoHistogramStruct](tao.md#taohistogramstruct) |  |
| `z_color` | [TaoCurveColorStruct](tao.md#taocurvecolorstruct) |  |
| `x_line` | 1D array of float | Coords for drawing a curve |
| `y_line` | 1D array of float |  |
| `y2_line` | 1D array of float | Second array needed for beam chamber curve. |
| `ix_line` | 1D array of int | Used by wave and aperture curves. |
| `x_symb` | 1D array of float | Coords for drawing the symbols |
| `y_symb` | 1D array of float |  |
| `z_symb` | 1D array of float | Symbol color |
| `err_symb` | 1D array of float | Error bars |
| `symb_size` | 1D array of float | Symbol size. Used with symbol_size_scale. |
| `ix_symb` | 1D array of int | Corresponding index in d1_data%d(:) array. |
| `y_axis_scale_factor` | float | y-axis conversion from internal to plotting units. |
| `line` | [QpLineStruct](sim_utils.md#qplinestruct) | Line attributes |
| `symbol` | [QpSymbolStruct](sim_utils.md#qpsymbolstruct) | Symbol attributes |
| `orbit` | [TaoCurveOrbitStruct](tao.md#taocurveorbitstruct) | Used for E/B field plotting. |
| `ix_universe` | int | Universe where data is. -1 => use s%global%default_universe |
| `symbol_every` | int | Symbol every how many points. |
| `ix_branch` | int |  |
| `ix_bunch` | int | Bunch to plot. |
| `n_turn` | int | Used for multi_turn_orbit plotting |
| `use_y2` | bool | Use y2 axis? |
| `draw_line` | bool | Draw a line through the data points? |
| `draw_symbols` | bool | Draw a symbol at the data points? |
| `draw_symbol_index` | bool | Draw the symbol index number curve%ix_symb? |
| `draw_error_bars` | bool | Draw error bars based upon data%error_rms if drawing data? !! logical :: draw_rms = .false.          ! Show mean and RMS values with legend? |
| `smooth_line_calc` | bool | Calculate data between element edge points? |
| `valid` | bool | valid data? |

::: pybmad.TaoD1DataStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoD1DataStruct

Fortran struct: `tao_d1_data_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L478))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Eg: 'x', etc. |
| `d2` | [TaoD2DataStruct](tao.md#taod2datastruct) | ptr to parent d2_data |
| `d` | [1D array of TaoDataStruct](tao.md#taodatastruct) | Pointer to the appropriate section in u%data |

::: pybmad.TaoD2DataStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoD2DataStruct

Fortran struct: `tao_d2_data_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L490))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Name to be used with commands. |
| `data_file_name` | str | Data file name . |
| `ref_file_name` | str | Reference file name. |
| `data_date` | str | Data measurement date. |
| `ref_date` | str | Reference data measurement date. |
| `descrip` | 1D array of str (shape: 10) | Array for descriptive information. |
| `d1` | [1D array of TaoD1DataStruct](tao.md#taod1datastruct) | Points to children |
| `ix_universe` | int | Index of universe this is in. |
| `ix_d2_data` | int | Index in u%d2_data(:) array. |
| `ix_ref` | int | Index of the reference data set. |
| `data_read_in` | bool | A data set has been read in? |
| `ref_read_in` | bool | A reference data set has been read in? |

::: pybmad.TaoDataStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoDataStruct

Fortran struct: `tao_data_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L424))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele_name` | str | Name of the lattice element where datum is evaluated. |
| `ele_start_name` | str | Name of starting lattice element when there is a range |
| `ele_ref_name` | str | Name of reference lattice element |
| `data_type` | str | Type of data: 'orbit.x', etc. |
| `merit_type` | str | Type of constraint: 'target', 'max', 'min', etc. |
| `id` | str | Used by Tao extension code. Not used by Tao directly. |
| `data_source` | str | 'lat', 'beam', 'data' or 'var'. Last two used for expressions. |
| `why_invalid` | str | Informational string if there is a problem. |
| `ix_uni` | int | Universe index of datum. |
| `ix_bunch` | int | Bunch number to get the data from. |
| `ix_branch` | int | Index of the associated lattice branch. |
| `ix_ele` | int | Index of the lattice element corresponding to ele_name |
| `ix_ele_start` | int | Index of lattice elment when there is a range |
| `ix_ele_ref` | int | Index of lattice elment when there is a reference. |
| `ix_ele_merit` | int | Index of lattice elment where merit is evaluated. |
| `ix_d1` | int | Index number in u%d2_data(i)%d1_data(j)%d(:) array. |
| `ix_data` | int | Index of this datum in the u%data(:) array of data_structs. |
| `ix_dModel` | int | Row number in the dModel_dVar derivative matrix. |
| `eval_point` | int | or anchor_center$, anchor_beginning$. Where to evaluate data relative to the element. |
| `meas_value` | float | Measured datum value. |
| `ref_value` | float | Measured datum value from the reference data set. |
| `model_value` | float | Datum value as calculated from the model. |
| `design_value` | float | What the datum value is in the design lattice. |
| `old_value` | float | The model_value at some previous time. |
| `base_value` | float | The value as calculated from the base model. |
| `error_rms` | float | Measurement error RMS. Used in plotting. |
| `delta_merit` | float | Diff used to calculate the merit function term. |
| `weight` | float | Weight for the merit function term. |
| `invalid_value` | float | Value used in merit calc if good_model = F (or possibly good_design & good_base). |
| `merit` | float | Merit function term value: weight * delta_merit^2 |
| `s` | float | longitudinal position of ele. |
| `s_offset` | float | Offset of the evaluation point. |
| `ref_s_offset` | float | Offset of the reference point. In development. |
| `err_message_printed` | bool | Used to prevent zillions of error messages being generated |
| `exists` | bool | See above |
| `good_model` | bool | See above |
| `good_base` | bool | See above |
| `good_design` | bool | See above |
| `good_meas` | bool | See above |
| `good_ref` | bool | See above |
| `good_user` | bool | See above |
| `good_opt` | bool | See above |
| `good_plot` | bool | See above |
| `useit_plot` | bool | See above |
| `useit_opt` | bool | See above |
| `spin_map` | [TaoSpinMapStruct](tao.md#taospinmapstruct) |  |
| `d1` | [TaoD1DataStruct](tao.md#taod1datastruct) | Pointer to the parent d1_data_struct |

::: pybmad.TaoDataVarComponentStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoDataVarComponentStruct

Fortran struct: `tao_data_var_component_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L160))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Eg: 'meas', 'ref', 'model', etc. |
| `sign` | float | +1 or -1 |

::: pybmad.TaoDrawingStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoDrawingStruct

Fortran struct: `tao_drawing_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L131))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele_shape` | [1D array of TaoEleShapeStruct](tao.md#taoeleshapestruct) |  |

::: pybmad.TaoDynamicApertureStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoDynamicApertureStruct

Fortran struct: `tao_dynamic_aperture_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1057))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `param` | [ApertureParamStruct](bmad.md#apertureparamstruct) |  |
| `scan` | [1D array of ApertureScanStruct](bmad.md#aperturescanstruct) | One scan for each pz. |
| `pz` | 1D array of float |  |
| `ellipse_scale` | float |  |
| `a_emit` | float |  |
| `b_emit` | float |  |

::: pybmad.TaoElePointerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoElePointerStruct

Fortran struct: `tao_ele_pointer_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L110))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `eles` | [1D array of ElePointerStruct](bmad.md#elepointerstruct) |  |
| `n_loc` | int |  |

::: pybmad.TaoEleShapeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoEleShapeStruct

Fortran struct: `tao_ele_shape_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L116))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele_id` | str | element "key::name" to match to. |
| `shape` | str | Shape to draw |
| `color` | str | Color of shape |
| `size` | float | plot vertical height |
| `label` | str | Can be: 'name', 's', 'none' |
| `draw` | bool | Draw the shape? |
| `multi` | bool | Can be part of a multi-shape. |
| `line_width` | int | Width of lines used to draw the shape. |
| `offset` | float | Vertical offset. |
| `ix_key` | int | Extracted from ele_id. 0 => all classes (quadrupole, etc.) |
| `name_ele` | str | Name of element. |
| `uni` | [1D array of TaoElePointerStruct](tao.md#taoelepointerstruct) |  |

::: pybmad.TaoEvalNodeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoEvalNodeStruct

Fortran struct: `tao_eval_node_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L98))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `type` | int |  |
| `name` | str |  |
| `scale` | float | Scale factor for ping data |
| `value` | 1D array of float |  |
| `info` | [1D array of TaoExpressionInfoStruct](tao.md#taoexpressioninfostruct) |  |
| `node` | [1D array of TaoEvalNodeStruct](tao.md#taoevalnodestruct) | Child nodes for tree construction. |

::: pybmad.TaoExpressionInfoStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoExpressionInfoStruct

Fortran struct: `tao_expression_info_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L92))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `good` | bool | Expression is valid. |
| `ele` | [EleStruct](bmad.md#elestruct) | Associated ele if it exists |
| `s` | float | Longitudinal position of expression. |

::: pybmad.TaoFloorPlanStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoFloorPlanStruct

Fortran struct: `tao_floor_plan_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L237))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `view` | str | or 'xz'. |
| `rotation` | float | Rotation of floor plan plot: 1.0 -> 360^deg |
| `correct_distortion` | bool | T -> Shrink one axis so x-scale = y-scale. |
| `flip_label_side` | bool | Draw element label on other side of element? |
| `size_is_absolute` | bool | Are shape sizes in meters or window pixels? |
| `draw_only_first_pass` | bool | Draw only first pass with multipass elements? |
| `draw_building_wall` | bool | Draw the building wall? |
| `orbit_scale` | float | Scale factor for drawing orbits. 0 -> Do not draw. |
| `orbit_color` | str |  |
| `orbit_pattern` | str |  |
| `orbit_lattice` | str | Or 'design' or 'base' |
| `orbit_width` | int |  |

::: pybmad.TaoGlobalStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoGlobalStruct

Fortran struct: `tao_global_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L641))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `beam_dead_cutoff` | float | Percentage of dead particles at which beam tracking is stopped. |
| `lm_opt_deriv_reinit` | float | Reinit derivative matrix cutoff |
| `de_lm_step_ratio` | float | Scaling for step sizes between DE and LM optimizers. |
| `de_var_to_population_factor` | float | DE population = max(n_var*factor, 20) |
| `lmdif_eps` | float | Tollerance for lmdif optimizer. |
| `lmdif_negligible_merit` | float |  |
| `svd_cutoff` | float | SVD singular value cutoff. |
| `unstable_penalty` | float | Used in unstable_ring datum merit calculation. |
| `merit_stop_value` | float | Merit value below which an optimizer will stop. |
| `dmerit_stop_value` | float | Fractional Merit change below which an optimizer will stop. |
| `random_sigma_cutoff` | float | Cut-off in sigmas. |
| `delta_e_chrom` | float | Delta E used from chrom calc. |
| `max_plot_time` | float | If plotting time (seconds) exceeds this than a message is generated. |
| `default_universe` | int | Default universe to work with. |
| `default_branch` | int | Default lattice branch to work with. |
| `n_opti_cycles` | int | Number of optimization cycles |
| `n_opti_loops` | int | Number of optimization loops |
| `n_threads` | int | Number of OpenMP threads for parallel calculations. |
| `phase_units` | int | Phase units on output. |
| `bunch_to_plot` | int | Which bunch to plot |
| `random_seed` | int | Use system clock by default |
| `n_top10_merit` | int | Number of top merit constraints to print. |
| `srdt_gen_n_slices` | int | Number times to slice elements for summation RDT calculation |
| `datum_err_messages_max` | int | Maximum number of error messages per call to lattice_calc. |
| `srdt_sxt_n_slices` | int | Number times to slice sextupoles for summation RDT calculation |
| `srdt_use_cache` | bool | Create cache for SRDT calculations.  Can use lots of memory if srdt_*_n_slices large. |
| `quiet` | str | Print I/O when running a command file? |
| `random_engine` | str | Non-beam random number engine |
| `random_gauss_converter` | str | Non-beam |
| `track_type` | str | or 'beam' |
| `lat_sigma_calc_uses_emit_from` | str | Lattice derived sigma matrix uses emit values from where? Other possibilities: "beam", "beam_init". |
| `prompt_string` | str |  |
| `prompt_color` | str | See read_a_line routine for possible settings. |
| `optimizer` | str | optimizer to use. |
| `print_command` | str |  |
| `var_out_file` | str |  |
| `history_file` | str |  |
| `beam_timer_on` | bool | For timing the beam tracking calculation. |
| `box_plots` | bool | For debugging plot layout issues. |
| `blank_line_between_commands` | bool | Add a blank line between command output? |
| `cmd_file_abort_on_error` | bool | Abort open command files if there is an error? |
| `concatenate_maps` | bool | False => tracking using DA. |
| `derivative_recalc` | bool | Recalc before each optimizer run? |
| `derivative_uses_design` | bool | Derivative calc uses design lattice instead of model? |
| `disable_smooth_line_calc` | bool | Global disable of the smooth line calculation. |
| `draw_curve_off_scale_warn` | bool | Display warning on graphs? |
| `external_plotting` | bool | Used with matplotlib and gui. |
| `label_lattice_elements` | bool | For lat_layout plots |
| `label_keys` | bool | For lat_layout plots |
| `lattice_calc_on` | bool | Turn on/off beam and single particle calculations. |
| `only_limit_opt_vars` | bool | Only apply limits to variables used in optimization. |
| `opt_with_ref` | bool | Use reference data in optimization? |
| `opt_with_base` | bool | Use base data in optimization? |
| `opt_match_auto_recalc` | bool | Set recalc = True for match elements before each cycle? |
| `opti_write_var_file` | bool | "run" command writes var_out_file |
| `optimizer_allow_user_abort` | bool | See Tao manual for more details. |
| `optimizer_var_limit_warn` | bool | Warn when vars reach a limit with optimization. |
| `plot_on` | bool | Do plotting? |
| `rad_int_user_calc_on` | bool | User set radiation integrals calculation on/off. |
| `rf_on` | bool | RFcavities on or off? Does not affect lcavities. |
| `single_step` | bool | For debugging and demonstrations: Single step through a command file? |
| `stop_on_error` | bool | For debugging: False prevents tao from exiting on an error. |
| `svd_retreat_on_merit_increase` | bool |  |
| `var_limits_on` | bool | Respect the variable limits? |
| `wait_for_CR_in_single_mode` | bool | For use with a python GUI. |
| `symbol_import` | bool | Import symbols from lattice file(s)? Internal stuff |
| `debug_on` | bool | For debugging. |
| `expression_tree_on` | bool | Use an expression tree instead of a stack? |
| `verbose_on` | bool | For verbose output. Used with debugging. |

::: pybmad.TaoGraphStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoGraphStruct

Fortran struct: `tao_graph_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L256))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Name identifying the graph |
| `type` | str | 'data', 'lat_layout', 'phase_space', 'histogram', 'dynamic_aperture' |
| `title` | str |  |
| `title_suffix` | str |  |
| `text_legend` | 1D array of str (shape: 10) | Array for holding descriptive info. |
| `text_legend_out` | 1D array of str (shape: 10) | Array for holding descriptive info. |
| `why_invalid` | str | Informative string to print. |
| `curve` | [1D array of TaoCurveStruct](tao.md#taocurvestruct) |  |
| `p` | [TaoPlotStruct](tao.md#taoplotstruct) | pointer to parent plot |
| `floor_plan` | [TaoFloorPlanStruct](tao.md#taofloorplanstruct) |  |
| `text_legend_origin` | [QpPointStruct](sim_utils.md#qppointstruct) |  |
| `curve_legend_origin` | [QpPointStruct](sim_utils.md#qppointstruct) |  |
| `curve_legend` | [QpLegendStruct](sim_utils.md#qplegendstruct) |  |
| `x` | [QpAxisStruct](sim_utils.md#qpaxisstruct) | X-axis parameters. |
| `y` | [QpAxisStruct](sim_utils.md#qpaxisstruct) | Y-axis attributes. |
| `x2` | [QpAxisStruct](sim_utils.md#qpaxisstruct) | X2-axis attributes (Not currently used). |
| `y2` | [QpAxisStruct](sim_utils.md#qpaxisstruct) | Y2-axis attributes. |
| `margin` | [QpRectStruct](sim_utils.md#qprectstruct) | Margin around the graph. |
| `scale_margin` | [QpRectStruct](sim_utils.md#qprectstruct) | Margin for scaling |
| `x_axis_scale_factor` | float | x-axis conversion from internal to plotting units. |
| `symbol_size_scale` | float | Symbol size scale factor for phase_space plots. |
| `box` | 1D array of int (shape: 4) | Defines which box the plot is put in. |
| `ix_branch` | int | Branch in lattice. Used when there are no associated curves. |
| `ix_universe` | int | Used for lat_layout plots. |
| `clip` | bool | Clip plot at graph boundary. |
| `y2_mirrors_y` | bool | Y2-axis same as Y-axis? |
| `limited` | bool | True if at least one data point past graph bounds. |
| `draw_axes` | bool | Draw axes, labels, etc? |
| `draw_curve_legend` | bool | Legend for displaying curve info. |
| `draw_grid` | bool | Draw a grid? |
| `draw_title` | bool |  |
| `draw_only_good_user_data_or_vars` | bool |  |
| `allow_wrap_around` | bool | "Wrap" curves to extend past lattice boundaries? |
| `is_valid` | bool | EG: Bad x_axis_type. |

::: pybmad.TaoHistogramStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoHistogramStruct

Fortran struct: `tao_histogram_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L165))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `density_normalized` | bool |  |
| `weight_by_charge` | bool |  |
| `minimum` | float | Computed by Tao. Not User settable. |
| `maximum` | float | Computed by Tao. Not User settable. |
| `width` | float |  |
| `center` | float |  |
| `number` | int |  |

::: pybmad.TaoInitStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoInitStruct

Fortran struct: `tao_init_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L802))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `parse_cmd_args` | bool | Used by custom programs to control Tao init |
| `debug_switch` | bool | Is the "-debug" switch present? |
| `external_plotting_switch` | bool | Is "-external_plotting" switch present? |
| `init_name` | str | label for initialization |
| `hook_init_file` | str |  |
| `hook_lat_file` | str | To be set by tao_hook_parse_command_args |
| `hook_beam_file` | str | To be set by tao_hook_parse_command_args |
| `hook_data_file` | str | To be set by tao_hook_parse_command_args |
| `hook_plot_file` | str | To be set by tao_hook_parse_command_args |
| `hook_startup_file` | str | To be set by tao_hook_parse_command_args |
| `hook_var_file` | str | To be set by tao_hook_parse_command_args |
| `hook_building_wall_file` | str | To be set by tao_hook_parse_command_args |
| `init_file_arg_path` | str | Path part of init_tao_file |
| `lattice_file_arg` | str | -lattice_file        command line argument. |
| `hook_init_file_arg` | str | -hook_init_file      command line argument |
| `init_file_arg` | str | -init_file           command line argument. |
| `beam_file_arg` | str | -beam_file           command line argument. |
| `beam_init_position_file_arg` | str | -beam_init_position_file command line argument. |
| `command_arg` | str | -command             command line argument. |
| `data_file_arg` | str | -data_file           command line argument. |
| `plot_file_arg` | str | -plot_file           command line argument. |
| `startup_file_arg` | str | -startup_file        command line argument. |
| `var_file_arg` | str | -var_file            command line argument. |
| `building_wall_file_arg` | str | -building_wall_file  command line argument. |
| `geometry_arg` | str | -geometry            command line argument. |
| `slice_lattice_arg` | str | -slice_lattice       command line argument. |
| `start_branch_at_arg` | str | -start_branch_at     command line argument. |
| `log_startup_arg` | str | -log_startup         command line argument |
| `no_stopping_arg` | str | -no_stopping         command line argument |
| `noplot_arg` | str | -noplot              command line argument |
| `no_rad_int_arg` | str | -no_rad_int          command line argument |
| `reverse_arg` | str | -reverse             command line argument |
| `debug_arg` | str | -debug               command line argument |
| `disable_smooth_line_calc_arg` | str | -disable_smooth_line_calc |
| `rf_on_arg` | str | -rf_on               command line argument |
| `prompt_color_arg` | str | -prompt_color        command line argument |
| `quiet_arg` | str | -quiet               command line argument |
| `noinit_arg` | str | -noinit              command line argument |
| `nostartup_arg` | str | -nostartup           command line argument |
| `symbol_import_arg` | str | -symbol_import       command line argument |
| `unique_name_suffix` | str |  |

::: pybmad.TaoLatSigmaStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoLatSigmaStruct

Fortran struct: `tao_lat_sigma_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L879))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `mat` | 2D array of float (shape: 6,6) |  |

::: pybmad.TaoLatticeBranchStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoLatticeBranchStruct

Fortran struct: `tao_lattice_branch_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L928))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `tao_lat` | [TaoLatticeStruct](tao.md#taolatticestruct) | Parent tao_lat |
| `lat_sigma` | [1D array of TaoLatSigmaStruct](tao.md#taolatsigmastruct) | Sigma matrix derived from lattice (not beam). |
| `spin_ele` | [1D array of TaoSpinEleStruct](tao.md#taospinelestruct) | Spin stuff |
| `bunch_params` | [1D array of BunchParamsStruct](bmad.md#bunchparamsstruct) | Per element |
| `bunch_params_comb` | [1D array of BunchTrackStruct](bmad.md#bunchtrackstruct) | A comb for each bunch in beam. |
| `orbit` | [1D array of CoordStruct](bmad.md#coordstruct) |  |
| `plot_cache` | [1D array of TaoPlotCacheStruct](tao.md#taoplotcachestruct) | Plotting data cache |
| `spin` | [TaoSpinPolarizationStruct](tao.md#taospinpolarizationstruct) |  |
| `srdt` | [SummationRdtStruct](bmad.md#summationrdtstruct) |  |
| `orb0` | [CoordStruct](bmad.md#coordstruct) | For saving beginning orbit in closed geometry branches. orb0 can then be used as an initial guess when closed_orbit is called again. |
| `modes_ri` | [NormalModesStruct](bmad.md#normalmodesstruct) | Synchrotron integrals stuff |
| `modes_6d` | [NormalModesStruct](bmad.md#normalmodesstruct) | 6D radiation matrices. |
| `ptc_normal_form` | [PtcNormalFormStruct](bmad.md#ptcnormalformstruct) | Collection of normal form structures defined in PTC |
| `bmad_normal_form` | [BmadNormalFormStruct](bmad.md#bmadnormalformstruct) | Collection of normal form structures defined in Bmad |
| `high_E_orb` | [1D array of CoordStruct](bmad.md#coordstruct) |  |
| `low_E_orb` | [1D array of CoordStruct](bmad.md#coordstruct) |  |
| `taylor_save` | [1D array of TaylorStruct (shape: 6)](bmad.md#taylorstruct) | Save to reduce computation time. |
| `cache_x_min` | float |  |
| `cache_x_max` | float |  |
| `comb_ds_save` | float | Master parameter for %bunch_params_comb(:)%ds_save |
| `ix_ref_taylor` | int |  |
| `ix_ele_taylor` | int |  |
| `track_state` | int |  |
| `cache_n_pts` | int |  |
| `ix_rad_int_cache` | int | Radiation integrals cache index. |
| `has_open_match_element` | bool |  |
| `plot_cache_valid` | bool | Valid plotting data cache? |
| `spin_map_valid` | bool |  |
| `twiss_valid` | bool | Invalid EG with unstable 1-turn matrix with a closed branch. With open branch: twiss_valid = T even if some Twiss (and orbit) is invalid. |
| `mode_flip_here` | bool | Twiss parameter mode flip seen? |
| `chrom_calc_ok` | bool |  |
| `rad_int_calc_ok` | bool |  |
| `emit_6d_calc_ok` | bool |  |
| `sigma_track_ok` | bool |  |

::: pybmad.TaoLatticeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoLatticeStruct

Fortran struct: `tao_lattice_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L967))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | "model", "base", or "design". |
| `lat` | [LatStruct](bmad.md#latstruct) | lattice structures |
| `high_E_lat` | [LatStruct](bmad.md#latstruct) | For chrom calc. |
| `low_E_lat` | [LatStruct](bmad.md#latstruct) | For chrom calc. |
| `rad_int_by_ele_ri` | [RadIntAllEleStruct](bmad.md#radintallelestruct) |  |
| `rad_int_by_ele_6d` | [RadIntAllEleStruct](bmad.md#radintallelestruct) |  |
| `tao_branch` | [1D array of TaoLatticeBranchStruct](tao.md#taolatticebranchstruct) |  |

::: pybmad.TaoModelBranchStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoModelBranchStruct

Fortran struct: `tao_model_branch_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1009))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele` | [1D array of TaoModelElementStruct](tao.md#taomodelelementstruct) | Per element information |
| `beam` | [TaoBeamBranchStruct](tao.md#taobeambranchstruct) |  |

::: pybmad.TaoModelElementStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoModelElementStruct

Fortran struct: `tao_model_element_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L984))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `beam` | [BeamStruct](bmad.md#beamstruct) | Beam distribution at element. |
| `save_beam_internally` | bool | Save beam here? Beam also saved at fork elements and at track ends. |
| `save_beam_to_file` | bool | Save beam to a file? Beam also saved at fork elements and at track ends. |

::: pybmad.TaoPingScaleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoPingScaleStruct

Fortran struct: `tao_ping_scale_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1099))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `a_mode_meas` | float |  |
| `a_mode_ref` | float |  |
| `b_mode_meas` | float |  |
| `b_mode_ref` | float |  |

::: pybmad.TaoPlotCacheStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoPlotCacheStruct

Fortran struct: `tao_plot_cache_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L918))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele_to_s` | [EleStruct](bmad.md#elestruct) | Integrated element from branch beginning. Will be marked as a hybrid element. |
| `orbit` | [CoordStruct](bmad.md#coordstruct) |  |
| `err` | bool |  |

::: pybmad.TaoPlotPageStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoPlotPageStruct

Fortran struct: `tao_plot_page_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L336))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `title` | [TaoTitleStruct](tao.md#taotitlestruct) | Title  at top of page. |
| `subtitle` | [TaoTitleStruct](tao.md#taotitlestruct) | Subtitle below title at top of page. |
| `border` | [QpRectStruct](sim_utils.md#qprectstruct) | Border around plots edge of page. |
| `floor_plan` | [TaoDrawingStruct](tao.md#taodrawingstruct) |  |
| `lat_layout` | [TaoDrawingStruct](tao.md#taodrawingstruct) |  |
| `pattern` | [1D array of TaoShapePatternStruct](tao.md#taoshapepatternstruct) |  |
| `template_` | [1D array of TaoPlotStruct](tao.md#taoplotstruct) | Templates for the plots. |
| `region` | [1D array of TaoPlotRegionStruct](tao.md#taoplotregionstruct) |  |
| `plot_display_type` | str | 'X' or 'TK' |
| `size` | 1D array of float (shape: 2) | width and height of plot window in pixels. |
| `text_height` | float | In points. Scales the height of all text |
| `main_title_text_scale` | float | Relative to text_height |
| `graph_title_text_scale` | float | Relative to text_height |
| `axis_number_text_scale` | float | Relative to text_height |
| `axis_label_text_scale` | float | Relative to text_height |
| `legend_text_scale` | float | Relative to text_height. For legends, plot_page, and lat_layout |
| `key_table_text_scale` | float | Relative to text_height |
| `floor_plan_shape_scale` | float |  |
| `floor_plan_text_scale` | float | Scale used = floor_plan_text_scale * legend_text_scale |
| `lat_layout_shape_scale` | float |  |
| `lat_layout_text_scale` | float | Scale used = lat_layout_text_scale * legend_text_scale |
| `n_curve_pts` | int | Default number of points for plotting a smooth curve. |
| `id_window` | int | X window id number. |
| `delete_overlapping_plots` | bool | Delete overlapping plots when a plot is placed? |
| `draw_graph_title_suffix` | bool | Draw the graph title suffix? |

::: pybmad.TaoPlotRegionStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoPlotRegionStruct

Fortran struct: `tao_plot_region_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L321))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Region name. Eg: 'r13', etc. |
| `plot` | [TaoPlotStruct](tao.md#taoplotstruct) | Plot associated with this region |
| `location` | 1D array of float (shape: 4) | [x1, x2, y1, y2] location on page. |
| `visible` | bool | To draw or not to draw. |
| `list_with_show_plot_command` | bool | False used for default plots to shorten the output of "show plot" |
| `setup_done` | bool | Used for plot bookkeeping. |

::: pybmad.TaoPlotStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoPlotStruct

Fortran struct: `tao_plot_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L298))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Identifying name. Rule: If name is blank, plot is not valid. |
| `description` | str | Descriptive string. |
| `graph` | [1D array of TaoGraphStruct](tao.md#taographstruct) | individual graphs of a plot |
| `r` | [TaoPlotRegionStruct](tao.md#taoplotregionstruct) | pointer to parent. |
| `ix_plot` | int | Index in s%plot_page%template(:) or %region(:) arrays. |
| `n_curve_pts` | int | Overrides s%plot_page%n_curve_pts. |
| `type` | str | or 'wave' |
| `x_axis_type` | str | 'index', 'ele_index', 's', 'none', 'floor', 'phase_space', etc. |
| `autoscale_x` | bool | Horizontal autoscale. |
| `autoscale_y` | bool | Vertical autoscale. |
| `autoscale_gang_x` | bool | scale cmd scales graphs together? |
| `autoscale_gang_y` | bool | scale cmd scales graphs together? |
| `list_with_show_plot_command` | bool | False used for default plots to shorten the output of "show plot" |
| `phantom` | bool | Used by tao_plot_init to add info lines to "show plot -templates" |
| `default_plot` | bool | One of Tao's default plots? |

::: pybmad.TaoShapePatternPointStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoShapePatternPointStruct

Fortran struct: `tao_shape_pattern_point_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L135))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `s` | float |  |
| `y` | float |  |
| `radius` | float |  |

::: pybmad.TaoShapePatternStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoShapePatternStruct

Fortran struct: `tao_shape_pattern_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L139))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `line` | [QpLineStruct](sim_utils.md#qplinestruct) | Line color and pattern set by shape using this pattern. |
| `pt` | [1D array of TaoShapePatternPointStruct](tao.md#taoshapepatternpointstruct) |  |

::: pybmad.TaoSpinDnDpzStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoSpinDnDpzStruct

Fortran struct: `tao_spin_dn_dpz_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L883))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `vec` | 1D array of float (shape: 3) | n0 derivative wrt pz. |
| `partial` | 2D array of float (shape: 3,3) | partial(i:) is spin n0 derivative wrt pz for i^th oscillation mode (1 => a-mode, etc.) |
| `partial2` | 2D array of float (shape: 3,3) | partial(i:) is spin n0 derivative wrt pz with i^th oscillation mode missing (1 => a-mode, etc.) |

::: pybmad.TaoSpinEleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoSpinEleStruct

Fortran struct: `tao_spin_ele_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L889))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `dn_dpz` | [TaoSpinDnDpzStruct](tao.md#taospindndpzstruct) |  |
| `orb_eigen_val` | 1D array of float (shape: 6) |  |
| `orb_eigen_vec` | 2D array of float (shape: 6,6) | (j,:) is j^th vector |
| `spin_eigen_vec` | 2D array of float (shape: 6,3) | (j,:) is j^th vector |
| `valid` | bool |  |

::: pybmad.TaoSpinMapStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoSpinMapStruct

Fortran struct: `tao_spin_map_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L384))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `valid` | bool |  |
| `map1` | [SpinOrbitMap1Struct](bmad.md#spinorbitmap1struct) |  |
| `axis_input` | [SpinAxisStruct](bmad.md#spinaxisstruct) | Input axes. |
| `axis0` | [SpinAxisStruct](bmad.md#spinaxisstruct) | Initial axes. |
| `axis1` | [SpinAxisStruct](bmad.md#spinaxisstruct) | Final axes. |
| `ix_ele` | int |  |
| `ix_ref` | int |  |
| `ix_uni` | int |  |
| `ix_branch` | int |  |
| `mat8` | 2D array of float (shape: 8,8) |  |

::: pybmad.TaoSpinPolarizationStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoSpinPolarizationStruct

Fortran struct: `tao_spin_polarization_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L897))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `tune` | float |  |
| `pol_limit_st` | float | Polarization calculated using Sokolov-Ternov formula. |
| `pol_limit_dk` | float | Equalibrium Polarization calculated via the Derbenev-Kondratenko-Mane formula. |
| `pol_limit_dk_partial` | 1D array of float (shape: 3) | Limit using only single mode to calc dn_dpz |
| `pol_limit_dk_partial2` | 1D array of float (shape: 3) | Limit using only single mode to calc dn_dpz |
| `pol_rate_bks` | float | BKS Polarization rate (1/sec). |
| `depol_rate` | float | Depolarization rate (1/sec). |
| `depol_rate_partial` | 1D array of float (shape: 3) | Depolarization rate (1/sec) using only single mode to calc dn_dpz. |
| `depol_rate_partial2` | 1D array of float (shape: 3) | Depolarization rate (1/sec) using only two modes to calc dn_dpz. |
| `integral_bn` | float | Integral of g^3 * b_hat * n_0 |
| `integral_bdn` | float | Integral of g^3 * b_hat * dn/ddelta |
| `integral_1ns` | float | Integral of g^3 (1 - 2(n * s_hat)/9) |
| `integral_dn2` | float | Integral of g^3 * 11 (dn/ddelta)^2 / 9 |
| `valid` | bool |  |
| `q_1turn` | [SpinOrbitMap1Struct](bmad.md#spinorbitmap1struct) | Save results from spin_concat_linear_maps in tao_spin_polarization. |
| `q_ele` | [1D array of SpinOrbitMap1Struct](bmad.md#spinorbitmap1struct) | Save results from spin_concat_linear_maps in tao_spin_polarization. |

::: pybmad.TaoSuperUniverseStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoSuperUniverseStruct

Fortran struct: `tao_super_universe_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1140))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `global_` | [TaoGlobalStruct](tao.md#taoglobalstruct) | User accessible global variables. |
| `init` | [TaoInitStruct](tao.md#taoinitstruct) | Initialization parameters |
| `com` | [TaoCommonStruct](tao.md#taocommonstruct) | Non-initialization common parameters |
| `plot_page` | [TaoPlotPageStruct](tao.md#taoplotpagestruct) | Defines the plot window. |
| `v1_var` | [1D array of TaoV1VarStruct](tao.md#taov1varstruct) | The variable types |
| `var` | [1D array of TaoVarStruct](tao.md#taovarstruct) | array of all variables. |
| `u` | [1D array of TaoUniverseStruct](tao.md#taouniversestruct) | array of universes. |
| `key` | 1D array of int |  |
| `building_wall` | [TaoBuildingWallStruct](tao.md#taobuildingwallstruct) |  |
| `wave` | [TaoWaveStruct](tao.md#taowavestruct) |  |
| `n_var_used` | int |  |
| `n_v1_var_used` | int |  |
| `history` | [1D array of TaoCmdHistoryStruct (shape: 1000)](tao.md#taocmdhistorystruct) | command history |
| `initialized` | bool | Does tao_init() need to be called? |

::: pybmad.TaoTitleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoTitleStruct

Fortran struct: `tao_title_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L152))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `string` | str | title character string. |
| `x` | float | x, y rwt lower left corner |
| `y` | float | x, y rwt lower left corner |
| `units` | str | %BOX, POINTS, etc... |
| `justify` | str | Left, Center, or Right justification. |
| `draw_it` | bool | draw the title? |

::: pybmad.TaoUniverseCalcStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoUniverseCalcStruct

Fortran struct: `tao_universe_calc_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1027))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `srdt_for_data` | int | 0 = false, 1 = 1st order, 2 = 1st & 2nd order |
| `rad_int_for_data` | bool | Do the radiation integrals need to be computed for |
| `rad_int_for_plotting` | bool | data or plotting? |
| `chrom_for_data` | bool | Does the chromaticity need to be computed for |
| `chrom_for_plotting` | bool | data or plotting? |
| `lat_sigma_for_data` | bool | Do the beam sigmas need to be computed for |
| `lat_sigma_for_plotting` | bool | data or plotting? |
| `dynamic_aperture` | bool | Do the dynamic_aperture calc? |
| `one_turn_map` | bool | Compute the one turn map? |
| `lattice` | bool | Used to indicate which lattices need tracking done. |
| `twiss` | bool | calc linear transfer matrix? |
| `track` | bool | tracking needs to be done? |
| `spin_matrices` | bool | Calculate G and D spin matrices? |

::: pybmad.TaoUniversePointerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoUniversePointerStruct

Fortran struct: `tao_universe_pointer_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1108))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `u` | [TaoUniverseStruct](tao.md#taouniversestruct) |  |

::: pybmad.TaoUniverseStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoUniverseStruct

Fortran struct: `tao_universe_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1115))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `model` | [TaoLatticeStruct](tao.md#taolatticestruct) |  |
| `design` | [TaoLatticeStruct](tao.md#taolatticestruct) |  |
| `base` | [TaoLatticeStruct](tao.md#taolatticestruct) |  |
| `beam` | [TaoBeamUniStruct](tao.md#taobeamunistruct) |  |
| `dynamic_aperture` | [TaoDynamicApertureStruct](tao.md#taodynamicaperturestruct) |  |
| `model_branch` | [1D array of TaoModelBranchStruct](tao.md#taomodelbranchstruct) | model specific information |
| `d2_data` | [1D array of TaoD2DataStruct](tao.md#taod2datastruct) | The data types |
| `data` | [1D array of TaoDataStruct](tao.md#taodatastruct) | Array of all data. |
| `ping_scale` | [TaoPingScaleStruct](tao.md#taopingscalestruct) |  |
| `scratch_lat` | [LatStruct](bmad.md#latstruct) | Scratch area. |
| `calc` | [TaoUniverseCalcStruct](tao.md#taouniversecalcstruct) | What needs to be calculated? |
| `ele_order` | [LatEleOrderStruct](bmad.md#lateleorderstruct) | Order of elements with same name. |
| `spin_map` | [TaoSpinMapStruct](tao.md#taospinmapstruct) |  |
| `dModel_dVar` | 2D array of float | Derivative matrix. |
| `ix_uni` | int | Universe index. |
| `n_d2_data_used` | int | Number of used %d2_data(:) components. |
| `n_data_used` | int | Number of used %data(:) components. |
| `is_on` | bool | universe turned on |
| `design_same_as_previous` | bool | Design lat same as the previous uni? |
| `picked_uni` | bool | Scratch logical. |

::: pybmad.TaoV1VarStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoV1VarStruct

Fortran struct: `tao_v1_var_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L592))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | V1 variable name. Eg: 'quad_k1'. |
| `ix_v1_var` | int | Index to s%v1_var(:) array |
| `v` | [1D array of TaoVarStruct](tao.md#taovarstruct) | Pointer to the appropriate section in s%var. |

::: pybmad.TaoVarSlaveStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoVarSlaveStruct

Fortran struct: `tao_var_slave_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L524))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ix_uni` | int | universe index. |
| `ix_branch` | int |  |
| `ix_ele` | int | Index of element in the u%lattice%ele(:) array. |
| `model_value` | float | Pointer to the variable in the model lat. |
| `base_value` | float | Pointer to the variable in the base lat. |

::: pybmad.TaoVarStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoVarStruct

Fortran struct: `tao_var_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L548))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele_name` | str | Associated lattice element name. |
| `attrib_name` | str | Name of the attribute to vary. |
| `id` | str | Used by Tao extension code. Not used by Tao directly. |
| `slave` | [1D array of TaoVarSlaveStruct](tao.md#taovarslavestruct) |  |
| `ix_v1` | int | Index of this var in the s%v1_var(i)%v(:) array. |
| `ix_var` | int | Index number of this var in the s%var(:) array. |
| `ix_dvar` | int | Column in the dData_dVar derivative matrix. |
| `ix_attrib` | int | Index in ele%value(:) array if appropriate. |
| `ix_key_table` | int | Has a key binding? |
| `model_value` | float | Model value. |
| `base_value` | float | Base value. |
| `design_value` | float | Design value from the design lattice. |
| `scratch_value` | float | Scratch space used by Tao. |
| `old_value` | float | Scratch space used by Tao. |
| `meas_value` | float | The value when the data measurement was taken. |
| `ref_value` | float | Value when the reference measurement was taken. |
| `correction_value` | float | Value determined by a fit to correct the lattice. |
| `high_lim` | float | High limit for the model_value. |
| `low_lim` | float | Low limit for the model_value. |
| `step` | float | Sets what is a small step for varying this var. |
| `weight` | float | Weight for the merit function term. |
| `delta_merit` | float | Diff used to calculate the merit function term. |
| `merit` | float | merit_term = weight * delta^2. |
| `dMerit_dVar` | float | Merit derivative. |
| `key_val0` | float | Key base value |
| `key_delta` | float | Change in value when a key is pressed. |
| `s` | float | longitudinal position of ele. |
| `extend_val` | float | For extension code. Not used by Tao. |
| `merit_type` | str | 'target' or 'limit' |
| `exists` | bool | See above |
| `good_var` | bool | See above |
| `good_user` | bool | See above |
| `good_opt` | bool | See above |
| `good_plot` | bool | See above |
| `useit_opt` | bool | See above |
| `useit_plot` | bool | See above |
| `key_bound` | bool | Variable bound to keyboard key? |
| `v1` | [TaoV1VarStruct](tao.md#taov1varstruct) | Pointer to the parent. |

::: pybmad.TaoWaveKickPtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoWaveKickPtStruct

Fortran struct: `tao_wave_kick_pt_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1068))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `phi_s` | float |  |
| `phi_r` | float |  |
| `phi` | float |  |
| `amp` | float |  |
| `s` | float | s-position of kick |
| `ix_dat_before_kick` | int | Index of datum in data array just before the kick. |
| `ele` | [EleStruct](bmad.md#elestruct) | lattice element at position of kick. |

::: pybmad.TaoWaveStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaoWaveStruct

Fortran struct: `tao_wave_struct` ([`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1075))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `data_type` | str |  |
| `rms_rel_a` | float |  |
| `rms_rel_b` | float |  |
| `rms_rel_as` | float |  |
| `rms_rel_bs` | float |  |
| `rms_rel_ar` | float |  |
| `rms_rel_br` | float |  |
| `rms_rel_k` | float |  |
| `rms_rel_ks` | float |  |
| `rms_rel_kr` | float |  |
| `rms_phi` | float |  |
| `rms_phi_s` | float |  |
| `rms_phi_r` | float |  |
| `amp_ba_s` | float |  |
| `amp_ba_r` | float |  |
| `chi_a` | float |  |
| `chi_c` | float |  |
| `chi_ba` | float |  |
| `amp_a` | 1D array of float (shape: 2) |  |
| `amp_b` | 1D array of float (shape: 2) |  |
| `amp_ba` | 1D array of float (shape: 2) |  |
| `coef_a` | 1D array of float (shape: 4) |  |
| `coef_b` | 1D array of float (shape: 4) |  |
| `coef_ba` | 1D array of float (shape: 4) |  |
| `n_func` | int | Number of functions used in the fit. |
| `ix_a1` | int |  |
| `ix_a2` | int |  |
| `ix_b1` | int |  |
| `ix_b2` | int |  |
| `i_a1` | int |  |
| `i_a2` | int |  |
| `i_b1` | int |  |
| `i_b2` | int |  |
| `n_a` | int |  |
| `n_b` | int |  |
| `i_curve_wrap_pt` | int | Index of last point before wrap in curve array. |
| `ix_data` | 1D array of int | Translates from plot point to datum index |
| `n_kick` | int |  |
| `kick` | [1D array of TaoWaveKickPtStruct](tao.md#taowavekickptstruct) |  |
| `base_graph` | [TaoGraphStruct](tao.md#taographstruct) | Graph before curves extended to 1.5 periods. |
| `region` | [TaoPlotRegionStruct](tao.md#taoplotregionstruct) | Where the wave plot is |
| `d1_dat` | [TaoD1DataStruct](tao.md#taod1datastruct) | D1 data for analysis |

## Procedures

### integrate_max

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L884)

::: pybmad.integrate_max
    options:
      show_root_heading: false
      show_root_toc_entry: false

### integrate_min

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L847)

::: pybmad.integrate_min
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_abort_command_file

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L47)

::: pybmad.tao_abort_command_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_add_to_normal_mode_h_array

Fortran source: [`tao/code/tao_init_data_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_data_mod.f90#L918)

::: pybmad.tao_add_to_normal_mode_h_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_alias_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L52)

::: pybmad.tao_alias_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_allocate_data_array

Fortran source: [`tao/code/tao_init_data_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_data_mod.f90#L701)

::: pybmad.tao_allocate_data_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_allocate_v1_var

Fortran source: [`tao/code/tao_init_variables_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_variables_mod.f90#L647)

::: pybmad.tao_allocate_v1_var
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_allocate_var_array

Fortran source: [`tao/code/tao_init_variables_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_variables_mod.f90#L950)

::: pybmad.tao_allocate_var_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_beam_emit_calc

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L58)

::: pybmad.tao_beam_emit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_beam_track

Fortran source: [`tao/code/tao_lattice_calc_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_lattice_calc_mod.f90#L435)

::: pybmad.tao_beam_track
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_beam_track_endpoint

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L67)

::: pybmad.tao_beam_track_endpoint
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_branch_index

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L76)

::: pybmad.tao_branch_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_calc_data_at_s_pts

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L2191)

::: pybmad.tao_calc_data_at_s_pts
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_cbar_wave_anal

Fortran source: [`tao/code/tao_wave_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_wave_mod.f90#L763)

::: pybmad.tao_cbar_wave_anal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_change_ele

Fortran source: [`tao/code/tao_change_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_change_mod.f90#L227)

::: pybmad.tao_change_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_change_tune

Fortran source: [`tao/code/tao_change_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_change_mod.f90#L25)

::: pybmad.tao_change_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_change_var

Fortran source: [`tao/code/tao_change_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_change_mod.f90#L89)

::: pybmad.tao_change_var
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_change_z_tune

Fortran source: [`tao/code/tao_change_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_change_mod.f90#L55)

::: pybmad.tao_change_z_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_chrom_calc_needed

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L96)

::: pybmad.tao_chrom_calc_needed
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_clear_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L103)

::: pybmad.tao_clear_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_clip_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L108)

::: pybmad.tao_clip_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_close_command_file

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L116)

::: pybmad.tao_close_command_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_cmd_history_record

Fortran source: [`tao/code/tao_command_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_command_mod.f90#L17)

::: pybmad.tao_cmd_history_record
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_command

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L119)

::: pybmad.tao_command
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_constraint_type_name

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L132)

::: pybmad.tao_constraint_type_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_control_tree_list

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L125)

::: pybmad.tao_control_tree_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_count_strings

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L139)

::: pybmad.tao_count_strings
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_create_plot_window

Fortran source: [`tao/code/tao_plot_window_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_window_mod.f90#L20)

::: pybmad.tao_create_plot_window
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_curve_beam_ellipse_setup

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L789)

::: pybmad.tao_curve_beam_ellipse_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_curve_check_universe

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L3077)

::: pybmad.tao_curve_check_universe
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_curve_data_setup

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L1241)

::: pybmad.tao_curve_data_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_curve_datum_calc

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L2886)

::: pybmad.tao_curve_datum_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_curve_ele_ref

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L146)

::: pybmad.tao_curve_ele_ref
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_curve_ix_uni

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L154)

::: pybmad.tao_curve_ix_uni
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_curve_name

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L161)

::: pybmad.tao_curve_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_curve_rms_calc

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L169)

::: pybmad.tao_curve_rms_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_d2_d1_name

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L178)

::: pybmad.tao_d2_d1_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_d2_data_stuffit

Fortran source: [`tao/code/tao_init_data_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_data_mod.f90#L769)

::: pybmad.tao_d2_data_stuffit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_data_check

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L186)

::: pybmad.tao_data_check
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_data_coupling_init

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L192)

::: pybmad.tao_data_coupling_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_data_sanity_check

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L198)

::: pybmad.tao_data_sanity_check
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_data_type_substitute

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L389)

::: pybmad.tao_data_type_substitute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_data_useit_plot_calc

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L2787)

::: pybmad.tao_data_useit_plot_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_datum_has_associated_ele

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L214)

::: pybmad.tao_datum_has_associated_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_datum_integrate

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L739)

::: pybmad.tao_datum_integrate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_datum_name

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L221)

::: pybmad.tao_datum_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_datum_s_position

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L702)

::: pybmad.tao_datum_s_position
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_de_optimizer

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L229)

::: pybmad.tao_de_optimizer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_deallocate_plot_cache

Fortran source: [`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1164)

::: pybmad.tao_deallocate_plot_cache
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_deallocate_tree

Fortran source: [`tao/code/tao_expression_tree_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_expression_tree_mod.f90#L208)

::: pybmad.tao_deallocate_tree
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_destroy_plot_window

Fortran source: [`tao/code/tao_plot_window_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_window_mod.f90#L62)

::: pybmad.tao_destroy_plot_window
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_dmerit_calc

Fortran source: [`tao/code/tao_dmerit_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_dmerit_mod.f90#L272)

::: pybmad.tao_dmerit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_dmodel_dvar_calc

Fortran source: [`tao/code/tao_dmerit_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_dmerit_mod.f90#L30)

::: pybmad.tao_dmodel_dvar_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_do_wire_scan

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L1032)

::: pybmad.tao_do_wire_scan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_draw_beam_chamber_wall

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L1560)

::: pybmad.tao_draw_beam_chamber_wall
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_draw_curve_data

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L1710)

::: pybmad.tao_draw_curve_data
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_draw_ele_for_floor_plan

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L649)

::: pybmad.tao_draw_ele_for_floor_plan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_draw_floor_plan

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L393)

::: pybmad.tao_draw_floor_plan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_draw_graph_axes

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L1871)

::: pybmad.tao_draw_graph_axes
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_draw_histogram_data

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L1831)

::: pybmad.tao_draw_histogram_data
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_draw_lat_layout

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L1190)

::: pybmad.tao_draw_lat_layout
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_draw_plots

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L22)

::: pybmad.tao_draw_plots
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_ele_geometry_with_misalignments

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L2040)

::: pybmad.tao_ele_geometry_with_misalignments
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_ele_shape_info

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L234)

::: pybmad.tao_ele_shape_info
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_eval_floor_orbit

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L2094)

::: pybmad.tao_eval_floor_orbit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_evaluate_a_datum

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L246)

::: pybmad.tao_evaluate_a_datum
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_evaluate_datum_at_s

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L1495)

::: pybmad.tao_evaluate_datum_at_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_evaluate_element_parameters

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L335)

::: pybmad.tao_evaluate_element_parameters
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_evaluate_expression

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L260)

::: pybmad.tao_evaluate_expression
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_evaluate_expression_new

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L282)

::: pybmad.tao_evaluate_expression_new
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_evaluate_expression_old

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L304)

::: pybmad.tao_evaluate_expression_old
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_evaluate_lat_or_beam_data

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L39)

::: pybmad.tao_evaluate_lat_or_beam_data
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_evaluate_stack_old

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L1619)

::: pybmad.tao_evaluate_stack_old
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_evaluate_tree

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L350)

::: pybmad.tao_evaluate_tree
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_evaluate_tune

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L326)

::: pybmad.tao_evaluate_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_expression_hash_substitute

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L426)

::: pybmad.tao_expression_hash_substitute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_expression_tree_to_string

Fortran source: [`tao/code/tao_expression_tree_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_expression_tree_mod.f90#L63)

::: pybmad.tao_expression_tree_to_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_find_plot_region

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L397)

::: pybmad.tao_find_plot_region
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_fixer

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L417)

::: pybmad.tao_fixer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_floor_to_screen

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L422)

::: pybmad.tao_floor_to_screen
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_floor_to_screen_coords

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L429)

::: pybmad.tao_floor_to_screen_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_geodesic_lm_optimizer

Fortran source: [`tao/code/tao_geodesic_lm_optimizer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_geodesic_lm_optimizer_mod.f90#L30)

::: pybmad.tao_geodesic_lm_optimizer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_get_data

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L368)

::: pybmad.tao_get_data
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_get_opt_vars

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L436)

::: pybmad.tao_get_opt_vars
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_get_user_input

Fortran source: [`tao/code/tao_get_user_input_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_get_user_input_mod.f90#L46)

::: pybmad.tao_get_user_input
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_graph_controller_setup

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L128)

::: pybmad.tao_graph_controller_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_graph_data_setup

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L1173)

::: pybmad.tao_graph_data_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_graph_data_slice_setup

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L257)

::: pybmad.tao_graph_data_slice_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_graph_dynamic_aperture_setup

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L674)

::: pybmad.tao_graph_dynamic_aperture_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_graph_histogram_setup

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L858)

::: pybmad.tao_graph_histogram_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_graph_name

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L447)

::: pybmad.tao_graph_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_graph_phase_space_setup

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L425)

::: pybmad.tao_graph_phase_space_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_graph_s_min_max_calc

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L3116)

::: pybmad.tao_graph_s_min_max_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_graph_setup

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L13)

::: pybmad.tao_graph_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L462)

::: pybmad.tao_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_beam_in_universe

Fortran source: [`tao/code/tao_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_mod.f90#L331)

::: pybmad.tao_init_beam_in_universe
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_beams

Fortran source: [`tao/code/tao_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_mod.f90#L153)

::: pybmad.tao_init_beams
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_data

Fortran source: [`tao/code/tao_init_data_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_data_mod.f90#L22)

::: pybmad.tao_init_data
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_data_end_stuff

Fortran source: [`tao/code/tao_init_data_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_data_mod.f90#L674)

::: pybmad.tao_init_data_end_stuff
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_data_in_universe

Fortran source: [`tao/code/tao_init_data_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_data_mod.f90#L807)

::: pybmad.tao_init_data_in_universe
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_dynamic_aperture

Fortran source: [`tao/code/tao_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_mod.f90#L417)

::: pybmad.tao_init_dynamic_aperture
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_find_elements

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L467)

::: pybmad.tao_init_find_elements
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_global

Fortran source: [`tao/code/tao_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_mod.f90#L22)

::: pybmad.tao_init_global
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_lattice

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L477)

::: pybmad.tao_init_lattice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_plotting

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L483)

::: pybmad.tao_init_plotting
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_init_variables

Fortran source: [`tao/code/tao_init_variables_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_variables_mod.f90#L23)

::: pybmad.tao_init_variables
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_inject_beam

Fortran source: [`tao/code/tao_lattice_calc_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_lattice_calc_mod.f90#L869)

::: pybmad.tao_inject_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_inject_particle

Fortran source: [`tao/code/tao_lattice_calc_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_lattice_calc_mod.f90#L775)

::: pybmad.tao_inject_particle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_is_valid_name

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L488)

::: pybmad.tao_is_valid_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_json_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L494)

::: pybmad.tao_json_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_key_info_to_str

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L500)

::: pybmad.tao_key_info_to_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_lat_bookkeeper

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L508)

::: pybmad.tao_lat_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_lat_emit_calc

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L516)

::: pybmad.tao_lat_emit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_lat_sigma_calc_needed

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L525)

::: pybmad.tao_lat_sigma_calc_needed
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_lat_sigma_track

Fortran source: [`tao/code/tao_lattice_calc_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_lattice_calc_mod.f90#L251)

::: pybmad.tao_lat_sigma_track
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_lattice_branches_equal_tao_lattice_branches

Fortran source: [`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1186)

::: pybmad.tao_lattice_branches_equal_tao_lattice_branches
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_lattice_calc

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L532)

::: pybmad.tao_lattice_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_lattice_equal_tao_lattice

Fortran source: [`tao/code/tao_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_struct.f90#L1209)

::: pybmad.tao_lattice_equal_tao_lattice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_limit_calc

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L538)

::: pybmad.tao_limit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_lm_optimizer

Fortran source: [`tao/code/tao_lm_optimizer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_lm_optimizer_mod.f90#L28)

::: pybmad.tao_lm_optimizer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_lmdif_optimizer

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L543)

::: pybmad.tao_lmdif_optimizer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_load_this_datum

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L461)

::: pybmad.tao_load_this_datum
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_locate_all_elements

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L548)

::: pybmad.tao_locate_all_elements
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_locate_elements

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L557)

::: pybmad.tao_locate_elements
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_mark_lattice_ele

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L570)

::: pybmad.tao_mark_lattice_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_merit

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L576)

::: pybmad.tao_merit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_next_word

Fortran source: [`tao/code/tao_command_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_command_mod.f90#L297)

::: pybmad.tao_next_word
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_one_turn_map_calc_needed

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L583)

::: pybmad.tao_one_turn_map_calc_needed
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_open_file

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L590)

::: pybmad.tao_open_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_open_scratch_file

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L598)

::: pybmad.tao_open_scratch_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_optimization_status

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L604)

::: pybmad.tao_optimization_status
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_orbit_beta_wave_anal

Fortran source: [`tao/code/tao_wave_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_wave_mod.f90#L393)

::: pybmad.tao_orbit_beta_wave_anal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_oreint_building_wall_pt

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L611)

::: pybmad.tao_oreint_building_wall_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_param_value_at_s

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L645)

::: pybmad.tao_param_value_at_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_param_value_routine

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L1158)

::: pybmad.tao_param_value_routine
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_parse_command_args

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L657)

::: pybmad.tao_parse_command_args
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_parse_element_param_str

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L664)

::: pybmad.tao_parse_element_param_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_particle_data_value

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L1105)

::: pybmad.tao_particle_data_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pause_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L672)

::: pybmad.tao_pause_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_phase_space_axis_index

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L1055)

::: pybmad.tao_phase_space_axis_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_phase_wave_anal

Fortran source: [`tao/code/tao_wave_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_wave_mod.f90#L584)

::: pybmad.tao_phase_wave_anal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pick_universe

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L678)

::: pybmad.tao_pick_universe
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pipe_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L688)

::: pybmad.tao_pipe_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_place_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L694)

::: pybmad.tao_place_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_plot_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L701)

::: pybmad.tao_plot_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_plot_data

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L1669)

::: pybmad.tao_plot_data
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_plot_histogram

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L197)

::: pybmad.tao_plot_histogram
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_plot_key_table

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L324)

::: pybmad.tao_plot_key_table
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_plot_setup

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L707)

::: pybmad.tao_plot_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_plot_struct_transfer

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L711)

::: pybmad.tao_plot_struct_transfer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_plot_wave

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L236)

::: pybmad.tao_plot_wave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pointer_to_building_wall_shape

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L717)

::: pybmad.tao_pointer_to_building_wall_shape
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pointer_to_datum

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L617)

::: pybmad.tao_pointer_to_datum
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pointer_to_datum_ele

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L1097)

::: pybmad.tao_pointer_to_datum_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pointer_to_ele_shape

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L724)

::: pybmad.tao_pointer_to_ele_shape
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pointer_to_tao_lat

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L736)

::: pybmad.tao_pointer_to_tao_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pointer_to_universe

Fortran sources (overloaded):

- `tao_pointer_to_universe_int`: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L1268)
- `tao_pointer_to_universe_str`: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L1300)

::: pybmad.tao_pointer_to_universe
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pointer_to_universes

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L634)

::: pybmad.tao_pointer_to_universes
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pointer_to_var_in_lattice

Fortran source: [`tao/code/tao_init_variables_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_variables_mod.f90#L776)

::: pybmad.tao_pointer_to_var_in_lattice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_pointer_to_var_in_lattice2

Fortran source: [`tao/code/tao_init_variables_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_init_variables_mod.f90#L866)

::: pybmad.tao_pointer_to_var_in_lattice2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_print_command_line_info

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L744)

::: pybmad.tao_print_command_line_info
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_ptc_normal_form

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L749)

::: pybmad.tao_ptc_normal_form
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_python_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L757)

::: pybmad.tao_python_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_quiet_set

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L876)

::: pybmad.tao_quiet_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_rad_int_calc_needed

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L763)

::: pybmad.tao_rad_int_calc_needed
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_re_allocate_expression_info

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L770)

::: pybmad.tao_re_allocate_expression_info
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_re_associate_node_array

Fortran source: [`tao/code/tao_expression_tree_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_expression_tree_mod.f90#L167)

::: pybmad.tao_re_associate_node_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_re_execute

Fortran source: [`tao/code/tao_command_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_command_mod.f90#L46)

::: pybmad.tao_re_execute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_read_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L788)

::: pybmad.tao_read_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_read_phase_space_index

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L794)

::: pybmad.tao_read_phase_space_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_regression_test

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L778)

::: pybmad.tao_regression_test
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_remove_blank_characters

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L783)

::: pybmad.tao_remove_blank_characters
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_run_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L802)

::: pybmad.tao_run_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_scale_cmd

Fortran source: [`tao/code/tao_scale_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_scale_mod.f90#L31)

::: pybmad.tao_scale_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_scale_graph

Fortran source: [`tao/code/tao_scale_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_scale_mod.f90#L311)

::: pybmad.tao_scale_graph
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_scale_ping_data

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L808)

::: pybmad.tao_scale_ping_data
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_scale_plot

Fortran source: [`tao/code/tao_scale_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_scale_mod.f90#L187)

::: pybmad.tao_scale_plot
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_scratch_values_calc

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L921)

::: pybmad.tao_scratch_values_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_beam_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L984)

::: pybmad.tao_set_beam_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_beam_init_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L1139)

::: pybmad.tao_set_beam_init_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_bmad_com_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L734)

::: pybmad.tao_set_bmad_com_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_branch_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L2296)

::: pybmad.tao_set_branch_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_calculate_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L196)

::: pybmad.tao_set_calculate_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_curve_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L1487)

::: pybmad.tao_set_curve_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_curve_invalid

Fortran source: [`tao/code/tao_graph_setup_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_graph_setup_mod.f90#L3045)

::: pybmad.tao_set_curve_invalid
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_data_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L2399)

::: pybmad.tao_set_data_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_data_useit_opt

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L814)

::: pybmad.tao_set_data_useit_opt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_default_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L2729)

::: pybmad.tao_set_default_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_drawing_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L3498)

::: pybmad.tao_set_drawing_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_dynamic_aperture_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L2782)

::: pybmad.tao_set_dynamic_aperture_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_elements_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L3022)

::: pybmad.tao_set_elements_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_floor_plan_axis_label

Fortran source: [`tao/code/tao_plot_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_plot_mod.f90#L591)

::: pybmad.tao_set_floor_plan_axis_label
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_geodesic_lm_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L823)

::: pybmad.tao_set_geodesic_lm_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_global_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L485)

::: pybmad.tao_set_global_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_graph_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L1921)

::: pybmad.tao_set_graph_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_integer_value

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L3315)

::: pybmad.tao_set_integer_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_invalid

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L830)

::: pybmad.tao_set_invalid
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_key_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L241)

::: pybmad.tao_set_key_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_lattice_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L345)

::: pybmad.tao_set_lattice_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_logical_value

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L3266)

::: pybmad.tao_set_logical_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_openmp_n_threads

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L156)

::: pybmad.tao_set_openmp_n_threads
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_opt_vars

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L851)

::: pybmad.tao_set_opt_vars
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_opti_de_param_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L875)

::: pybmad.tao_set_opti_de_param_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_particle_start_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L1311)

::: pybmad.tao_set_particle_start_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_plot_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L1742)

::: pybmad.tao_set_plot_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_plot_page_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L1400)

::: pybmad.tao_set_plot_page_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_ptc_com_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L785)

::: pybmad.tao_set_ptc_com_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_qp_axis_struct

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L3749)

::: pybmad.tao_set_qp_axis_struct
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_qp_point_struct

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L3868)

::: pybmad.tao_set_qp_point_struct
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_qp_rect_struct

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L3688)

::: pybmad.tao_set_qp_rect_struct
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_ran_state_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L282)

::: pybmad.tao_set_ran_state_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_real_value

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L3440)

::: pybmad.tao_set_real_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_region_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L1864)

::: pybmad.tao_set_region_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_space_charge_com_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L662)

::: pybmad.tao_set_space_charge_com_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_symbolic_number_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L3597)

::: pybmad.tao_set_symbolic_number_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_tune_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L26)

::: pybmad.tao_set_tune_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_universe_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L2869)

::: pybmad.tao_set_universe_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_var_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L2150)

::: pybmad.tao_set_var_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_var_model_value

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L840)

::: pybmad.tao_set_var_model_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_var_useit_opt

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L848)

::: pybmad.tao_set_var_useit_opt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_wave_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L923)

::: pybmad.tao_set_wave_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_set_z_tune_cmd

Fortran source: [`tao/code/tao_set_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_set_mod.f90#L109)

::: pybmad.tao_set_z_tune_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_setup_key_table

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L858)

::: pybmad.tao_setup_key_table
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_shape_init

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L923)

::: pybmad.tao_shape_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_show_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L931)

::: pybmad.tao_show_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_show_constraints

Fortran source: [`tao/code/tao_top10_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_top10_mod.f90#L292)

::: pybmad.tao_show_constraints
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_single_mode

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L882)

::: pybmad.tao_single_mode
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_single_track

Fortran source: [`tao/code/tao_lattice_calc_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_lattice_calc_mod.f90#L34)

::: pybmad.tao_single_track
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_spin_matrices_calc_needed

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L913)

::: pybmad.tao_spin_matrices_calc_needed
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_spin_tracking_turn_on

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L920)

::: pybmad.tao_spin_tracking_turn_on
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_split_component

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L887)

::: pybmad.tao_split_component
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_srdt_calc_needed

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L863)

::: pybmad.tao_srdt_calc_needed
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_subin_uni_number

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L944)

::: pybmad.tao_subin_uni_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_svd_optimizer

Fortran source: [`tao/code/tao_svd_optimizer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_svd_optimizer_mod.f90#L20)

::: pybmad.tao_svd_optimizer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_symbol_import_from_lat

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L870)

::: pybmad.tao_symbol_import_from_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_taper_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L952)

::: pybmad.tao_taper_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_to_change_number

Fortran source: [`tao/code/tao_change_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_change_mod.f90#L398)

::: pybmad.tao_to_change_number
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_to_int

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L1995)

::: pybmad.tao_to_int
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_to_phase_and_coupling_reading

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L312)

::: pybmad.tao_to_phase_and_coupling_reading
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_to_real

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L958)

::: pybmad.tao_to_real
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_too_many_particles_lost

Fortran source: [`tao/code/tao_lattice_calc_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_lattice_calc_mod.f90#L738)

::: pybmad.tao_too_many_particles_lost
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_top10_derivative_print

Fortran source: [`tao/code/tao_top10_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_top10_mod.f90#L130)

::: pybmad.tao_top10_derivative_print
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_top10_merit_categories_print

Fortran source: [`tao/code/tao_top10_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_top10_mod.f90#L31)

::: pybmad.tao_top10_merit_categories_print
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_top_level

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L966)

::: pybmad.tao_top_level
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_tracking_ele_index

Fortran source: [`tao/code/tao_data_and_eval_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_data_and_eval_mod.f90#L821)

::: pybmad.tao_tracking_ele_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_turn_on_special_calcs_if_needed_for_plotting

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L972)

::: pybmad.tao_turn_on_special_calcs_if_needed_for_plotting
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_type_expression_tree

Fortran source: [`tao/code/tao_expression_tree_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_expression_tree_mod.f90#L21)

::: pybmad.tao_type_expression_tree
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_uni_atsign_index

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L1364)

::: pybmad.tao_uni_atsign_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_universe_index

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L977)

::: pybmad.tao_universe_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_use_data

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L984)

::: pybmad.tao_use_data
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_use_var

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L990)

::: pybmad.tao_use_var
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_user_is_terminating_optimization

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L996)

::: pybmad.tao_user_is_terminating_optimization
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_var1_name

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L1001)

::: pybmad.tao_var1_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_var_attrib_name

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L1008)

::: pybmad.tao_var_attrib_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_var_check

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L88)

::: pybmad.tao_var_check
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_var_repoint

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L1015)

::: pybmad.tao_var_repoint
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_var_target_calc

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L1026)

::: pybmad.tao_var_target_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_var_useit_plot_calc

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L1031)

::: pybmad.tao_var_useit_plot_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_var_write

Fortran source: [`tao/code/tao_top10_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_top10_mod.f90#L535)

::: pybmad.tao_var_write
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_veto_vars_with_zero_dmodel

Fortran source: [`tao/code/tao_dmerit_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_dmerit_mod.f90#L236)

::: pybmad.tao_veto_vars_with_zero_dmodel
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_wave_analysis

Fortran source: [`tao/code/tao_wave_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_wave_mod.f90#L202)

::: pybmad.tao_wave_analysis
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_wave_cmd

Fortran source: [`tao/code/tao_wave_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_wave_mod.f90#L23)

::: pybmad.tao_wave_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_wave_fit

Fortran source: [`tao/code/tao_wave_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_wave_mod.f90#L988)

::: pybmad.tao_wave_fit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_write_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L1038)

::: pybmad.tao_write_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_x_axis_cmd

Fortran source: [`tao/code/tao_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_interface.f90#L1043)

::: pybmad.tao_x_axis_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_x_scale_cmd

Fortran source: [`tao/code/tao_x_scale_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_x_scale_mod.f90#L33)

::: pybmad.tao_x_scale_cmd
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_x_scale_graph

Fortran source: [`tao/code/tao_x_scale_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_x_scale_mod.f90#L266)

::: pybmad.tao_x_scale_graph
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tao_x_scale_plot

Fortran source: [`tao/code/tao_x_scale_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/tao/code/tao_x_scale_mod.f90#L177)

::: pybmad.tao_x_scale_plot
    options:
      show_root_heading: false
      show_root_toc_entry: false
