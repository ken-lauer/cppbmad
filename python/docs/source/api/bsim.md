# Bsim

Bsim application-level structs and routines (BBU, etc.).

## Classes (Fortran Structures)

::: pybmad.BbuBeamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BbuBeamStruct

Fortran struct: `bbu_beam_struct` ([`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L19))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `bunch` | [1D array of BunchStruct](bmad.md#bunchstruct) | Bunches in the lattice |
| `stage` | [1D array of BbuStageStruct](bsim.md#bbustagestruct) |  |
| `ix_ele_bunch` | 1D array of int | element where bunch is |
| `ix_bunch_head` | int | Index to head bunch(:) |
| `ix_bunch_end` | int | Index of the end bunch(:). -1 -> no bunches. |
| `n_bunch_in_lat` | int | Number of bunches transversing the lattice. |
| `ix_stage_voltage_max` | int |  |
| `hom_voltage_max` | float |  |
| `time_now` | float |  |
| `one_turn_time` | float |  |
| `rf_wavelength_max` | float |  |

::: pybmad.BbuParamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BbuParamStruct

Fortran struct: `bbu_param_struct` ([`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L33))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `lat_filename` | str | Bmad lattice file name |
| `lat2_filename` | str | Bmad lattice2 file name for secondary parser |
| `bunch_by_bunch_info_file` | str | For outputting bunch-by-bunch info. |
| `hybridize` | bool | Combine non-hom elements to speed up simulation? |
| `write_digested_hybrid_lat` | bool | For debugging purposes. |
| `write_voltage_vs_time_dat` | bool | For debugging purposes. |
| `keep_overlays_and_groups` | bool | Keep when hybridizing? |
| `keep_all_lcavities` | bool | Keep when hybridizing? |
| `use_taylor_for_hybrids` | bool | Use taylor map for hybrids when true. Otherwise tracking method is linear. |
| `stable_orbit_anal` | bool | Write stable_orbit.out and hom_voltage.out? |
| `limit_factor` | float | Init_hom_amp * limit_factor = simulation unstable limit |
| `simulation_turns_max` | float | Sets the duration of the simulation. |
| `bunch_freq` | float | Freq in Hz. |
| `init_particle_offset` | float | Initial particle offset for particles born in the first turn period. |
| `current` | float | Starting current (amps) |
| `rel_tol` | float | Final threshold current accuracy. |
| `drscan` | bool | If true, scan DR variable as in PRSTAB 7 (2004) Fig. 3. |
| `use_interpolated_threshold` | bool |  |
| `write_hom_info` | bool | Write HOM parameters to main output file? |
| `elindex` | int |  |
| `elname` | str | Element to step length for DRSCAN |
| `nstep` | int | Number of steps for DRSCAN. |
| `begdr` | float | Beginning DR value for DRSCAN. |
| `enddr` | float | End DR value for DRSCAN. |
| `nrep` | int | Number of times to repeat threshold calculation |
| `ran_seed` | int | If set to 0, the output results will vary from run to run. |
| `hom_order_cutoff` | int | If positive -> ignore HOM's with order greater than this. |
| `ran_gauss_sigma_cut` | float |  |
| `ele_track_end` | str |  |
| `ix_ele_track_end` | int | Default: set to last element with a wake |
| `regression` | bool | Do regression test? |
| `normalize_z_to_rf` | bool | make starting z = mod(z, rf_wavelength)? Ramp parameters |
| `ramp_on` | bool |  |
| `ramp_pattern` | 1D array of float (shape: 1000) |  |
| `ramp_n_start` | int | Index of start of ramp Internal parameters |
| `n_ramp_pattern` | int | Number of valid ramp_pattern |

::: pybmad.BbuStageStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BbuStageStruct

Fortran struct: `bbu_stage_struct` ([`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L6))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ix_ele_lr_wake` | int | Element index of element with the wake |
| `ix_ele_stage_end` | int | Element at end of stage. |
| `ix_pass` | int | Pass index when in multipass section |
| `ix_stage_pass1` | int | Index of corresponding stage on first pass |
| `ix_head_bunch` | int |  |
| `ix_hom_max` | int |  |
| `hom_voltage_max` | float |  |
| `time_at_wake_ele` | float |  |
| `ave_orb` | 1D array of float (shape: 6) |  |
| `rms_orb` | 1D array of float (shape: 6) |  |
| `min_orb` | 1D array of float (shape: 6) |  |
| `max_orb` | 1D array of float (shape: 6) |  |
| `n_orb` | int |  |

## Procedures

### bbu_add_a_bunch

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L559)

::: pybmad.bbu_add_a_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bbu_hom_voltage_calc

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L657)

::: pybmad.bbu_hom_voltage_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bbu_remove_head_bunch

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L630)

::: pybmad.bbu_remove_head_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bbu_setup

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L82)

::: pybmad.bbu_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bbu_track_a_stage

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L453)

::: pybmad.bbu_track_a_stage
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bbu_track_all

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L213)

::: pybmad.bbu_track_all
    options:
      show_root_heading: false
      show_root_toc_entry: false

### check_rf_freq

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L757)

::: pybmad.check_rf_freq
    options:
      show_root_heading: false
      show_root_toc_entry: false

### count_lines_in_file

Fortran source: [`bsim/code/count_lines_in_file.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/count_lines_in_file.f90#L5)

::: pybmad.count_lines_in_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hom_voltage

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L819)

::: pybmad.hom_voltage
    options:
      show_root_heading: false
      show_root_toc_entry: false

### insert_phase_trombone

Fortran source: [`bsim/code/bsim_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bsim_interface.f90#L18)

::: pybmad.insert_phase_trombone
    options:
      show_root_heading: false
      show_root_toc_entry: false

### logical_to_python

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L831)

::: pybmad.logical_to_python
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rf_cav_names

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L732)

::: pybmad.rf_cav_names
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_tune_3d

Fortran source: [`bsim/code/bsim_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bsim_interface.f90#L7)

::: pybmad.set_tune_3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_bunch_by_bunch_info

Fortran source: [`bsim/code/bbu_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/fe4ddb3f52a29e32330bd95dc4b7effa37c36092/bsim/code/bbu_track_mod.f90#L785)

::: pybmad.write_bunch_by_bunch_info
    options:
      show_root_heading: false
      show_root_toc_entry: false
