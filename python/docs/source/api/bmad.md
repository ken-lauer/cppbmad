# Bmad

Core Bmad particle accelerator simulation library.

## Classes (Fortran Structures)

::: pybmad.AcKickerFreqStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### AcKickerFreqStruct

Fortran struct: `ac_kicker_freq_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L698))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `f` | float |  |
| `amp` | float |  |
| `phi` | float |  |

::: pybmad.AcKickerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### AcKickerStruct

Fortran struct: `ac_kicker_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L704))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `amp_vs_time` | [1D array of AcKickerTimeStruct](bmad.md#ackickertimestruct) |  |
| `frequency` | [1D array of AcKickerFreqStruct](bmad.md#ackickerfreqstruct) |  |

::: pybmad.AcKickerTimeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### AcKickerTimeStruct

Fortran struct: `ac_kicker_time_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L692))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `amp` | float |  |
| `time` | float |  |
| `spline` | [SplineStruct](sim_utils.md#splinestruct) |  |

::: pybmad.AnormalModeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### AnormalModeStruct

Fortran struct: `anormal_mode_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2001))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `emittance` | float | Beam emittance (unnormalized). Includes vertical photon opening angle. |
| `emittance_no_vert` | float | Unnormalized beam emittance without the vertical photon opening angle taken into account. |
| `synch_int` | 1D array of float (shape: 4:6) | Synchrotron integrals |
| `j_damp` | float | damping partition number |
| `alpha_damp` | float | damping per turn |
| `chrom` | float | Chromaticity |
| `tune` | float | "Fractional" tune in radians |

::: pybmad.ApertureParamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ApertureParamStruct

Fortran struct: `aperture_param_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2153))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `min_angle` | float |  |
| `max_angle` | float |  |
| `n_angle` | int |  |
| `n_turn` | int | Number of turns a particle must survive. |
| `x_init` | float | Initial x coordinate to start with for theta_xy = 0. |
| `y_init` | float | Initial y coordinate to start with for theta_xy = pi/2. |
| `rel_accuracy` | float | Relative resolution of bracketed aperture. |
| `abs_accuracy` | float | Absolute resolution of bracketed aperture (meters). |
| `start_ele` | str | Element to start tracking at. |

::: pybmad.AperturePointStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### AperturePointStruct

Fortran struct: `aperture_point_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2144))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x` | float | (x,y) aperture point with respect to the reference orbit. |
| `y` | float | (x,y) aperture point with respect to the reference orbit. |
| `plane` | int | plane determining loss |
| `ix_ele` | int | ele index particle lost at |
| `i_turn` | int | turn particle lost at |

::: pybmad.ApertureScanStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ApertureScanStruct

Fortran struct: `aperture_scan_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2167))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `point` | [1D array of AperturePointStruct](bmad.md#aperturepointstruct) | Set of aperture points at different angles. |
| `ref_orb` | [CoordStruct](bmad.md#coordstruct) | Ref orbit around which the scan is made. |
| `pz_start` | float | Starting pz. |

::: pybmad.AstraLatticeParamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### AstraLatticeParamStruct

Fortran struct: `astra_lattice_param_struct` ([`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/astra_interface_mod.f90#L7))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `fieldmap_dimension` | int | Dimensions for field map. 1 or 3 |

::: pybmad.BaseLineEleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BaseLineEleStruct

Fortran struct: `base_line_ele_struct` ([`bmad/parsing/bmad_parser_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_struct.f90#L29))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Name of sequence or element |
| `tag` | str | Tag name. |
| `ix_multi` | int | Multipass indentifier |
| `orientation` | int | Element reversed? |
| `ix_ele_in_in_lat` | int |  |
| `ele_order_reflect` | bool | Part of reflection or reversed line? |

::: pybmad.BeamInitStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BeamInitStruct

Fortran struct: `beam_init_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1162))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `position_file` | str | File with particle positions. |
| `distribution_type` | 1D array of str (shape: 3) | distribution type (in x-px, y-py, and z-pz planes) "ELLIPSE", "KV", "GRID", "FILE", "RAN_GAUSS" or "" = "RAN_GAUSS" |
| `spin` | 1D array of float (shape: 3) | Spin (x, y, z) |
| `ellipse` | [1D array of EllipseBeamInitStruct (shape: 3)](bmad.md#ellipsebeaminitstruct) | Ellipse beam distribution |
| `KV` | [KvBeamInitStruct](bmad.md#kvbeaminitstruct) | KV beam distribution |
| `grid` | [1D array of GridBeamInitStruct (shape: 3)](bmad.md#gridbeaminitstruct) | Grid beam distribution |
| `center_jitter` | 1D array of float (shape: 6) | Bunch center rms jitter |
| `emit_jitter` | 1D array of float (shape: 2) | a and b bunch emittance rms jitter normalized to emittance |
| `sig_z_jitter` | float | bunch length RMS jitter |
| `sig_pz_jitter` | float | RMS pz spread jitter |
| `n_particle` | int | Number of particles per bunch. |
| `renorm_center` | bool | Renormalize centroid? |
| `renorm_sigma` | bool | Renormalize sigma? |
| `random_engine` | str | Or 'quasi'. Random number engine to use. |
| `random_gauss_converter` | str | Or 'quick' or 'exact'. Uniform to gauss conversion method. |
| `random_sigma_cutoff` | float | Cut-off in sigmas. |
| `a_norm_emit` | float | a-mode normalized emittance (emit * beta * gamma) |
| `b_norm_emit` | float | b-mode normalized emittance (emit * beta * gamma) |
| `a_emit` | float | a-mode emittance |
| `b_emit` | float | b-mode emittance |
| `dPz_dz` | float | Correlation of Pz with long position. |
| `center` | 1D array of float (shape: 6) | Bench phase space center offset relative to reference. |
| `t_offset` | float | Time center offset |
| `dt_bunch` | float | Time between bunches. |
| `sig_z` | float | Z sigma in m. |
| `sig_pz` | float | pz sigma |
| `bunch_charge` | float | charge (Coul) in a bunch. |
| `n_bunch` | int | Number of bunches. |
| `ix_turn` | int | Turn index used to adjust particles time if needed. |
| `species` | str | "positron", etc. "" => use referece particle. |
| `full_6D_coupling_calc` | bool | Use V from 6x6 1-turn mat to match distribution? Else use 4x4 1-turn mat used. |
| `use_particle_start` | bool | Use lat%particle_start instead of beam_init%center, %t_offset, and %spin? |
| `use_t_coords` | bool | If true, the distributions will be taken as in t-coordinates |
| `use_z_as_t` | bool | Only used if  use_t_coords = .true. If true,  z describes the t distribution If false, z describes the s distribution |
| `file_name` | str | OLD!! DO NOT USE!! |

::: pybmad.BeamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BeamStruct

Fortran struct: `beam_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1137))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `bunch` | [1D array of BunchStruct](bmad.md#bunchstruct) |  |

::: pybmad.BmadCommonStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BmadCommonStruct

Fortran struct: `bmad_common_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2313))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `max_aperture_limit` | float | Max Aperture. |
| `d_orb` | 1D array of float (shape: 6) | Orbit deltas for the mat6 via tracking calc. |
| `default_ds_step` | float | Default integration step for eles without an explicit step calc. |
| `significant_length` | float | meter |
| `rel_tol_tracking` | float | Closed orbit relative tolerance. |
| `abs_tol_tracking` | float | Closed orbit absolute tolerance. |
| `rel_tol_adaptive_tracking` | float | Runge-Kutta tracking relative tolerance. |
| `abs_tol_adaptive_tracking` | float | Runge-Kutta tracking absolute tolerance. |
| `init_ds_adaptive_tracking` | float | Initial step size |
| `min_ds_adaptive_tracking` | float | Min step size to take. |
| `fatal_ds_adaptive_tracking` | float | If actual step size is below this particle is lost. |
| `autoscale_amp_abs_tol` | float | Autoscale absolute amplitude tolerance (eV). |
| `autoscale_amp_rel_tol` | float | Autoscale relative amplitude tolerance |
| `autoscale_phase_tol` | float | Autoscale phase tolerance. |
| `electric_dipole_moment` | float | Particle's EDM. Call set_ptc to transfer value to PTC. |
| `synch_rad_scale` | float | Synch radiation kick scale. 1 => normal, 0 => no kicks. |
| `sad_eps_scale` | float | Used in sad_mult step length calc. |
| `sad_amp_max` | float | Used in sad_mult step length calc. |
| `sad_n_div_max` | int | Used in sad_mult step length calc. |
| `taylor_order` | int | Taylor order to use. 0 -> default = ptc_private%taylor_order_saved. |
| `runge_kutta_order` | int | Runge Kutta order. |
| `default_integ_order` | int | PTC integration order. |
| `max_num_runge_kutta_step` | int | Maximum number of RK steps before particle is considered lost. |
| `rf_phase_below_transition_ref` | bool | Autoscale uses below transition stable point for RFCavities? |
| `sr_wakes_on` | bool | Short range wakefields? |
| `lr_wakes_on` | bool | Long range wakefields |
| `auto_bookkeeper` | bool | Deprecated and no longer used. |
| `high_energy_space_charge_on` | bool | High energy space charge effect switch. |
| `high_energy_space_charge_linear` | bool | High energy space charge effect switch. |
| `csr_and_space_charge_on` | bool | Space charge switch. |
| `spin_tracking_on` | bool | spin tracking? |
| `spin_sokolov_ternov_flipping_on` | bool | Spin flipping during synchrotron radiation emission? |
| `radiation_damping_on` | bool | Radiation damping toggle. |
| `radiation_zero_average` | bool | Shift damping to be zero on the zero orbit to get rid of sawtooth? |
| `radiation_fluctuations_on` | bool | Radiation fluctuations toggle. |
| `conserve_taylor_maps` | bool | Enable bookkeeper to set ele%taylor_map_includes_offsets = F? |
| `absolute_time_tracking` | bool | Absolute or relative time tracking? |
| `absolute_time_ref_shift` | bool | Apply reference time shift when using absolute time tracking? |
| `convert_to_kinetic_momentum` | bool | Cancel kicks due to finite vector potential when doing symplectic tracking? Set to True to test symp_lie_bmad against runge_kutta. |
| `normalize_twiss` | bool | Normalize matrix when computing Twiss for off-energy ref? |
| `aperture_limit_on` | bool | Use apertures in tracking? |
| `spin_n0_direction_user_set` | bool | User sets direction of n0 for closed geometry branches? |
| `debug` | bool | Used for code debugging. |

::: pybmad.BmadNormalFormStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BmadNormalFormStruct

Fortran struct: `bmad_normal_form_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1589))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele_origin` | [EleStruct](bmad.md#elestruct) | Element at which the on-turn map was created. |
| `M` | [1D array of TaylorStruct (shape: 6)](bmad.md#taylorstruct) | One-turn taylor map: M = A o N o A_inv, N = exp(:h:) |
| `A` | [1D array of TaylorStruct (shape: 6)](bmad.md#taylorstruct) | Map from Floquet -> Lab coordinates |
| `A_inv` | [1D array of TaylorStruct (shape: 6)](bmad.md#taylorstruct) | Map from Lab -> Floquet coordinates |
| `dhdj` | [1D array of TaylorStruct (shape: 6)](bmad.md#taylorstruct) | Nonlinear tune function operating on Floquet coordinates |
| `F` | [1D array of ComplexTaylorStruct (shape: 6)](bmad.md#complextaylorstruct) | Vector field factorization in phasor basis: |
| `L` | [1D array of ComplexTaylorStruct (shape: 6)](bmad.md#complextaylorstruct) | L component |
| `h` | [1D array of ResonanceHStruct](bmad.md#resonancehstruct) |  |

::: pybmad.BookkeepingStateStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BookkeepingStateStruct

Fortran struct: `bookkeeping_state_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L958))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `attributes` | int | Element dependent attributes: super_ok$, ok$ or stale$ |
| `control` | int | Lord/slave bookkeeping status: super_ok$, ok$ or stale$ |
| `floor_position` | int | Global (floor) geometry: super_ok$, ok$ or stale$ |
| `s_position` | int | Longitudinal position & element length: super_ok$, ok$ or stale$ |
| `ref_energy` | int | Reference energy and ref time: super_ok$, ok$ or stale$ |
| `mat6` | int | Linear transfer map status: super_ok$, ok$ or stale$ |
| `rad_int` | int | Radiation integrals cache status |
| `ptc` | int | Associated PTC fibre (or layout) status. |
| `has_misalign` | bool | Used to avoid unnecessary calls to offset_particle. |

::: pybmad.BpmPhaseCouplingStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BpmPhaseCouplingStruct

Fortran struct: `bpm_phase_coupling_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L589))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `K_22a` | float | In-phase y/x for a-mode oscillations. |
| `K_12a` | float | Out-of-phase y/x for a-mode oscillations. |
| `K_11b` | float | In-phase x/y for b-mode oscillations. |
| `K_12b` | float | Out-of-phase x/y for b-mode oscillations. |
| `Cbar22_a` | float | Cbar22 as calculated from K_22a. |
| `Cbar12_a` | float | Cbar12 as calculated from K_12a. |
| `Cbar11_b` | float | Cbar11 as calculated from K_11b. |
| `Cbar12_b` | float | Cbar12 as calculated from K_12b. |
| `phi_a` | float | a-mode betatron phase. |
| `phi_b` | float | b-mode betatron phase. |

::: pybmad.BranchPointerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BranchPointerStruct

Fortran struct: `branch_pointer_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L931))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `branch` | [BranchStruct](bmad.md#branchstruct) |  |

::: pybmad.BranchStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BranchStruct

Fortran struct: `branch_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1618))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Name of line that defines the branch. |
| `ix_branch` | int | Index of this branch. 0 => Main branch |
| `ix_from_branch` | int | -1 => No creating fork element to this branch. |
| `ix_from_ele` | int | Index of creating fork element which forks to this branch. |
| `ix_to_ele` | int | Index of element in this branch that creating fork element forks to. |
| `ix_fixer` | int | Index of active fixer or beginning_ele element. |
| `n_ele_track` | int |  |
| `n_ele_max` | int |  |
| `lat` | [LatStruct](bmad.md#latstruct) |  |
| `a` | [ModeInfoStruct](bmad.md#modeinfostruct) | Note: Tunes are the fractional part. |
| `b` | [ModeInfoStruct](bmad.md#modeinfostruct) | Note: Tunes are the fractional part. |
| `z` | [ModeInfoStruct](bmad.md#modeinfostruct) | Note: Tunes are the fractional part. |
| `ele` | [1D array of EleStruct](bmad.md#elestruct) |  |
| `param` | [LatParamStruct](bmad.md#latparamstruct) |  |
| `particle_start` | [CoordStruct](bmad.md#coordstruct) |  |
| `wall3d` | 1D array of Wall3dStruct |  |
| `ptc` | [PtcBranch1Struct](bmad.md#ptcbranch1struct) | Pointer to layout. Note: ptc info not transferred with "branch1 = branch2" set. |
| `b_logic` | bool | For Bmad internal use only. |

::: pybmad.BunchParamsStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BunchParamsStruct

Fortran struct: `bunch_params_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1212))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `centroid` | [CoordStruct](bmad.md#coordstruct) | Lab frame |
| `x` | [TwissStruct](bmad.md#twissstruct) | Projected Twiss parameters |
| `y` | [TwissStruct](bmad.md#twissstruct) | Projected Twiss parameters |
| `z` | [TwissStruct](bmad.md#twissstruct) | Projected Twiss parameters |
| `a` | [TwissStruct](bmad.md#twissstruct) | Normal mode twiss parameters |
| `b` | [TwissStruct](bmad.md#twissstruct) | Normal mode twiss parameters |
| `c` | [TwissStruct](bmad.md#twissstruct) | Normal mode twiss parameters |
| `sigma` | 2D array of float (shape: 6,6) | beam size matrix |
| `rel_max` | 1D array of float (shape: 7) | Max orbit relative to centroid. 7 -> time. |
| `rel_min` | 1D array of float (shape: 7) | Min orbit relative to_centroid. 7 -> time. |
| `s` | float | Longitudinal position. |
| `t` | float | Time. |
| `sigma_t` | float | RMS of time spread. |
| `charge_live` | float | Charge of all non-lost particle |
| `charge_tot` | float | Charge of all particles. |
| `n_particle_tot` | int | Total number of particles |
| `n_particle_live` | int | Number of non-lost particles |
| `n_particle_lost_in_ele` | int | Number lost in element (not calculated by Bmad) |
| `n_good_steps` | int | Number of good steps (set when tracking with space charge) |
| `n_bad_steps` | int | Number of bad steps (set when tracking with space charge) |
| `ix_ele` | int | Lattice element where params evaluated at. |
| `location` | int | Location in element: upstream_end$, inside$, or downstream_end$ |
| `twiss_valid` | bool | Is the data here valid? Note: IF there is no energy variation (RF off) twiss_valid may be true but in this case the z-twiss will not be valid. |

::: pybmad.BunchStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BunchStruct

Fortran struct: `bunch_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1117))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `particle` | [1D array of CoordStruct](bmad.md#coordstruct) |  |
| `ix_z` | 1D array of int | bunch%ix_z(1) is index of head particle, etc. |
| `charge_tot` | float | Total charge in a bunch (Coul). |
| `charge_live` | float | Charge of live particles (Coul). |
| `z_center` | float | Longitudinal center of bunch at creation time. Note: Generally, z_center of bunch #1 is 0 and z_center of the other bunches is negative. |
| `t_center` | float | Center of bunch at creation time relative to head bunch. |
| `t0` | float | Used by track1_bunch_space_charge for tracking so particles have constant t. |
| `drift_between_t_and_s` | bool | Drift (ignore any fields) instead of tracking to speed up the calculation? This can only be done under certain circumstances. |
| `ix_ele` | int | Nominal element bunch is at. But, EG, dead particles can be someplace else. |
| `ix_bunch` | int | Bunch index. Head bunch = 1, etc. |
| `ix_turn` | int | Turn index for long term tracking. ix_turn = 0 before end of first turn, etc. |
| `n_live` | int |  |
| `n_good` | int | Number of accepted steps when using adaptive step size control. |
| `n_bad` | int | Number of rejected steps when using adaptive step size control. |

::: pybmad.BunchTrackStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BunchTrackStruct

Fortran struct: `bunch_track_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1238))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `pt` | [1D array of BunchParamsStruct](bmad.md#bunchparamsstruct) | Array indexed from 0 |
| `ds_save` | float | Min distance between points. |
| `n_pt` | int | Track upper bound |

::: pybmad.CartesianMapStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CartesianMapStruct

Fortran struct: `cartesian_map_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L732))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `field_scale` | float | Factor to scale the fields by |
| `r0` | 1D array of float (shape: 3) | Field origin offset. |
| `master_parameter` | int | Master parameter in ele%value(:) array to use for scaling the field. |
| `ele_anchor_pt` | int | anchor_beginning$, anchor_center$, or anchor_end$ |
| `field_type` | int | or electric$ |
| `ptr` | [CartesianMapTermStruct](bmad.md#cartesianmaptermstruct) |  |

::: pybmad.CartesianMapTerm1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CartesianMapTerm1Struct

Fortran struct: `cartesian_map_term1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L718))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | float |  |
| `kx` | float |  |
| `ky` | float |  |
| `kz` | float |  |
| `x0` | float |  |
| `y0` | float |  |
| `phi_z` | float |  |
| `family` | int | family_x$, etc. |
| `form` | int | hyper_y$, etc. |

::: pybmad.CartesianMapTermStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CartesianMapTermStruct

Fortran struct: `cartesian_map_term_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L726))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `file` | str | Input file name. Used also as ID for instances. |
| `n_link` | int | For memory management of %term |
| `term` | [1D array of CartesianMapTerm1Struct](bmad.md#cartesianmapterm1struct) |  |

::: pybmad.ComplexTaylorStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ComplexTaylorStruct

Fortran struct: `complex_taylor_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L509))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ref` | complex |  |
| `term` | [1D array of ComplexTaylorTermStruct](bmad.md#complextaylortermstruct) |  |

::: pybmad.ComplexTaylorTermStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ComplexTaylorTermStruct

Fortran struct: `complex_taylor_term_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L491))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | complex |  |
| `expn` | 1D array of int (shape: 6) |  |

::: pybmad.ControlRamp1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ControlRamp1Struct

Fortran struct: `control_ramp1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1380))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `y_knot` | 1D array of float |  |
| `stack` | [1D array of ExpressionAtomStruct](bmad.md#expressionatomstruct) | Evaluation stack |
| `attribute` | str | Name of attribute controlled. Set to "FIELD_OVERLAPS" for field overlaps. |
| `slave_name` | str | Name of slave. |
| `is_controller` | bool | Is the slave a controller? If so bookkeeping is different. |

::: pybmad.ControlStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ControlStruct

Fortran struct: `control_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1362))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `value` | float | Used by group, and overlay elements. |
| `y_knot` | 1D array of float |  |
| `stack` | [1D array of ExpressionAtomStruct](bmad.md#expressionatomstruct) | Evaluation stack |
| `slave` | [LatEleLocStruct](bmad.md#latelelocstruct) |  |
| `lord` | [LatEleLocStruct](bmad.md#latelelocstruct) |  |
| `slave_name` | str | Name of slave. |
| `attribute` | str | Name of attribute controlled. Set to "FIELD_OVERLAPS" for field overlaps. Set to "INPUT" or "OUTPUT" for feedback slaves. |
| `ix_attrib` | int | Index of attribute controlled. See note above! |

::: pybmad.ControlVar1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ControlVar1Struct

Fortran struct: `control_var1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1374))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `value` | float |  |
| `old_value` | float |  |

::: pybmad.ControllerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ControllerStruct

Fortran struct: `controller_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1397))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `var` | [1D array of ControlVar1Struct](bmad.md#controlvar1struct) |  |
| `ramp` | [1D array of ControlRamp1Struct](bmad.md#controlramp1struct) | For ramper lord elements |
| `ramper_lord` | [1D array of RamperLordStruct](bmad.md#ramperlordstruct) | Ramper lord info for this slave |
| `x_knot` | 1D array of float |  |

::: pybmad.ConverterDir1DStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ConverterDir1DStruct

Fortran struct: `converter_dir_1D_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1266))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `pc_out` | float | pc_out value at fit |
| `poly` | 1D array of float (shape: 0:4) | param(r) = Sum: poly(i) * r^i |

::: pybmad.ConverterDir2DStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ConverterDir2DStruct

Fortran struct: `converter_dir_2D_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1271))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `k` | float |  |
| `poly` | 1D array of float (shape: 0:3) |  |

::: pybmad.ConverterDirCoefStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ConverterDirCoefStruct

Fortran struct: `converter_dir_coef_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1276))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `fit_1d_r` | 1D array of ConverterDir1dStruct |  |
| `fit_2d_r` | ConverterDir2dStruct |  |
| `fit_2d_pc` | ConverterDir2dStruct |  |
| `c0` | float |  |

::: pybmad.ConverterDirectionOutStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ConverterDirectionOutStruct

Fortran struct: `converter_direction_out_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1283))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `beta` | [ConverterDirCoefStruct](bmad.md#converterdircoefstruct) |  |
| `alpha_x` | [ConverterDirCoefStruct](bmad.md#converterdircoefstruct) |  |
| `alpha_y` | [ConverterDirCoefStruct](bmad.md#converterdircoefstruct) |  |
| `dxds_min` | [ConverterDirCoefStruct](bmad.md#converterdircoefstruct) |  |
| `dxds_max` | [ConverterDirCoefStruct](bmad.md#converterdircoefstruct) |  |
| `dyds_max` | [ConverterDirCoefStruct](bmad.md#converterdircoefstruct) |  |
| `c_x` | [ConverterDirCoefStruct](bmad.md#converterdircoefstruct) |  |

::: pybmad.ConverterDistributionStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ConverterDistributionStruct

Fortran struct: `converter_distribution_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1344))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `thickness` | float |  |
| `sub_dist` | [1D array of ConverterSubDistributionStruct](bmad.md#convertersubdistributionstruct) | Distribution at various pc_in values. |

::: pybmad.ConverterProbPcRStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ConverterProbPcRStruct

Fortran struct: `converter_prob_pc_r_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1249))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `pc_out` | 1D array of float | Grid pc_out values. |
| `r` | 1D array of float | Grid r_out values. |
| `prob` | 2D array of float | Probability grid. |
| `spin_z` | 2D array of float | Z polarization grid. Stuff below is calculated rather than read in from the lattice file. |
| `pc_out_min` | float |  |
| `pc_out_max` | float |  |
| `integrated_prob` | float | Integrated probability over (pc_out, r) with restrictions factered in. |
| `p_norm` | 2D array of float | Normalized probability taking into account. angle_out_max, pc_out_min, and pc_out_max restrictions. |
| `integ_pc_out` | 1D array of float | Normalized probability integrated from min pc_out up. |
| `integ_r` | 2D array of float |  |
| `integ_r_ave` | 1D array of float |  |

::: pybmad.ConverterStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ConverterStruct

Fortran struct: `converter_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1352))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `species_out` | int | Output species |
| `material_type` | str |  |
| `dist` | [1D array of ConverterDistributionStruct](bmad.md#converterdistributionstruct) | Distribution at various thicknesses |

::: pybmad.ConverterSubDistributionStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ConverterSubDistributionStruct

Fortran struct: `converter_sub_distribution_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1292))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `pc_in` | float |  |
| `spin_in` | 1D array of float (shape: 3) |  |
| `prob_pc_r` | [ConverterProbPcRStruct](bmad.md#converterprobpcrstruct) |  |
| `dir_out` | [ConverterDirectionOutStruct](bmad.md#converterdirectionoutstruct) |  |

::: pybmad.CoordArrayStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CoordArrayStruct

Fortran struct: `coord_array_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L581))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `orbit` | [1D array of CoordStruct](bmad.md#coordstruct) |  |

::: pybmad.CoordStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CoordStruct

Fortran struct: `coord_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L551))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `vec` | 1D array of float (shape: 6) | (x, px, y, py, z, pz). Generally phase space for charged particles. See Bmad manual. |
| `s` | float | Longitudinal position |
| `t` | float | Absolute time (not relative to reference). Note: Quad precision! |
| `spin` | 1D array of float (shape: 3) | Spin. |
| `field` | 1D array of float (shape: 2) | Photon E-field intensity (x,y). |
| `phase` | 1D array of float (shape: 2) | Photon E-field phase (x,y). For charged particles, phase(1) is RF phase. |
| `charge` | float | Macroparticle weight (which is different from particle species charge). For some space charge calcs the weight is in Coulombs. |
| `dt_ref` | float | Used in: * time tracking for computing z. * by coherent photons = path_length/c_light. |
| `r` | float | For general use. Not used by Bmad. |
| `p0c` | float | For non-photons: Reference momentum. For photons: Photon momentum (not reference). |
| `E_potential` | float | Potential energy. |
| `beta` | float | Velocity / c_light. |
| `ix_ele` | int | Index of the lattice element the particle is in. May be -1 if element is not associated with a lattice. |
| `ix_branch` | int | Index of the lattice branch the particle is in. |
| `ix_turn` | int | Turn index for multiturn tracking. |
| `ix_user` | int | For general use, not used by Bmad. |
| `state` | int | alive$, lost$, lost_neg_x_aperture$, lost_pz$, etc. |
| `direction` | int | +1 or -1. Sign of longitudinal direction of motion (ds/dt). This is independent of the element orientation. |
| `time_dir` | int | +1 or -1. Time direction. -1 => Traveling backwards in time. |
| `species` | int | positron$, proton$, etc. |
| `location` | int | upstream_end$, inside$, or downstream_end$ |

::: pybmad.CrystalParamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CrystalParamStruct

Fortran struct: `crystal_param_struct` ([`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_utils_mod.f90#L9))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `cap_gamma` | float |  |
| `dtheta_sin_2theta` | float |  |
| `b_eff` | float |  |
| `wavelength` | float |  |
| `old_vvec` | 1D array of float (shape: 3) |  |
| `new_vvec` | 1D array of float (shape: 3) |  |

::: pybmad.CsrBunchSliceStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CsrBunchSliceStruct

Fortran struct: `csr_bunch_slice_struct` ([`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L34))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x0` | float | Transverse center of the particle distrubution |
| `y0` | float | Transverse center of the particle distrubution |
| `z0_edge` | float | Left (min z) edge of bin |
| `z1_edge` | float | Right (max z) edge of bin |
| `z_center` | float | z at center of bin. |
| `sig_x` | float | particle's RMS width |
| `sig_y` | float | particle's RMS width |
| `charge` | float | charge of the particles |
| `dcharge_density_dz` | float | Charge density gradient |
| `edge_dcharge_density_dz` | float | gradient between this and preceeding bin. [Evaluated at bin edge.] |
| `kick_csr` | float | CSR kick |
| `coef_lsc_plus` | 2D array of float (shape: 0:2,0:2) | LSC Kick coefs. |
| `coef_lsc_minus` | 2D array of float (shape: 0:2,0:2) | LSC Kick coefs. |
| `kick_lsc` | float |  |
| `n_particle` | float | Number of particles in slice can be a fraction since particles span multiple bins. |

::: pybmad.CsrEleInfoStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CsrEleInfoStruct

Fortran struct: `csr_ele_info_struct` ([`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L22))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele` | [EleStruct](bmad.md#elestruct) | lattice element |
| `orbit0` | [CoordStruct](bmad.md#coordstruct) | centroid orbit at entrance/exit ends |
| `orbit1` | [CoordStruct](bmad.md#coordstruct) | centroid orbit at entrance/exit ends |
| `floor0` | [FloorPositionStruct](bmad.md#floorpositionstruct) | Floor position of centroid at entrance/exit ends |
| `floor1` | [FloorPositionStruct](bmad.md#floorpositionstruct) | Floor position of centroid at entrance/exit ends |
| `ref_floor0` | [FloorPositionStruct](bmad.md#floorpositionstruct) | Floor position of element ref coords at entrance/exit ends |
| `ref_floor1` | [FloorPositionStruct](bmad.md#floorpositionstruct) | Floor position of element ref coords at entrance/exit ends |
| `spline` | [SplineStruct](sim_utils.md#splinestruct) | Spline for centroid orbit. spline%x = distance along chord. The spline is zero at the ends by construction. |
| `theta_chord` | float | Reference angle of chord in z-x plane |
| `L_chord` | float | Chord Length. Negative if bunch moves backwards in element. |
| `dL_s` | float | L_s(of element) - L_chord |

::: pybmad.CsrKick1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CsrKick1Struct

Fortran struct: `csr_kick1_struct` ([`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L57))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `I_csr` | float | Kick integral. |
| `I_int_csr` | float | Integrated Kick integral. |
| `image_kick_csr` | float | kick. |
| `L_vec` | 1D array of float (shape: 3) | L vector in global coordinates. |
| `L` | float | Distance between source and kick points. |
| `dL` | float | = epsilon_L = Ls - L |
| `dz_particles` | float | Kicked particle - source particle position at constant time. |
| `s_chord_source` | float | Source point coordinate along chord. |
| `theta_L` | float | Angle of L vector |
| `theta_sl` | float | Angle between velocity of particle at source pt and L |
| `theta_lk` | float | Angle between L and velocity of kicked particle |
| `ix_ele_source` | int | Source element index. |
| `floor_s` | [FloorPositionStruct](bmad.md#floorpositionstruct) | Floor position of source pt |

::: pybmad.CsrParticlePositionStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CsrParticlePositionStruct

Fortran struct: `csr_particle_position_struct` ([`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L73))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `r` | 1D array of float (shape: 3) | particle position |
| `charge` | float | particle charge |

::: pybmad.CsrStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CsrStruct

Fortran struct: `csr_struct` ([`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L78))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `gamma` | float | Relativistic gamma factor. |
| `gamma2` | float | Relativistic gamma factor. |
| `rel_mass` | float | m_particle / m_electron |
| `beta` | float | Relativistic beta factor. |
| `dz_slice` | float | Bin width |
| `ds_track_step` | float | True step size |
| `s_kick` | float | Kick point longitudinal location (element ref coords) from entrance end |
| `s_chord_kick` | float | Kick point along beam centroid line |
| `y_source` | float | Height of source particle. |
| `kick_factor` | float | Coefficient to scale the kick |
| `actual_track_step` | float | ds_track_step scalled by Length_centroid_chord / Length_element ratio |
| `x0_bunch` | float | Bunch centroid |
| `y0_bunch` | float | Bunch centroid |
| `floor_k` | [FloorPositionStruct](bmad.md#floorpositionstruct) | Floor coords at kick point |
| `species` | int | Particle type |
| `ix_ele_kick` | int | Same as element being tracked through. |
| `slice` | [1D array of CsrBunchSliceStruct](bmad.md#csrbunchslicestruct) | slice(i) refers to the i^th bunch slice. |
| `kick1` | [1D array of CsrKick1Struct](bmad.md#csrkick1struct) | kick1(i) referes to the kick between two slices i bins apart. |
| `eleinfo` | [1D array of CsrEleInfoStruct](bmad.md#csreleinfostruct) | Element-by-element information. |
| `kick_ele` | [EleStruct](bmad.md#elestruct) | Element where the kick pt is == ele tracked through. |
| `mesh3d` | Mesh3dStruct |  |
| `position` | [1D array of CsrParticlePositionStruct](bmad.md#csrparticlepositionstruct) |  |

::: pybmad.CylindricalMapStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CylindricalMapStruct

Fortran struct: `cylindrical_map_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L754))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `m` | int | Azimuthal Mode: varies as cos(m*phi - theta0_azimuth) |
| `harmonic` | int | Harmonic of fundamental |
| `phi0_fieldmap` | float | Mode oscillates as: twopi * (f * t + phi0_fieldmap) |
| `theta0_azimuth` | float | Azimuthal ((x, y) plane) orientation of mode. |
| `field_scale` | float | Factor to scale the fields by |
| `master_parameter` | int | Master parameter in ele%value(:) array to use for scaling the field. |
| `ele_anchor_pt` | int | anchor_beginning$, anchor_center$, or anchor_end$ |
| `dz` | float | Distance between sampled field points. |
| `r0` | 1D array of float (shape: 3) | Field origin offset. |
| `ptr` | [CylindricalMapTermStruct](bmad.md#cylindricalmaptermstruct) |  |

::: pybmad.CylindricalMapTerm1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CylindricalMapTerm1Struct

Fortran struct: `cylindrical_map_term1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L743))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `e_coef` | complex |  |
| `b_coef` | complex |  |

::: pybmad.CylindricalMapTermStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CylindricalMapTermStruct

Fortran struct: `cylindrical_map_term_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L748))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `file` | str | Input file name. Used also as ID for instances. |
| `n_link` | int | For memory management of this structure |
| `term` | [1D array of CylindricalMapTerm1Struct](bmad.md#cylindricalmapterm1struct) |  |

::: pybmad.DiffuseParamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### DiffuseParamStruct

Fortran struct: `diffuse_param_struct` ([`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L13))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x` | float |  |
| `y` | float |  |
| `lambda_` | float |  |
| `c_norm` | float |  |
| `chx_norm` | float |  |
| `prob_spline` | [1D array of SplineStruct (shape: 0:50)](sim_utils.md#splinestruct) |  |
| `n_pt_spline` | int |  |

::: pybmad.EleAttributeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### EleAttributeStruct

Fortran struct: `ele_attribute_struct` ([`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L25))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `state` | int | See above. |
| `kind` | int | Is_switch$, is_real$, etc. See attribute_type routine. |
| `units` | str | EG: 'T*m'. |
| `ix_attrib` | int | Attribute index. Frequently will be where in the ele%value(:) array the attribute is. |
| `value` | float | Used by type_ele. |

::: pybmad.ElePointerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ElePointerStruct

Fortran struct: `ele_pointer_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L923))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele` | [EleStruct](bmad.md#elestruct) |  |
| `loc` | [LatEleLocStruct](bmad.md#latelelocstruct) |  |
| `id` | int | For general use. Not used by Bmad. In particular, used by Tao to designate universe ele is in. |

::: pybmad.EleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### EleStruct

Fortran struct: `ele_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1422))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | name of element. |
| `type` | str | type name. |
| `alias` | str | Another name. |
| `component_name` | str | Used by overlays, multipass patch, etc. |
| `descrip` | str | Description string. |
| `a` | [TwissStruct](bmad.md#twissstruct) | Twiss parameters at end of element |
| `b` | [TwissStruct](bmad.md#twissstruct) | Twiss parameters at end of element |
| `z` | [TwissStruct](bmad.md#twissstruct) | Twiss parameters at end of element |
| `x` | [XyDispStruct](bmad.md#xydispstruct) | Projected dispersions. |
| `y` | [XyDispStruct](bmad.md#xydispstruct) | Projected dispersions. |
| `ac_kick` | [AcKickerStruct](bmad.md#ackickerstruct) | ac_kicker element parameters. |
| `bookkeeping_state` | [BookkeepingStateStruct](bmad.md#bookkeepingstatestruct) | Attribute bookkeeping |
| `branch` | [BranchStruct](bmad.md#branchstruct) | Pointer to branch containing element. |
| `control` | [ControllerStruct](bmad.md#controllerstruct) | group & overlay variables. |
| `converter` | [ConverterStruct](bmad.md#converterstruct) | EG: Positron converter in linac. |
| `rf` | [RfEleStruct](bmad.md#rfelestruct) | RF parameters. |
| `foil` | [FoilStruct](bmad.md#foilstruct) |  |
| `lord` | [EleStruct](bmad.md#elestruct) | Pointer to a slice lord. |
| `floor` | [FloorPositionStruct](bmad.md#floorpositionstruct) |  |
| `high_energy_space_charge` | [HighEnergySpaceChargeStruct](bmad.md#highenergyspacechargestruct) |  |
| `mode3` | [Mode3Struct](bmad.md#mode3struct) | 6D normal mode structure. |
| `photon` | [PhotonElementStruct](bmad.md#photonelementstruct) |  |
| `multipole_cache` | [MultipoleCacheStruct](bmad.md#multipolecachestruct) |  |
| `rad_map` | [RadMapEleStruct](bmad.md#radmapelestruct) | Radiation kick parameters Note: The reference orbits for spin and orbit Taylor maps are not necessarily the same. For example, Sprint spin Taylor maps can be with respect to the zero orbit independent of the orbital map. |
| `taylor` | [1D array of TaylorStruct (shape: 6)](bmad.md#taylorstruct) | Phase space Taylor map. |
| `spin_taylor_ref_orb_in` | 1D array of float (shape: 6) |  |
| `spin_taylor` | [1D array of TaylorStruct (shape: 0:3)](bmad.md#taylorstruct) | Quaternion Spin Taylor map. |
| `wake` | [WakeStruct](bmad.md#wakestruct) | Wakes |
| `wall3d` | 1D array of Wall3dStruct | Chamber or capillary wall E/M field structs. |
| `cartesian_map` | [1D array of CartesianMapStruct](bmad.md#cartesianmapstruct) | Used to define E/M fields |
| `cylindrical_map` | [1D array of CylindricalMapStruct](bmad.md#cylindricalmapstruct) | Used to define E/M fields |
| `gen_gradients` | [1D array of GenGradientsStruct](bmad.md#gengradientsstruct) | Used to define E/M fields. |
| `grid_field` | [1D array of GridFieldStruct](bmad.md#gridfieldstruct) | Used to define E/M fields. The difference between map_ref_orb and time_ref_orb is that map_ref_orb is the reference orbit for the 1st order spin/orbit map which, in general, is non-zero while time_ref_orb follows the reference particle which is generally the zero orbit (non-zero, for example, in the second slice of a sliced wiggler). |
| `map_ref_orb_in` | [CoordStruct](bmad.md#coordstruct) | Entrance end transfer map ref orbit |
| `map_ref_orb_out` | [CoordStruct](bmad.md#coordstruct) | Exit end transfer map ref orbit |
| `time_ref_orb_in` | [CoordStruct](bmad.md#coordstruct) | Reference orbit at entrance end for ref_time calc. |
| `time_ref_orb_out` | [CoordStruct](bmad.md#coordstruct) | Reference orbit at exit end for ref_time calc. |
| `value` | 1D array of float (shape: num_ele_attrib$) | attribute values. |
| `old_value` | 1D array of float (shape: num_ele_attrib$) | Used to see if %value(:) array has changed. Note: The reference orbit for spin/orbit matrices is %map_ref_orb_in/out |
| `spin_q` | 2D array of float (shape: 0:3,0:6) | 0th and 1st order Spin transport quaternion. |
| `vec0` | 1D array of float (shape: 6) | 0th order transport vector. |
| `mat6` | 2D array of float (shape: 6,6) | 1st order transport matrix. |
| `c_mat` | 2D array of float (shape: 2,2) | 2x2 C coupling matrix |
| `dc_mat_dpz` | 2D array of float (shape: 2,2) | d(c_mat)/dpz variation. |
| `gamma_c` | float | gamma associated with C matrix |
| `s_start` | float | longitudinal ref position at entrance_end |
| `s` | float | longitudinal ref position at the exit end. |
| `ref_time` | float | Time ref particle passes exit end. |
| `a_pole` | 1D array of float | knl for multipole elements. |
| `b_pole` | 1D array of float | tilt for multipole elements. |
| `a_pole_elec` | 1D array of float | Electrostatic multipoles. ksnl for multipole elements. |
| `b_pole_elec` | 1D array of float | Electrostatic multipoles. |
| `custom` | 1D array of float | Custom attributes. |
| `r` | 3D array of float | For general use. Not used by Bmad. |
| `key` | int | Element class (quadrupole, etc.). |
| `sub_key` | int | Records bend input type. |
| `ix_ele` | int | Index in branch ele(0:) array. Set to ix_slice_slave$ = -2 for slice_slave$ elements. |
| `ix_branch` | int | Index in lat%branch(:) array. Note: lat%ele => lat%branch(0). |
| `lord_status` | int | Type of lord element this is. overlay_lord$, etc. |
| `n_slave` | int | Number of slaves (except field overlap slaves) of this element. |
| `n_slave_field` | int | Number of field slaves of this element. |
| `ix1_slave` | int | Pointer index to this element's slaves. |
| `slave_status` | int | Type of slave element this is. multipass_slave$, slice_slave$, etc. |
| `n_lord` | int | Number of lords (except field overlap and ramper lords). |
| `n_lord_field` | int | Number of field lords of this element. |
| `n_lord_ramper` | int | Number of ramper lords. |
| `ic1_lord` | int | Pointer index to this element's lords. |
| `ix_pointer` | int | For general use. Not used by Bmad. |
| `ixx` | int | Index for Bmad internal use. |
| `iyy` | int | Index for Bmad internal use. |
| `izz` | int | Index for Bmad internal use. |
| `mat6_calc_method` | int | taylor$, symp_lie_ptc$, etc. |
| `tracking_method` | int | taylor$, linear$, etc. |
| `spin_tracking_method` | int | symp_lie_ptc$, etc. |
| `csr_method` | int | or one_dim$ ("1_dim"), steady_state_3d$ |
| `space_charge_method` | int | slice$, slice_longitudinal$, slice_transverse$, fft_3D$, cathode_fft_3d$ |
| `ptc_integration_type` | int | drift_kick$, matrix_kick$, or ripken_kick$ |
| `field_calc` | int | no_field$, fieldmap$, refer_to_lords$, or custom$ |
| `aperture_at` | int | Aperture location: entrance_end$, ... |
| `aperture_type` | int | rectangular$, elliptical$, auto_aperture$, ... |
| `ref_species` | int | Reference species |
| `orientation` | int | -1 -> Element is longitudinally reversed. +1 -> Normal. |
| `symplectify` | bool | Symplectify mat6 matrices. |
| `mode_flip` | bool | Have the normal modes traded places? |
| `multipoles_on` | bool | For turning multipoles on/off |
| `scale_multipoles` | bool | Are ab_multipoles within other elements (EG: quads, etc.) scaled by the strength of the element? |
| `taylor_map_includes_offsets` | bool | Taylor map calculated with element misalignments? |
| `field_master` | bool | Calculate strength from the field value? |
| `is_on` | bool | For turning element on/off. |
| `logic` | bool | For general use. Not used by Bmad (except during lattice parsing). |
| `bmad_logic` | bool | For Bmad internal use only. |
| `select` | bool | For Bmad internal use only. |
| `offset_moves_aperture` | bool | element offsets affects aperture? ! final :: ele_finalizer |

::: pybmad.EllipseBeamInitStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### EllipseBeamInitStruct

Fortran struct: `ellipse_beam_init_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1141))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `part_per_ellipse` | int | number of particles per ellipse |
| `n_ellipse` | int | number of ellipses (>= 1) |
| `sigma_cutoff` | float | sigma cutoff of the representation |

::: pybmad.EmFieldStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### EmFieldStruct

Fortran struct: `em_field_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2052))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `E` | 1D array of float (shape: 3) | electric field. |
| `B` | 1D array of float (shape: 3) | magnetic field. |
| `dE` | 2D array of float (shape: 3,3) | electric field gradient. |
| `dB` | 2D array of float (shape: 3,3) | magnetic field gradient. |
| `phi` | float | Electric scalar potential. |
| `phi_B` | float | Magnetic scalar potential. |
| `A` | 1D array of float (shape: 3) | Magnetic vector potential. |

::: pybmad.ExpressionAtomStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ExpressionAtomStruct

Fortran struct: `expression_atom_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L57))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `type` | int | plus$, minum$, sin$, cos$, etc. To convert to string use: expression_op_name |
| `value` | float |  |

::: pybmad.ExpressionTreeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ExpressionTreeStruct

Fortran struct: `expression_tree_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L65))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `type` | int | plus$, minum$, sin$, cos$, etc. |
| `value` | float |  |
| `reverse_polish` | bool | Can the children node(:) array be in reverse Polish order? |
| `node` | [1D array of ExpressionTreeStruct](bmad.md#expressiontreestruct) | Child nodes. Note: Pointer used here since Ifort does not support allocatable. |

::: pybmad.ExtraParsingInfoStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ExtraParsingInfoStruct

Fortran struct: `extra_parsing_info_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2225))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ran_state` | [RandomStateStruct](sim_utils.md#randomstatestruct) |  |
| `ran_seed` | int |  |
| `undeterministic_ran_function_called` | bool | Used with bmad_com |
| `d_orb_set` | bool |  |
| `max_aperture_limit_set` | bool |  |
| `default_ds_step_set` | bool |  |
| `significant_length_set` | bool |  |
| `rel_tol_tracking_set` | bool |  |
| `abs_tol_tracking_set` | bool |  |
| `rel_tol_adaptive_tracking_set` | bool |  |
| `abs_tol_adaptive_tracking_set` | bool |  |
| `init_ds_adaptive_tracking_set` | bool |  |
| `min_ds_adaptive_tracking_set` | bool |  |
| `fatal_ds_adaptive_tracking_set` | bool |  |
| `synch_rad_scale_set` | bool |  |
| `autoscale_amp_abs_tol_set` | bool |  |
| `autoscale_amp_rel_tol_set` | bool |  |
| `autoscale_phase_tol_set` | bool |  |
| `rf_phase_below_transition_ref_set` | bool |  |
| `electric_dipole_moment_set` | bool |  |
| `taylor_order_set` | bool |  |
| `runge_kutta_order_set` | bool |  |
| `default_integ_order_set` | bool |  |
| `sr_wakes_on_set` | bool |  |
| `lr_wakes_on_set` | bool |  |
| `high_energy_space_charge_on_set` | bool |  |
| `high_energy_space_charge_linear_set` | bool |  |
| `csr_and_space_charge_on_set` | bool |  |
| `spin_tracking_on_set` | bool |  |
| `spin_sokolov_ternov_flipping_on_set` | bool |  |
| `radiation_damping_on_set` | bool |  |
| `radiation_zero_average_set` | bool |  |
| `radiation_fluctuations_on_set` | bool |  |
| `conserve_taylor_maps_set` | bool |  |
| `absolute_time_tracking_set` | bool |  |
| `absolute_time_ref_shift_set` | bool |  |
| `convert_to_kinetic_momentum_set` | bool |  |
| `aperture_limit_on_set` | bool |  |
| `normalize_twiss_set` | bool |  |
| `sad_eps_scale_set` | bool |  |
| `sad_amp_max_set` | bool |  |
| `sad_n_div_max_set` | bool |  |
| `max_num_runge_kutta_step_set` | bool |  |
| `spin_n0_direction_user_set_set` | bool |  |
| `debug_set` | bool | Used with space_charge_com |
| `ds_track_step_set` | bool |  |
| `dt_track_step_set` | bool |  |
| `cathode_strength_cutoff_set` | bool |  |
| `sc_rel_tol_tracking_set` | bool | For: space_charge_com%rel_tol_tracking |
| `sc_abs_tol_tracking_set` | bool | For: space_charge_com%abs_tol_tracking |
| `beam_chamber_height_set` | bool |  |
| `lsc_sigma_cutoff_set` | bool |  |
| `particle_sigma_cutoff_set` | bool |  |
| `space_charge_mesh_size_set` | bool |  |
| `csr3d_mesh_size_set` | bool |  |
| `n_bin_set` | bool |  |
| `particle_bin_span_set` | bool |  |
| `n_shield_images_set` | bool |  |
| `sc_min_in_bin_set` | bool |  |
| `lsc_kick_transverse_dependence_set` | bool |  |
| `sc_debug_set` | bool |  |
| `diagnostic_output_file_set` | bool | Used with ptc_com |
| `old_integrator_set` | bool |  |
| `use_orientation_patches_set` | bool |  |
| `print_info_messages_set` | bool |  |
| `max_fringe_order_set` | bool |  |
| `exact_model_set` | bool |  |
| `exact_misalign_set` | bool |  |
| `vertical_kick_set` | bool |  |
| `cut_factor_set` | bool |  |
| `translate_patch_drift_time_set` | bool |  |
| `pancake_symplectic_set` | bool |  |
| `pancake_canonical_set` | bool |  |

::: pybmad.FloorPositionStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### FloorPositionStruct

Fortran struct: `floor_position_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L855))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `r` | 1D array of float (shape: 3) | (x, y, z) offset from origin |
| `w` | 2D array of float (shape: 3,3) | W matrix. Columns are unit vectors of the frame axes. |
| `theta` | float | angular orientation consistent with W matrix |
| `phi` | float | angular orientation consistent with W matrix |
| `psi` | float | angular orientation consistent with W matrix |

::: pybmad.FoilStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### FoilStruct

Fortran struct: `foil_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1309))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `material` | [1D array of MaterialStruct](bmad.md#materialstruct) |  |

::: pybmad.FringeFieldInfoStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### FringeFieldInfoStruct

Fortran struct: `fringe_field_info_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L35))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `hard_ele` | [EleStruct](bmad.md#elestruct) |  |
| `s_edge_hard` | float |  |
| `ds_edge` | float | Distance from particle to edge in hard_ele frame. |
| `particle_at` | int | first_track_edge$, second_track_edge$, or none$ |
| `hard_location` | int | Particle location wrt hard_ele. Points to element in location(:). |
| `location` | 1D array of int | Particle location in an element. entrance_end$, inside$, or exit_end$ Elements in list are the tracking element or its lords. |
| `has_fringe` | bool | Has a fringe to worry about? |

::: pybmad.GenGradCurveStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GenGradCurveStruct

Fortran struct: `gen_grad_curve_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L775))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `kind` | int | gg_a$ (skew), gg_b$ (normal), or gg_bs$ (solenoid). |
| `n` | int | Azimuthal harmonic index (n = 0 for gg_bs$). |
| `m_max` | int | Max GG derivative order stored. deriv(iz, 0:m_max) are the GG derivatives d^m/ds^m; columns m_max+1:2*m_max+1 hold the interpolating spline extension (see n_spline_create). |
| `deriv` | 2D array of float | Range: (iz0:iz1, 0:2*m_max+1) |

::: pybmad.GenGradientsStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GenGradientsStruct

Fortran struct: `gen_gradients_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L784))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `file` | str | Input file name. Used also as ID for instances. |
| `curve` | [1D array of GenGradCurveStruct](bmad.md#gengradcurvestruct) |  |
| `ele_anchor_pt` | int | anchor_beginning$, anchor_center$, or anchor_end$ |
| `field_type` | int | or electric$ |
| `iz0` | int | curve%deriv(iz0:iz1, :) lower bound. |
| `iz1` | int | curve%deriv(iz0:iz1, :) upper bound. |
| `dz` | float | Point spacing between base planes. |
| `g_ref` | float | Reference-frame curvature 1/rho (0 => straight frame). Must be equal to g for a bend and zero for all else. |
| `r0` | 1D array of float (shape: 3) | Field origin relative to ele_anchor_pt. r0(1:2) = transverse expansion axis (GGCoefs origin), r0(3) = longitudinal offset. |
| `field_scale` | float | Factor to scale the fields by. |
| `master_parameter` | int | Master parameter in ele%value(:) array to use for scaling the field. |

::: pybmad.GgTaylorStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GgTaylorStruct

Fortran struct: `gg_taylor_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L839))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ref` | float |  |
| `term` | [1D array of GgTaylorTermStruct](bmad.md#ggtaylortermstruct) |  |

::: pybmad.GgTaylorTermStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GgTaylorTermStruct

Fortran struct: `gg_taylor_term_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L834))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | float |  |
| `expn` | 1D array of int (shape: 2) |  |

::: pybmad.GptLatParamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GptLatParamStruct

Fortran struct: `gpt_lat_param_struct` ([`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L7))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `fieldmap_dimension` | int | Dimensions for field map. 1 or 3 |
| `only_write_autophase_parameters` | bool | Option to only write phasing info |
| `gpt_filename` | str | Blank => Append '.gpt' to Bmad lattice file name. |
| `header_file_name` | str | Header file to include in gpt file. |
| `tracking_end_element` | str | Bmad lattice element name or index. |

::: pybmad.GridBeamInitStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GridBeamInitStruct

Fortran struct: `grid_beam_init_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1153))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_x` | int | Number of columns. |
| `n_px` | int | Number of rows. |
| `x_min` | float | Lower x limit. |
| `x_max` | float | Upper x limit. |
| `px_min` | float | Lower px limit. |
| `px_max` | float | Upper px limit. |

::: pybmad.GridFieldPt1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GridFieldPt1Struct

Fortran struct: `grid_field_pt1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L802))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `E` | 1D array of complex (shape: 3) |  |
| `B` | 1D array of complex (shape: 3) |  |

::: pybmad.GridFieldPtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GridFieldPtStruct

Fortran struct: `grid_field_pt_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L807))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `file` | str | Input file name. Used also as ID for instances. |
| `n_link` | int | For memory management of this structure |
| `pt` | [3D array of GridFieldPt1Struct](bmad.md#gridfieldpt1struct) |  |

::: pybmad.GridFieldStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GridFieldStruct

Fortran struct: `grid_field_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L813))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `geometry` | int | Type of grid: xyz$, or rotationally_symmetric_rz$ |
| `harmonic` | int | Harmonic of fundamental for AC fields. |
| `phi0_fieldmap` | float | Mode oscillates as: twopi * (f * t + phi0_fieldmap) |
| `field_scale` | float | Factor to scale the fields by |
| `field_type` | int | or magnetic$ or electric$ |
| `master_parameter` | int | Master parameter in ele%value(:) array to use for scaling the field. |
| `ele_anchor_pt` | int | anchor_beginning$, anchor_center$, or anchor_end$ |
| `interpolation_order` | int | Possibilities are 1 or 3. |
| `dr` | 1D array of float (shape: 3) | Grid spacing. |
| `r0` | 1D array of float (shape: 3) | Field origin relative to ele_anchor_pt. |
| `curved_ref_frame` | bool |  |
| `ptr` | [GridFieldPtStruct](bmad.md#gridfieldptstruct) |  |
| `bi_coef` | [3D array of BicubicCmplxCoefStruct (shape: 4, 2, 3)](sim_utils.md#bicubiccmplxcoefstruct) | Save computed coefs for faster tracking |
| `tri_coef` | [3D array of TricubicCmplxCoefStruct (shape: 4, 2, 3)](sim_utils.md#tricubiccmplxcoefstruct) | Save computed coefs for faster tracking |

::: pybmad.HighEnergySpaceChargeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### HighEnergySpaceChargeStruct

Fortran struct: `high_energy_space_charge_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L863))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `closed_orb` | [CoordStruct](bmad.md#coordstruct) | beam orbit |
| `kick_const` | float |  |
| `sig_x` | float |  |
| `sig_y` | float |  |
| `phi` | float | Rotation angle to go from lab frame to rotated frame. |
| `sin_phi` | float |  |
| `cos_phi` | float |  |
| `sig_z` | float |  |

::: pybmad.IbsLifetimeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### IbsLifetimeStruct

Fortran struct: `ibs_lifetime_struct` ([`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L23))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `Tlx` | float |  |
| `Tly` | float |  |
| `Tlp` | float |  |

::: pybmad.IbsMaxratioStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### IbsMaxratioStruct

Fortran struct: `ibs_maxratio_struct` ([`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L30))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `rx` | float |  |
| `ry` | float |  |
| `r_p` | float |  |

::: pybmad.IbsSimParamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### IbsSimParamStruct

Fortran struct: `ibs_sim_param_struct` ([`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L10))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `tau_a` | float | horizontal damping rate (needed for coulomb log tail cut) |
| `clog_to_use` | int | see multi_coulomb_log subroutine for valid settings.  Set to 1 to disable tail-cut.  Set to 1 for linacs. |
| `set_dispersion` | bool | True: add vertical dispersion to transfer matrix.  Valid for kubo method. |
| `eta_set` | float | If set_dispersion, then this value is used to add y-z coupling to the transfer matrix. |
| `etap_set` | float | If set_dispersion, then this value is used to add y-z coupling to the transfer matrix. |
| `do_pwd` | bool | If true, then use potential well distortion to calculate bunch lengths.  If false, bunch length is proportional to energy spread. |
| `inductance` | float | Inductive part of impedance for pwd calc. |
| `formula` | str | Which IBS formulation to use.  See subroutine ibs1 for a list. real(rp) :: fake_3HC = -1   ! If greater than zero, divide growth rates by this factor. |

::: pybmad.IbsStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### IbsStruct

Fortran struct: `ibs_struct` ([`bmad/multiparticle/ibs_rates_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_rates_mod.f90#L7))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `inv_Ta` | float |  |
| `inv_Tb` | float |  |
| `inv_Tz` | float |  |

::: pybmad.Interval1CoefStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Interval1CoefStruct

Fortran struct: `interval1_coef_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L217))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `c0` | float |  |
| `c1` | float |  |
| `n_exp` | float |  |

::: pybmad.KvBeamInitStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### KvBeamInitStruct

Fortran struct: `kv_beam_init_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1147))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `part_per_phi` | 1D array of int (shape: 2) | number of particles per angle variable. |
| `n_I2` | int | number of I2 |
| `A` | float | A = I1/e |

::: pybmad.LatEleLocStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LatEleLocStruct

Fortran struct: `lat_ele_loc_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L881))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ix_ele` | int |  |
| `ix_branch` | int |  |

::: pybmad.LatEleOrder1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LatEleOrder1Struct

Fortran struct: `lat_ele_order1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L888))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ix_branch` | int | Branch index |
| `ix_order` | int | Order index. -1 -> Unique in lattice, 0 -> unique in branch. |

::: pybmad.LatEleOrderArrayStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LatEleOrderArrayStruct

Fortran struct: `lat_ele_order_array_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L893))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele` | [1D array of LatEleOrder1Struct](bmad.md#lateleorder1struct) |  |

::: pybmad.LatEleOrderStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LatEleOrderStruct

Fortran struct: `lat_ele_order_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L912))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `branch` | [1D array of LatEleOrderArrayStruct](bmad.md#lateleorderarraystruct) |  |

::: pybmad.LatParamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LatParamStruct

Fortran struct: `lat_param_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1532))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_part` | float | Particles/bunch (for BeamBeam elements). |
| `total_length` | float | total_length of branch. Warning: branch may not start at s = 0. |
| `unstable_factor` | float | If positive: Growth rate/turn if unstable in closed branches or \|orbit-aperture\|/aperture if particle hits wall. Zero otherwise. |
| `t1_with_RF` | 2D array of float (shape: 6,6) | Full 1-turn matrix with RF on. |
| `t1_no_RF` | 2D array of float (shape: 6,6) | Full 1-turn matrix with RF off. |
| `spin_tune` | float | Closed orbit spin tune. |
| `particle` | int | Reference particle: positron$, electron$, etc. Call lattice_bookkeeper if this is changed. |
| `default_tracking_species` | int | Default particle type to use in tracking. |
| `geometry` | int | open$ or closed$ |
| `ixx` | int | Integer for general use |
| `stable` | bool | is closed lat stable? |
| `live_branch` | bool | Should tracking be done on the branch? |
| `g1_integral` | float | Approximate \|g\| (bending strength) integral of branch. |
| `g2_integral` | float | Approximate g^2 integral of branch. |
| `g3_integral` | float | Approximate g^2 integral of branch. |
| `bookkeeping_state` | [BookkeepingStateStruct](bmad.md#bookkeepingstatestruct) | Overall status for the branch. |
| `beam_init` | [BeamInitStruct](bmad.md#beaminitstruct) | For beam initialization. |

::: pybmad.LatPointerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LatPointerStruct

Fortran struct: `lat_pointer_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L937))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `lat` | [LatStruct](bmad.md#latstruct) |  |

::: pybmad.LatStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LatStruct

Fortran struct: `lat_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1656))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `use_name` | str | Name of lat given by USE statement |
| `lattice` | str | Lattice |
| `machine` | str | Name of the machine the lattice is for ("LHC", etc). |
| `input_file_name` | str | Name of the lattice input file |
| `title` | str | General title |
| `print_str` | 1D array of str | Saved print statements. |
| `constant` | [1D array of ExpressionAtomStruct](bmad.md#expressionatomstruct) | Constants defined in the lattice |
| `a` | [ModeInfoStruct](bmad.md#modeinfostruct) | Tunes (fractional part), etc. |
| `b` | [ModeInfoStruct](bmad.md#modeinfostruct) | Tunes (fractional part), etc. |
| `z` | [ModeInfoStruct](bmad.md#modeinfostruct) | Tunes (fractional part), etc. |
| `param` | [LatParamStruct](bmad.md#latparamstruct) | Parameters |
| `lord_state` | [BookkeepingStateStruct](bmad.md#bookkeepingstatestruct) | lord bookkeeping status. |
| `ele_init` | [EleStruct](bmad.md#elestruct) | For use by any program |
| `ele` | [1D array of EleStruct](bmad.md#elestruct) | Array of elements [=> branch(0)]. |
| `branch` | [1D array of BranchStruct](bmad.md#branchstruct) | Branch(0:) array |
| `control` | [1D array of ControlStruct](bmad.md#controlstruct) | Control list |
| `particle_start` | [CoordStruct](bmad.md#coordstruct) | Starting particle_coords. |
| `beam_init` | [BeamInitStruct](bmad.md#beaminitstruct) | Beam initialization. |
| `pre_tracker` | [PreTrackerStruct](bmad.md#pretrackerstruct) | For OPAL/IMPACT-T |
| `nametable` | [NametableStruct](sim_utils.md#nametablestruct) | For quick searching by element name. |
| `custom` | 1D array of float | Custom attributes. |
| `version` | int | Version number |
| `n_ele_track` | int | Number of lat elements to track through. |
| `n_ele_max` | int | Index of last valid element in %ele(:) array |
| `n_control_max` | int | Last index used in control_array |
| `n_ic_max` | int | Last index used in ic_array |
| `input_taylor_order` | int | As set in the input file |
| `ic` | 1D array of int | Index to %control(:) from slaves. |
| `photon_type` | int | Or coherent$. For X-ray simulations. |
| `creation_hash` | int | Set by bmad_parser. creation_hash will vary if any of the lattice files are modified. |
| `ramper_slave_bookkeeping` | int |  |
| `parser_make_xfer_mats` | bool | Is Bmad parser to make element transfer matrices? |

::: pybmad.LinacNormalModeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LinacNormalModeStruct

Fortran struct: `linac_normal_mode_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2011))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `i2_E4` | float | Integral: g^2 * gamma^4 |
| `i3_E7` | float | Integral: g^3 * gamma^7 |
| `i5a_E6` | float | Integral: (g^3 * H_a) * gamma^6 |
| `i5b_E6` | float | Integral: (g^3 * H_b) * gamma^6 |
| `sig_E1` | float | Energy spread after 1 pass (eV) |
| `a_emittance_end` | float | a mode emittance at end of linac |
| `b_emittance_end` | float | b mode emittance at end of linac |

::: pybmad.LinearEleIsfStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LinearEleIsfStruct

Fortran struct: `linear_ele_isf_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L339))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `node` | [1D array of LinearIsf1Struct](bmad.md#linearisf1struct) | Array per PTC integration node. |

::: pybmad.LinearIsf1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LinearIsf1Struct

Fortran struct: `linear_isf1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L330))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `orb0` | 1D array of float (shape: 6) | Closed orbit. |
| `isf` | 2D array of float (shape: 0:3, 0:6) | Linear ISF map at a given point. |
| `s` | float | Offset from beginning of element. !! real(rp) :: m_1turn(6,6) = 0   ! Orbital 1-turn matrix. |

::: pybmad.MadEnergyStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MadEnergyStruct

Fortran struct: `mad_energy_struct` ([`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L14))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `total` | float |  |
| `beta` | float | normalized velocity: v/c |
| `gamma` | float | relativistic factor: 1/sqrt(1-beta^2) |
| `kinetic` | float | kinetic energy |
| `p0c` | float | particle momentum |
| `particle` | int | particle species |

::: pybmad.MadMapStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MadMapStruct

Fortran struct: `mad_map_struct` ([`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L23))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `k` | 1D array of float (shape: 6) | 0th order map. |
| `r` | 2D array of float (shape: 6,6) | 1st order map. |
| `t` | 3D array of float (shape: 6,6,6) | 2nd order map. |

::: pybmad.MaterialStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MaterialStruct

Fortran struct: `material_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1301))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `species` | int |  |
| `number` | int | Relative number |
| `density` | float |  |
| `density_used` | float |  |
| `area_density` | float |  |
| `area_density_used` | float |  |
| `radiation_length` | float |  |
| `radiation_length_used` | float |  |

::: pybmad.Mesh3DStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Mesh3DStruct

Fortran struct: `mesh3d_struct` ([`bmad/space_charge/open_spacecharge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_mod.f90#L14))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `nlo` | 1D array of int (shape: 3) | Lowest  grid index in x, y, z (m) of rho and the quantity being computed (phi or E) |
| `nhi` | 1D array of int (shape: 3) | Highest grid index in x, y, z (m) of rho and the quantity being computed (phi or E) |
| `npad` | 1D array of int (shape: 3) | Array padding for cyclic convolution |
| `min` | 1D array of float (shape: 3) | Minimim in each dimension |
| `max` | 1D array of float (shape: 3) | Maximum in each dimension |
| `delta` | 1D array of float (shape: 3) | Grid spacing |
| `gamma` | float | Relativistic gamma |
| `charge` | float | Total charge on mesh |
| `rho` | 3D array of float | Charge density grid |
| `phi` | 3D array of float | electric potential grid |
| `efield` | 4D array of float | electric field grid |
| `bfield` | 4D array of float | magnetic field grid |

::: pybmad.Mode3Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Mode3Struct

Fortran struct: `mode3_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L943))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `v` | 2D array of float (shape: 6,6) |  |
| `a` | [TwissStruct](bmad.md#twissstruct) |  |
| `b` | [TwissStruct](bmad.md#twissstruct) |  |
| `c` | [TwissStruct](bmad.md#twissstruct) |  |
| `x` | [TwissStruct](bmad.md#twissstruct) |  |
| `y` | [TwissStruct](bmad.md#twissstruct) |  |

::: pybmad.ModeInfoStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ModeInfoStruct

Fortran struct: `mode_info_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1568))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `stable` | bool | Is the mode stable? |
| `tune` | float | "fractional" tune in radians |
| `emit` | float | Emittance (unnormalized). |
| `chrom` | float | Chromaticity. |
| `sigma` | float | Beam size. |
| `sigmap` | float | Beam divergence. |

::: pybmad.MomentumApertureStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MomentumApertureStruct

Fortran struct: `momentum_aperture_struct` ([`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/touschek_mod.f90#L23))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `s` | float |  |
| `pos` | float |  |
| `neg` | float |  |

::: pybmad.MultipassAllInfoStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MultipassAllInfoStruct

Fortran struct: `multipass_all_info_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2135))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `lord` | [1D array of MultipassLordInfoStruct](bmad.md#multipasslordinfostruct) | Array of lords |
| `branch` | [1D array of MultipassBranchInfoStruct](bmad.md#multipassbranchinfostruct) |  |

::: pybmad.MultipassBranchInfoStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MultipassBranchInfoStruct

Fortran struct: `multipass_branch_info_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2128))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele` | [1D array of MultipassEleInfoStruct](bmad.md#multipasseleinfostruct) |  |

::: pybmad.MultipassEleInfoStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MultipassEleInfoStruct

Fortran struct: `multipass_ele_info_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2121))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `multipass` | bool | True if involved in multipass. False otherwise |
| `ix_pass` | int | Pass number |
| `ix_lord` | 1D array of int | Pointers to lord(:) array |
| `ix_super` | 1D array of int | Indexes to slave(ix_pass, super_slave%ix_ele) matrix |

::: pybmad.MultipassLordInfoStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MultipassLordInfoStruct

Fortran struct: `multipass_lord_info_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2110))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `lord` | [EleStruct](bmad.md#elestruct) | Lord element |
| `n_pass` | int | Number of passes (= number of slaves) |
| `n_super_slave` | int | Number of super_slaves per super_lord. |
| `super_lord` | [1D array of ElePointerStruct](bmad.md#elepointerstruct) | Super_lord list if they exist. |
| `slave` | [2D array of ElePointerStruct](bmad.md#elepointerstruct) | Slaves list in tracking part. |

::: pybmad.MultipassRegionBranchStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MultipassRegionBranchStruct

Fortran struct: `multipass_region_branch_struct` ([`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L12))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele` | [1D array of MultipassRegionEleStruct](bmad.md#multipassregionelestruct) |  |

::: pybmad.MultipassRegionEleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MultipassRegionEleStruct

Fortran struct: `multipass_region_ele_struct` ([`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L6))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ix_region` | int |  |
| `region_start_pt` | bool |  |
| `region_stop_pt` | bool |  |

::: pybmad.MultipassRegionLatStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MultipassRegionLatStruct

Fortran struct: `multipass_region_lat_struct` ([`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L16))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `branch` | [1D array of MultipassRegionBranchStruct](bmad.md#multipassregionbranchstruct) |  |

::: pybmad.MultipoleCacheStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MultipoleCacheStruct

Fortran struct: `multipole_cache_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L973))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `a_pole_mag` | 1D array of float |  |
| `b_pole_mag` | 1D array of float |  |
| `a_kick_mag` | 1D array of float |  |
| `b_kick_mag` | 1D array of float |  |
| `ix_pole_mag_max` | int |  |
| `ix_kick_mag_max` | int |  |
| `mag_valid` | bool | From elseparator hkick and vkick. |
| `a_pole_elec` | 1D array of float |  |
| `b_pole_elec` | 1D array of float |  |
| `a_kick_elec` | 1D array of float |  |
| `b_kick_elec` | 1D array of float |  |
| `ix_pole_elec_max` | int |  |
| `ix_kick_elec_max` | int |  |
| `elec_valid` | bool |  |

::: pybmad.NormalModesStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### NormalModesStruct

Fortran struct: `normal_modes_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2021))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `synch_int` | 1D array of float (shape: 0:3) | Synchrotron integrals I0, I1, I2, and I3 |
| `sigE_E` | float | SigmaE/E |
| `sig_z` | float | Sigma_Z |
| `e_loss` | float | Energy loss / turn (eV) |
| `rf_voltage` | float | Total rfcavity voltage (eV) |
| `pz_aperture` | float | pz aperture limit. Used with Touschek calculations. |
| `pz_average` | float | Average over branch due to damping. |
| `momentum_compaction` | float |  |
| `dpz_damp` | float | Change in pz without RF |
| `a` | [AnormalModeStruct](bmad.md#anormalmodestruct) |  |
| `b` | [AnormalModeStruct](bmad.md#anormalmodestruct) |  |
| `z` | [AnormalModeStruct](bmad.md#anormalmodestruct) |  |
| `lin` | [LinacNormalModeStruct](bmad.md#linacnormalmodestruct) |  |

::: pybmad.ParserControllerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ParserControllerStruct

Fortran struct: `parser_controller_struct` ([`bmad/parsing/bmad_parser_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_struct.f90#L78))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `attrib_name` | str |  |
| `stack` | [1D array of ExpressionAtomStruct](bmad.md#expressionatomstruct) | Arithmetic expression stack |
| `y_knot` | 1D array of float |  |
| `n_stk` | int |  |

::: pybmad.ParserEleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ParserEleStruct

Fortran struct: `parser_ele_struct` ([`bmad/parsing/bmad_parser_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_struct.f90#L86))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `control` | [1D array of ParserControllerStruct](bmad.md#parsercontrollerstruct) |  |
| `field_overlaps` | 1D array of str |  |
| `ref_name` | str |  |
| `ix_super_ref_multipass` | int | Multipass index for superimpose reference element. |
| `ele_name` | str | For fork element or superimpose statement. |
| `names1` | 1D array of str | Currently just used by feedback element. |
| `names2` | 1D array of str | Currently just used by feedback element. |
| `lat_file` | str | File where element was defined. |
| `offset` | float |  |
| `ix_line_in_file` | int | Line in file where element was defined. |
| `ix_count` | int |  |
| `ele_pt` | int |  |
| `ref_pt` | int |  |
| `index` | int |  |
| `superposition_command_here` | bool |  |
| `superposition_has_been_set` | bool |  |
| `wrap_superimpose` | bool |  |
| `create_jumbo_slave` | bool |  |
| `is_range` | bool | For girders |
| `default_attrib` | str | For group/overlay elements: slave attribute |

::: pybmad.ParserLatStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ParserLatStruct

Fortran struct: `parser_lat_struct` ([`bmad/parsing/bmad_parser_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_struct.f90#L108))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele` | [1D array of ParserEleStruct](bmad.md#parserelestruct) |  |

::: pybmad.PhotonCoordStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonCoordStruct

Fortran struct: `photon_coord_struct` ([`bmad/photon/capillary_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/capillary_mod.f90#L5))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `orb` | [CoordStruct](bmad.md#coordstruct) | Phase space: orb%vec = (x, vx/c, y, vy/c, s, vs/c) |
| `track_len` | float | Total track length from the start of the element. |
| `ix_section` | int | Cross section index |

::: pybmad.PhotonElementStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonElementStruct

Fortran struct: `photon_element_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1099))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `curvature` | [SurfaceCurvatureStruct](bmad.md#surfacecurvaturestruct) |  |
| `target` | [PhotonTargetStruct](bmad.md#photontargetstruct) |  |
| `material` | [PhotonMaterialStruct](bmad.md#photonmaterialstruct) |  |
| `segmented` | [SurfaceSegmentedStruct](bmad.md#surfacesegmentedstruct) |  |
| `h_misalign` | [SurfaceHMisalignStruct](bmad.md#surfacehmisalignstruct) |  |
| `displacement` | [SurfaceDisplacementStruct](bmad.md#surfacedisplacementstruct) |  |
| `pixel` | [PixelDetecStruct](bmad.md#pixeldetecstruct) |  |
| `reflectivity_table_type` | int |  |
| `reflectivity_table_sigma` | [PhotonReflectTableStruct](bmad.md#photonreflecttablestruct) | If polarization is ignored use sigma table. |
| `reflectivity_table_pi` | [PhotonReflectTableStruct](bmad.md#photonreflecttablestruct) |  |
| `init_energy_prob` | [1D array of SplineStruct](sim_utils.md#splinestruct) | Initial energy probability density |
| `integrated_init_energy_prob` | 1D array of float |  |

::: pybmad.PhotonInitSplinesStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonInitSplinesStruct

Fortran struct: `photon_init_splines_struct` ([`bmad/photon/photon_init_spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_spline_mod.f90#L24))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `source_type` | str | 'bend', 'wiggler', 'undulator' |
| `spline_space_dimensions` | int | Dimensions: [energy, y_angle, x_angle, x, y] |
| `energy_prob` | [1D array of SplineStruct](sim_utils.md#splinestruct) |  |
| `y_angle` | [1D array of PhotonInitYAngleSplineStruct](bmad.md#photoninityanglesplinestruct) |  |

::: pybmad.PhotonInitXAngleSplineStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonInitXAngleSplineStruct

Fortran struct: `photon_init_x_angle_spline_struct` ([`bmad/photon/photon_init_spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_spline_mod.f90#L15))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `prob` | [1D array of SplineStruct](sim_utils.md#splinestruct) |  |
| `pl` | [1D array of SplineStruct](sim_utils.md#splinestruct) |  |
| `pc` | [1D array of SplineStruct](sim_utils.md#splinestruct) |  |
| `pl45` | [1D array of SplineStruct](sim_utils.md#splinestruct) |  |

::: pybmad.PhotonInitYAngleSplineStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonInitYAngleSplineStruct

Fortran struct: `photon_init_y_angle_spline_struct` ([`bmad/photon/photon_init_spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_spline_mod.f90#L19))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `prob` | [1D array of SplineStruct](sim_utils.md#splinestruct) |  |
| `pl` | [1D array of SplineStruct](sim_utils.md#splinestruct) |  |
| `pc` | [1D array of SplineStruct](sim_utils.md#splinestruct) |  |
| `pl45` | [1D array of SplineStruct](sim_utils.md#splinestruct) |  |
| `x_angle` | [1D array of PhotonInitXAngleSplineStruct](bmad.md#photoninitxanglesplinestruct) |  |

::: pybmad.PhotonMaterialStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonMaterialStruct

Fortran struct: `photon_material_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1084))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `f0_m1` | complex | For multilayer_mirror only. |
| `f0_m2` | complex | For multilayer_mirror only. |
| `f_0` | complex |  |
| `f_h` | complex | Structure factor for H direction. |
| `f_hbar` | complex | Structure factor for -H direction. |
| `f_hkl` | complex | = sqrt(f_h * f_hbar) |
| `h_norm` | 1D array of float (shape: 3) | Normalized H vector for crystals. |
| `l_ref` | 1D array of float (shape: 3) | Crystal reference orbit displacement vector in element coords. |

::: pybmad.PhotonReflectSurfaceStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonReflectSurfaceStruct

Fortran struct: `photon_reflect_surface_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L235))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `description` | str | Descriptive name |
| `reflectivity_file` | str |  |
| `table` | [1D array of PhotonReflectTableStruct](bmad.md#photonreflecttablestruct) |  |
| `surface_roughness_rms` | float | sigma in Dugan's notation |
| `roughness_correlation_len` | float | T in Dugan's notation |
| `ix_surface` | int |  |

::: pybmad.PhotonReflectTableStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonReflectTableStruct

Fortran struct: `photon_reflect_table_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L221))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `angle` | 1D array of float | Vector of angle values for %p_reflect |
| `energy` | 1D array of float | Vector of energy values for %p_reflect |
| `int1` | [1D array of Interval1CoefStruct](bmad.md#interval1coefstruct) |  |
| `p_reflect` | 2D array of float | (angle, ev) probability. Log used for smooth surface reflection |
| `max_energy` | float | maximum energy for this table |
| `p_reflect_scratch` | 1D array of float | Scratch space |
| `bragg_angle` | 1D array of float | Bragg angle at energy values. |

::: pybmad.PhotonTargetStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonTargetStruct

Fortran struct: `photon_target_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1076))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `type` | int | or rectangular$ |
| `n_corner` | int |  |
| `ele_loc` | [LatEleLocStruct](bmad.md#latelelocstruct) |  |
| `corner` | [1D array of TargetPointStruct (shape: 8)](bmad.md#targetpointstruct) |  |
| `center` | [TargetPointStruct](bmad.md#targetpointstruct) |  |

::: pybmad.PhotonTrackStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonTrackStruct

Fortran struct: `photon_track_struct` ([`bmad/photon/capillary_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/capillary_mod.f90#L11))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `old` | [PhotonCoordStruct](bmad.md#photoncoordstruct) |  |
| `now` | [PhotonCoordStruct](bmad.md#photoncoordstruct) |  |

::: pybmad.PixelDetecStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PixelDetecStruct

Fortran struct: `pixel_detec_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1053))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `dr` | 1D array of float (shape: 2) |  |
| `r0` | 1D array of float (shape: 2) |  |
| `n_track_tot` | int | How many photons were launched from source element. |
| `n_hit_detec` | int | How many photons hit the detector. |
| `n_hit_pixel` | int | How many photons hit the pixel grid of the detector. |
| `pt` | [2D array of PixelPtStruct](bmad.md#pixelptstruct) | Grid of pixels |

::: pybmad.PixelPtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PixelPtStruct

Fortran struct: `pixel_pt_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1043))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_photon` | int |  |
| `E_x` | complex |  |
| `E_y` | complex |  |
| `intensity_x` | float |  |
| `intensity_y` | float |  |
| `intensity` | float |  |
| `orbit` | 1D array of float (shape: 6) | x, Vx/c, y, Vy/c, dummy, E - E_ref. |
| `orbit_rms` | 1D array of float (shape: 6) | RMS statistics. |
| `init_orbit` | 1D array of float (shape: 6) | Initial orbit at start of lattice statistics. |
| `init_orbit_rms` | 1D array of float (shape: 6) | Initial orbit at start of lattice RMS statistics. |

::: pybmad.PmdHeaderStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PmdHeaderStruct

Fortran struct: `pmd_header_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2453))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `openPMD` | str |  |
| `openPMDextension` | str |  |
| `basePath` | str |  |
| `particlesPath` | str |  |
| `meshesPath` | str |  |
| `author` | str |  |
| `software` | str |  |
| `softwareVersion` | str |  |
| `date` | str |  |
| `latticeFile` | str |  |
| `latticeName` | str |  |

::: pybmad.PreTrackerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PreTrackerStruct

Fortran struct: `pre_tracker_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1640))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `who` | int | Can be opal$, or impactt$ |
| `ix_ele_start` | int |  |
| `ix_ele_end` | int |  |
| `input_file` | str |  |

::: pybmad.PtcBranch1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PtcBranch1Struct

Fortran struct: `ptc_branch1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1561))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `m_u_layout` | [1D array of PtcLayoutPointerStruct](bmad.md#ptclayoutpointerstruct) |  |

::: pybmad.PtcLayoutPointerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PtcLayoutPointerStruct

Fortran struct: `ptc_layout_pointer_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1557))

All attributes may be passed to the initializer as arguments:

::: pybmad.PtcNormalFormStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PtcNormalFormStruct

Fortran struct: `ptc_normal_form_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1600))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele_origin` | [EleStruct](bmad.md#elestruct) | Element at which the on-turn map was created. |
| `orb0` | 1D array of float (shape: 6) | Closed orbit at element. |
| `valid_map` | bool |  |

::: pybmad.PtcRadMapStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PtcRadMapStruct

Fortran struct: `ptc_rad_map_struct` ([`bmad/ptc/ptc_map_with_radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_map_with_radiation_mod.f90#L9))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `lattice_file` | str | Name of the lattice file |
| `dref_time` | float | Time ref particle takes. |
| `p0c_start` | float | ref momentum at start |
| `p0c_end` | float | ref momentum at end |
| `s_end` | float | Ending s-position |
| `map_order` | int |  |
| `radiation_damping_on` | bool |  |
| `ix_branch` | int |  |
| `ix_ele_start` | int | Start point for making the map |
| `ix_ele_end` | int | End point for making the map |
| `nodamp_mat` | 2D array of float (shape: 6,6) | Nondamped orbital matrix. M_orbit = M_damp * M_nodamp |
| `damp_mat` | 2D array of float (shape: 6,6) | Damping "correction" to orbital matrix. Stoc_mat is referenced to the start of the map. That is, it is applied before the transport matrix. |
| `stoc_mat` | 2D array of float (shape: 6,6) | Stochatic matrix for the orbit. |
| `ref0` | 1D array of float (shape: 6) | Reference orbit at start. |
| `ref1` | 1D array of float (shape: 6) | Reference orbit at end. |

::: pybmad.RadInt1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RadInt1Struct

Fortran struct: `rad_int1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2420))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `i0` | float |  |
| `i1` | float |  |
| `i2` | float |  |
| `i3` | float |  |
| `i4a` | float |  |
| `i4b` | float |  |
| `i4z` | float |  |
| `i5a` | float |  |
| `i5b` | float |  |
| `i6b` | float |  |
| `lin_i2_E4` | float |  |
| `lin_i3_E7` | float |  |
| `lin_i5a_E6` | float |  |
| `lin_i5b_E6` | float |  |
| `lin_norm_emit_a` | float | Running sum |
| `lin_norm_emit_b` | float | Running sum |
| `lin_sig_E` | float | Running sum |
| `n_steps` | float | number of qromb steps needed |

::: pybmad.RadIntAllEleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RadIntAllEleStruct

Fortran struct: `rad_int_all_ele_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2447))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `branch` | [1D array of RadIntBranchStruct](bmad.md#radintbranchstruct) | Array is indexed from 0 |

::: pybmad.RadIntBranchStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RadIntBranchStruct

Fortran struct: `rad_int_branch_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2443))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele` | [1D array of RadInt1Struct](bmad.md#radint1struct) | Array is indexed from 0 |

::: pybmad.RadIntCache1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RadIntCache1Struct

Fortran struct: `rad_int_cache1_struct` ([`bmad/modules/rad_int_common.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/rad_int_common.f90#L28))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `pt` | [1D array of RadIntTrackPointStruct](bmad.md#radinttrackpointstruct) | pt(0:n_pt) |
| `n_pt` | int | Upper bound of pt(0:n_pt) |
| `cache_type` | int |  |

::: pybmad.RadIntInfoStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RadIntInfoStruct

Fortran struct: `rad_int_info_struct` ([`bmad/modules/rad_int_common.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/rad_int_common.f90#L43))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `branch` | [BranchStruct](bmad.md#branchstruct) |  |
| `ele` | [EleStruct](bmad.md#elestruct) |  |
| `orbit` | [1D array of CoordStruct](bmad.md#coordstruct) |  |
| `a` | [TwissStruct](bmad.md#twissstruct) |  |
| `b` | [TwissStruct](bmad.md#twissstruct) |  |
| `cache_ele` | [RadIntCache1Struct](bmad.md#radintcache1struct) | pointer to cache in use |
| `eta_a` | 1D array of float (shape: 4) |  |
| `eta_b` | 1D array of float (shape: 4) |  |
| `g` | float | bending strength (1/bending_radius) |
| `g2` | float | bending strength (1/bending_radius) |
| `g_x` | float | components in x-y plane |
| `g_y` | float | components in x-y plane |
| `dg2_x` | float |  |
| `dg2_y` | float |  |

::: pybmad.RadIntTrackPointStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RadIntTrackPointStruct

Fortran struct: `rad_int_track_point_struct` ([`bmad/modules/rad_int_common.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/rad_int_common.f90#L15))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `s_body` | float |  |
| `mat6` | 2D array of float (shape: 6,6) |  |
| `vec0` | 1D array of float (shape: 6) |  |
| `ref_orb_in` | [CoordStruct](bmad.md#coordstruct) |  |
| `ref_orb_out` | [CoordStruct](bmad.md#coordstruct) |  |
| `g_x0` | float | Additional g factors for bends. |
| `g_y0` | float | Additional g factors for bends. |
| `dgx_dx` | float | bending strength gradient |
| `dgx_dy` | float | bending strength gradient |
| `dgy_dx` | float | bending strength gradient |
| `dgy_dy` | float | bending strength gradient |

::: pybmad.RadMapEleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RadMapEleStruct

Fortran struct: `rad_map_ele_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L996))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `rm0` | [RadMapStruct](bmad.md#radmapstruct) | Upstream half and downstream half matrices for an element. |
| `rm1` | [RadMapStruct](bmad.md#radmapstruct) | Upstream half and downstream half matrices for an element. |
| `stale` | bool |  |

::: pybmad.RadMapStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RadMapStruct

Fortran struct: `rad_map_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L988))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ref_orb` | 1D array of float (shape: 6) | Reference point around which damp_mat is calculated. |
| `damp_dmat` | 2D array of float (shape: 6,6) | damp_correction = xfer_mat_with_damping - xfer_mat_without_damping. |
| `xfer_damp_vec` | 1D array of float (shape: 6) | Transfer map with damping 0th order vector. |
| `xfer_damp_mat` | 2D array of float (shape: 6,6) | 1st order matrix: xfer_no_damp_mat + xfer_damp_correction. |
| `stoc_mat` | 2D array of float (shape: 6,6) | Stochastic variance or "kick" (Cholesky decomposed) matrix. |

::: pybmad.RamperLordStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RamperLordStruct

Fortran struct: `ramper_lord_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1391))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ix_ele` | int | Lord index |
| `ix_con` | int | Index in lord%control%ramp(:) array |
| `attrib_ptr` | float | Pointer to attribute in this element. |

::: pybmad.ResonanceHStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ResonanceHStruct

Fortran struct: `resonance_h_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1579))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `id` | str | 6 digit ID. EG: '003100' |
| `c_val` | complex | Resonance value |

::: pybmad.RfEleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RfEleStruct

Fortran struct: `rf_ele_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1337))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `steps` | [1D array of RfStairStepStruct](bmad.md#rfstairstepstruct) | Energy stair step array indexed from zero. |
| `ds_step` | float | length of a stair step. |

::: pybmad.RfStairStepStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RfStairStepStruct

Fortran struct: `rf_stair_step_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1316))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `E_tot0` | float | Reference energy in the drift region (before the kick point). |
| `E_tot1` | float | Reference energy after the kick point. |
| `p0c` | float | Reference momentum in the drift region (before the kick point). |
| `p1c` | float | Reference momentum after the kick point. |
| `scale` | float | Scale for multipole kick at the kick point. Sum over all steps will be 1. |
| `time` | float | Reference particle time at the kick point with respect to beginning of element. |
| `s0` | float | S-position at beginning of drift region relative to the beginning of the element. |
| `s` | float | S-position at the kick point relative to the beginning of the element. |
| `ix_step` | int | Step index in ele%rf%steps(:) array |

::: pybmad.SeqEleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SeqEleStruct

Fortran struct: `seq_ele_struct` ([`bmad/parsing/bmad_parser_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_struct.f90#L14))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | name of element, subline, or sublist |
| `actual_arg` | 1D array of str |  |
| `tag` | str | tag name. |
| `slice_start` | str | For "my_line[start:end]" slice constructs. |
| `slice_end` | str | For "my_line[start:end]" slice constructs. |
| `type` | int | LINE$, REPLACEMENT_LINE$, LIST$, ELEMENT$ |
| `ix_ele` | int | if an element: pointer to ELE array if a line or list: pointer to SEQ array |
| `ix_arg` | int | index in arg list (for replacement lines) |
| `rep_count` | int | how many copies of an element |
| `ele_order_reflect` | bool | Travel through ele sequence in reverse order |
| `ele_orientation` | int | element has reverse orientation. |

::: pybmad.SeqStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SeqStruct

Fortran struct: `seq_struct` ([`bmad/parsing/bmad_parser_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_struct.f90#L38))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | name of sequence |
| `ele` | [1D array of SeqEleStruct](bmad.md#seqelestruct) | Elements in the sequence |
| `dummy_arg` | 1D array of str |  |
| `corresponding_actual_arg` | 1D array of str |  |
| `type` | int | LINE$, REPLACEMENT_LINE$ or LIST$ |
| `ix_list` | int | Current index for lists |
| `list_upcount` | int |  |
| `index` | int | Alphabetical order sorted index |
| `file_name` | str | File where sequence is defined |
| `ix_file_line` | int | Line number in file where sequence is defined |
| `multipass` | bool |  |
| `ptc_layout` | bool | Put in separate PTC layout |
| `active` | bool | Used to prevent infinite loops. |

::: pybmad.SpaceChargeCommonStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SpaceChargeCommonStruct

Fortran struct: `space_charge_common_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2176))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ds_track_step` | float | CSR tracking step size |
| `dt_track_step` | float | Time Runge kutta initial step. |
| `cathode_strength_cutoff` | float | Cutoff for the cathode field calc. |
| `rel_tol_tracking` | float | Relative tolerance for tracking. |
| `abs_tol_tracking` | float | Absolute tolerance for tracking. |
| `beam_chamber_height` | float | Used in shielding calculation. |
| `lsc_sigma_cutoff` | float | Cutoff for the 1-dim longitudinal SC calc. If a bin sigma is < cutoff * sigma_ave then ignore. |
| `particle_sigma_cutoff` | float | 3D SC calc cutoff for particles with (x,y,z) position far from the center. Negative or zero means ignore. |
| `mesh_growth_factor` | float | Fractional padding when growing SC mesh (default: 10%). Set to 0 for tight-fit (no caching speedup). |
| `mesh_shrink_factor` | float | Fractional threshold for shrinking SC mesh (default: 10%). Mesh shrinks when bunch fills < (1-this) of the mesh range. |
| `space_charge_mesh_size` | 1D array of int (shape: 3) | Gird size for fft_3d space charge calc. |
| `csr3d_mesh_size` | 1D array of int (shape: 3) | Gird size for CSR. |
| `n_bin` | int | Number of bins used |
| `particle_bin_span` | int | Longitudinal particle length / dz_bin |
| `n_shield_images` | int | Chamber wall shielding. 0 = no shielding. |
| `sc_min_in_bin` | int | Minimum number of particles in a bin for sigmas to be valid. |
| `lsc_kick_transverse_dependence` | bool |  |
| `debug` | bool |  |
| `diagnostic_output_file` | str | If non-blank write a diagnostic (EG wake) file |

::: pybmad.SpinAxisStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SpinAxisStruct

Fortran struct: `spin_axis_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L291))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `l` | 1D array of float (shape: 3) | Transverse axis. |
| `n0` | 1D array of float (shape: 3) | Invariant spin axis on closed orbit. |
| `m` | 1D array of float (shape: 3) | Transverse axis. |

::: pybmad.SpinEigenStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SpinEigenStruct

Fortran struct: `spin_eigen_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L286))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `vec` | 1D array of complex (shape: 8) |  |
| `val` | complex |  |

::: pybmad.SpinMatchingStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SpinMatchingStruct

Fortran struct: `spin_matching_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L297))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `axis` | [SpinAxisStruct](bmad.md#spinaxisstruct) |  |
| `eigen` | [1D array of SpinEigenStruct (shape: 8)](bmad.md#spineigenstruct) |  |
| `dn_dpz` | 1D array of float (shape: 3) | Invariant spin derivative |
| `alpha` | 1D array of float (shape: 6) | Alpha vector |
| `beta` | 1D array of float (shape: 6) | Beta vector |
| `orb0` | 1D array of float (shape: 6) | Closed orbit |
| `M_1turn` | 2D array of float (shape: 8,8) | 1-turn matrix |
| `M_ele` | 2D array of float (shape: 8,8) | Transfer matrix through element. |
| `sq_ele` | 1D array of float (shape: 0:3) |  |
| `sq_1turn` | 1D array of float (shape: 0:3) |  |
| `valid` | bool |  |

::: pybmad.SpinOrbitMap1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SpinOrbitMap1Struct

Fortran struct: `spin_orbit_map1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L322))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `orb_mat` | 2D array of float (shape: 6,6) | Orbital matrix |
| `vec0` | 1D array of float (shape: 6) | Orbital 0th order map: r_out = mat6 * r_in + vec0 |
| `spin_q` | 2D array of float (shape: 0:3,0:6) | 0th and 1st order quaternion spin map |

::: pybmad.SpinPolarStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SpinPolarStruct

Fortran struct: `spin_polar_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L313))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `polarization` | float |  |
| `theta` | float | Spherical coords: Angle from z-axis. |
| `phi` | float | Spherical coords: Angle in (x,y) plane. |
| `xi` | float | Spinor phase angle (See Bmad manual). |

::: pybmad.StrongBeamStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### StrongBeamStruct

Fortran struct: `strong_beam_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2070))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ix_slice` | int | 0 -> at element center and not at slice. |
| `x_center` | float | Strong beam slice center. |
| `y_center` | float | Strong beam slice center. |
| `x_sigma` | float | Strong beam slice sigma. |
| `y_sigma` | float | Strong beam slice sigma. |
| `dx` | float | Particle - beam slice distance. |
| `dy` | float | Particle - beam slice distance. |

::: pybmad.SummationRdtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SummationRdtStruct

Fortran struct: `summation_rdt_struct` ([`bmad/modules/srdt_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/srdt_mod.f90#L11))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `h11001` | complex |  |
| `h00111` | complex |  |
| `h20001` | complex |  |
| `h00201` | complex |  |
| `h10002` | complex |  |
| `h21000` | complex |  |
| `h30000` | complex |  |
| `h10110` | complex |  |
| `h10020` | complex |  |
| `h10200` | complex | 2nd order in K2 moments |
| `h31000` | complex |  |
| `h40000` | complex |  |
| `h20110` | complex |  |
| `h11200` | complex |  |
| `h20020` | complex |  |
| `h20200` | complex |  |
| `h00310` | complex |  |
| `h00400` | complex |  |
| `h22000` | complex |  |
| `h00220` | complex |  |
| `h11110` | complex |  |

::: pybmad.SurfaceCurvatureStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SurfaceCurvatureStruct

Fortran struct: `surface_curvature_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1063))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `xy` | 2D array of float (shape: 0:6,0:6) |  |
| `spherical` | float |  |
| `elliptical` | 1D array of float (shape: 3) | Total curvature = elliptical + spherical |
| `has_curvature` | bool | Dependent var. Will be set by Bmad |

::: pybmad.SurfaceDisplacementPtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SurfaceDisplacementPtStruct

Fortran struct: `surface_displacement_pt_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1030))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x0` | float | Position at center |
| `y0` | float | Position at center |
| `z0` | float |  |
| `dz_dx` | float |  |
| `dz_dy` | float |  |
| `d2z_dxdy` | float |  |

::: pybmad.SurfaceDisplacementStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SurfaceDisplacementStruct

Fortran struct: `surface_displacement_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1035))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `active` | bool |  |
| `dr` | 1D array of float (shape: 2) |  |
| `r0` | 1D array of float (shape: 2) |  |
| `pt` | [2D array of SurfaceDisplacementPtStruct](bmad.md#surfacedisplacementptstruct) |  |

::: pybmad.SurfaceHMisalignPtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SurfaceHMisalignPtStruct

Fortran struct: `surface_h_misalign_pt_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1017))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x0` | float | Position at center |
| `y0` | float | Position at center |
| `rot_y` | float | rot_t = x-rotation for Bragg and z-rotation for Laue. |
| `rot_t` | float | rot_t = x-rotation for Bragg and z-rotation for Laue. |
| `rot_y_rms` | float | rot_t = x-rotation for Bragg and z-rotation for Laue. |
| `rot_t_rms` | float | rot_t = x-rotation for Bragg and z-rotation for Laue. |

::: pybmad.SurfaceHMisalignStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SurfaceHMisalignStruct

Fortran struct: `surface_h_misalign_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1022))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `active` | bool |  |
| `dr` | 1D array of float (shape: 2) |  |
| `r0` | 1D array of float (shape: 2) |  |
| `pt` | [2D array of SurfaceHMisalignPtStruct](bmad.md#surfacehmisalignptstruct) |  |

::: pybmad.SurfaceSegmentedPtStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SurfaceSegmentedPtStruct

Fortran struct: `surface_segmented_pt_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1004))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x0` | float | Position at center |
| `y0` | float | Position at center |
| `z0` | float | Position at center |
| `dz_dx` | float | Slope at center |
| `dz_dy` | float | Slope at center |

::: pybmad.SurfaceSegmentedStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SurfaceSegmentedStruct

Fortran struct: `surface_segmented_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1009))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `active` | bool |  |
| `dr` | 1D array of float (shape: 2) |  |
| `r0` | 1D array of float (shape: 2) |  |
| `pt` | [2D array of SurfaceSegmentedPtStruct](bmad.md#surfacesegmentedptstruct) |  |

::: pybmad.TargetPointStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TargetPointStruct

Fortran struct: `target_point_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L1072))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `r` | 1D array of float (shape: 3) | (x, y, z) |

::: pybmad.TaylorStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaylorStruct

Fortran struct: `taylor_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L500))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ref` | float |  |
| `term` | [1D array of TaylorTermStruct](bmad.md#taylortermstruct) |  |

::: pybmad.TaylorTermStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TaylorTermStruct

Fortran struct: `taylor_term_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L486))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | float |  |
| `expn` | 1D array of int (shape: 6) |  |

::: pybmad.TrackPointStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TrackPointStruct

Fortran struct: `track_point_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2079))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `s_lab` | float | Longitudinal lab coord with respect to the upstream end. |
| `s_body` | float | Longitudinal body coord with respect to the entrance end. |
| `orb` | [CoordStruct](bmad.md#coordstruct) | Particle position in lab coords. |
| `field` | [EmFieldStruct](bmad.md#emfieldstruct) | E&M fields in lab coordinates. |
| `strong_beam` | [StrongBeamStruct](bmad.md#strongbeamstruct) | Strong beam info for beambeam element. |
| `vec0` | 1D array of float (shape: 6) | 0th order part of xfer map from the beginning. |
| `mat6` | 2D array of float (shape: 6,6) | 1st order part of xfer map (transfer matrix). |

::: pybmad.TrackStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TrackStruct

Fortran struct: `track_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2091))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `pt` | [1D array of TrackPointStruct](bmad.md#trackpointstruct) | Array of track points indexed from 0. |
| `ds_save` | float | Min distance between points. Not positive => Save at all points. |
| `n_pt` | int | Track upper bound for %pt(0:) array. n_bad and n_ok are used by adaptive trackers to record the number of times the step length had to be shortened. |
| `n_bad` | int | Number of "bad" steps where the step length was shortened. |
| `n_ok` | int | Number of "good" steps where the step length was not shortened. |

::: pybmad.TwissStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### TwissStruct

Fortran struct: `twiss_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L194))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `beta` | float |  |
| `alpha` | float |  |
| `gamma` | float |  |
| `phi` | float |  |
| `eta` | float |  |
| `etap` | float |  |
| `deta_ds` | float |  |
| `sigma` | float |  |
| `sigma_p` | float |  |
| `emit` | float |  |
| `norm_emit` | float |  |
| `chrom` | float |  |
| `dbeta_dpz` | float |  |
| `dalpha_dpz` | float |  |
| `deta_dpz` | float |  |
| `detap_dpz` | float |  |

::: pybmad.WakeLrModeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### WakeLrModeStruct

Fortran struct: `wake_lr_mode_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L656))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `freq` | float | Actual Frequency in Hz. |
| `freq_in` | float | Input frequency in Hz. |
| `R_over_Q` | float | Strength in V/C/m^(2*m_mode). |
| `Q` | float | Used for backwards compatability. |
| `damp` | float | Damping factor = omega / 2 * Q = pi * freq / Q |
| `phi` | float | Phase in radians/2pi. |
| `angle` | float | polarization angle (radians/2pi). |
| `b_sin` | float | non-skew sin-like component of the wake. |
| `b_cos` | float | non-skew cos-like component of the wake. |
| `a_sin` | float | skew sin-like component of the wake. |
| `a_cos` | float | skew cos-like component of the wake. |
| `m` | int | Mode order (1 = dipole, 2 = quad, etc.) |
| `polarized` | bool | Polaraized mode? |

::: pybmad.WakeLrStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### WakeLrStruct

Fortran struct: `wake_lr_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L672))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `file` | str |  |
| `mode` | [1D array of WakeLrModeStruct](bmad.md#wakelrmodestruct) |  |
| `t_ref` | float | time reference value for computing the wake amplitude. This is used to prevent value overflow with long trains. |
| `freq_spread` | float | Random frequency spread of long range modes. |
| `amp_scale` | float | Wake amplitude scale factor. |
| `time_scale` | float | time scale factor. |
| `self_wake_on` | bool | Long range self-wake used in tracking? |

::: pybmad.WakeSrModeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### WakeSrModeStruct

Fortran struct: `wake_sr_mode_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L626))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `amp` | float | Amplitude |
| `damp` | float | Dampling factor (1/m). |
| `k` | float | k factor |
| `phi` | float | Phase in radians/2pi |
| `b_sin` | float | non-skew (x) sin-like component of the wake |
| `b_cos` | float | non-skew (x) cos-like component of the wake |
| `a_sin` | float | skew (y) sin-like component of the wake |
| `a_cos` | float | skew (y) cos-like component of the wake |
| `polarization` | int | Transverse: none$, x_axis$, y_axis$. Not used for longitudinal. |
| `position_dependence` | int | Transverse: leading$, trailing$, none$ Longitudinal: x_leading$, ..., y_trailing$, none$ |

::: pybmad.WakeSrStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### WakeSrStruct

Fortran struct: `wake_sr_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L640))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `file` | str |  |
| `z_long` | [WakeSrZLongStruct](bmad.md#wakesrzlongstruct) |  |
| `long_wake` | [1D array of WakeSrModeStruct](bmad.md#wakesrmodestruct) |  |
| `trans_wake` | [1D array of WakeSrModeStruct](bmad.md#wakesrmodestruct) |  |
| `z_ref_long` | float | z reference value for computing the wake amplitude. |
| `z_ref_trans` | float | This is used to prevent value overflow with long bunches. |
| `z_max` | float | Max allowable z value. 0-> ignore |
| `amp_scale` | float | Wake amplitude scale factor. |
| `z_scale` | float | z-distance scale factor. |
| `scale_with_length` | bool | Scale wake with element length? |

::: pybmad.WakeSrZLongStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### WakeSrZLongStruct

Fortran struct: `wake_sr_z_long_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L614))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `w` | 1D array of float | Input single particle Wake. Indexed from 1. |
| `fw` | 1D array of complex | Fourier transform of w. |
| `fbunch` | 1D array of complex | Scratch space. |
| `w_out` | 1D array of complex | Scratch space. |
| `dz` | float | Distance between points. If zero there is no wake. |
| `z0` | float | Wake extent is [-z0, z0]. |
| `smoothing_sigma` | float | 0 => No smoothing. |
| `position_dependence` | int | Transverse: leading$, trailing$, none$ Longitudinal: x_leading$, ..., y_trailing$, none$ |
| `time_based` | bool | Was input time based? |

::: pybmad.WakeStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### WakeStruct

Fortran struct: `wake_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L685))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `sr` | [WakeSrStruct](bmad.md#wakesrstruct) | Short-range wake |
| `lr` | [WakeLrStruct](bmad.md#wakelrstruct) | Long-range wake |

::: pybmad.Wall3DSectionStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Wall3DSectionStruct

Fortran struct: `wall3d_section_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L439))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str | Identifying name |
| `material` | str | Material. |
| `v` | 1D array of Wall3dVertexStruct | Array of vertices. Always stored relative. |
| `surface` | [PhotonReflectSurfaceStruct](bmad.md#photonreflectsurfacestruct) | Surface reflectivity tables. |
| `type` | int | normal$, clear$, opaque$, wall_start$, wall_end$ |
| `n_vertex_input` | int | Number of vertices specified by the user. |
| `ix_ele` | int | index of lattice element containing section |
| `ix_branch` | int | Index of branch lattice element is in. |
| `vertices_state` | int | absolute$, or shifted_to_relative$. If set to absolute$ on input, will be changed to shifted_to_relative$ by section initalizer. |
| `patch_in_region` | bool | Patch element exists between this section and previous one? |
| `thickness` | float | Material thickness. |
| `s` | float | Longitudinal position |
| `r0` | 1D array of float (shape: 2) | Center of section Section-to-section spline interpolation of the center of the section |
| `dx0_ds` | float | Center of wall derivative |
| `dy0_ds` | float | Center of wall derivative |
| `x0_coef` | 1D array of float (shape: 0:3) | Spline coefs for x-center |
| `y0_coef` | 1D array of float (shape: 0:3) | Spline coefs for y-center Section-to_section spline interpolation of the wall. |
| `dr_ds` | float | derivative of wall radius |
| `p1_coef` | 1D array of float (shape: 3) | Spline coefs for p0 function |
| `p2_coef` | 1D array of float (shape: 3) | Spline coefs for p1 function |

::: pybmad.Wall3DStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Wall3DStruct

Fortran struct: `wall3d_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L471))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `type` | int | or mask_plate$ |
| `ix_wall3d` | int | Index in branch%wall3d(:) array. |
| `n_link` | int | For memory management of ele%wall3d |
| `thickness` | float | For diffraction_plate elements |
| `clear_material` | str |  |
| `opaque_material` | str |  |
| `superimpose` | bool | Can overlap another wall |
| `ele_anchor_pt` | int | anchor_beginning$, anchor_center$, or anchor_end$ |
| `section` | 1D array of Wall3dSectionStruct | Indexed from 1. |

::: pybmad.Wall3DVertexStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Wall3DVertexStruct

Fortran struct: `wall3d_vertex_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L423))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `x` | float | Coordinates of the vertex. |
| `y` | float | Coordinates of the vertex. |
| `radius_x` | float | Radius of arc or ellipse x-axis half width. 0 => Straight line. |
| `radius_y` | float | Ellipse y-axis half height. |
| `tilt` | float | Tilt of ellipse |
| `angle` | float | Angle of (x, y) point. |
| `x0` | float | Center of ellipse |
| `y0` | float | Center of ellipse |
| `type` | int | No longer used. |

::: pybmad.XyDispStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### XyDispStruct

Fortran struct: `xy_disp_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L874))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `eta` | float |  |
| `etap` | float |  |
| `deta_ds` | float |  |
| `sigma` | float |  |
| `deta_dpz` | float |  |
| `detap_dpz` | float |  |

## Procedures

### ab_multipole_kick

Fortran source: [`bmad/modules/multipole_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/multipole_mod.f90#L314)

::: pybmad.bmad.ab_multipole_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ab_multipole_kicks

Fortran source: [`bmad/modules/multipole_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/multipole_mod.f90#L82)

::: pybmad.bmad.ab_multipole_kicks
    options:
      show_root_heading: false
      show_root_toc_entry: false

### absolute_photon_position

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L70)

::: pybmad.bmad.absolute_photon_position
    options:
      show_root_heading: false
      show_root_toc_entry: false

### absolute_time_tracking

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L374)

::: pybmad.bmad.absolute_time_tracking
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ac_kicker_amp

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L381)

::: pybmad.bmad.ac_kicker_amp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### action_to_xyz

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L398)

::: pybmad.bmad.action_to_xyz
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_lattice_control_structs

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L390)

::: pybmad.bmad.add_lattice_control_structs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_ptc_layout_to_list

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L195)

::: pybmad.bmad.add_ptc_layout_to_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_superimpose

Fortran source: [`bmad/modules/superimpose_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/superimpose_mod.f90#L59)

::: pybmad.bmad.add_superimpose
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_this_multipass

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L2884)

::: pybmad.bmad.add_this_multipass
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_this_name_to_list

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L489)

::: pybmad.bmad.add_this_name_to_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_this_taylor_term

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L145)

::: pybmad.bmad.add_this_taylor_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### adjust_super_slave_names

Fortran source: [`bmad/modules/superimpose_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/superimpose_mod.f90#L813)

::: pybmad.bmad.adjust_super_slave_names
    options:
      show_root_heading: false
      show_root_toc_entry: false

### allocate_branch_array

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L398)

::: pybmad.bmad.allocate_branch_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### allocate_grid_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L412)

::: pybmad.bmad.allocate_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### allocate_lat_ele_array

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L418)

::: pybmad.bmad.allocate_lat_ele_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### allocate_plat

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L4296)

::: pybmad.bmad.allocate_plat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### angle_between_polars

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L437)

::: pybmad.bmad.angle_between_polars
    options:
      show_root_heading: false
      show_root_toc_entry: false

### angle_to_canonical_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L444)

::: pybmad.bmad.angle_to_canonical_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### aperture_bookkeeper

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bookkeeper_mod.f90#L1850)

::: pybmad.bmad.aperture_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apply_all_rampers

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L451)

::: pybmad.bmad.apply_all_rampers
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apply_element_edge_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L458)

::: pybmad.bmad.apply_element_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apply_energy_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L470)

::: pybmad.bmad.apply_energy_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apply_fft_3d_kicks

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L1666)

::: pybmad.bmad.apply_fft_3d_kicks
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apply_patch_to_ptc_fibre

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L3302)

::: pybmad.bmad.apply_patch_to_ptc_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apply_rampers_to_slave

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L479)

::: pybmad.bmad.apply_rampers_to_slave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### array_re_str

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L226)

::: pybmad.bmad.array_re_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### astra_max_field_reference

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/astra_interface_mod.f90#L979)

::: pybmad.bmad.astra_max_field_reference
    options:
      show_root_heading: false
      show_root_toc_entry: false

### at_this_ele_end

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L486)

::: pybmad.bmad.at_this_ele_end
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_bookkeeper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L493)

::: pybmad.bmad.attribute_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_free

Fortran sources (overloaded):

- `attribute_free1`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L3072)
- `attribute_free2`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L3101)
- `attribute_free3`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L3141)

::: pybmad.bmad.attribute_free
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_index

Fortran sources (overloaded):

- `attribute_index1`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L207)
- `attribute_index2`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L253)

::: pybmad.bmad.attribute_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_info

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L510)

::: pybmad.bmad.attribute_info
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_name

Fortran sources (overloaded):

- `attribute_name1`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L389)
- `attribute_name2`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L412)

::: pybmad.bmad.attribute_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_set_bookkeeping

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L500)

::: pybmad.bmad.attribute_set_bookkeeping
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_type

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L2037)

::: pybmad.bmad.attribute_type
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_units

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L2140)

::: pybmad.bmad.attribute_units
    options:
      show_root_heading: false
      show_root_toc_entry: false

### autoscale_phase_and_amp

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L509)

::: pybmad.bmad.autoscale_phase_and_amp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### average_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L518)

::: pybmad.bmad.average_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bane1

Fortran source: [`bmad/multiparticle/ibs_rates_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_rates_mod.f90#L226)

::: pybmad.bmad.bane1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bbi_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L525)

::: pybmad.bmad.bbi_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bbi_slice_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L532)

::: pybmad.bmad.bbi_slice_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### beam_envelope_ibs

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/envelope_mod.f90#L770)

::: pybmad.bmad.beam_envelope_ibs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### beam_equal_beam

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L5323)

::: pybmad.bmad.beam_equal_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### beam_init_setup

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L566)

::: pybmad.bmad.beam_init_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### beam_tilts

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L1020)

::: pybmad.bmad.beam_tilts
    options:
      show_root_heading: false
      show_root_toc_entry: false

### beambeam_fibre_setup

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L2711)

::: pybmad.bmad.beambeam_fibre_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/fringe_mod.f90#L39)

::: pybmad.bmad.bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_exact_multipole_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L540)

::: pybmad.bmad.bend_exact_multipole_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_length_has_been_set

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L551)

::: pybmad.bmad.bend_length_has_been_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_e_rel_init

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L915)

::: pybmad.bmad.bend_photon_e_rel_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_energy_integ_prob

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L241)

::: pybmad.bmad.bend_photon_energy_integ_prob
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_energy_normalized_probability

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L1088)

::: pybmad.bmad.bend_photon_energy_normalized_probability
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_init

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L143)

::: pybmad.bmad.bend_photon_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_polarization_init

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L388)

::: pybmad.bmad.bend_photon_polarization_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_vert_angle_init

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L446)

::: pybmad.bmad.bend_photon_vert_angle_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_shift

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L558)

::: pybmad.bmad.bend_shift
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_vert_angle_integ_prob

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L313)

::: pybmad.bmad.bend_vert_angle_integ_prob
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bjmt1

Fortran source: [`bmad/multiparticle/ibs_rates_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_rates_mod.f90#L43)

::: pybmad.bmad.bjmt1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bl_via_mat

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L882)

::: pybmad.bmad.bl_via_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bl_via_vlassov

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L844)

::: pybmad.bmad.bl_via_vlassov
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bmad_parser

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L586)

::: pybmad.bmad.bmad_parser
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bmad_parser2

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L597)

::: pybmad.bmad.bmad_parser2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bmad_parser_string_attribute_set

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L1356)

::: pybmad.bmad.bmad_parser_string_attribute_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bmad_patch_parameters_to_ptc

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L3043)

::: pybmad.bmad.bmad_patch_parameters_to_ptc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bp_set_ran_status

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L6053)

::: pybmad.bmad.bp_set_ran_status
    options:
      show_root_heading: false
      show_root_toc_entry: false

### branch_equal_branch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4839)

::: pybmad.bmad.branch_equal_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### branch_name

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L607)

::: pybmad.bmad.branch_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### branch_to_ptc_m_u

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L75)

::: pybmad.bmad.branch_to_ptc_m_u
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bunch_equal_bunch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L5270)

::: pybmad.bmad.bunch_equal_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### c_to_cbar

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L645)

::: pybmad.bmad.c_to_cbar
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_bunch_params

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_utils.f90#L1326)

::: pybmad.bmad.calc_bunch_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_bunch_params_slice

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_utils.f90#L1198)

::: pybmad.bmad.calc_bunch_params_slice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_bunch_params_z_slice

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_utils.f90#L1268)

::: pybmad.bmad.calc_bunch_params_z_slice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_bunch_sigma_matrix_etc

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_utils.f90#L1757)

::: pybmad.bmad.calc_bunch_sigma_matrix_etc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_emittances_and_twiss_from_sigma_matrix

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_utils.f90#L1454)

::: pybmad.bmad.calc_emittances_and_twiss_from_sigma_matrix
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_next_fringe_edge

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L652)

::: pybmad.bmad.calc_next_fringe_edge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_spin_params

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_utils.f90#L1708)

::: pybmad.bmad.calc_spin_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_super_slave_key

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L661)

::: pybmad.bmad.calc_super_slave_key
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_wall_radius

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/wall3d_mod.f90#L505)

::: pybmad.bmad.calc_wall_radius
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_wiggler_g_params

Fortran source: [`bmad/modules/rad_int_common.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/rad_int_common.f90#L475)

::: pybmad.bmad.calc_wiggler_g_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_z_tune

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L668)

::: pybmad.bmad.calc_z_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### canonical_to_angle_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L674)

::: pybmad.bmad.canonical_to_angle_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### capillary_photon_hit_spot_calc

Fortran source: [`bmad/photon/capillary_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/capillary_mod.f90#L253)

::: pybmad.bmad.capillary_photon_hit_spot_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### capillary_propagate_photon_a_step

Fortran source: [`bmad/photon/capillary_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/capillary_mod.f90#L178)

::: pybmad.bmad.capillary_propagate_photon_a_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### capillary_reflect_photon

Fortran source: [`bmad/photon/capillary_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/capillary_mod.f90#L364)

::: pybmad.bmad.capillary_reflect_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### capillary_track_photon_to_wall

Fortran source: [`bmad/photon/capillary_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/capillary_mod.f90#L92)

::: pybmad.bmad.capillary_track_photon_to_wall
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cbar_to_c

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L681)

::: pybmad.bmad.cbar_to_c
    options:
      show_root_heading: false
      show_root_toc_entry: false

### check_aperture_limit

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L688)

::: pybmad.bmad.check_aperture_limit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### check_controller_controls

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L787)

::: pybmad.bmad.check_controller_controls
    options:
      show_root_heading: false
      show_root_toc_entry: false

### check_for_superimpose_problem

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L3936)

::: pybmad.bmad.check_for_superimpose_problem
    options:
      show_root_heading: false
      show_root_toc_entry: false

### check_if_s_in_bounds

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L796)

::: pybmad.bmad.check_if_s_in_bounds
    options:
      show_root_heading: false
      show_root_toc_entry: false

### choose_quads_for_set_tune

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L806)

::: pybmad.bmad.choose_quads_for_set_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### chrom_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L816)

::: pybmad.bmad.chrom_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### chrom_tune

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L832)

::: pybmad.bmad.chrom_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cimp1

Fortran source: [`bmad/multiparticle/ibs_rates_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_rates_mod.f90#L633)

::: pybmad.bmad.cimp1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### classical_radius

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L699)

::: pybmad.bmad.classical_radius
    options:
      show_root_heading: false
      show_root_toc_entry: false

### clear_lat_1turn_mats

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L843)

::: pybmad.bmad.clear_lat_1turn_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### clear_taylor_maps_from_elements

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L849)

::: pybmad.bmad.clear_taylor_maps_from_elements
    options:
      show_root_heading: false
      show_root_toc_entry: false

### closed_orbit_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L855)

::: pybmad.bmad.closed_orbit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### closed_orbit_from_tracking

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L865)

::: pybmad.bmad.closed_orbit_from_tracking
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cmplx_re_str

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L251)

::: pybmad.bmad.cmplx_re_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### combine_consecutive_elements

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L876)

::: pybmad.bmad.combine_consecutive_elements
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_clean

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L132)

::: pybmad.bmad.complex_taylor_clean
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_coef

Fortran sources (overloaded):

- `complex_taylor_coef1`: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L179)
- `complex_taylor_coef2`: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L214)

::: pybmad.bmad.complex_taylor_coef
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_equal_complex_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L5136)

::: pybmad.bmad.complex_taylor_equal_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_exponent_index

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L661)

::: pybmad.bmad.complex_taylor_exponent_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_make_unit

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L440)

::: pybmad.bmad.complex_taylor_make_unit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_to_mat6

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L689)

::: pybmad.bmad.complex_taylor_to_mat6
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylors_equal_complex_taylors

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L5172)

::: pybmad.bmad.complex_taylors_equal_complex_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### compute_slave_coupler

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bookkeeper_mod.f90#L1520)

::: pybmad.bmad.compute_slave_coupler
    options:
      show_root_heading: false
      show_root_toc_entry: false

### compute_super_lord_s

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L3710)

::: pybmad.bmad.compute_super_lord_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### concat_ele_taylor

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L2147)

::: pybmad.bmad.concat_ele_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### concat_taylor

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L2084)

::: pybmad.bmad.concat_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### concat_transfer_mat

Fortran source: [`bmad/modules/transfer_map_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/transfer_map_mod.f90#L560)

::: pybmad.bmad.concat_transfer_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### control_bookkeeper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L883)

::: pybmad.bmad.control_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_bend_exact_multipole

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L926)

::: pybmad.bmad.convert_bend_exact_multipole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1023)

::: pybmad.bmad.convert_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_field_ele_to_lab

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/em_field_mod.f90#L752)

::: pybmad.bmad.convert_field_ele_to_lab
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_local_cartesian_to_local_curvilinear

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L1629)

::: pybmad.bmad.convert_local_cartesian_to_local_curvilinear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_local_curvilinear_to_local_cartesian

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L1612)

::: pybmad.bmad.convert_local_curvilinear_to_local_cartesian
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_particle_coordinates_s_to_t

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L891)

::: pybmad.bmad.convert_particle_coordinates_s_to_t
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_particle_coordinates_t_to_s

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L899)

::: pybmad.bmad.convert_particle_coordinates_t_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_pc_to

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L917)

::: pybmad.bmad.convert_pc_to
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_total_energy_to

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L908)

::: pybmad.bmad.convert_total_energy_to
    options:
      show_root_heading: false
      show_root_toc_entry: false

### converter_distribution_parser

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2064)

::: pybmad.bmad.converter_distribution_parser
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coord_equal_coord

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4888)

::: pybmad.bmad.coord_equal_coord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coord_state_name

Fortran source: [`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2635)

::: pybmad.bmad.coord_state_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_body_to_local

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L714)

::: pybmad.bmad.coords_body_to_local
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_body_to_rel_exit

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L705)

::: pybmad.bmad.coords_body_to_rel_exit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_curvilinear_to_floor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L778)

::: pybmad.bmad.coords_curvilinear_to_floor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_floor_to_curvilinear

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L748)

::: pybmad.bmad.coords_floor_to_curvilinear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_floor_to_local_curvilinear

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L738)

::: pybmad.bmad.coords_floor_to_local_curvilinear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_floor_to_relative

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L731)

::: pybmad.bmad.coords_floor_to_relative
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_local_curvilinear_to_body

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L758)

::: pybmad.bmad.coords_local_curvilinear_to_body
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_local_curvilinear_to_floor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L767)

::: pybmad.bmad.coords_local_curvilinear_to_floor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_relative_to_floor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L723)

::: pybmad.bmad.coords_relative_to_floor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cos_phi

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L1388)

::: pybmad.bmad.cos_phi
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coulombfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L189)

::: pybmad.bmad.coulombfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_concatenated_wall3d

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/wall3d_mod.f90#L1204)

::: pybmad.bmad.create_concatenated_wall3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_element_slice

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L941)

::: pybmad.bmad.create_element_slice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_feedback

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L933)

::: pybmad.bmad.create_feedback
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_field_overlap

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L953)

::: pybmad.bmad.create_field_overlap
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_girder

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L961)

::: pybmad.bmad.create_girder
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_group

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L971)

::: pybmad.bmad.create_group
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_lat_ele_nametable

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L979)

::: pybmad.bmad.create_lat_ele_nametable
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_overlay

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L986)

::: pybmad.bmad.create_overlay
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_planar_wiggler_model

Fortran source: [`bmad/modules/element_modeling_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/element_modeling_mod.f90#L114)

::: pybmad.bmad.create_planar_wiggler_model
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_ramper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L994)

::: pybmad.bmad.create_ramper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_sol_quad_model

Fortran source: [`bmad/modules/element_modeling_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/element_modeling_mod.f90#L62)

::: pybmad.bmad.create_sol_quad_model
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_unique_ele_names

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1002)

::: pybmad.bmad.create_unique_ele_names
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_wiggler_cartesian_map

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1010)

::: pybmad.bmad.create_wiggler_cartesian_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### crystal_attribute_bookkeeper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1017)

::: pybmad.bmad.crystal_attribute_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### crystal_diffraction_field_calc

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_utils_mod.f90#L256)

::: pybmad.bmad.crystal_diffraction_field_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### crystal_h_misalign

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L1054)

::: pybmad.bmad.crystal_h_misalign
    options:
      show_root_heading: false
      show_root_toc_entry: false

### crystal_type_to_crystal_params

Fortran source: [`bmad/interface/xraylib_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/xraylib_interface.f90#L314)

::: pybmad.bmad.crystal_type_to_crystal_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### csr_and_sc_apply_kicks

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L1504)

::: pybmad.bmad.csr_and_sc_apply_kicks
    options:
      show_root_heading: false
      show_root_toc_entry: false

### csr_bin_kicks

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L692)

::: pybmad.bmad.csr_bin_kicks
    options:
      show_root_heading: false
      show_root_toc_entry: false

### csr_bin_particles

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L425)

::: pybmad.bmad.csr_bin_particles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cumulr

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L1113)

::: pybmad.bmad.cumulr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### custom_attribute_ubound_index

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L2832)

::: pybmad.bmad.custom_attribute_ubound_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### custom_ele_attrib_name_list

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L3004)

::: pybmad.bmad.custom_ele_attrib_name_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### d_integral

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L1088)

::: pybmad.bmad.d_integral
    options:
      show_root_heading: false
      show_root_toc_entry: false

### damping_matrix_d

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/envelope_mod.f90#L148)

::: pybmad.bmad.damping_matrix_d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ddz_calc_csr

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L1014)

::: pybmad.bmad.ddz_calc_csr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### deallocate_ele_pointers

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1040)

::: pybmad.bmad.deallocate_ele_pointers
    options:
      show_root_heading: false
      show_root_toc_entry: false

### deallocate_expression_tree

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L664)

::: pybmad.bmad.deallocate_expression_tree
    options:
      show_root_heading: false
      show_root_toc_entry: false

### deallocate_lat_pointers

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1047)

::: pybmad.bmad.deallocate_lat_pointers
    options:
      show_root_heading: false
      show_root_toc_entry: false

### default_tracking_species

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1053)

::: pybmad.bmad.default_tracking_species
    options:
      show_root_heading: false
      show_root_toc_entry: false

### deposit_particles

Fortran source: [`bmad/space_charge/open_spacecharge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_mod.f90#L159)

::: pybmad.bmad.deposit_particles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### detector_pixel_pt

Fortran source: [`bmad/photon/photon_target_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_target_mod.f90#L320)

::: pybmad.bmad.detector_pixel_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### diffraction_plate_or_mask_hit_spot

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1060)

::: pybmad.bmad.diffraction_plate_or_mask_hit_spot
    options:
      show_root_heading: false
      show_root_toc_entry: false

### diffusion_matrix_b

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/envelope_mod.f90#L127)

::: pybmad.bmad.diffusion_matrix_b
    options:
      show_root_heading: false
      show_root_toc_entry: false

### distance_to_aperture

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1068)

::: pybmad.bmad.distance_to_aperture
    options:
      show_root_heading: false
      show_root_toc_entry: false

### do_mode_flip

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1078)

::: pybmad.bmad.do_mode_flip
    options:
      show_root_heading: false
      show_root_toc_entry: false

### dpc_given_de

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1085)

::: pybmad.bmad.dpc_given_de
    options:
      show_root_heading: false
      show_root_toc_entry: false

### drift_and_pipe_track_methods_adjustment

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L4930)

::: pybmad.bmad.drift_and_pipe_track_methods_adjustment
    options:
      show_root_heading: false
      show_root_toc_entry: false

### drift_multipass_name_correction

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L3064)

::: pybmad.bmad.drift_multipass_name_correction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### drift_orbit_time

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/time_tracker_mod.f90#L925)

::: pybmad.bmad.drift_orbit_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### drift_particle_to_s

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/space_charge_mod.f90#L638)

::: pybmad.bmad.drift_particle_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### drift_particle_to_t

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/space_charge_mod.f90#L685)

::: pybmad.bmad.drift_particle_to_t
    options:
      show_root_heading: false
      show_root_toc_entry: false

### dspline_len

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L1765)

::: pybmad.bmad.dspline_len
    options:
      show_root_heading: false
      show_root_toc_entry: false

### dynamic_aperture_point

Fortran source: [`bmad/modules/dynamic_aperture_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/dynamic_aperture_mod.f90#L284)

::: pybmad.bmad.dynamic_aperture_point
    options:
      show_root_heading: false
      show_root_toc_entry: false

### dynamic_aperture_scan

Fortran source: [`bmad/modules/dynamic_aperture_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/dynamic_aperture_mod.f90#L27)

::: pybmad.bmad.dynamic_aperture_scan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### e_accel_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1091)

::: pybmad.bmad.e_accel_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### e_crit_photon

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L1375)

::: pybmad.bmad.e_crit_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eigen_decomp_6mat

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L449)

::: pybmad.bmad.eigen_decomp_6mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_compute_ref_energy_and_time

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1100)

::: pybmad.bmad.ele_compute_ref_energy_and_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_equal_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4262)

::: pybmad.bmad.ele_equal_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_equals_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4293)

::: pybmad.bmad.ele_equals_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_finalizer

Fortran source: [`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2815)

::: pybmad.bmad.ele_finalizer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_full_name

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1109)

::: pybmad.bmad.ele_full_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_geometry

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1117)

::: pybmad.bmad.ele_geometry
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_geometry_with_misalignments

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1127)

::: pybmad.bmad.ele_geometry_with_misalignments
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_has_constant_ds_dt_ref

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1135)

::: pybmad.bmad.ele_has_constant_ds_dt_ref
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_has_nonzero_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1142)

::: pybmad.bmad.ele_has_nonzero_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_has_nonzero_offset

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1149)

::: pybmad.bmad.ele_has_nonzero_offset
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_is_monitor

Fortran source: [`bmad/modules/measurement_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/measurement_mod.f90#L27)

::: pybmad.bmad.ele_is_monitor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_loc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1172)

::: pybmad.bmad.ele_loc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_loc_name

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1156)

::: pybmad.bmad.ele_loc_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_misalignment_l_s_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1165)

::: pybmad.bmad.ele_misalignment_l_s_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_nametable_index

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1179)

::: pybmad.bmad.ele_nametable_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_order_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1186)

::: pybmad.bmad.ele_order_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_reference_energy_correction

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1193)

::: pybmad.bmad.ele_reference_energy_correction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_rf_step_index

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1203)

::: pybmad.bmad.ele_rf_step_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_to_fibre

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1211)

::: pybmad.bmad.ele_to_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_to_ptc_magnetic_bn_an

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L3098)

::: pybmad.bmad.ele_to_ptc_magnetic_bn_an
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_to_spin_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1222)

::: pybmad.bmad.ele_to_spin_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_to_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1230)

::: pybmad.bmad.ele_to_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_unique_name

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1239)

::: pybmad.bmad.ele_unique_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_value_has_changed

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1247)

::: pybmad.bmad.ele_value_has_changed
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_vec_equal_ele_vec

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4673)

::: pybmad.bmad.ele_vec_equal_ele_vec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### elec_multipole_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1256)

::: pybmad.bmad.elec_multipole_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### element_at_s

Fortran sources (overloaded):

- `element_at_s_branch`: [`bmad/modules/element_at_s_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/element_at_s_mod.f90#L75)
- `element_at_s_lat`: [`bmad/modules/element_at_s_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/element_at_s_mod.f90#L198)

::: pybmad.bmad.element_at_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### element_slice_iterator

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1267)

::: pybmad.bmad.element_slice_iterator
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ellipinc_test

Fortran source: [`bmad/space_charge/csr3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr3d_mod.f90#L628)

::: pybmad.bmad.ellipinc_test
    options:
      show_root_heading: false
      show_root_toc_entry: false

### em_field_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1276)

::: pybmad.bmad.em_field_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### em_field_derivatives

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/em_field_mod.f90#L655)

::: pybmad.bmad.em_field_derivatives
    options:
      show_root_heading: false
      show_root_toc_entry: false

### em_field_kick_vector_time

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/time_tracker_mod.f90#L692)

::: pybmad.bmad.em_field_kick_vector_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### em_field_plus_em_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4228)

::: pybmad.bmad.em_field_plus_em_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### emit_6d

Fortran source: [`bmad/modules/rad_6d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/rad_6d_mod.f90#L50)

::: pybmad.bmad.emit_6d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### energy_func

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L279)

::: pybmad.bmad.energy_func
    options:
      show_root_heading: false
      show_root_toc_entry: false

### entering_element

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1292)

::: pybmad.bmad.entering_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### envelope_radints

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/envelope_mod.f90#L614)

::: pybmad.bmad.envelope_radints
    options:
      show_root_heading: false
      show_root_toc_entry: false

### envelope_radints_ibs

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/envelope_mod.f90#L510)

::: pybmad.bmad.envelope_radints_ibs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ac_kicker

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L135)

::: pybmad.bmad.eq_ac_kicker
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ac_kicker_freq

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L113)

::: pybmad.bmad.eq_ac_kicker_freq
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ac_kicker_time

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L91)

::: pybmad.bmad.eq_ac_kicker_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_anormal_mode

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2264)

::: pybmad.bmad.eq_anormal_mode
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_aperture_param

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L3383)

::: pybmad.bmad.eq_aperture_param
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_aperture_point

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L3357)

::: pybmad.bmad.eq_aperture_point
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_aperture_scan

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L3417)

::: pybmad.bmad.eq_aperture_scan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_beam

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L3335)

::: pybmad.bmad.eq_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_beam_init

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2076)

::: pybmad.bmad.eq_beam_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_bmad_common

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2540)

::: pybmad.bmad.eq_bmad_common
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_bookkeeping_state

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1172)

::: pybmad.bmad.eq_bookkeeping_state
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_bpm_phase_coupling

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L353)

::: pybmad.bmad.eq_bpm_phase_coupling
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_branch

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L3040)

::: pybmad.bmad.eq_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_bunch

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L3221)

::: pybmad.bmad.eq_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_bunch_params

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L3273)

::: pybmad.bmad.eq_bunch_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cartesian_map

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L805)

::: pybmad.bmad.eq_cartesian_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cartesian_map_term

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L779)

::: pybmad.bmad.eq_cartesian_map_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cartesian_map_term1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L745)

::: pybmad.bmad.eq_cartesian_map_term1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_complex_taylor

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L3016)

::: pybmad.bmad.eq_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_complex_taylor_term

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2996)

::: pybmad.bmad.eq_complex_taylor_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_control

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1868)

::: pybmad.bmad.eq_control
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_control_ramp1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1930)

::: pybmad.bmad.eq_control_ramp1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_control_var1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1908)

::: pybmad.bmad.eq_control_var1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_controller

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1964)

::: pybmad.bmad.eq_controller
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_coord

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L273)

::: pybmad.bmad.eq_coord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_coord_array

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L331)

::: pybmad.bmad.eq_coord_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cylindrical_map

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L882)

::: pybmad.bmad.eq_cylindrical_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cylindrical_map_term

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L856)

::: pybmad.bmad.eq_cylindrical_map_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cylindrical_map_term1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L836)

::: pybmad.bmad.eq_cylindrical_map_term1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ele

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2738)

::: pybmad.bmad.eq_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ellipse_beam_init

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2004)

::: pybmad.bmad.eq_ellipse_beam_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_em_field

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2366)

::: pybmad.bmad.eq_em_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_expression_atom

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L389)

::: pybmad.bmad.eq_expression_atom
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_floor_position

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1010)

::: pybmad.bmad.eq_floor_position
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_gen_grad_curve

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1254)

::: pybmad.bmad.eq_gen_grad_curve
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_gen_gradients

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1282)

::: pybmad.bmad.eq_gen_gradients
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_gg_taylor

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L721)

::: pybmad.bmad.eq_gg_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_gg_taylor_term

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L701)

::: pybmad.bmad.eq_gg_taylor_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_grid_beam_init

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2048)

::: pybmad.bmad.eq_grid_beam_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_grid_field

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L967)

::: pybmad.bmad.eq_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_grid_field_pt

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L941)

::: pybmad.bmad.eq_grid_field_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_grid_field_pt1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L921)

::: pybmad.bmad.eq_grid_field_pt1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_high_energy_space_charge

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1036)

::: pybmad.bmad.eq_high_energy_space_charge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_interval1_coef

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L163)

::: pybmad.bmad.eq_interval1_coef
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_kv_beam_init

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2026)

::: pybmad.bmad.eq_kv_beam_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_lat

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L3094)

::: pybmad.bmad.eq_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_lat_ele_loc

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L617)

::: pybmad.bmad.eq_lat_ele_loc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_lat_param

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2162)

::: pybmad.bmad.eq_lat_param
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_linac_normal_mode

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2294)

::: pybmad.bmad.eq_linac_normal_mode
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_mode3

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1144)

::: pybmad.bmad.eq_mode3
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_mode_info

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2212)

::: pybmad.bmad.eq_mode_info
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_normal_modes

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2324)

::: pybmad.bmad.eq_normal_modes
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_photon_element

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1658)

::: pybmad.bmad.eq_photon_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_photon_material

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1558)

::: pybmad.bmad.eq_photon_material
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_photon_reflect_surface

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L239)

::: pybmad.bmad.eq_photon_reflect_surface
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_photon_reflect_table

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L185)

::: pybmad.bmad.eq_photon_reflect_table
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_photon_target

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1532)

::: pybmad.bmad.eq_photon_target
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_pixel_detec

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1626)

::: pybmad.bmad.eq_pixel_detec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_pixel_pt

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1590)

::: pybmad.bmad.eq_pixel_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_pre_tracker

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2240)

::: pybmad.bmad.eq_pre_tracker
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_rad_int1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2642)

::: pybmad.bmad.eq_rad_int1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_rad_int_all_ele

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2716)

::: pybmad.bmad.eq_rad_int_all_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_rad_int_branch

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2694)

::: pybmad.bmad.eq_rad_int_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_rad_map

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1206)

::: pybmad.bmad.eq_rad_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_rad_map_ele

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1232)

::: pybmad.bmad.eq_rad_map_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ramper_lord

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1843)

::: pybmad.bmad.eq_ramper_lord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_space_charge_common

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2486)

::: pybmad.bmad.eq_space_charge_common
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_spin_polar

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L67)

::: pybmad.bmad.eq_spin_polar
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_spline

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L43)

::: pybmad.bmad.eq_spline
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_strong_beam

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2396)

::: pybmad.bmad.eq_strong_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_curvature

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1508)

::: pybmad.bmad.eq_surface_curvature
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_displacement

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1462)

::: pybmad.bmad.eq_surface_displacement
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_displacement_pt

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1434)

::: pybmad.bmad.eq_surface_displacement_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_h_misalign

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1406)

::: pybmad.bmad.eq_surface_h_misalign
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_h_misalign_pt

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1378)

::: pybmad.bmad.eq_surface_h_misalign_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_segmented

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1350)

::: pybmad.bmad.eq_surface_segmented
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_segmented_pt

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1324)

::: pybmad.bmad.eq_surface_segmented_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_target_point

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1490)

::: pybmad.bmad.eq_target_point
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_taylor

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L677)

::: pybmad.bmad.eq_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_taylor_term

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L657)

::: pybmad.bmad.eq_taylor_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_track

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2456)

::: pybmad.bmad.eq_track
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_track_point

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L2426)

::: pybmad.bmad.eq_track_point
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_twiss

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1096)

::: pybmad.bmad.eq_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L637)

::: pybmad.bmad.eq_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake_lr

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L583)

::: pybmad.bmad.eq_wake_lr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake_lr_mode

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L541)

::: pybmad.bmad.eq_wake_lr_mode
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake_sr

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L497)

::: pybmad.bmad.eq_wake_sr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake_sr_mode

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L461)

::: pybmad.bmad.eq_wake_sr_mode
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake_sr_z_long

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L411)

::: pybmad.bmad.eq_wake_sr_z_long
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wall3d

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1803)

::: pybmad.bmad.eq_wall3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wall3d_section

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1740)

::: pybmad.bmad.eq_wall3d_section
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wall3d_vertex

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1706)

::: pybmad.bmad.eq_wall3d_vertex
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_xy_disp

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/equality_mod.f90#L1068)

::: pybmad.bmad.eq_xy_disp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### equal_sign_here

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L8016)

::: pybmad.bmad.equal_sign_here
    options:
      show_root_heading: false
      show_root_toc_entry: false

### equivalent_taylor_attributes

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1300)

::: pybmad.bmad.equivalent_taylor_attributes
    options:
      show_root_heading: false
      show_root_toc_entry: false

### etdiv

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/envelope_mod.f90#L1736)

::: pybmad.bmad.etdiv
    options:
      show_root_heading: false
      show_root_toc_entry: false

### evaluate_array_index

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L808)

::: pybmad.bmad.evaluate_array_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### evaluate_logical

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L859)

::: pybmad.bmad.evaluate_logical
    options:
      show_root_heading: false
      show_root_toc_entry: false

### exact_bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/fringe_mod.f90#L1467)

::: pybmad.bmad.exact_bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### exp_bessi0

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/touschek_mod.f90#L662)

::: pybmad.bmad.exp_bessi0
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expect_one_of

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L7958)

::: pybmad.bmad.expect_one_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expect_this

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L7817)

::: pybmad.bmad.expect_this
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_stack_to_string

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L1773)

::: pybmad.bmad.expression_stack_to_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_stack_value

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L1489)

::: pybmad.bmad.expression_stack_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_string_to_stack

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L903)

::: pybmad.bmad.expression_string_to_stack
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_string_to_tree

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L84)

::: pybmad.bmad.expression_string_to_tree
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_tree_to_string

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L738)

::: pybmad.bmad.expression_tree_to_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_value

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L1421)

::: pybmad.bmad.expression_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fft1

Fortran source: [`bmad/space_charge/fast_fourier_am.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/fast_fourier_am.f90#L37)

::: pybmad.bmad.fft1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fibre_to_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1307)

::: pybmad.bmad.fibre_to_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### field_attribute_free

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L3474)

::: pybmad.bmad.field_attribute_free
    options:
      show_root_heading: false
      show_root_toc_entry: false

### finalize_reflectivity_table

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L365)

::: pybmad.bmad.finalize_reflectivity_table
    options:
      show_root_heading: false
      show_root_toc_entry: false

### find_element_ends

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1317)

::: pybmad.bmad.find_element_ends
    options:
      show_root_heading: false
      show_root_toc_entry: false

### find_fwhm

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/longitudinal_profile_mod.f90#L401)

::: pybmad.bmad.find_fwhm
    options:
      show_root_heading: false
      show_root_toc_entry: false

### find_matching_fieldmap

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1325)

::: pybmad.bmad.find_matching_fieldmap
    options:
      show_root_heading: false
      show_root_toc_entry: false

### find_normalization

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/longitudinal_profile_mod.f90#L355)

::: pybmad.bmad.find_normalization
    options:
      show_root_heading: false
      show_root_toc_entry: false

### floor_angles_to_w_mat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1335)

::: pybmad.bmad.floor_angles_to_w_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### floor_w_mat_to_angles

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1342)

::: pybmad.bmad.floor_w_mat_to_angles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### form_complex_taylor

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L1593)

::: pybmad.bmad.form_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### form_digested_bmad_file_name

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L5436)

::: pybmad.bmad.form_digested_bmad_file_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fringe_here

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1349)

::: pybmad.bmad.fringe_here
    options:
      show_root_heading: false
      show_root_toc_entry: false

### g_bend_from_em_field

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/em_field_mod.f90#L63)

::: pybmad.bmad.g_bend_from_em_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### g_bending_strength_from_em_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1358)

::: pybmad.bmad.g_bending_strength_from_em_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### g_integrals_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1453)

::: pybmad.bmad.g_integrals_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gamma_ref

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1370)

::: pybmad.bmad.gamma_ref
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gen_grad_at_s_to_gg_a_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1386)

::: pybmad.bmad.gen_grad_at_s_to_gg_a_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gen_grad_at_s_to_gg_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1377)

::: pybmad.bmad.gen_grad_at_s_to_gg_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_astra_fieldgrid_name_and_scaling

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/astra_interface_mod.f90#L400)

::: pybmad.bmad.get_astra_fieldgrid_name_and_scaling
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_bl_from_fwhm

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/longitudinal_profile_mod.f90#L518)

::: pybmad.bmad.get_bl_from_fwhm
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_called_file

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L68)

::: pybmad.bmad.get_called_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_emit_from_sigma_mat

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L946)

::: pybmad.bmad.get_emit_from_sigma_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_gpt_fieldgrid_name_and_scaling

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L661)

::: pybmad.bmad.get_gpt_fieldgrid_name_and_scaling
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_list_of_names

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L2271)

::: pybmad.bmad.get_list_of_names
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_next_word

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L257)

::: pybmad.bmad.get_next_word
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_opal_fieldgrid_name_and_scaling

Fortran source: [`bmad/interface/opal_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/opal_interface_mod.f90#L359)

::: pybmad.bmad.get_opal_fieldgrid_name_and_scaling
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_overlay_group_names

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L2344)

::: pybmad.bmad.get_overlay_group_names
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_sequence_args

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L4026)

::: pybmad.bmad.get_sequence_args
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_slave_list

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1396)

::: pybmad.bmad.get_slave_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_switch

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L7858)

::: pybmad.bmad.get_switch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_coef_table_init

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L50)

::: pybmad.bmad.gg_coef_table_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_001

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L69)

::: pybmad.bmad.gg_set_block_001
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_002

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L672)

::: pybmad.bmad.gg_set_block_002
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_003

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L1275)

::: pybmad.bmad.gg_set_block_003
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_004

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L1878)

::: pybmad.bmad.gg_set_block_004
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_005

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L2481)

::: pybmad.bmad.gg_set_block_005
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_006

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L3084)

::: pybmad.bmad.gg_set_block_006
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_007

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L3687)

::: pybmad.bmad.gg_set_block_007
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_008

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L4290)

::: pybmad.bmad.gg_set_block_008
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_009

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L4893)

::: pybmad.bmad.gg_set_block_009
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_010

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L5496)

::: pybmad.bmad.gg_set_block_010
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_011

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L6099)

::: pybmad.bmad.gg_set_block_011
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_012

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L6702)

::: pybmad.bmad.gg_set_block_012
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_013

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L7305)

::: pybmad.bmad.gg_set_block_013
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_set_block_014

Fortran source: [`bmad/modules/gg_coef_table_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/gg_coef_table_mod.f90#L7908)

::: pybmad.bmad.gg_set_block_014
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_taylor_equal_gg_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4993)

::: pybmad.bmad.gg_taylor_equal_gg_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_taylors_equal_gg_taylors

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L5039)

::: pybmad.bmad.gg_taylors_equal_gg_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gpt_field_grid_scaling

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L747)

::: pybmad.bmad.gpt_field_grid_scaling
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gpt_max_field_reference

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L1555)

::: pybmad.bmad.gpt_max_field_reference
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gpt_to_particle_bunch

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L35)

::: pybmad.bmad.gpt_to_particle_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gradient_shift_sr_wake

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1404)

::: pybmad.bmad.gradient_shift_sr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### grid_field_interpolate

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/em_field_mod.f90#L222)

::: pybmad.bmad.grid_field_interpolate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hard_multipole_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/fringe_mod.f90#L691)

::: pybmad.bmad.hard_multipole_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### has_attribute

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L2804)

::: pybmad.bmad.has_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### has_curvature

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_utils_mod.f90#L31)

::: pybmad.bmad.has_curvature
    options:
      show_root_heading: false
      show_root_toc_entry: false

### has_orientation_attributes

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L1988)

::: pybmad.bmad.has_orientation_attributes
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hdf5_read_beam

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1412)

::: pybmad.bmad.hdf5_read_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hdf5_read_grid_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1423)

::: pybmad.bmad.hdf5_read_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hdf5_write_beam

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1434)

::: pybmad.bmad.hdf5_write_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hdf5_write_grid_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1444)

::: pybmad.bmad.hdf5_write_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hwang_bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/fringe_mod.f90#L240)

::: pybmad.bmad.hwang_bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### i_csr

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L1359)

::: pybmad.bmad.i_csr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ibs1

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L607)

::: pybmad.bmad.ibs1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ibs_blowup1turn

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L536)

::: pybmad.bmad.ibs_blowup1turn
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ibs_delta_calc

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L410)

::: pybmad.bmad.ibs_delta_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ibs_equib_der

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L246)

::: pybmad.bmad.ibs_equib_der
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ibs_equib_rlx

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L73)

::: pybmad.bmad.ibs_equib_rlx
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ibs_lifetime

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L364)

::: pybmad.bmad.ibs_lifetime
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ibs_matrix_c

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/envelope_mod.f90#L733)

::: pybmad.bmad.ibs_matrix_c
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ibs_rates1turn

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L453)

::: pybmad.bmad.ibs_rates1turn
    options:
      show_root_heading: false
      show_root_toc_entry: false

### igfcoulombfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L205)

::: pybmad.bmad.igfcoulombfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### igfexfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L243)

::: pybmad.bmad.igfexfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### igfeyfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L288)

::: pybmad.bmad.igfeyfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### igfezfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L333)

::: pybmad.bmad.igfezfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### image_charge_kick_calc

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L1442)

::: pybmad.bmad.image_charge_kick_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_attribute_name1

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L1942)

::: pybmad.bmad.init_attribute_name1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_attribute_name_array

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L572)

::: pybmad.bmad.init_attribute_name_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_beam_distribution

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_utils.f90#L212)

::: pybmad.bmad.init_beam_distribution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_bmad

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1468)

::: pybmad.bmad.init_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_bmad_parser_common

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1473)

::: pybmad.bmad.init_bmad_parser_common
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_bunch_distribution

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_utils.f90#L328)

::: pybmad.bmad.init_bunch_distribution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_complex_taylor_series

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L5210)

::: pybmad.bmad.init_complex_taylor_series
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_coord

Fortran sources (overloaded):

- `init_coord1`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L247)
- `init_coord2`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L258)
- `init_coord3`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L268)

::: pybmad.bmad.init_coord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_custom

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1479)

::: pybmad.bmad.init_custom
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1485)

::: pybmad.bmad.init_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_fringe_info

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1494)

::: pybmad.bmad.init_fringe_info
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_gg_taylor_series

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L5077)

::: pybmad.bmad.init_gg_taylor_series
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1503)

::: pybmad.bmad.init_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_multipole_cache

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1511)

::: pybmad.bmad.init_multipole_cache
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_photon_from_a_photon_init_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1459)

::: pybmad.bmad.init_photon_from_a_photon_init_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_photon_integ_prob

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L1415)

::: pybmad.bmad.init_photon_integ_prob
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_spin_distribution

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_utils.f90#L1141)

::: pybmad.bmad.init_spin_distribution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_surface_segment

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L8222)

::: pybmad.bmad.init_surface_segment
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_taylor_series

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1517)

::: pybmad.bmad.init_taylor_series
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_wake

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1525)

::: pybmad.bmad.init_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### insert_element

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1533)

::: pybmad.bmad.insert_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### integrand_base

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/touschek_mod.f90#L585)

::: pybmad.bmad.integrand_base
    options:
      show_root_heading: false
      show_root_toc_entry: false

### integrate_psi

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/longitudinal_profile_mod.f90#L290)

::: pybmad.bmad.integrate_psi
    options:
      show_root_heading: false
      show_root_toc_entry: false

### integrated_mats

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/envelope_mod.f90#L409)

::: pybmad.bmad.integrated_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### integration_timer

Fortran sources (overloaded):

- `integration_timer_ele`: [`bmad/modules/integration_timer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/integration_timer_mod.f90#L45)
- `integration_timer_fibre`: [`bmad/modules/integration_timer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/integration_timer_mod.f90#L74)

::: pybmad.bmad.integration_timer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### interpolate_field

Fortran source: [`bmad/space_charge/open_spacecharge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_mod.f90#L417)

::: pybmad.bmad.interpolate_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ion_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1543)

::: pybmad.bmad.ion_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_attribute

Fortran source: [`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2677)

::: pybmad.bmad.is_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### key_name_to_key_index

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1551)

::: pybmad.bmad.key_name_to_key_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### kick_vector_calc

Fortran source: [`bmad/modules/runge_kutta_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/runge_kutta_mod.f90#L663)

::: pybmad.bmad.kick_vector_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### kill_complex_taylor

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L565)

::: pybmad.bmad.kill_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### kill_ptc_layouts

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1559)

::: pybmad.bmad.kill_ptc_layouts
    options:
      show_root_heading: false
      show_root_toc_entry: false

### kill_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1565)

::: pybmad.bmad.kill_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### kind_name

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L944)

::: pybmad.bmad.kind_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### knot_interpolate

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1571)

::: pybmad.bmad.knot_interpolate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### knots_to_string

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1579)

::: pybmad.bmad.knots_to_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lafun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L229)

::: pybmad.bmad.lafun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_compute_ref_energy_and_time

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1586)

::: pybmad.bmad.lat_compute_ref_energy_and_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_ele_locator

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1593)

::: pybmad.bmad.lat_ele_locator
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_equal_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4716)

::: pybmad.bmad.lat_equal_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_geometry

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1605)

::: pybmad.bmad.lat_geometry
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_make_mat6

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1611)

::: pybmad.bmad.lat_make_mat6
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_sanity_check

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1620)

::: pybmad.bmad.lat_sanity_check
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_to_ptc_layout

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1627)

::: pybmad.bmad.lat_to_ptc_layout
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_vec_equal_lat_vec

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4797)

::: pybmad.bmad.lat_vec_equal_lat_vec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lattice_bookkeeper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1633)

::: pybmad.bmad.lattice_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lcavity_rf_step_setup

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1640)

::: pybmad.bmad.lcavity_rf_step_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### linear_bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/fringe_mod.f90#L151)

::: pybmad.bmad.linear_bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### linear_coef

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L2014)

::: pybmad.bmad.linear_coef
    options:
      show_root_heading: false
      show_root_toc_entry: false

### linear_to_spin_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1646)

::: pybmad.bmad.linear_to_spin_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### load_parse_line

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L595)

::: pybmad.bmad.load_parse_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lord_edge_aligned

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1652)

::: pybmad.bmad.lord_edge_aligned
    options:
      show_root_heading: false
      show_root_toc_entry: false

### low_energy_z_correction

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1660)

::: pybmad.bmad.low_energy_z_correction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lsc_kick_params_calc

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L1098)

::: pybmad.bmad.lsc_kick_params_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_add_offsets_and_multipoles

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L185)

::: pybmad.bmad.mad_add_offsets_and_multipoles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_concat_map2

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1473)

::: pybmad.bmad.mad_concat_map2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_drift

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L336)

::: pybmad.bmad.mad_drift
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_elsep

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L394)

::: pybmad.bmad.mad_elsep
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_map_to_taylor

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1651)

::: pybmad.bmad.mad_map_to_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_quadrupole

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1049)

::: pybmad.bmad.mad_quadrupole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_rfcavity

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1165)

::: pybmad.bmad.mad_rfcavity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_sbend

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L614)

::: pybmad.bmad.mad_sbend
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_sbend_body

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L762)

::: pybmad.bmad.mad_sbend_body
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_sbend_fringe

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L678)

::: pybmad.bmad.mad_sbend_fringe
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_sextupole

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L523)

::: pybmad.bmad.mad_sextupole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_solenoid

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1215)

::: pybmad.bmad.mad_solenoid
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_tmfoc

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1000)

::: pybmad.bmad.mad_tmfoc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_tmsymm

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1327)

::: pybmad.bmad.mad_tmsymm
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_tmtilt

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1367)

::: pybmad.bmad.mad_tmtilt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_track1

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1554)

::: pybmad.bmad.mad_track1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_g2_mats

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1678)

::: pybmad.bmad.make_g2_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_g_mats

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1670)

::: pybmad.bmad.make_g_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_hvbp

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L198)

::: pybmad.bmad.make_hvbp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_hybrid_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1685)

::: pybmad.bmad.make_hybrid_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mad_map

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L98)

::: pybmad.bmad.make_mad_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1706)

::: pybmad.bmad.make_mat6
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_bmad

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1723)

::: pybmad.bmad.make_mat6_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_bmad_photon

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1732)

::: pybmad.bmad.make_mat6_bmad_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_high_energy_space_charge

Fortran source: [`bmad/space_charge/high_energy_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/high_energy_space_charge_mod.f90#L232)

::: pybmad.bmad.make_mat6_high_energy_space_charge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_mad

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L56)

::: pybmad.bmad.make_mat6_mad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_symp_lie_ptc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1741)

::: pybmad.bmad.make_mat6_symp_lie_ptc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1715)

::: pybmad.bmad.make_mat6_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_tracking

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1748)

::: pybmad.bmad.make_mat6_tracking
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_n

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L825)

::: pybmad.bmad.make_n
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_pbrh

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/envelope_mod.f90#L654)

::: pybmad.bmad.make_pbrh
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_smat_from_abc

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L1066)

::: pybmad.bmad.make_smat_from_abc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_unit_mad_map

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1888)

::: pybmad.bmad.make_unit_mad_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_v

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/envelope_mod.f90#L275)

::: pybmad.bmad.make_v
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_v_mats

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1758)

::: pybmad.bmad.make_v_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### makeup_control_slave

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bookkeeper_mod.f90#L1557)

::: pybmad.bmad.makeup_control_slave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### makeup_group_lord

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bookkeeper_mod.f90#L22)

::: pybmad.bmad.makeup_group_lord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### makeup_multipass_slave

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bookkeeper_mod.f90#L343)

::: pybmad.bmad.makeup_multipass_slave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### makeup_super_slave

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bookkeeper_mod.f90#L511)

::: pybmad.bmad.makeup_super_slave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### makeup_super_slave1

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bookkeeper_mod.f90#L1172)

::: pybmad.bmad.makeup_super_slave1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### map1_inverse

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1694)

::: pybmad.bmad.map1_inverse
    options:
      show_root_heading: false
      show_root_toc_entry: false

### map1_make_unit

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1700)

::: pybmad.bmad.map1_make_unit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### map1_times_map1

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4188)

::: pybmad.bmad.map1_times_map1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### map_to_angle_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1766)

::: pybmad.bmad.map_to_angle_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mark_patch_regions

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/wall3d_mod.f90#L1538)

::: pybmad.bmad.mark_patch_regions
    options:
      show_root_heading: false
      show_root_toc_entry: false

### master_parameter_value

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1772)

::: pybmad.bmad.master_parameter_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mat4_multipole

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1805)

::: pybmad.bmad.mat4_multipole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mat6_add_offsets

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1780)

::: pybmad.bmad.mat6_add_offsets
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mat6_add_pitch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1787)

::: pybmad.bmad.mat6_add_pitch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mat6_to_complex_taylor

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L772)

::: pybmad.bmad.mat6_to_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mat_symp_decouple

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1794)

::: pybmad.bmad.mat_symp_decouple
    options:
      show_root_heading: false
      show_root_toc_entry: false

### match_ele_to_mat6

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1814)

::: pybmad.bmad.match_ele_to_mat6
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mexp

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1824)

::: pybmad.bmad.mexp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mfft1

Fortran source: [`bmad/space_charge/fast_fourier_am.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/fast_fourier_am.f90#L53)

::: pybmad.bmad.mfft1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### misalign_ptc_fibre

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L2822)

::: pybmad.bmad.misalign_ptc_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### momentum_compaction

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1831)

::: pybmad.bmad.momentum_compaction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mpxx1

Fortran source: [`bmad/multiparticle/ibs_rates_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_rates_mod.f90#L337)

::: pybmad.bmad.mpxx1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mpzt1

Fortran source: [`bmad/multiparticle/ibs_rates_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_rates_mod.f90#L475)

::: pybmad.bmad.mpzt1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multi_coulomb_log

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L694)

::: pybmad.bmad.multi_coulomb_log
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multi_turn_tracking_analysis

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1838)

::: pybmad.bmad.multi_turn_tracking_analysis
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multilayer_type_to_multilayer_params

Fortran source: [`bmad/interface/xraylib_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/xraylib_interface.f90#L152)

::: pybmad.bmad.multilayer_type_to_multilayer_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipass_all_info

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1859)

::: pybmad.bmad.multipass_all_info
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipass_chain

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1866)

::: pybmad.bmad.multipass_chain
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipass_region_info

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L34)

::: pybmad.bmad.multipass_region_info
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole1_ab_to_kt

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1875)

::: pybmad.bmad.multipole1_ab_to_kt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole1_kt_to_ab

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1883)

::: pybmad.bmad.multipole1_kt_to_ab
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_ab_to_kt

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1891)

::: pybmad.bmad.multipole_ab_to_kt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_ele_to_ab

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1898)

::: pybmad.bmad.multipole_ele_to_ab
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_ele_to_kt

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1911)

::: pybmad.bmad.multipole_ele_to_kt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_init

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1928)

::: pybmad.bmad.multipole_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_kick

Fortran source: [`bmad/modules/multipole_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/multipole_mod.f90#L203)

::: pybmad.bmad.multipole_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_kick_mat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1936)

::: pybmad.bmad.multipole_kick_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_kicks

Fortran source: [`bmad/modules/multipole_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/multipole_mod.f90#L32)

::: pybmad.bmad.multipole_kicks
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_kt_to_ab

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1921)

::: pybmad.bmad.multipole_kt_to_ab
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_spin_tracking

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1946)

::: pybmad.bmad.multipole_spin_tracking
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mytan

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L533)

::: pybmad.bmad.mytan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### n_attrib_string_max_len

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L2774)

::: pybmad.bmad.n_attrib_string_max_len
    options:
      show_root_heading: false
      show_root_toc_entry: false

### new_control

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1954)

::: pybmad.bmad.new_control
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nint_chk

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L42)

::: pybmad.bmad.nint_chk
    options:
      show_root_heading: false
      show_root_toc_entry: false

### normal_form_complex_taylors

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L976)

::: pybmad.bmad.normal_form_complex_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### normal_form_taylors

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L911)

::: pybmad.bmad.normal_form_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### normal_mode3_calc

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L130)

::: pybmad.bmad.normal_mode3_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### normal_mode_dispersion

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1962)

::: pybmad.bmad.normal_mode_dispersion
    options:
      show_root_heading: false
      show_root_toc_entry: false

### normalize_evecs

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L1129)

::: pybmad.bmad.normalize_evecs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### num_field_eles

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1968)

::: pybmad.bmad.num_field_eles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### num_lords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1975)

::: pybmad.bmad.num_lords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### odeint_bmad

Fortran source: [`bmad/modules/runge_kutta_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/runge_kutta_mod.f90#L62)

::: pybmad.bmad.odeint_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### odeint_bmad_time

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/time_tracker_mod.f90#L57)

::: pybmad.bmad.odeint_bmad_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### offset_particle

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1982)

::: pybmad.bmad.offset_particle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### offset_photon

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L1995)

::: pybmad.bmad.offset_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### one_turn_mat_at_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2005)

::: pybmad.bmad.one_turn_mat_at_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### open_binary_file

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/binary_parser_mod.f90#L406)

::: pybmad.bmad.open_binary_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### orbit_amplitude_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2014)

::: pybmad.bmad.orbit_amplitude_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### orbit_reference_energy_correction

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2022)

::: pybmad.bmad.orbit_reference_energy_correction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### orbit_to_floor_phase_space

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2031)

::: pybmad.bmad.orbit_to_floor_phase_space
    options:
      show_root_heading: false
      show_root_toc_entry: false

### orbit_to_local_curvilinear

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2039)

::: pybmad.bmad.orbit_to_local_curvilinear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### orbit_too_large

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2048)

::: pybmad.bmad.orbit_too_large
    options:
      show_root_heading: false
      show_root_toc_entry: false

### order_evecs_by_n_similarity

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L568)

::: pybmad.bmad.order_evecs_by_n_similarity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### order_evecs_by_plane_dominance

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L661)

::: pybmad.bmad.order_evecs_by_plane_dominance
    options:
      show_root_heading: false
      show_root_toc_entry: false

### order_evecs_by_tune

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L710)

::: pybmad.bmad.order_evecs_by_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### order_particles_in_z

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/wake_mod.f90#L637)

::: pybmad.bmad.order_particles_in_z
    options:
      show_root_heading: false
      show_root_toc_entry: false

### order_super_lord_slaves

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2057)

::: pybmad.bmad.order_super_lord_slaves
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_alloc_freespace_array

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L151)

::: pybmad.bmad.osc_alloc_freespace_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_alloc_image_array

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L1095)

::: pybmad.bmad.osc_alloc_image_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_alloc_rectpipe_arrays

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L876)

::: pybmad.bmad.osc_alloc_rectpipe_arrays
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_getgrnpipe

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L632)

::: pybmad.bmad.osc_getgrnpipe
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_read_rectpipe_grn

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L825)

::: pybmad.bmad.osc_read_rectpipe_grn
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_write_rectpipe_grn

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L851)

::: pybmad.bmad.osc_write_rectpipe_grn
    options:
      show_root_heading: false
      show_root_toc_entry: false

### p_func

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L1066)

::: pybmad.bmad.p_func
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_cartesian_map

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L6244)

::: pybmad.bmad.parse_cartesian_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_cylindrical_map

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L6443)

::: pybmad.bmad.parse_cylindrical_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_gen_gradients

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L7048)

::: pybmad.bmad.parse_gen_gradients
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_grid_field

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L6666)

::: pybmad.bmad.parse_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_integer_list

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L7286)

::: pybmad.bmad.parse_integer_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_integer_list2

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L7352)

::: pybmad.bmad.parse_integer_list2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_line_or_list

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L4073)

::: pybmad.bmad.parse_line_or_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_real_list

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L7463)

::: pybmad.bmad.parse_real_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_real_list2

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L7629)

::: pybmad.bmad.parse_real_list2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_superimpose_command

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L8156)

::: pybmad.bmad.parse_superimpose_command
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser2_add_superimpose

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L3420)

::: pybmad.bmad.parser2_add_superimpose
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_add_branch

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L5483)

::: pybmad.bmad.parser_add_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_add_constant

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L1284)

::: pybmad.bmad.parser_add_constant
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_add_lords

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L4356)

::: pybmad.bmad.parser_add_lords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_add_superimpose

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L3151)

::: pybmad.bmad.parser_add_superimpose
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_call_check

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L182)

::: pybmad.bmad.parser_call_check
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_debug_print_info

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L6072)

::: pybmad.bmad.parser_debug_print_info
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_error

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L2762)

::: pybmad.bmad.parser_error
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_expand_line

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L5662)

::: pybmad.bmad.parser_expand_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_fast_complex_read

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L8493)

::: pybmad.bmad.parser_fast_complex_read
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_fast_integer_read

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L8426)

::: pybmad.bmad.parser_fast_integer_read
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_fast_real_read

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L8595)

::: pybmad.bmad.parser_fast_real_read
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_file_stack

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L364)

::: pybmad.bmad.parser_file_stack
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_get_integer

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L7741)

::: pybmad.bmad.parser_get_integer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_get_logical

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L7771)

::: pybmad.bmad.parser_get_logical
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_identify_fork_to_element

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L5559)

::: pybmad.bmad.parser_identify_fork_to_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_init_custom_elements

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L8095)

::: pybmad.bmad.parser_init_custom_elements
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_print_line

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L8046)

::: pybmad.bmad.parser_print_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_read_lr_wake

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L1688)

::: pybmad.bmad.parser_read_lr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_read_old_format_lr_wake

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L1826)

::: pybmad.bmad.parser_read_old_format_lr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_read_old_format_sr_wake

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L1957)

::: pybmad.bmad.parser_read_old_format_sr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_read_sr_wake

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L1450)

::: pybmad.bmad.parser_read_sr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_set_attribute

Fortran source: [`bmad/parsing/parser_set_attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/parser_set_attribute_mod.f90#L41)

::: pybmad.bmad.parser_set_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_transfer_control_struct

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L8322)

::: pybmad.bmad.parser_transfer_control_struct
    options:
      show_root_heading: false
      show_root_toc_entry: false

### particle_in_global_frame

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/time_tracker_mod.f90#L857)

::: pybmad.bmad.particle_in_global_frame
    options:
      show_root_heading: false
      show_root_toc_entry: false

### particle_is_moving_backwards

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2073)

::: pybmad.bmad.particle_is_moving_backwards
    options:
      show_root_heading: false
      show_root_toc_entry: false

### particle_is_moving_forward

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2080)

::: pybmad.bmad.particle_is_moving_forward
    options:
      show_root_heading: false
      show_root_toc_entry: false

### particle_rf_time

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2088)

::: pybmad.bmad.particle_rf_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### patch_flips_propagation_direction

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2098)

::: pybmad.bmad.patch_flips_propagation_direction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### patch_length

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2105)

::: pybmad.bmad.patch_length
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_absorption_and_phase_shift

Fortran source: [`bmad/interface/xraylib_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/xraylib_interface.f90#L34)

::: pybmad.bmad.photon_absorption_and_phase_shift
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_add_to_detector_statistics

Fortran source: [`bmad/photon/photon_target_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_target_mod.f90#L239)

::: pybmad.bmad.photon_add_to_detector_statistics
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_diffuse_scattering

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L826)

::: pybmad.bmad.photon_diffuse_scattering
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_hit_func

Fortran source: [`bmad/photon/capillary_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/capillary_mod.f90#L317)

::: pybmad.bmad.photon_hit_func
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_read_spline

Fortran source: [`bmad/photon/photon_init_spline_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_spline_mod.f90#L48)

::: pybmad.bmad.photon_read_spline
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_reflection

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L776)

::: pybmad.bmad.photon_reflection
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_reflection_std_surface_init

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L66)

::: pybmad.bmad.photon_reflection_std_surface_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_reflectivity

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L635)

::: pybmad.bmad.photon_reflectivity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_target_corner_calc

Fortran source: [`bmad/photon/photon_target_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_target_mod.f90#L173)

::: pybmad.bmad.photon_target_corner_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_target_setup

Fortran source: [`bmad/photon/photon_target_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_target_mod.f90#L29)

::: pybmad.bmad.photon_target_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_type

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_utils_mod.f90#L57)

::: pybmad.bmad.photon_type
    options:
      show_root_heading: false
      show_root_toc_entry: false

### physical_ele_end

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2113)

::: pybmad.bmad.physical_ele_end
    options:
      show_root_heading: false
      show_root_toc_entry: false

### point_photon_emission

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L304)

::: pybmad.bmad.point_photon_emission
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_attribute

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2121)

::: pybmad.bmad.pointer_to_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_branch

Fortran sources (overloaded):

- `pointer_to_branch_given_name`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L75)
- `pointer_to_branch_given_ele`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L85)

::: pybmad.bmad.pointer_to_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_ele

Fortran sources (overloaded):

- `pointer_to_ele1`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4049)
- `pointer_to_ele2`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4097)
- `pointer_to_ele3`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4125)
- `pointer_to_ele4`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4158)

::: pybmad.bmad.pointer_to_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_element_at_s

Fortran source: [`bmad/modules/element_at_s_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/element_at_s_mod.f90#L267)

::: pybmad.bmad.pointer_to_element_at_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_fibre

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2141)

::: pybmad.bmad.pointer_to_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_field_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2148)

::: pybmad.bmad.pointer_to_field_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_girder

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2157)

::: pybmad.bmad.pointer_to_girder
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_indexed_attribute

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2165)

::: pybmad.bmad.pointer_to_indexed_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_lord

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2175)

::: pybmad.bmad.pointer_to_lord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_multipass_lord

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2185)

::: pybmad.bmad.pointer_to_multipass_lord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_next_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2194)

::: pybmad.bmad.pointer_to_next_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_slave

Fortran source: [`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_struct.f90#L2743)

::: pybmad.bmad.pointer_to_slave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_super_lord

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2203)

::: pybmad.bmad.pointer_to_super_lord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_surface_displacement_pt

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_utils_mod.f90#L586)

::: pybmad.bmad.pointer_to_surface_displacement_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_surface_segmented_pt

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_utils_mod.f90#L484)

::: pybmad.bmad.pointer_to_surface_segmented_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_wake_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2212)

::: pybmad.bmad.pointer_to_wake_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_wall3d

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/wall3d_mod.f90#L1089)

::: pybmad.bmad.pointer_to_wall3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointers_to_attribute

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2220)

::: pybmad.bmad.pointers_to_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### polar_to_spinor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2234)

::: pybmad.bmad.polar_to_spinor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### polar_to_vec

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2241)

::: pybmad.bmad.polar_to_vec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### print_mesh3d

Fortran source: [`bmad/space_charge/open_spacecharge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_mod.f90#L37)

::: pybmad.bmad.print_mesh3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### prob_x_diffuse

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L1147)

::: pybmad.bmad.prob_x_diffuse
    options:
      show_root_heading: false
      show_root_toc_entry: false

### project_emit_to_xyz

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L1191)

::: pybmad.bmad.project_emit_to_xyz
    options:
      show_root_heading: false
      show_root_toc_entry: false

### propagate_part_way

Fortran source: [`bmad/modules/rad_int_common.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/rad_int_common.f90#L259)

::: pybmad.bmad.propagate_part_way
    options:
      show_root_heading: false
      show_root_toc_entry: false

### psi_prime_sca

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/longitudinal_profile_mod.f90#L86)

::: pybmad.bmad.psi_prime_sca
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_bookkeeper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2248)

::: pybmad.bmad.ptc_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_calculate_tracking_step_size

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L1431)

::: pybmad.bmad.ptc_calculate_tracking_step_size
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_check_for_lost_particle

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L1550)

::: pybmad.bmad.ptc_check_for_lost_particle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_closed_orbit_calc

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L657)

::: pybmad.bmad.ptc_closed_orbit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_emit_calc

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L347)

::: pybmad.bmad.ptc_emit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_kill_map_with_radiation

Fortran source: [`bmad/ptc/ptc_map_with_radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_map_with_radiation_mod.f90#L570)

::: pybmad.bmad.ptc_kill_map_with_radiation
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_layouts_resplit

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L1507)

::: pybmad.bmad.ptc_layouts_resplit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_linear_isf_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2254)

::: pybmad.bmad.ptc_linear_isf_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_one_turn_mat_and_closed_orbit_calc

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L305)

::: pybmad.bmad.ptc_one_turn_mat_and_closed_orbit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_ran_seed_put

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2261)

::: pybmad.bmad.ptc_ran_seed_put
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_read_flat_file

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2266)

::: pybmad.bmad.ptc_read_flat_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_read_map_with_radiation

Fortran source: [`bmad/ptc/ptc_map_with_radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_map_with_radiation_mod.f90#L430)

::: pybmad.bmad.ptc_read_map_with_radiation
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_set_rf_state_for_c_normal

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2275)

::: pybmad.bmad.ptc_set_rf_state_for_c_normal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_set_taylor_order_if_needed

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L47)

::: pybmad.bmad.ptc_set_taylor_order_if_needed
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_setup_map_with_radiation

Fortran source: [`bmad/ptc/ptc_map_with_radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_map_with_radiation_mod.f90#L68)

::: pybmad.bmad.ptc_setup_map_with_radiation
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_spin_calc

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L457)

::: pybmad.bmad.ptc_spin_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_spin_matching_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2280)

::: pybmad.bmad.ptc_spin_matching_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_track_all

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L587)

::: pybmad.bmad.ptc_track_all
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_track_map_with_radiation

Fortran source: [`bmad/ptc/ptc_map_with_radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_map_with_radiation_mod.f90#L253)

::: pybmad.bmad.ptc_track_map_with_radiation
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_transfer_map_with_spin

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2287)

::: pybmad.bmad.ptc_transfer_map_with_spin
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_write_map_with_radiation

Fortran source: [`bmad/ptc/ptc_map_with_radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_map_with_radiation_mod.f90#L320)

::: pybmad.bmad.ptc_write_map_with_radiation
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptwo

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L1259)

::: pybmad.bmad.ptwo
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pwd_mat

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/longitudinal_profile_mod.f90#L679)

::: pybmad.bmad.pwd_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rad1_damp_and_stoc_mats

Fortran source: [`bmad/modules/rad_6d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/rad_6d_mod.f90#L406)

::: pybmad.bmad.rad1_damp_and_stoc_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rad_damp_and_stoc_mats

Fortran source: [`bmad/modules/rad_6d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/rad_6d_mod.f90#L261)

::: pybmad.bmad.rad_damp_and_stoc_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rad_g_integrals

Fortran source: [`bmad/modules/rad_6d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/rad_6d_mod.f90#L744)

::: pybmad.bmad.rad_g_integrals
    options:
      show_root_heading: false
      show_root_toc_entry: false

### radiation_integrals

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2306)

::: pybmad.bmad.radiation_integrals
    options:
      show_root_heading: false
      show_root_toc_entry: false

### radiation_map_setup

Fortran source: [`bmad/modules/radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/radiation_mod.f90#L166)

::: pybmad.bmad.radiation_map_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ramper_slave_setup

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L614)

::: pybmad.bmad.ramper_slave_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ramper_value

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L621)

::: pybmad.bmad.ramper_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### randomize_lr_wake_frequencies

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/wake_mod.f90#L29)

::: pybmad.bmad.randomize_lr_wake_frequencies
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rchomp

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L270)

::: pybmad.bmad.rchomp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### re_allocate

Fortran sources (overloaded):

- `re_allocate_wall3d_vertex_array`: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/wall3d_mod.f90#L36)
- `re_allocate_wall3d_section_array`: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/wall3d_mod.f90#L84)

::: pybmad.bmad.re_allocate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### re_allocate_eles

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2316)

::: pybmad.bmad.re_allocate_eles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### re_associate_node_array

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L849)

::: pybmad.bmad.re_associate_node_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### re_str

Fortran sources (overloaded):

- `re_str_rp`: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L184)
- `re_str_qp`: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L205)

::: pybmad.bmad.re_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_beam_ascii

Fortran source: [`bmad/multiparticle/beam_file_io.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_file_io.f90#L834)

::: pybmad.bmad.read_beam_ascii
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_beam_file

Fortran source: [`bmad/multiparticle/beam_file_io.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_file_io.f90#L380)

::: pybmad.bmad.read_beam_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_binary_cartesian_map

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/binary_parser_mod.f90#L85)

::: pybmad.bmad.read_binary_cartesian_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_binary_cylindrical_map

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/binary_parser_mod.f90#L211)

::: pybmad.bmad.read_binary_cylindrical_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_binary_grid_field

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/binary_parser_mod.f90#L338)

::: pybmad.bmad.read_binary_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_digested_bmad_file

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2360)

::: pybmad.bmad.read_digested_bmad_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_surface_reflection_file

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_reflection_mod.f90#L412)

::: pybmad.bmad.read_surface_reflection_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_beam

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2370)

::: pybmad.bmad.reallocate_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_bp_com_const

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L3126)

::: pybmad.bmad.reallocate_bp_com_const
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_bunch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2379)

::: pybmad.bmad.reallocate_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_control

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2324)

::: pybmad.bmad.reallocate_control
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_coord

Fortran sources (overloaded):

- `reallocate_coord_n`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L168)
- `reallocate_coord_lat`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L175)
- `reallocate_coord_array`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L183)

::: pybmad.bmad.reallocate_coord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_expression_stack

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2331)

::: pybmad.bmad.reallocate_expression_stack
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_sequence

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L8133)

::: pybmad.bmad.reallocate_sequence
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rel_tracking_charge_to_mass

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2339)

::: pybmad.bmad.rel_tracking_charge_to_mass
    options:
      show_root_heading: false
      show_root_toc_entry: false

### relative_mode_flip

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2387)

::: pybmad.bmad.relative_mode_flip
    options:
      show_root_heading: false
      show_root_toc_entry: false

### release_rad_int_cache

Fortran source: [`bmad/modules/radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/radiation_mod.f90#L23)

::: pybmad.bmad.release_rad_int_cache
    options:
      show_root_heading: false
      show_root_toc_entry: false

### remove_constant_taylor

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L1915)

::: pybmad.bmad.remove_constant_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### remove_dead_from_bunch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L630)

::: pybmad.bmad.remove_dead_from_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### remove_eles_from_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2347)

::: pybmad.bmad.remove_eles_from_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### remove_lord_slave_link

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2354)

::: pybmad.bmad.remove_lord_slave_link
    options:
      show_root_heading: false
      show_root_toc_entry: false

### residual_pwd_sig_z

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/ibs_mod.f90#L914)

::: pybmad.bmad.residual_pwd_sig_z
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reverse_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2395)

::: pybmad.bmad.reverse_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rf_coupler_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2411)

::: pybmad.bmad.rf_coupler_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rf_is_on

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2423)

::: pybmad.bmad.rf_is_on
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rf_ref_time_offset

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2431)

::: pybmad.bmad.rf_ref_time_offset
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L723)

::: pybmad.bmad.rfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rk_adaptive_time_step

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/time_tracker_mod.f90#L434)

::: pybmad.bmad.rk_adaptive_time_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rk_time_step1

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/time_tracker_mod.f90#L560)

::: pybmad.bmad.rk_time_step1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate3

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/astra_interface_mod.f90#L362)

::: pybmad.bmad.rotate3
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_em_field

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/em_field_mod.f90#L173)

::: pybmad.bmad.rotate_em_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_field_zx

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L1588)

::: pybmad.bmad.rotate_field_zx
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_for_curved_surface

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2439)

::: pybmad.bmad.rotate_for_curved_surface
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_spin

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2448)

::: pybmad.bmad.rotate_spin
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_spin_a_step

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2455)

::: pybmad.bmad.rotate_spin_a_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_spin_given_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2464)

::: pybmad.bmad.rotate_spin_given_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### s_body_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2472)

::: pybmad.bmad.s_body_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### s_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2480)

::: pybmad.bmad.s_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### s_ref_to_s_chord

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L1812)

::: pybmad.bmad.s_ref_to_s_chord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### s_source_calc

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L848)

::: pybmad.bmad.s_source_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sad_mult_hard_bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/fringe_mod.f90#L447)

::: pybmad.bmad.sad_mult_hard_bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sad_soft_bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/fringe_mod.f90#L892)

::: pybmad.bmad.sad_soft_bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### save_a_beam_step

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2486)

::: pybmad.bmad.save_a_beam_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### save_a_bunch_step

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2496)

::: pybmad.bmad.save_a_bunch_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### save_a_step

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2506)

::: pybmad.bmad.save_a_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sbend_body_with_k1_map

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2520)

::: pybmad.bmad.sbend_body_with_k1_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sc_adaptive_step

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/space_charge_mod.f90#L297)

::: pybmad.bmad.sc_adaptive_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sc_step

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/space_charge_mod.f90#L227)

::: pybmad.bmad.sc_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_active_fixer

Fortran source: [`bmad/modules/fixer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/fixer_mod.f90#L32)

::: pybmad.bmad.set_active_fixer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_branch_and_ele_for_omp

Fortran source: [`bmad/modules/dynamic_aperture_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/dynamic_aperture_mod.f90#L232)

::: pybmad.bmad.set_branch_and_ele_for_omp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_custom_attribute_name

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L2863)

::: pybmad.bmad.set_custom_attribute_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_attribute

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2532)

::: pybmad.bmad.set_ele_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_defaults

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2644)

::: pybmad.bmad.set_ele_defaults
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_misalignments

Not exposed at the top level — import as `pybmad.extra.set_ele_misalignments`.

Fortran source: `ele_misalignments.f90:68`

::: pybmad.extra.set_ele_misalignments
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_name

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2542)

::: pybmad.bmad.set_ele_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_real_attribute

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2549)

::: pybmad.bmad.set_ele_real_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_status_stale

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2559)

::: pybmad.bmad.set_ele_status_stale
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_flags_for_changed_attribute

Fortran sources (overloaded):

- `set_flags_for_changed_all_attribute`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L327)
- `set_flags_for_changed_integer_attribute`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L335)
- `set_flags_for_changed_logical_attribute`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L343)
- `set_flags_for_changed_lat_attribute`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L351)
- `set_flags_for_changed_real_attribute`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L358)

::: pybmad.bmad.set_flags_for_changed_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_fringe_on_off

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2568)

::: pybmad.bmad.set_fringe_on_off
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_lords_status_stale

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2575)

::: pybmad.bmad.set_lords_status_stale
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_on_off

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2584)

::: pybmad.bmad.set_on_off
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_orbit_to_zero

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2596)

::: pybmad.bmad.set_orbit_to_zero
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ptc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2604)

::: pybmad.bmad.set_ptc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ptc_base_state

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2612)

::: pybmad.bmad.set_ptc_base_state
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ptc_com_pointers

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L1005)

::: pybmad.bmad.set_ptc_com_pointers
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ptc_quiet

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L3365)

::: pybmad.bmad.set_ptc_quiet
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ptc_verbose

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L1176)

::: pybmad.bmad.set_ptc_verbose
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_pwd_ele

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/longitudinal_profile_mod.f90#L547)

::: pybmad.bmad.set_pwd_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_status_flags

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2619)

::: pybmad.bmad.set_status_flags
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_tune

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2651)

::: pybmad.bmad.set_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_tune_via_group_knobs

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2664)

::: pybmad.bmad.set_tune_via_group_knobs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2626)

::: pybmad.bmad.set_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_z_tune

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2636)

::: pybmad.bmad.set_z_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### settable_dep_var_bookkeeping

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L4993)

::: pybmad.bmad.settable_dep_var_bookkeeping
    options:
      show_root_heading: false
      show_root_toc_entry: false

### setup_high_energy_space_charge_calc

Fortran source: [`bmad/space_charge/high_energy_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/high_energy_space_charge_mod.f90#L32)

::: pybmad.bmad.setup_high_energy_space_charge_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sigma_mat_ptc_to_bmad

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L1327)

::: pybmad.bmad.sigma_mat_ptc_to_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### significant_difference

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2675)

::: pybmad.bmad.significant_difference
    options:
      show_root_heading: false
      show_root_toc_entry: false

### skip_ele_blender

Fortran source: [`bmad/interface/blender_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/blender_interface_mod.f90#L87)

::: pybmad.bmad.skip_ele_blender
    options:
      show_root_heading: false
      show_root_toc_entry: false

### slice_lattice

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2683)

::: pybmad.bmad.slice_lattice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### soft_quadrupole_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/fringe_mod.f90#L555)

::: pybmad.bmad.soft_quadrupole_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sol_quad_mat6_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2692)

::: pybmad.bmad.sol_quad_mat6_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### solve_psi_adaptive

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/longitudinal_profile_mod.f90#L168)

::: pybmad.bmad.solve_psi_adaptive
    options:
      show_root_heading: false
      show_root_toc_entry: false

### solve_psi_fixed_steps

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/longitudinal_profile_mod.f90#L227)

::: pybmad.bmad.solve_psi_fixed_steps
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sort_complex_taylor_terms

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L604)

::: pybmad.bmad.sort_complex_taylor_terms
    options:
      show_root_heading: false
      show_root_toc_entry: false

### space_charge_cathodeimages

Fortran source: [`bmad/space_charge/open_spacecharge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_mod.f90#L84)

::: pybmad.bmad.space_charge_cathodeimages
    options:
      show_root_heading: false
      show_root_toc_entry: false

### space_charge_freespace

Fortran source: [`bmad/space_charge/open_spacecharge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_mod.f90#L61)

::: pybmad.bmad.space_charge_freespace
    options:
      show_root_heading: false
      show_root_toc_entry: false

### space_charge_rectpipe

Fortran source: [`bmad/space_charge/open_spacecharge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_mod.f90#L111)

::: pybmad.bmad.space_charge_rectpipe
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_concat_linear_maps

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2713)

::: pybmad.bmad.spin_concat_linear_maps
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_depolarization_rate

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2725)

::: pybmad.bmad.spin_depolarization_rate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_dn_dpz_from_mat8

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2733)

::: pybmad.bmad.spin_dn_dpz_from_mat8
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_dn_dpz_from_qmap

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2741)

::: pybmad.bmad.spin_dn_dpz_from_qmap
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_map1_normalize

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2750)

::: pybmad.bmad.spin_map1_normalize
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_mat8_resonance_strengths

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2764)

::: pybmad.bmad.spin_mat8_resonance_strengths
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_mat_to_eigen

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2756)

::: pybmad.bmad.spin_mat_to_eigen
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_omega

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2771)

::: pybmad.bmad.spin_omega
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_quat_resonance_strengths

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2781)

::: pybmad.bmad.spin_quat_resonance_strengths
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_taylor_to_linear

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2788)

::: pybmad.bmad.spin_taylor_to_linear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spinor_to_polar

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2796)

::: pybmad.bmad.spinor_to_polar
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spinor_to_vec

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2803)

::: pybmad.bmad.spinor_to_vec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spline_fit_orbit

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2810)

::: pybmad.bmad.spline_fit_orbit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### split_expression_string

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L1920)

::: pybmad.bmad.split_expression_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### split_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2817)

::: pybmad.bmad.split_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sprint_spin_taylor_map

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2830)

::: pybmad.bmad.sprint_spin_taylor_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sr_longitudinal_wake_particle

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/wake_mod.f90#L309)

::: pybmad.bmad.sr_longitudinal_wake_particle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sr_transverse_wake_particle

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/wake_mod.f90#L405)

::: pybmad.bmad.sr_transverse_wake_particle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sr_z_long_wake

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/wake_mod.f90#L512)

::: pybmad.bmad.sr_z_long_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### srdt_calc

Fortran source: [`bmad/modules/srdt_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/srdt_mod.f90#L81)

::: pybmad.bmad.srdt_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### srdt_lsq_solution

Fortran source: [`bmad/modules/srdt_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/srdt_mod.f90#L656)

::: pybmad.bmad.srdt_lsq_solution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### start_branch_at

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2837)

::: pybmad.bmad.start_branch_at
    options:
      show_root_heading: false
      show_root_toc_entry: false

### stream_ele_end

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2845)

::: pybmad.bmad.stream_ele_end
    options:
      show_root_heading: false
      show_root_toc_entry: false

### string_attrib

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L2350)

::: pybmad.bmad.string_attrib
    options:
      show_root_heading: false
      show_root_toc_entry: false

### strong_beam_sigma_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2851)

::: pybmad.bmad.strong_beam_sigma_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### strong_beam_strength

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2858)

::: pybmad.bmad.strong_beam_strength
    options:
      show_root_heading: false
      show_root_toc_entry: false

### surface_grid_displacement

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_utils_mod.f90#L186)

::: pybmad.bmad.surface_grid_displacement
    options:
      show_root_heading: false
      show_root_toc_entry: false

### switch_attrib_value_name

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/attribute_mod.f90#L2431)

::: pybmad.bmad.switch_attrib_value_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### symp_lie_bmad

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2865)

::: pybmad.bmad.symp_lie_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### t6_to_b123

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L49)

::: pybmad.bmad.t6_to_b123
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taper_mag_strengths

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2876)

::: pybmad.bmad.taper_mag_strengths
    options:
      show_root_heading: false
      show_root_toc_entry: false

### target_min_max_calc

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L1190)

::: pybmad.bmad.target_min_max_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### target_rot_mats

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L1141)

::: pybmad.bmad.target_rot_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taylor_equal_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4921)

::: pybmad.bmad.taylor_equal_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taylor_inverse

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L1981)

::: pybmad.bmad.taylor_inverse
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taylor_propagate1

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L2490)

::: pybmad.bmad.taylor_propagate1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taylor_to_mad_map

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1785)

::: pybmad.bmad.taylor_to_mad_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taylors_equal_taylors

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L4958)

::: pybmad.bmad.taylors_equal_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tilt_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2885)

::: pybmad.bmad.tilt_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tilt_coords_photon

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2894)

::: pybmad.bmad.tilt_coords_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tilt_mat6

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2901)

::: pybmad.bmad.tilt_mat6
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_eta_reading

Fortran source: [`bmad/modules/measurement_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/measurement_mod.f90#L209)

::: pybmad.bmad.to_eta_reading
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_fieldmap_coords

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/em_field_mod.f90#L101)

::: pybmad.bmad.to_fieldmap_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_orbit_reading

Fortran source: [`bmad/modules/measurement_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/measurement_mod.f90#L127)

::: pybmad.bmad.to_orbit_reading
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_phase_and_coupling_reading

Fortran source: [`bmad/modules/measurement_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/measurement_mod.f90#L288)

::: pybmad.bmad.to_phase_and_coupling_reading
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_photon_angle_coords

Fortran source: [`bmad/photon/photon_target_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_target_mod.f90#L355)

::: pybmad.bmad.to_photon_angle_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_surface_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2907)

::: pybmad.bmad.to_surface_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### touschek_lifetime

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/touschek_mod.f90#L82)

::: pybmad.bmad.touschek_lifetime
    options:
      show_root_heading: false
      show_root_toc_entry: false

### touschek_lifetime_ele_by_ele

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/touschek_mod.f90#L148)

::: pybmad.bmad.touschek_lifetime_ele_by_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### touschek_lifetime_with_aperture

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/touschek_mod.f90#L226)

::: pybmad.bmad.touschek_lifetime_with_aperture
    options:
      show_root_heading: false
      show_root_toc_entry: false

### touschek_rate1

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/touschek_mod.f90#L427)

::: pybmad.bmad.touschek_rate1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### touschek_rate1_zap

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/touschek_mod.f90#L269)

::: pybmad.bmad.touschek_rate1_zap
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3172)

::: pybmad.bmad.track1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_beam

Fortran source: [`bmad/multiparticle/beam_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_mod.f90#L193)

::: pybmad.bmad.track1_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bmad

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3185)

::: pybmad.bmad.track1_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bmad_photon

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3197)

::: pybmad.bmad.track1_bmad_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bunch

Fortran source: [`bmad/multiparticle/beam_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_mod.f90#L238)

::: pybmad.bmad.track1_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bunch_csr

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L137)

::: pybmad.bmad.track1_bunch_csr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bunch_csr3d

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/csr_and_space_charge_mod.f90#L1881)

::: pybmad.bmad.track1_bunch_csr3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bunch_hom

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_utils.f90#L32)

::: pybmad.bmad.track1_bunch_hom
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bunch_space_charge

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3206)

::: pybmad.bmad.track1_bunch_space_charge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_crystal

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L845)

::: pybmad.bmad.track1_crystal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_diffraction_plate_or_mask

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L145)

::: pybmad.bmad.track1_diffraction_plate_or_mask
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_high_energy_space_charge

Fortran source: [`bmad/space_charge/high_energy_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/high_energy_space_charge_mod.f90#L169)

::: pybmad.bmad.track1_high_energy_space_charge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_lens

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L27)

::: pybmad.bmad.track1_lens
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_linear

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3216)

::: pybmad.bmad.track1_linear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_lr_wake

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/wake_mod.f90#L110)

::: pybmad.bmad.track1_lr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_mad

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mad_mod.f90#L1597)

::: pybmad.bmad.track1_mad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_mirror

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L435)

::: pybmad.bmad.track1_mirror
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_mosaic_crystal

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L618)

::: pybmad.bmad.track1_mosaic_crystal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_multilayer_mirror

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L482)

::: pybmad.bmad.track1_multilayer_mirror
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_radiation

Fortran source: [`bmad/modules/radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/radiation_mod.f90#L61)

::: pybmad.bmad.track1_radiation
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_radiation_center

Fortran source: [`bmad/modules/radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/radiation_mod.f90#L275)

::: pybmad.bmad.track1_radiation_center
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_runge_kutta

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3224)

::: pybmad.bmad.track1_runge_kutta
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_sample

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L235)

::: pybmad.bmad.track1_sample
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_spin

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3236)

::: pybmad.bmad.track1_spin
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_spin_integration

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3245)

::: pybmad.bmad.track1_spin_integration
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_spin_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3254)

::: pybmad.bmad.track1_spin_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_sr_wake

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/wake_mod.f90#L702)

::: pybmad.bmad.track1_sr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_symp_lie_ptc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3262)

::: pybmad.bmad.track1_symp_lie_ptc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3271)

::: pybmad.bmad.track1_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_time_runge_kutta

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3281)

::: pybmad.bmad.track1_time_runge_kutta
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_beambeam

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2914)

::: pybmad.bmad.track_a_beambeam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_bend

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2925)

::: pybmad.bmad.track_a_bend
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_bend_photon

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L1246)

::: pybmad.bmad.track_a_bend_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_capillary

Fortran source: [`bmad/photon/capillary_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/capillary_mod.f90#L39)

::: pybmad.bmad.track_a_capillary
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_converter

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2935)

::: pybmad.bmad.track_a_converter
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_crab_cavity

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2945)

::: pybmad.bmad.track_a_crab_cavity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_drift

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2955)

::: pybmad.bmad.track_a_drift
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_drift_photon

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2965)

::: pybmad.bmad.track_a_drift_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_foil

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3082)

::: pybmad.bmad.track_a_foil
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_gkicker

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2973)

::: pybmad.bmad.track_a_gkicker
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_lcavity

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2983)

::: pybmad.bmad.track_a_lcavity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_lcavity_old

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L2993)

::: pybmad.bmad.track_a_lcavity_old
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_mask

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3003)

::: pybmad.bmad.track_a_mask
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_match

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3013)

::: pybmad.bmad.track_a_match
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_patch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3033)

::: pybmad.bmad.track_a_patch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_patch_photon

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/track1_photon_mod.f90#L65)

::: pybmad.bmad.track_a_patch_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_pickup

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3023)

::: pybmad.bmad.track_a_pickup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_quadrupole

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3042)

::: pybmad.bmad.track_a_quadrupole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_rfcavity

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3052)

::: pybmad.bmad.track_a_rfcavity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_sad_mult

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3062)

::: pybmad.bmad.track_a_sad_mult
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_sol_quad

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3072)

::: pybmad.bmad.track_a_sol_quad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_thick_multipole

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3092)

::: pybmad.bmad.track_a_thick_multipole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_wiggler

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3102)

::: pybmad.bmad.track_a_wiggler
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_zero_length_element

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3112)

::: pybmad.bmad.track_a_zero_length_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_all

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3122)

::: pybmad.bmad.track_all
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_beam

Fortran source: [`bmad/multiparticle/beam_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_mod.f90#L40)

::: pybmad.bmad.track_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_bunch

Fortran source: [`bmad/multiparticle/beam_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_mod.f90#L103)

::: pybmad.bmad.track_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_bunch_time

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3132)

::: pybmad.bmad.track_bunch_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_bunch_to_s

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/space_charge_mod.f90#L436)

::: pybmad.bmad.track_bunch_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_bunch_to_t

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/space_charge_mod.f90#L504)

::: pybmad.bmad.track_bunch_to_t
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_complex_taylor

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L828)

::: pybmad.bmad.track_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_from_s_to_s

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3142)

::: pybmad.bmad.track_from_s_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_func

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/space_charge_mod.f90#L586)

::: pybmad.bmad.track_func
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_many

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3152)

::: pybmad.bmad.track_many
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_to_surface

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3163)

::: pybmad.bmad.track_to_surface
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_until_dead

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/time_tracker_mod.f90#L1183)

::: pybmad.bmad.track_until_dead
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tracking_rad_map_setup

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3292)

::: pybmad.bmad.tracking_rad_map_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_ac_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3302)

::: pybmad.bmad.transfer_ac_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_branch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3308)

::: pybmad.bmad.transfer_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_branch_parameters

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3315)

::: pybmad.bmad.transfer_branch_parameters
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_branches

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3322)

::: pybmad.bmad.transfer_branches
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3329)

::: pybmad.bmad.transfer_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_ele_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3337)

::: pybmad.bmad.transfer_ele_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_eles

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3344)

::: pybmad.bmad.transfer_eles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_fieldmap

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3351)

::: pybmad.bmad.transfer_fieldmap
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_fixer_params

Fortran source: [`bmad/modules/fixer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/fixer_mod.f90#L113)

::: pybmad.bmad.transfer_fixer_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3358)

::: pybmad.bmad.transfer_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_lat_parameters

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3365)

::: pybmad.bmad.transfer_lat_parameters
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_map_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3372)

::: pybmad.bmad.transfer_map_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_map_from_s_to_s

Fortran source: [`bmad/modules/transfer_map_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/transfer_map_mod.f90#L59)

::: pybmad.bmad.transfer_map_from_s_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_mat2_from_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3393)

::: pybmad.bmad.transfer_mat2_from_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_mat_from_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3385)

::: pybmad.bmad.transfer_mat_from_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_matrix_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3400)

::: pybmad.bmad.transfer_matrix_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3410)

::: pybmad.bmad.transfer_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_wake

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3417)

::: pybmad.bmad.transfer_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### truncate_complex_taylor_to_order

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L900)

::: pybmad.bmad.truncate_complex_taylor_to_order
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss1_propagate

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3472)

::: pybmad.bmad.twiss1_propagate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss3_at_start

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L1384)

::: pybmad.bmad.twiss3_at_start
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss3_from_twiss2

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L1332)

::: pybmad.bmad.twiss3_from_twiss2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss3_propagate1

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L1265)

::: pybmad.bmad.twiss3_propagate1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss3_propagate_all

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/mode3_mod.f90#L1236)

::: pybmad.bmad.twiss3_propagate_all
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_and_track

Fortran sources (overloaded):

- `twiss_and_track_branch`: [`bmad/modules/twiss_and_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/twiss_and_track_mod.f90#L88)
- `twiss_and_track_all`: [`bmad/modules/twiss_and_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/twiss_and_track_mod.f90#L133)

::: pybmad.bmad.twiss_and_track
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_and_track_at_s

Fortran source: [`bmad/modules/twiss_and_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/twiss_and_track_mod.f90#L347)

::: pybmad.bmad.twiss_and_track_at_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_and_track_from_s_to_s

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3429)

::: pybmad.bmad.twiss_and_track_from_s_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_and_track_intra_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3441)

::: pybmad.bmad.twiss_and_track_intra_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_at_element

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3454)

::: pybmad.bmad.twiss_at_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_at_start

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3463)

::: pybmad.bmad.twiss_at_start
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_from_tracking

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3500)

::: pybmad.bmad.twiss_from_tracking
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_propagate1

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3510)

::: pybmad.bmad.twiss_propagate1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_propagate_all

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3517)

::: pybmad.bmad.twiss_propagate_all
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_to_1_turn_mat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3525)

::: pybmad.bmad.twiss_to_1_turn_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_complex_taylors

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/complex_taylor_mod.f90#L281)

::: pybmad.bmad.type_complex_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_coord

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3532)

::: pybmad.bmad.type_coord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3538)

::: pybmad.bmad.type_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_end_stuff

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L909)

::: pybmad.bmad.type_end_stuff
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_expression_tree

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/expression_mod.f90#L695)

::: pybmad.bmad.type_expression_tree
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_ptc_fibre

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_interface_mod.f90#L385)

::: pybmad.bmad.type_ptc_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_ptc_layout

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L28)

::: pybmad.bmad.type_ptc_layout
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_taylors

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3551)

::: pybmad.bmad.type_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3562)

::: pybmad.bmad.type_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### update_ele_from_fibre

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/ptc/ptc_layout_mod.f90#L1204)

::: pybmad.bmad.update_ele_from_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### update_fibre_from_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3594)

::: pybmad.bmad.update_fibre_from_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### update_floor_angles

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3587)

::: pybmad.bmad.update_floor_angles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### valid_field_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3601)

::: pybmad.bmad.valid_field_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### valid_fringe_type

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3609)

::: pybmad.bmad.valid_fringe_type
    options:
      show_root_heading: false
      show_root_toc_entry: false

### valid_mat6_calc_method

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3617)

::: pybmad.bmad.valid_mat6_calc_method
    options:
      show_root_heading: false
      show_root_toc_entry: false

### valid_spin_tracking_method

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3625)

::: pybmad.bmad.valid_spin_tracking_method
    options:
      show_root_heading: false
      show_root_toc_entry: false

### valid_tracking_method

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3633)

::: pybmad.bmad.valid_tracking_method
    options:
      show_root_heading: false
      show_root_toc_entry: false

### value_of_attribute

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3641)

::: pybmad.bmad.value_of_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### value_to_line

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L434)

::: pybmad.bmad.value_to_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### vec_to_polar

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3652)

::: pybmad.bmad.vec_to_polar
    options:
      show_root_heading: false
      show_root_toc_entry: false

### vec_to_spinor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3660)

::: pybmad.bmad.vec_to_spinor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### verify_valid_name

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L2596)

::: pybmad.bmad.verify_valid_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### vert_angle_func

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_init_mod.f90#L352)

::: pybmad.bmad.vert_angle_func
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_for_bend_angle

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3668)

::: pybmad.bmad.w_mat_for_bend_angle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_for_tilt

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3691)

::: pybmad.bmad.w_mat_for_tilt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_for_x_pitch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3675)

::: pybmad.bmad.w_mat_for_x_pitch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_for_y_pitch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3683)

::: pybmad.bmad.w_mat_for_y_pitch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### wall3d_d_radius

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/wall3d_mod.f90#L659)

::: pybmad.bmad.wall3d_d_radius
    options:
      show_root_heading: false
      show_root_toc_entry: false

### wall3d_initializer

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/wall3d_mod.f90#L135)

::: pybmad.bmad.wall3d_initializer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### wall3d_section_initializer

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/wall3d_mod.f90#L207)

::: pybmad.bmad.wall3d_section_initializer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### wall3d_to_position

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/wall3d_mod.f90#L1153)

::: pybmad.bmad.wall3d_to_position
    options:
      show_root_heading: false
      show_root_toc_entry: false

### word_to_value

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/bmad_parser_mod.f90#L1076)

::: pybmad.bmad.word_to_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_ascii_beam_file

Fortran source: [`bmad/multiparticle/beam_file_io.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_file_io.f90#L123)

::: pybmad.bmad.write_ascii_beam_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_astra_bend

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/astra_interface_mod.f90#L344)

::: pybmad.bmad.write_astra_bend
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_astra_ele

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/astra_interface_mod.f90#L146)

::: pybmad.bmad.write_astra_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_astra_field_grid_file

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/astra_interface_mod.f90#L504)

::: pybmad.bmad.write_astra_field_grid_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_astra_field_grid_file_3d

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/astra_interface_mod.f90#L718)

::: pybmad.bmad.write_astra_field_grid_file_3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_astra_lattice_file

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/astra_interface_mod.f90#L29)

::: pybmad.bmad.write_astra_lattice_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_beam_file

Fortran source: [`bmad/multiparticle/beam_file_io.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/beam_file_io.f90#L29)

::: pybmad.bmad.write_beam_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_beam_floor_positions

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3699)

::: pybmad.bmad.write_beam_floor_positions
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_binary_cartesian_map

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/binary_parser_mod.f90#L25)

::: pybmad.bmad.write_binary_cartesian_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_binary_cylindrical_map

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/binary_parser_mod.f90#L149)

::: pybmad.bmad.write_binary_cylindrical_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_binary_grid_field

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/parsing/binary_parser_mod.f90#L275)

::: pybmad.bmad.write_binary_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_blender_ele

Fortran source: [`bmad/interface/blender_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/blender_interface_mod.f90#L112)

::: pybmad.bmad.write_blender_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_blender_lat_layout

Fortran source: [`bmad/interface/blender_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/blender_interface_mod.f90#L15)

::: pybmad.bmad.write_blender_lat_layout
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_bmad_lattice_file

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3719)

::: pybmad.bmad.write_bmad_lattice_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_digested_bmad_file

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3708)

::: pybmad.bmad.write_digested_bmad_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_gpt_ele

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L338)

::: pybmad.bmad.write_gpt_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_gpt_field_grid_file_1d

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L797)

::: pybmad.bmad.write_gpt_field_grid_file_1d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_gpt_field_grid_file_2d

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L1008)

::: pybmad.bmad.write_gpt_field_grid_file_2d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_gpt_field_grid_file_3d

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L1264)

::: pybmad.bmad.write_gpt_field_grid_file_3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_gpt_lattice_file

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/gpt_interface_mod.f90#L163)

::: pybmad.bmad.write_gpt_lattice_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lat_line

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L331)

::: pybmad.bmad.write_lat_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_elegant_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3729)

::: pybmad.bmad.write_lattice_elegant_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_foreign_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3741)

::: pybmad.bmad.write_lattice_foreign_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_mad_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3753)

::: pybmad.bmad.write_lattice_mad_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_pals_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3765)

::: pybmad.bmad.write_lattice_pals_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_sad_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3773)

::: pybmad.bmad.write_lattice_sad_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_scibmad_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3782)

::: pybmad.bmad.write_lattice_scibmad_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_line_element

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/output/write_lattice_file_mod.f90#L142)

::: pybmad.bmad.write_line_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_opal_field_grid_file

Fortran source: [`bmad/interface/opal_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/opal_interface_mod.f90#L419)

::: pybmad.bmad.write_opal_field_grid_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_opal_lattice_file

Fortran source: [`bmad/interface/opal_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/opal_interface_mod.f90#L26)

::: pybmad.bmad.write_opal_lattice_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_time_particle_distribution

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/time_tracker_mod.f90#L1005)

::: pybmad.bmad.write_time_particle_distribution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### xlafun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L377)

::: pybmad.bmad.xlafun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### xraylib_nist_compound

Fortran source: [`bmad/interface/xraylib_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/interface/xraylib_interface.f90#L457)

::: pybmad.bmad.xraylib_nist_compound
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ylafun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L397)

::: pybmad.bmad.ylafun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### z_at_surface

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/photon/photon_utils_mod.f90#L96)

::: pybmad.bmad.z_at_surface
    options:
      show_root_heading: false
      show_root_toc_entry: false

### zero_ele_kicks

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3801)

::: pybmad.bmad.zero_ele_kicks
    options:
      show_root_heading: false
      show_root_toc_entry: false

### zero_ele_offsets

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/modules/bmad_routine_interface.f90#L3807)

::: pybmad.bmad.zero_ele_offsets
    options:
      show_root_heading: false
      show_root_toc_entry: false

### zero_lr_wakes_in_lat

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/multiparticle/wake_mod.f90#L75)

::: pybmad.bmad.zero_lr_wakes_in_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### zlafun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/b0e0f0853fdde1ed1efb6c17caea95ad55521d76/bmad/space_charge/open_spacecharge_core_mod.f90#L412)

::: pybmad.bmad.zlafun
    options:
      show_root_heading: false
      show_root_toc_entry: false
