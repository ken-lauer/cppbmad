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

Fortran struct: `ac_kicker_freq_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L693))

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

Fortran struct: `ac_kicker_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L699))

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

Fortran struct: `ac_kicker_time_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L687))

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

Fortran struct: `anormal_mode_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1986))

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

Fortran struct: `aperture_param_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2138))

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

Fortran struct: `aperture_point_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2129))

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

Fortran struct: `aperture_scan_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2152))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `point` | [1D array of AperturePointStruct](bmad.md#aperturepointstruct) | Set of aperture points at different angles. |
| `ref_orb` | [CoordStruct](bmad.md#coordstruct) | Ref orbit around which the scan is made. |
| `pz_start` | float | Starting pz. |

::: pybmad.BeamInitStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BeamInitStruct

Fortran struct: `beam_init_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1148))

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

Fortran struct: `beam_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1123))

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

Fortran struct: `bmad_common_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2296))

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

Fortran struct: `bmad_normal_form_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1575))

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

Fortran struct: `bookkeeping_state_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L944))

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

Fortran struct: `bpm_phase_coupling_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L584))

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

::: pybmad.BranchStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BranchStruct

Fortran struct: `branch_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1604))

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

::: pybmad.BunchParamsStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### BunchParamsStruct

Fortran struct: `bunch_params_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1198))

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

Fortran struct: `bunch_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1103))

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

Fortran struct: `bunch_track_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1224))

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

Fortran struct: `cartesian_map_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L727))

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

Fortran struct: `cartesian_map_term1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L713))

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

Fortran struct: `cartesian_map_term_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L721))

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

Fortran struct: `complex_taylor_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L504))

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

Fortran struct: `complex_taylor_term_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L486))

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

Fortran struct: `control_ramp1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1366))

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

Fortran struct: `control_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1348))

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

Fortran struct: `control_var1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1360))

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

Fortran struct: `controller_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1383))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `var` | [1D array of ControlVar1Struct](bmad.md#controlvar1struct) |  |
| `ramp` | [1D array of ControlRamp1Struct](bmad.md#controlramp1struct) | For ramper lord elements |
| `ramper_lord` | [1D array of RamperLordStruct](bmad.md#ramperlordstruct) | Ramper lord info for this slave |
| `x_knot` | 1D array of float |  |

::: pybmad.CoordArrayStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CoordArrayStruct

Fortran struct: `coord_array_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L576))

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

Fortran struct: `coord_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L546))

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

::: pybmad.CylindricalMapStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### CylindricalMapStruct

Fortran struct: `cylindrical_map_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L749))

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

Fortran struct: `cylindrical_map_term1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L738))

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

Fortran struct: `cylindrical_map_term_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L743))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `file` | str | Input file name. Used also as ID for instances. |
| `n_link` | int | For memory management of this structure |
| `term` | [1D array of CylindricalMapTerm1Struct](bmad.md#cylindricalmapterm1struct) |  |

::: pybmad.ElePointerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### ElePointerStruct

Fortran struct: `ele_pointer_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L909))

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

Fortran struct: `ele_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1408))

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
| `rf` | [RfEleStruct](bmad.md#rfelestruct) | RF parameters. |
| `lord` | [EleStruct](bmad.md#elestruct) | Pointer to a slice lord. |
| `ptc_fibre` | [Fibre](forest.md#fibre) | PTC track corresponding to this ele. %floor is floor coord of lab coordinates at the downstream end. Notice that if ele%direction = -1, the lab coords have the z-axis antiparallel to the +s-direction. |
| `floor` | [FloorPositionStruct](bmad.md#floorpositionstruct) |  |
| `high_energy_space_charge` | [HighEnergySpaceChargeStruct](bmad.md#highenergyspacechargestruct) |  |
| `mode3` | [Mode3Struct](bmad.md#mode3struct) | 6D normal mode structure. |
| `photon` | [PhotonElementStruct](bmad.md#photonelementstruct) |  |
| `rad_map` | [RadMapEleStruct](bmad.md#radmapelestruct) | Radiation kick parameters Note: The reference orbits for spin and orbit Taylor maps are not necessarily the same. For example, Sprint spin Taylor maps can be with respect to the zero orbit independent of the orbital map. |
| `taylor` | [1D array of TaylorStruct (shape: 6)](bmad.md#taylorstruct) | Phase space Taylor map. |
| `spin_taylor_ref_orb_in` | 1D array of float (shape: 6) |  |
| `spin_taylor` | [1D array of TaylorStruct (shape: 0:3)](bmad.md#taylorstruct) | Quaternion Spin Taylor map. |
| `wake` | [WakeStruct](bmad.md#wakestruct) | Wakes |
| `wall3d` | 1D array of Wall3dStruct | Chamber or capillary wall E/M field structs. |
| `cartesian_map` | [1D array of CartesianMapStruct](bmad.md#cartesianmapstruct) | Used to define E/M fields |
| `cylindrical_map` | [1D array of CylindricalMapStruct](bmad.md#cylindricalmapstruct) | Used to define E/M fields |
| `gen_grad_map` | [1D array of GenGradMapStruct](bmad.md#gengradmapstruct) | Used to define E/M fields. |
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

Fortran struct: `ellipse_beam_init_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1127))

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

Fortran struct: `em_field_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2037))

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

Fortran struct: `expression_atom_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L57))

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

Fortran struct: `expression_tree_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L65))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | str |  |
| `type` | int | plus$, minum$, sin$, cos$, etc. |
| `value` | float |  |
| `node` | [1D array of ExpressionTreeStruct](bmad.md#expressiontreestruct) | Child nodes. Note: Pointer used here since Ifort does not support allocatable. |

::: pybmad.FloorPositionStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### FloorPositionStruct

Fortran struct: `floor_position_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L841))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `r` | 1D array of float (shape: 3) | (x, y, z) offset from origin |
| `w` | 2D array of float (shape: 3,3) | W matrix. Columns are unit vectors of the frame axes. |
| `theta` | float | angular orientation consistent with W matrix |
| `phi` | float | angular orientation consistent with W matrix |
| `psi` | float | angular orientation consistent with W matrix |

::: pybmad.GenGrad1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GenGrad1Struct

Fortran struct: `gen_grad1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L764))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `m` | int | Azimuthal index |
| `sincos` | int | sin$ or cos$ |
| `n_deriv_max` | int | Max GG derivative The derivative matrix is extended to include the interpolating spline polynomial. |
| `deriv` | 2D array of float | Range: (iz0:iz1, 0:2*n_deriv_max+1) |

::: pybmad.GenGradMapStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GenGradMapStruct

Fortran struct: `gen_grad_map_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L772))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `file` | str | Input file name. Used also as ID for instances. |
| `gg` | [1D array of GenGrad1Struct](bmad.md#gengrad1struct) |  |
| `ele_anchor_pt` | int | anchor_beginning$, anchor_center$, or anchor_end$ |
| `field_type` | int | or electric$ |
| `iz0` | int | gg%deriv(iz0:iz1, :) lower bound. |
| `iz1` | int | gg%deriv(iz0:iz1, :) upper bound. |
| `dz` | float | Point spacing. |
| `r0` | 1D array of float (shape: 3) | field origin relative to ele_anchor_pt. |
| `field_scale` | float | Factor to scale the fields by |
| `master_parameter` | int | Master parameter in ele%value(:) array to use for scaling the field. |
| `curved_ref_frame` | bool |  |

::: pybmad.GgTaylorStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GgTaylorStruct

Fortran struct: `gg_taylor_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L825))

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

Fortran struct: `gg_taylor_term_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L820))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `coef` | float |  |
| `expn` | 1D array of int (shape: 2) |  |

::: pybmad.GridBeamInitStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### GridBeamInitStruct

Fortran struct: `grid_beam_init_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1139))

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

Fortran struct: `grid_field_pt1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L788))

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

Fortran struct: `grid_field_pt_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L793))

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

Fortran struct: `grid_field_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L799))

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

Fortran struct: `high_energy_space_charge_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L849))

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

::: pybmad.Interval1CoefStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Interval1CoefStruct

Fortran struct: `interval1_coef_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L216))

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

Fortran struct: `kv_beam_init_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1133))

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

Fortran struct: `lat_ele_loc_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L867))

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

Fortran struct: `lat_ele_order1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L874))

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

Fortran struct: `lat_ele_order_array_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L879))

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

Fortran struct: `lat_ele_order_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L898))

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

Fortran struct: `lat_param_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1518))

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

::: pybmad.LatStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### LatStruct

Fortran struct: `lat_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1641))

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

Fortran struct: `linac_normal_mode_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1996))

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

::: pybmad.MadEnergyStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### MadEnergyStruct

Fortran struct: `mad_energy_struct` ([`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L14))

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

Fortran struct: `mad_map_struct` ([`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L23))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `k` | 1D array of float (shape: 6) | 0th order map. |
| `r` | 2D array of float (shape: 6,6) | 1st order map. |
| `t` | 3D array of float (shape: 6,6,6) | 2nd order map. |

::: pybmad.Mode3Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Mode3Struct

Fortran struct: `mode3_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L929))

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

Fortran struct: `mode_info_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1554))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `stable` | bool | Is the mode stable? |
| `tune` | float | "fractional" tune in radians |
| `emit` | float | Emittance (unnormalized). |
| `chrom` | float | Chromaticity. |
| `sigma` | float | Beam size. |
| `sigmap` | float | Beam divergence. |

::: pybmad.NormalModesStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### NormalModesStruct

Fortran struct: `normal_modes_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2006))

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

::: pybmad.PhotonElementStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonElementStruct

Fortran struct: `photon_element_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1085))

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

::: pybmad.PhotonMaterialStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PhotonMaterialStruct

Fortran struct: `photon_material_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1070))

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

Fortran struct: `photon_reflect_surface_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L234))

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

Fortran struct: `photon_reflect_table_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L220))

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

Fortran struct: `photon_target_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1062))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `type` | int | or rectangular$ |
| `n_corner` | int |  |
| `ele_loc` | [LatEleLocStruct](bmad.md#latelelocstruct) |  |
| `corner` | [1D array of TargetPointStruct (shape: 8)](bmad.md#targetpointstruct) |  |
| `center` | [TargetPointStruct](bmad.md#targetpointstruct) |  |

::: pybmad.PixelDetecStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PixelDetecStruct

Fortran struct: `pixel_detec_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1039))

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

Fortran struct: `pixel_pt_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1029))

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

::: pybmad.PreTrackerStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PreTrackerStruct

Fortran struct: `pre_tracker_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1625))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `who` | int | Can be opal$, or impactt$ |
| `ix_ele_start` | int |  |
| `ix_ele_end` | int |  |
| `input_file` | str |  |

::: pybmad.PtcNormalFormStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### PtcNormalFormStruct

Fortran struct: `ptc_normal_form_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1586))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele_origin` | [EleStruct](bmad.md#elestruct) | Element at which the on-turn map was created. |
| `orb0` | 1D array of float (shape: 6) | Closed orbit at element. |
| `valid_map` | bool |  |

::: pybmad.RadInt1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RadInt1Struct

Fortran struct: `rad_int1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2401))

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

Fortran struct: `rad_int_all_ele_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2428))

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

Fortran struct: `rad_int_branch_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2424))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `ele` | [1D array of RadInt1Struct](bmad.md#radint1struct) | Array is indexed from 0 |

::: pybmad.RadMapEleStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### RadMapEleStruct

Fortran struct: `rad_map_ele_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L982))

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

Fortran struct: `rad_map_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L974))

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

Fortran struct: `ramper_lord_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1377))

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

Fortran struct: `resonance_h_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1565))

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

Fortran struct: `rf_ele_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1323))

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

Fortran struct: `rf_stair_step_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1302))

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

::: pybmad.SpaceChargeCommonStruct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SpaceChargeCommonStruct

Fortran struct: `space_charge_common_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2161))

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

Fortran struct: `spin_axis_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L290))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `l` | 1D array of float (shape: 3) | Transverse axis. |
| `n0` | 1D array of float (shape: 3) | Invariant spin axis on closed orbit. |
| `m` | 1D array of float (shape: 3) | Transverse axis. |

::: pybmad.SpinOrbitMap1Struct
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### SpinOrbitMap1Struct

Fortran struct: `spin_orbit_map1_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L321))

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

Fortran struct: `spin_polar_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L312))

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

Fortran struct: `strong_beam_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2055))

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

Fortran struct: `summation_rdt_struct` ([`bmad/modules/srdt_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/srdt_mod.f90#L11))

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

Fortran struct: `surface_curvature_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1049))

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

Fortran struct: `surface_displacement_pt_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1016))

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

Fortran struct: `surface_displacement_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1021))

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

Fortran struct: `surface_h_misalign_pt_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1003))

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

Fortran struct: `surface_h_misalign_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1008))

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

Fortran struct: `surface_segmented_pt_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L990))

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

Fortran struct: `surface_segmented_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L995))

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

Fortran struct: `target_point_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L1058))

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

Fortran struct: `taylor_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L495))

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

Fortran struct: `taylor_term_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L481))

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

Fortran struct: `track_point_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2064))

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

Fortran struct: `track_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2076))

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

Fortran struct: `twiss_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L193))

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

Fortran struct: `wake_lr_mode_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L651))

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

Fortran struct: `wake_lr_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L667))

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

Fortran struct: `wake_sr_mode_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L621))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `amp` | float | Amplitude |
| `damp` | float | Dampling factor. |
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

Fortran struct: `wake_sr_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L635))

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

Fortran struct: `wake_sr_z_long_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L609))

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

Fortran struct: `wake_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L680))

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

Fortran struct: `wall3d_section_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L434))

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

Fortran struct: `wall3d_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L466))

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

Fortran struct: `wall3d_vertex_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L418))

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

Fortran struct: `xy_disp_struct` ([`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L860))

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

Fortran source: [`bmad/modules/multipole_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/multipole_mod.f90#L314)

::: pybmad.ab_multipole_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ab_multipole_kicks

Fortran source: [`bmad/modules/multipole_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/multipole_mod.f90#L82)

::: pybmad.ab_multipole_kicks
    options:
      show_root_heading: false
      show_root_toc_entry: false

### absolute_photon_position

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_init_mod.f90#L61)

::: pybmad.absolute_photon_position
    options:
      show_root_heading: false
      show_root_toc_entry: false

### absolute_time_tracking

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L374)

::: pybmad.absolute_time_tracking
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ac_kicker_amp

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L381)

::: pybmad.ac_kicker_amp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### action_to_xyz

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L398)

::: pybmad.action_to_xyz
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_lattice_control_structs

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L390)

::: pybmad.add_lattice_control_structs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_superimpose

Fortran source: [`bmad/modules/superimpose_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/superimpose_mod.f90#L59)

::: pybmad.add_superimpose
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_this_multipass

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L2870)

::: pybmad.add_this_multipass
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_this_name_to_list

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/output/write_lattice_file_mod.f90#L489)

::: pybmad.add_this_name_to_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### add_this_taylor_term

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L145)

::: pybmad.add_this_taylor_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### adjust_super_slave_names

Fortran source: [`bmad/modules/superimpose_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/superimpose_mod.f90#L813)

::: pybmad.adjust_super_slave_names
    options:
      show_root_heading: false
      show_root_toc_entry: false

### allocate_branch_array

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L398)

::: pybmad.allocate_branch_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### allocate_grid_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L412)

::: pybmad.allocate_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### allocate_lat_ele_array

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L418)

::: pybmad.allocate_lat_ele_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### angle_between_polars

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L437)

::: pybmad.angle_between_polars
    options:
      show_root_heading: false
      show_root_toc_entry: false

### angle_to_canonical_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L444)

::: pybmad.angle_to_canonical_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### aperture_bookkeeper

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bookkeeper_mod.f90#L1840)

::: pybmad.aperture_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apply_all_rampers

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L451)

::: pybmad.apply_all_rampers
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apply_energy_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L470)

::: pybmad.apply_energy_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apply_patch_to_ptc_fibre

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L3302)

::: pybmad.apply_patch_to_ptc_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### apply_rampers_to_slave

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L479)

::: pybmad.apply_rampers_to_slave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### array_re_str

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/output/write_lattice_file_mod.f90#L226)

::: pybmad.array_re_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### astra_max_field_reference

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/astra_interface_mod.f90#L979)

::: pybmad.astra_max_field_reference
    options:
      show_root_heading: false
      show_root_toc_entry: false

### at_this_ele_end

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L486)

::: pybmad.at_this_ele_end
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_bookkeeper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L493)

::: pybmad.attribute_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_free

Fortran sources (overloaded):

- `attribute_free1`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L3072)
- `attribute_free2`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L3101)
- `attribute_free3`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L3140)

::: pybmad.attribute_free
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_index

Fortran sources (overloaded):

- `attribute_index1`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L207)
- `attribute_index2`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L253)

::: pybmad.attribute_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_name

Fortran sources (overloaded):

- `attribute_name1`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L389)
- `attribute_name2`: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L412)

::: pybmad.attribute_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_type

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L2037)

::: pybmad.attribute_type
    options:
      show_root_heading: false
      show_root_toc_entry: false

### attribute_units

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L2140)

::: pybmad.attribute_units
    options:
      show_root_heading: false
      show_root_toc_entry: false

### autoscale_phase_and_amp

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L509)

::: pybmad.autoscale_phase_and_amp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### average_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L518)

::: pybmad.average_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bbi_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L525)

::: pybmad.bbi_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bbi_slice_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L532)

::: pybmad.bbi_slice_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### beam_envelope_ibs

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/envelope_mod.f90#L770)

::: pybmad.beam_envelope_ibs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### beam_equal_beam

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L5321)

::: pybmad.beam_equal_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### beam_init_setup

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L566)

::: pybmad.beam_init_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### beam_tilts

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L1020)

::: pybmad.beam_tilts
    options:
      show_root_heading: false
      show_root_toc_entry: false

### beambeam_fibre_setup

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L2711)

::: pybmad.beambeam_fibre_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/fringe_mod.f90#L39)

::: pybmad.bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_exact_multipole_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L540)

::: pybmad.bend_exact_multipole_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_length_has_been_set

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L551)

::: pybmad.bend_length_has_been_set
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_e_rel_init

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_init_mod.f90#L897)

::: pybmad.bend_photon_e_rel_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_energy_integ_prob

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_init_mod.f90#L232)

::: pybmad.bend_photon_energy_integ_prob
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_energy_normalized_probability

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_init_mod.f90#L1067)

::: pybmad.bend_photon_energy_normalized_probability
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_init

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_init_mod.f90#L134)

::: pybmad.bend_photon_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_polarization_init

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_init_mod.f90#L370)

::: pybmad.bend_photon_polarization_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_photon_vert_angle_init

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_init_mod.f90#L428)

::: pybmad.bend_photon_vert_angle_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_shift

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L558)

::: pybmad.bend_shift
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bend_vert_angle_integ_prob

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_init_mod.f90#L301)

::: pybmad.bend_vert_angle_integ_prob
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bl_via_vlassov

Fortran source: [`bmad/multiparticle/ibs_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/ibs_mod.f90#L837)

::: pybmad.bl_via_vlassov
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bmad_parser

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L586)

::: pybmad.bmad_parser
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bmad_parser2

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L597)

::: pybmad.bmad_parser2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bmad_patch_parameters_to_ptc

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L3043)

::: pybmad.bmad_patch_parameters_to_ptc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bp_set_ran_status

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L6018)

::: pybmad.bp_set_ran_status
    options:
      show_root_heading: false
      show_root_toc_entry: false

### branch_equal_branch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4837)

::: pybmad.branch_equal_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### branch_name

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L607)

::: pybmad.branch_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### branch_to_ptc_m_u

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L75)

::: pybmad.branch_to_ptc_m_u
    options:
      show_root_heading: false
      show_root_toc_entry: false

### bunch_equal_bunch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L5268)

::: pybmad.bunch_equal_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### c_to_cbar

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L645)

::: pybmad.c_to_cbar
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_bunch_params

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_utils.f90#L1314)

::: pybmad.calc_bunch_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_bunch_params_slice

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_utils.f90#L1186)

::: pybmad.calc_bunch_params_slice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_bunch_params_z_slice

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_utils.f90#L1256)

::: pybmad.calc_bunch_params_z_slice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_bunch_sigma_matrix_etc

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_utils.f90#L1745)

::: pybmad.calc_bunch_sigma_matrix_etc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_emittances_and_twiss_from_sigma_matrix

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_utils.f90#L1442)

::: pybmad.calc_emittances_and_twiss_from_sigma_matrix
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_spin_params

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_utils.f90#L1696)

::: pybmad.calc_spin_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_super_slave_key

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L661)

::: pybmad.calc_super_slave_key
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_wall_radius

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/wall3d_mod.f90#L505)

::: pybmad.calc_wall_radius
    options:
      show_root_heading: false
      show_root_toc_entry: false

### calc_z_tune

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L668)

::: pybmad.calc_z_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### canonical_to_angle_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L674)

::: pybmad.canonical_to_angle_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cbar_to_c

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L681)

::: pybmad.cbar_to_c
    options:
      show_root_heading: false
      show_root_toc_entry: false

### check_aperture_limit

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L688)

::: pybmad.check_aperture_limit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### check_controller_controls

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L787)

::: pybmad.check_controller_controls
    options:
      show_root_heading: false
      show_root_toc_entry: false

### check_for_superimpose_problem

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L3922)

::: pybmad.check_for_superimpose_problem
    options:
      show_root_heading: false
      show_root_toc_entry: false

### check_if_s_in_bounds

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L796)

::: pybmad.check_if_s_in_bounds
    options:
      show_root_heading: false
      show_root_toc_entry: false

### choose_quads_for_set_tune

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L806)

::: pybmad.choose_quads_for_set_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### chrom_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L816)

::: pybmad.chrom_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### chrom_tune

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L832)

::: pybmad.chrom_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### classical_radius

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L699)

::: pybmad.classical_radius
    options:
      show_root_heading: false
      show_root_toc_entry: false

### clear_lat_1turn_mats

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L843)

::: pybmad.clear_lat_1turn_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### clear_taylor_maps_from_elements

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L849)

::: pybmad.clear_taylor_maps_from_elements
    options:
      show_root_heading: false
      show_root_toc_entry: false

### closed_orbit_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L855)

::: pybmad.closed_orbit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### closed_orbit_from_tracking

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L865)

::: pybmad.closed_orbit_from_tracking
    options:
      show_root_heading: false
      show_root_toc_entry: false

### cmplx_re_str

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/output/write_lattice_file_mod.f90#L251)

::: pybmad.cmplx_re_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### combine_consecutive_elements

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L876)

::: pybmad.combine_consecutive_elements
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_clean

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L132)

::: pybmad.complex_taylor_clean
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_coef

Fortran sources (overloaded):

- `complex_taylor_coef1`: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L179)
- `complex_taylor_coef2`: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L214)

::: pybmad.complex_taylor_coef
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_equal_complex_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L5134)

::: pybmad.complex_taylor_equal_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_exponent_index

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L661)

::: pybmad.complex_taylor_exponent_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_make_unit

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L440)

::: pybmad.complex_taylor_make_unit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylor_to_mat6

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L689)

::: pybmad.complex_taylor_to_mat6
    options:
      show_root_heading: false
      show_root_toc_entry: false

### complex_taylors_equal_complex_taylors

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L5170)

::: pybmad.complex_taylors_equal_complex_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### compute_slave_coupler

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bookkeeper_mod.f90#L1520)

::: pybmad.compute_slave_coupler
    options:
      show_root_heading: false
      show_root_toc_entry: false

### concat_ele_taylor

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L2147)

::: pybmad.concat_ele_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### concat_taylor

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L2084)

::: pybmad.concat_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### concat_transfer_mat

Fortran source: [`bmad/modules/transfer_map_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/transfer_map_mod.f90#L560)

::: pybmad.concat_transfer_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### control_bookkeeper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L883)

::: pybmad.control_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_bend_exact_multipole

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L926)

::: pybmad.convert_bend_exact_multipole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1023)

::: pybmad.convert_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_field_ele_to_lab

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/em_field_mod.f90#L827)

::: pybmad.convert_field_ele_to_lab
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_local_cartesian_to_local_curvilinear

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/gpt_interface_mod.f90#L1629)

::: pybmad.convert_local_cartesian_to_local_curvilinear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_local_curvilinear_to_local_cartesian

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/gpt_interface_mod.f90#L1612)

::: pybmad.convert_local_curvilinear_to_local_cartesian
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_particle_coordinates_s_to_t

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L891)

::: pybmad.convert_particle_coordinates_s_to_t
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_particle_coordinates_t_to_s

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L899)

::: pybmad.convert_particle_coordinates_t_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_pc_to

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L917)

::: pybmad.convert_pc_to
    options:
      show_root_heading: false
      show_root_toc_entry: false

### convert_total_energy_to

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L908)

::: pybmad.convert_total_energy_to
    options:
      show_root_heading: false
      show_root_toc_entry: false

### converter_distribution_parser

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2063)

::: pybmad.converter_distribution_parser
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coord_equal_coord

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4886)

::: pybmad.coord_equal_coord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coord_state_name

Fortran source: [`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2616)

::: pybmad.coord_state_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_body_to_local

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L714)

::: pybmad.coords_body_to_local
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_body_to_rel_exit

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L705)

::: pybmad.coords_body_to_rel_exit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_curvilinear_to_floor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L778)

::: pybmad.coords_curvilinear_to_floor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_floor_to_curvilinear

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L748)

::: pybmad.coords_floor_to_curvilinear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_floor_to_local_curvilinear

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L738)

::: pybmad.coords_floor_to_local_curvilinear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_floor_to_relative

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L731)

::: pybmad.coords_floor_to_relative
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_local_curvilinear_to_body

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L758)

::: pybmad.coords_local_curvilinear_to_body
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_local_curvilinear_to_floor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L767)

::: pybmad.coords_local_curvilinear_to_floor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coords_relative_to_floor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L723)

::: pybmad.coords_relative_to_floor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### coulombfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L189)

::: pybmad.coulombfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_concatenated_wall3d

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/wall3d_mod.f90#L1204)

::: pybmad.create_concatenated_wall3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_element_slice

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L941)

::: pybmad.create_element_slice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_feedback

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L933)

::: pybmad.create_feedback
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_field_overlap

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L953)

::: pybmad.create_field_overlap
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_girder

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L961)

::: pybmad.create_girder
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_group

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L971)

::: pybmad.create_group
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_lat_ele_nametable

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L979)

::: pybmad.create_lat_ele_nametable
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_overlay

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L986)

::: pybmad.create_overlay
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_planar_wiggler_model

Fortran source: [`bmad/modules/element_modeling_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/element_modeling_mod.f90#L114)

::: pybmad.create_planar_wiggler_model
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_ramper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L994)

::: pybmad.create_ramper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_sol_quad_model

Fortran source: [`bmad/modules/element_modeling_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/element_modeling_mod.f90#L62)

::: pybmad.create_sol_quad_model
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_unique_ele_names

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1002)

::: pybmad.create_unique_ele_names
    options:
      show_root_heading: false
      show_root_toc_entry: false

### create_wiggler_cartesian_map

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1010)

::: pybmad.create_wiggler_cartesian_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### crystal_attribute_bookkeeper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1017)

::: pybmad.crystal_attribute_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### crystal_h_misalign

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L1054)

::: pybmad.crystal_h_misalign
    options:
      show_root_heading: false
      show_root_toc_entry: false

### crystal_type_to_crystal_params

Fortran source: [`bmad/interface/xraylib_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/xraylib_interface.f90#L314)

::: pybmad.crystal_type_to_crystal_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### custom_attribute_ubound_index

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L2832)

::: pybmad.custom_attribute_ubound_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### custom_ele_attrib_name_list

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L3004)

::: pybmad.custom_ele_attrib_name_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### damping_matrix_d

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/envelope_mod.f90#L148)

::: pybmad.damping_matrix_d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### deallocate_ele_pointers

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1040)

::: pybmad.deallocate_ele_pointers
    options:
      show_root_heading: false
      show_root_toc_entry: false

### deallocate_expression_tree

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L652)

::: pybmad.deallocate_expression_tree
    options:
      show_root_heading: false
      show_root_toc_entry: false

### deallocate_lat_pointers

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1047)

::: pybmad.deallocate_lat_pointers
    options:
      show_root_heading: false
      show_root_toc_entry: false

### default_tracking_species

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1053)

::: pybmad.default_tracking_species
    options:
      show_root_heading: false
      show_root_toc_entry: false

### detector_pixel_pt

Fortran source: [`bmad/photon/photon_target_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_target_mod.f90#L320)

::: pybmad.detector_pixel_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### diffraction_plate_or_mask_hit_spot

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1060)

::: pybmad.diffraction_plate_or_mask_hit_spot
    options:
      show_root_heading: false
      show_root_toc_entry: false

### diffusion_matrix_b

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/envelope_mod.f90#L127)

::: pybmad.diffusion_matrix_b
    options:
      show_root_heading: false
      show_root_toc_entry: false

### distance_to_aperture

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1068)

::: pybmad.distance_to_aperture
    options:
      show_root_heading: false
      show_root_toc_entry: false

### do_mode_flip

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1078)

::: pybmad.do_mode_flip
    options:
      show_root_heading: false
      show_root_toc_entry: false

### dpc_given_de

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1085)

::: pybmad.dpc_given_de
    options:
      show_root_heading: false
      show_root_toc_entry: false

### drift_and_pipe_track_methods_adjustment

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L4916)

::: pybmad.drift_and_pipe_track_methods_adjustment
    options:
      show_root_heading: false
      show_root_toc_entry: false

### drift_multipass_name_correction

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L3050)

::: pybmad.drift_multipass_name_correction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### drift_orbit_time

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/time_tracker_mod.f90#L883)

::: pybmad.drift_orbit_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### drift_particle_to_s

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/space_charge_mod.f90#L621)

::: pybmad.drift_particle_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### drift_particle_to_t

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/space_charge_mod.f90#L668)

::: pybmad.drift_particle_to_t
    options:
      show_root_heading: false
      show_root_toc_entry: false

### dspline_len

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/csr_and_space_charge_mod.f90#L1746)

::: pybmad.dspline_len
    options:
      show_root_heading: false
      show_root_toc_entry: false

### dynamic_aperture_point

Fortran source: [`bmad/modules/dynamic_aperture_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/dynamic_aperture_mod.f90#L284)

::: pybmad.dynamic_aperture_point
    options:
      show_root_heading: false
      show_root_toc_entry: false

### dynamic_aperture_scan

Fortran source: [`bmad/modules/dynamic_aperture_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/dynamic_aperture_mod.f90#L27)

::: pybmad.dynamic_aperture_scan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### e_accel_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1091)

::: pybmad.e_accel_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### e_crit_photon

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_init_mod.f90#L1354)

::: pybmad.e_crit_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eigen_decomp_6mat

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L449)

::: pybmad.eigen_decomp_6mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_compute_ref_energy_and_time

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1100)

::: pybmad.ele_compute_ref_energy_and_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_equal_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4260)

::: pybmad.ele_equal_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_equals_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4291)

::: pybmad.ele_equals_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_finalizer

Fortran source: [`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2796)

::: pybmad.ele_finalizer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_full_name

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1109)

::: pybmad.ele_full_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_geometry

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1117)

::: pybmad.ele_geometry
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_geometry_with_misalignments

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1127)

::: pybmad.ele_geometry_with_misalignments
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_has_constant_ds_dt_ref

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1135)

::: pybmad.ele_has_constant_ds_dt_ref
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_has_nonzero_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1142)

::: pybmad.ele_has_nonzero_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_has_nonzero_offset

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1149)

::: pybmad.ele_has_nonzero_offset
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_is_monitor

Fortran source: [`bmad/modules/measurement_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/measurement_mod.f90#L27)

::: pybmad.ele_is_monitor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_loc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1172)

::: pybmad.ele_loc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_loc_name

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1156)

::: pybmad.ele_loc_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_misalignment_l_s_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1165)

::: pybmad.ele_misalignment_l_s_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_nametable_index

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1179)

::: pybmad.ele_nametable_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_order_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1186)

::: pybmad.ele_order_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_reference_energy_correction

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1193)

::: pybmad.ele_reference_energy_correction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_rf_step_index

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1203)

::: pybmad.ele_rf_step_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_to_fibre

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1211)

::: pybmad.ele_to_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_to_ptc_magnetic_bn_an

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L3098)

::: pybmad.ele_to_ptc_magnetic_bn_an
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_to_spin_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1222)

::: pybmad.ele_to_spin_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_to_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1230)

::: pybmad.ele_to_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_unique_name

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1239)

::: pybmad.ele_unique_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_value_has_changed

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1247)

::: pybmad.ele_value_has_changed
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ele_vec_equal_ele_vec

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4671)

::: pybmad.ele_vec_equal_ele_vec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### elec_multipole_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1256)

::: pybmad.elec_multipole_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### element_at_s

Fortran sources (overloaded):

- `element_at_s_branch`: [`bmad/modules/element_at_s_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/element_at_s_mod.f90#L75)
- `element_at_s_lat`: [`bmad/modules/element_at_s_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/element_at_s_mod.f90#L198)

::: pybmad.element_at_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### element_slice_iterator

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1267)

::: pybmad.element_slice_iterator
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ellipinc_test

Fortran source: [`bmad/space_charge/csr3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/csr3d_mod.f90#L628)

::: pybmad.ellipinc_test
    options:
      show_root_heading: false
      show_root_toc_entry: false

### em_field_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1276)

::: pybmad.em_field_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### em_field_derivatives

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/em_field_mod.f90#L655)

::: pybmad.em_field_derivatives
    options:
      show_root_heading: false
      show_root_toc_entry: false

### em_field_kick_vector_time

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/time_tracker_mod.f90#L650)

::: pybmad.em_field_kick_vector_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### em_field_plus_em_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4226)

::: pybmad.em_field_plus_em_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### emit_6d

Fortran source: [`bmad/modules/rad_6d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/rad_6d_mod.f90#L50)

::: pybmad.emit_6d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### entering_element

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1292)

::: pybmad.entering_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### envelope_radints

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/envelope_mod.f90#L614)

::: pybmad.envelope_radints
    options:
      show_root_heading: false
      show_root_toc_entry: false

### envelope_radints_ibs

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/envelope_mod.f90#L510)

::: pybmad.envelope_radints_ibs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ac_kicker

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L135)

::: pybmad.eq_ac_kicker
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ac_kicker_freq

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L113)

::: pybmad.eq_ac_kicker_freq
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ac_kicker_time

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L91)

::: pybmad.eq_ac_kicker_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_anormal_mode

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2264)

::: pybmad.eq_anormal_mode
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_aperture_param

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L3373)

::: pybmad.eq_aperture_param
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_aperture_point

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L3347)

::: pybmad.eq_aperture_point
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_aperture_scan

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L3407)

::: pybmad.eq_aperture_scan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_beam

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L3325)

::: pybmad.eq_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_beam_init

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2076)

::: pybmad.eq_beam_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_bmad_common

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2536)

::: pybmad.eq_bmad_common
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_bookkeeping_state

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1172)

::: pybmad.eq_bookkeeping_state
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_bpm_phase_coupling

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L353)

::: pybmad.eq_bpm_phase_coupling
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_branch

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L3032)

::: pybmad.eq_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_bunch

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L3211)

::: pybmad.eq_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_bunch_params

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L3263)

::: pybmad.eq_bunch_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cartesian_map

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L805)

::: pybmad.eq_cartesian_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cartesian_map_term

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L779)

::: pybmad.eq_cartesian_map_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cartesian_map_term1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L745)

::: pybmad.eq_cartesian_map_term1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_complex_taylor

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L3008)

::: pybmad.eq_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_complex_taylor_term

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2988)

::: pybmad.eq_complex_taylor_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_control

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1868)

::: pybmad.eq_control
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_control_ramp1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1930)

::: pybmad.eq_control_ramp1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_control_var1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1908)

::: pybmad.eq_control_var1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_controller

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1964)

::: pybmad.eq_controller
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_coord

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L273)

::: pybmad.eq_coord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_coord_array

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L331)

::: pybmad.eq_coord_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cylindrical_map

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L882)

::: pybmad.eq_cylindrical_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cylindrical_map_term

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L856)

::: pybmad.eq_cylindrical_map_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_cylindrical_map_term1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L836)

::: pybmad.eq_cylindrical_map_term1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ele

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2732)

::: pybmad.eq_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ellipse_beam_init

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2004)

::: pybmad.eq_ellipse_beam_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_em_field

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2366)

::: pybmad.eq_em_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_expression_atom

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L389)

::: pybmad.eq_expression_atom
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_floor_position

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1010)

::: pybmad.eq_floor_position
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_gen_grad1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1254)

::: pybmad.eq_gen_grad1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_gen_grad_map

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1282)

::: pybmad.eq_gen_grad_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_gg_taylor

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L721)

::: pybmad.eq_gg_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_gg_taylor_term

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L701)

::: pybmad.eq_gg_taylor_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_grid_beam_init

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2048)

::: pybmad.eq_grid_beam_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_grid_field

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L967)

::: pybmad.eq_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_grid_field_pt

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L941)

::: pybmad.eq_grid_field_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_grid_field_pt1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L921)

::: pybmad.eq_grid_field_pt1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_high_energy_space_charge

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1036)

::: pybmad.eq_high_energy_space_charge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_interval1_coef

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L163)

::: pybmad.eq_interval1_coef
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_kv_beam_init

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2026)

::: pybmad.eq_kv_beam_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_lat

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L3086)

::: pybmad.eq_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_lat_ele_loc

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L617)

::: pybmad.eq_lat_ele_loc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_lat_param

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2162)

::: pybmad.eq_lat_param
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_linac_normal_mode

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2294)

::: pybmad.eq_linac_normal_mode
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_mode3

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1144)

::: pybmad.eq_mode3
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_mode_info

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2212)

::: pybmad.eq_mode_info
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_normal_modes

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2324)

::: pybmad.eq_normal_modes
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_photon_element

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1658)

::: pybmad.eq_photon_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_photon_material

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1558)

::: pybmad.eq_photon_material
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_photon_reflect_surface

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L239)

::: pybmad.eq_photon_reflect_surface
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_photon_reflect_table

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L185)

::: pybmad.eq_photon_reflect_table
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_photon_target

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1532)

::: pybmad.eq_photon_target
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_pixel_detec

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1626)

::: pybmad.eq_pixel_detec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_pixel_pt

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1590)

::: pybmad.eq_pixel_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_pre_tracker

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2240)

::: pybmad.eq_pre_tracker
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_rad_int1

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2636)

::: pybmad.eq_rad_int1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_rad_int_all_ele

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2710)

::: pybmad.eq_rad_int_all_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_rad_int_branch

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2688)

::: pybmad.eq_rad_int_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_rad_map

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1206)

::: pybmad.eq_rad_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_rad_map_ele

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1232)

::: pybmad.eq_rad_map_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_ramper_lord

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1843)

::: pybmad.eq_ramper_lord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_space_charge_common

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2486)

::: pybmad.eq_space_charge_common
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_spin_polar

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L67)

::: pybmad.eq_spin_polar
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_spline

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L43)

::: pybmad.eq_spline
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_strong_beam

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2396)

::: pybmad.eq_strong_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_curvature

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1508)

::: pybmad.eq_surface_curvature
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_displacement

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1462)

::: pybmad.eq_surface_displacement
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_displacement_pt

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1434)

::: pybmad.eq_surface_displacement_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_h_misalign

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1406)

::: pybmad.eq_surface_h_misalign
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_h_misalign_pt

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1378)

::: pybmad.eq_surface_h_misalign_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_segmented

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1350)

::: pybmad.eq_surface_segmented
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_surface_segmented_pt

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1324)

::: pybmad.eq_surface_segmented_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_target_point

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1490)

::: pybmad.eq_target_point
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_taylor

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L677)

::: pybmad.eq_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_taylor_term

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L657)

::: pybmad.eq_taylor_term
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_track

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2456)

::: pybmad.eq_track
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_track_point

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L2426)

::: pybmad.eq_track_point
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_twiss

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1096)

::: pybmad.eq_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L637)

::: pybmad.eq_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake_lr

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L583)

::: pybmad.eq_wake_lr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake_lr_mode

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L541)

::: pybmad.eq_wake_lr_mode
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake_sr

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L497)

::: pybmad.eq_wake_sr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake_sr_mode

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L461)

::: pybmad.eq_wake_sr_mode
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wake_sr_z_long

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L411)

::: pybmad.eq_wake_sr_z_long
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wall3d

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1803)

::: pybmad.eq_wall3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wall3d_section

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1740)

::: pybmad.eq_wall3d_section
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_wall3d_vertex

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1706)

::: pybmad.eq_wall3d_vertex
    options:
      show_root_heading: false
      show_root_toc_entry: false

### eq_xy_disp

Fortran source: [`bmad/modules/equality_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/equality_mod.f90#L1068)

::: pybmad.eq_xy_disp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### equal_sign_here

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7980)

::: pybmad.equal_sign_here
    options:
      show_root_heading: false
      show_root_toc_entry: false

### equivalent_taylor_attributes

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1300)

::: pybmad.equivalent_taylor_attributes
    options:
      show_root_heading: false
      show_root_toc_entry: false

### etdiv

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/envelope_mod.f90#L1736)

::: pybmad.etdiv
    options:
      show_root_heading: false
      show_root_toc_entry: false

### evaluate_array_index

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L808)

::: pybmad.evaluate_array_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### evaluate_logical

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L859)

::: pybmad.evaluate_logical
    options:
      show_root_heading: false
      show_root_toc_entry: false

### exact_bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/fringe_mod.f90#L1465)

::: pybmad.exact_bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### exp_bessi0

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/touschek_mod.f90#L662)

::: pybmad.exp_bessi0
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expect_one_of

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7934)

::: pybmad.expect_one_of
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expect_this

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7793)

::: pybmad.expect_this
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_stack_to_string

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L1760)

::: pybmad.expression_stack_to_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_stack_value

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L1476)

::: pybmad.expression_stack_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_string_to_stack

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L890)

::: pybmad.expression_string_to_stack
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_string_to_tree

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L81)

::: pybmad.expression_string_to_tree
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_tree_to_string

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L726)

::: pybmad.expression_tree_to_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### expression_value

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L1408)

::: pybmad.expression_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fft1

Fortran source: [`bmad/space_charge/fast_fourier_am.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/fast_fourier_am.f90#L37)

::: pybmad.fft1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fibre_to_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1307)

::: pybmad.fibre_to_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### field_attribute_free

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L3468)

::: pybmad.field_attribute_free
    options:
      show_root_heading: false
      show_root_toc_entry: false

### finalize_reflectivity_table

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_reflection_mod.f90#L356)

::: pybmad.finalize_reflectivity_table
    options:
      show_root_heading: false
      show_root_toc_entry: false

### find_element_ends

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1317)

::: pybmad.find_element_ends
    options:
      show_root_heading: false
      show_root_toc_entry: false

### find_fwhm

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/longitudinal_profile_mod.f90#L401)

::: pybmad.find_fwhm
    options:
      show_root_heading: false
      show_root_toc_entry: false

### find_matching_fieldmap

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1325)

::: pybmad.find_matching_fieldmap
    options:
      show_root_heading: false
      show_root_toc_entry: false

### find_normalization

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/longitudinal_profile_mod.f90#L355)

::: pybmad.find_normalization
    options:
      show_root_heading: false
      show_root_toc_entry: false

### floor_angles_to_w_mat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1335)

::: pybmad.floor_angles_to_w_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### floor_w_mat_to_angles

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1342)

::: pybmad.floor_w_mat_to_angles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### form_complex_taylor

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L1593)

::: pybmad.form_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### form_digested_bmad_file_name

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L5401)

::: pybmad.form_digested_bmad_file_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### fringe_here

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1349)

::: pybmad.fringe_here
    options:
      show_root_heading: false
      show_root_toc_entry: false

### g_bend_from_em_field

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/em_field_mod.f90#L63)

::: pybmad.g_bend_from_em_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### g_bending_strength_from_em_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1358)

::: pybmad.g_bending_strength_from_em_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### g_integrals_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1452)

::: pybmad.g_integrals_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gamma_ref

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1370)

::: pybmad.gamma_ref
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gen_grad1_to_gg_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1386)

::: pybmad.gen_grad1_to_gg_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gen_grad_at_s_to_gg_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1377)

::: pybmad.gen_grad_at_s_to_gg_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gen_grad_field

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/em_field_mod.f90#L749)

::: pybmad.gen_grad_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_bl_from_fwhm

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/longitudinal_profile_mod.f90#L518)

::: pybmad.get_bl_from_fwhm
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_called_file

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L68)

::: pybmad.get_called_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_emit_from_sigma_mat

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L946)

::: pybmad.get_emit_from_sigma_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_list_of_names

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L2257)

::: pybmad.get_list_of_names
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_next_word

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L257)

::: pybmad.get_next_word
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_sequence_args

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L4012)

::: pybmad.get_sequence_args
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_slave_list

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1395)

::: pybmad.get_slave_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### get_switch

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7834)

::: pybmad.get_switch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_taylor_equal_gg_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4991)

::: pybmad.gg_taylor_equal_gg_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gg_taylors_equal_gg_taylors

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L5037)

::: pybmad.gg_taylors_equal_gg_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gpt_field_grid_scaling

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/gpt_interface_mod.f90#L747)

::: pybmad.gpt_field_grid_scaling
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gpt_max_field_reference

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/gpt_interface_mod.f90#L1555)

::: pybmad.gpt_max_field_reference
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gpt_to_particle_bunch

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/gpt_interface_mod.f90#L35)

::: pybmad.gpt_to_particle_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### gradient_shift_sr_wake

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1403)

::: pybmad.gradient_shift_sr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### grid_field_interpolate

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/em_field_mod.f90#L222)

::: pybmad.grid_field_interpolate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hard_multipole_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/fringe_mod.f90#L691)

::: pybmad.hard_multipole_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### has_attribute

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L2804)

::: pybmad.has_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### has_curvature

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_utils_mod.f90#L31)

::: pybmad.has_curvature
    options:
      show_root_heading: false
      show_root_toc_entry: false

### has_orientation_attributes

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L1988)

::: pybmad.has_orientation_attributes
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hdf5_write_beam

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1433)

::: pybmad.hdf5_write_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hdf5_write_grid_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1443)

::: pybmad.hdf5_write_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### hwang_bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/fringe_mod.f90#L240)

::: pybmad.hwang_bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ibs_matrix_c

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/envelope_mod.f90#L733)

::: pybmad.ibs_matrix_c
    options:
      show_root_heading: false
      show_root_toc_entry: false

### igfcoulombfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L205)

::: pybmad.igfcoulombfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### igfexfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L243)

::: pybmad.igfexfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### igfeyfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L288)

::: pybmad.igfeyfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### igfezfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L333)

::: pybmad.igfezfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_attribute_name1

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L1942)

::: pybmad.init_attribute_name1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_attribute_name_array

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L572)

::: pybmad.init_attribute_name_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_beam_distribution

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_utils.f90#L200)

::: pybmad.init_beam_distribution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_bmad

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1467)

::: pybmad.init_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_bmad_parser_common

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1472)

::: pybmad.init_bmad_parser_common
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_bunch_distribution

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_utils.f90#L316)

::: pybmad.init_bunch_distribution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_complex_taylor_series

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L5208)

::: pybmad.init_complex_taylor_series
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_coord

Fortran sources (overloaded):

- `init_coord1`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L247)
- `init_coord2`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L258)
- `init_coord3`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L268)

::: pybmad.init_coord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_custom

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1478)

::: pybmad.init_custom
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1484)

::: pybmad.init_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_gg_taylor_series

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L5075)

::: pybmad.init_gg_taylor_series
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1502)

::: pybmad.init_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_multipole_cache

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1510)

::: pybmad.init_multipole_cache
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_photon_from_a_photon_init_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1458)

::: pybmad.init_photon_from_a_photon_init_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_photon_integ_prob

Fortran source: [`bmad/photon/photon_init_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_init_mod.f90#L1394)

::: pybmad.init_photon_integ_prob
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_spin_distribution

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_utils.f90#L1129)

::: pybmad.init_spin_distribution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_surface_segment

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L8186)

::: pybmad.init_surface_segment
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_taylor_series

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1516)

::: pybmad.init_taylor_series
    options:
      show_root_heading: false
      show_root_toc_entry: false

### init_wake

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1524)

::: pybmad.init_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### insert_element

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1532)

::: pybmad.insert_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### integrand_base

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/touschek_mod.f90#L585)

::: pybmad.integrand_base
    options:
      show_root_heading: false
      show_root_toc_entry: false

### integrate_psi

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/longitudinal_profile_mod.f90#L290)

::: pybmad.integrate_psi
    options:
      show_root_heading: false
      show_root_toc_entry: false

### integrated_mats

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/envelope_mod.f90#L409)

::: pybmad.integrated_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### integration_timer

Fortran sources (overloaded):

- `integration_timer_ele`: [`bmad/modules/integration_timer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/integration_timer_mod.f90#L45)
- `integration_timer_fibre`: [`bmad/modules/integration_timer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/integration_timer_mod.f90#L74)

::: pybmad.integration_timer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ion_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1542)

::: pybmad.ion_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### is_attribute

Fortran source: [`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2658)

::: pybmad.is_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### key_name_to_key_index

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1550)

::: pybmad.key_name_to_key_index
    options:
      show_root_heading: false
      show_root_toc_entry: false

### kick_vector_calc

Fortran source: [`bmad/modules/runge_kutta_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/runge_kutta_mod.f90#L638)

::: pybmad.kick_vector_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### kill_complex_taylor

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L565)

::: pybmad.kill_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### kill_ptc_layouts

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1558)

::: pybmad.kill_ptc_layouts
    options:
      show_root_heading: false
      show_root_toc_entry: false

### kill_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1564)

::: pybmad.kill_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### kind_name

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L944)

::: pybmad.kind_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### knot_interpolate

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1570)

::: pybmad.knot_interpolate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### knots_to_string

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1578)

::: pybmad.knots_to_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lafun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L229)

::: pybmad.lafun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_compute_ref_energy_and_time

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1585)

::: pybmad.lat_compute_ref_energy_and_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_ele_locator

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1592)

::: pybmad.lat_ele_locator
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_equal_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4714)

::: pybmad.lat_equal_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_geometry

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1604)

::: pybmad.lat_geometry
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_make_mat6

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1610)

::: pybmad.lat_make_mat6
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_sanity_check

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1619)

::: pybmad.lat_sanity_check
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_to_ptc_layout

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1626)

::: pybmad.lat_to_ptc_layout
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lat_vec_equal_lat_vec

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4795)

::: pybmad.lat_vec_equal_lat_vec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lattice_bookkeeper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1632)

::: pybmad.lattice_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lcavity_rf_step_setup

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1639)

::: pybmad.lcavity_rf_step_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### linear_bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/fringe_mod.f90#L151)

::: pybmad.linear_bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### linear_coef

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L2001)

::: pybmad.linear_coef
    options:
      show_root_heading: false
      show_root_toc_entry: false

### linear_to_spin_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1645)

::: pybmad.linear_to_spin_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### load_parse_line

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L595)

::: pybmad.load_parse_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### lord_edge_aligned

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1651)

::: pybmad.lord_edge_aligned
    options:
      show_root_heading: false
      show_root_toc_entry: false

### low_energy_z_correction

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1659)

::: pybmad.low_energy_z_correction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_add_offsets_and_multipoles

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L185)

::: pybmad.mad_add_offsets_and_multipoles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_concat_map2

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1473)

::: pybmad.mad_concat_map2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_drift

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L336)

::: pybmad.mad_drift
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_elsep

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L394)

::: pybmad.mad_elsep
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_map_to_taylor

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1651)

::: pybmad.mad_map_to_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_quadrupole

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1049)

::: pybmad.mad_quadrupole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_rfcavity

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1165)

::: pybmad.mad_rfcavity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_sbend

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L614)

::: pybmad.mad_sbend
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_sbend_body

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L762)

::: pybmad.mad_sbend_body
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_sbend_fringe

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L678)

::: pybmad.mad_sbend_fringe
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_sextupole

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L523)

::: pybmad.mad_sextupole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_solenoid

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1215)

::: pybmad.mad_solenoid
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_tmfoc

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1000)

::: pybmad.mad_tmfoc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_tmsymm

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1327)

::: pybmad.mad_tmsymm
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_tmtilt

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1367)

::: pybmad.mad_tmtilt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mad_track1

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1554)

::: pybmad.mad_track1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_g2_mats

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1677)

::: pybmad.make_g2_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_g_mats

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1669)

::: pybmad.make_g_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_hvbp

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L198)

::: pybmad.make_hvbp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_hybrid_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1684)

::: pybmad.make_hybrid_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mad_map

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L98)

::: pybmad.make_mad_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1705)

::: pybmad.make_mat6
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_bmad

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1722)

::: pybmad.make_mat6_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_bmad_photon

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1731)

::: pybmad.make_mat6_bmad_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_high_energy_space_charge

Fortran source: [`bmad/space_charge/high_energy_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/high_energy_space_charge_mod.f90#L232)

::: pybmad.make_mat6_high_energy_space_charge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_mad

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L56)

::: pybmad.make_mat6_mad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_symp_lie_ptc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1740)

::: pybmad.make_mat6_symp_lie_ptc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1714)

::: pybmad.make_mat6_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_mat6_tracking

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1747)

::: pybmad.make_mat6_tracking
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_n

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L825)

::: pybmad.make_n
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_pbrh

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/envelope_mod.f90#L654)

::: pybmad.make_pbrh
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_smat_from_abc

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L1066)

::: pybmad.make_smat_from_abc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_unit_mad_map

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1888)

::: pybmad.make_unit_mad_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_v

Fortran source: [`bmad/multiparticle/envelope_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/envelope_mod.f90#L275)

::: pybmad.make_v
    options:
      show_root_heading: false
      show_root_toc_entry: false

### make_v_mats

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1757)

::: pybmad.make_v_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### makeup_control_slave

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bookkeeper_mod.f90#L1557)

::: pybmad.makeup_control_slave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### makeup_group_lord

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bookkeeper_mod.f90#L22)

::: pybmad.makeup_group_lord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### makeup_multipass_slave

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bookkeeper_mod.f90#L343)

::: pybmad.makeup_multipass_slave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### makeup_super_slave

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bookkeeper_mod.f90#L511)

::: pybmad.makeup_super_slave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### makeup_super_slave1

Fortran source: [`bmad/modules/bookkeeper_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bookkeeper_mod.f90#L1172)

::: pybmad.makeup_super_slave1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### map1_inverse

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1693)

::: pybmad.map1_inverse
    options:
      show_root_heading: false
      show_root_toc_entry: false

### map1_make_unit

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1699)

::: pybmad.map1_make_unit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### map1_times_map1

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4186)

::: pybmad.map1_times_map1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### map_to_angle_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1765)

::: pybmad.map_to_angle_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mark_patch_regions

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/wall3d_mod.f90#L1538)

::: pybmad.mark_patch_regions
    options:
      show_root_heading: false
      show_root_toc_entry: false

### master_parameter_value

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1771)

::: pybmad.master_parameter_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mat4_multipole

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1804)

::: pybmad.mat4_multipole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mat6_add_offsets

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1779)

::: pybmad.mat6_add_offsets
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mat6_add_pitch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1786)

::: pybmad.mat6_add_pitch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mat6_to_complex_taylor

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L772)

::: pybmad.mat6_to_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mat_symp_decouple

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1793)

::: pybmad.mat_symp_decouple
    options:
      show_root_heading: false
      show_root_toc_entry: false

### match_ele_to_mat6

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1813)

::: pybmad.match_ele_to_mat6
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mexp

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1823)

::: pybmad.mexp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mfft1

Fortran source: [`bmad/space_charge/fast_fourier_am.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/fast_fourier_am.f90#L53)

::: pybmad.mfft1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### misalign_ptc_fibre

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L2822)

::: pybmad.misalign_ptc_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### momentum_compaction

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1830)

::: pybmad.momentum_compaction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multi_turn_tracking_analysis

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1837)

::: pybmad.multi_turn_tracking_analysis
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multilayer_type_to_multilayer_params

Fortran source: [`bmad/interface/xraylib_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/xraylib_interface.f90#L152)

::: pybmad.multilayer_type_to_multilayer_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipass_chain

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1865)

::: pybmad.multipass_chain
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole1_ab_to_kt

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1874)

::: pybmad.multipole1_ab_to_kt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole1_kt_to_ab

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1882)

::: pybmad.multipole1_kt_to_ab
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_ab_to_kt

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1890)

::: pybmad.multipole_ab_to_kt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_ele_to_ab

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1897)

::: pybmad.multipole_ele_to_ab
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_ele_to_kt

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1910)

::: pybmad.multipole_ele_to_kt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_init

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1927)

::: pybmad.multipole_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_kick

Fortran source: [`bmad/modules/multipole_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/multipole_mod.f90#L203)

::: pybmad.multipole_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_kick_mat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1935)

::: pybmad.multipole_kick_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_kicks

Fortran source: [`bmad/modules/multipole_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/multipole_mod.f90#L32)

::: pybmad.multipole_kicks
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_kt_to_ab

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1920)

::: pybmad.multipole_kt_to_ab
    options:
      show_root_heading: false
      show_root_toc_entry: false

### multipole_spin_tracking

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1945)

::: pybmad.multipole_spin_tracking
    options:
      show_root_heading: false
      show_root_toc_entry: false

### mytan

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L533)

::: pybmad.mytan
    options:
      show_root_heading: false
      show_root_toc_entry: false

### n_attrib_string_max_len

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L2774)

::: pybmad.n_attrib_string_max_len
    options:
      show_root_heading: false
      show_root_toc_entry: false

### new_control

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1953)

::: pybmad.new_control
    options:
      show_root_heading: false
      show_root_toc_entry: false

### nint_chk

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L42)

::: pybmad.nint_chk
    options:
      show_root_heading: false
      show_root_toc_entry: false

### normal_form_complex_taylors

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L976)

::: pybmad.normal_form_complex_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### normal_form_taylors

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L911)

::: pybmad.normal_form_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### normal_mode3_calc

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L130)

::: pybmad.normal_mode3_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### normal_mode_dispersion

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1961)

::: pybmad.normal_mode_dispersion
    options:
      show_root_heading: false
      show_root_toc_entry: false

### normalize_evecs

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L1129)

::: pybmad.normalize_evecs
    options:
      show_root_heading: false
      show_root_toc_entry: false

### num_field_eles

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1967)

::: pybmad.num_field_eles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### num_lords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1974)

::: pybmad.num_lords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### odeint_bmad

Fortran source: [`bmad/modules/runge_kutta_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/runge_kutta_mod.f90#L50)

::: pybmad.odeint_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### odeint_bmad_time

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/time_tracker_mod.f90#L46)

::: pybmad.odeint_bmad_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### offset_particle

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1981)

::: pybmad.offset_particle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### offset_photon

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L1994)

::: pybmad.offset_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### one_turn_mat_at_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2004)

::: pybmad.one_turn_mat_at_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### open_binary_file

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/binary_parser_mod.f90#L406)

::: pybmad.open_binary_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### orbit_amplitude_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2013)

::: pybmad.orbit_amplitude_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### orbit_reference_energy_correction

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2021)

::: pybmad.orbit_reference_energy_correction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### orbit_to_floor_phase_space

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2030)

::: pybmad.orbit_to_floor_phase_space
    options:
      show_root_heading: false
      show_root_toc_entry: false

### orbit_to_local_curvilinear

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2038)

::: pybmad.orbit_to_local_curvilinear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### orbit_too_large

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2047)

::: pybmad.orbit_too_large
    options:
      show_root_heading: false
      show_root_toc_entry: false

### order_evecs_by_n_similarity

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L568)

::: pybmad.order_evecs_by_n_similarity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### order_evecs_by_plane_dominance

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L661)

::: pybmad.order_evecs_by_plane_dominance
    options:
      show_root_heading: false
      show_root_toc_entry: false

### order_evecs_by_tune

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L710)

::: pybmad.order_evecs_by_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### order_particles_in_z

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/wake_mod.f90#L625)

::: pybmad.order_particles_in_z
    options:
      show_root_heading: false
      show_root_toc_entry: false

### order_super_lord_slaves

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2056)

::: pybmad.order_super_lord_slaves
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_alloc_freespace_array

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L151)

::: pybmad.osc_alloc_freespace_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_alloc_image_array

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L1095)

::: pybmad.osc_alloc_image_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_alloc_rectpipe_arrays

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L876)

::: pybmad.osc_alloc_rectpipe_arrays
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_getgrnpipe

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L632)

::: pybmad.osc_getgrnpipe
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_read_rectpipe_grn

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L825)

::: pybmad.osc_read_rectpipe_grn
    options:
      show_root_heading: false
      show_root_toc_entry: false

### osc_write_rectpipe_grn

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L851)

::: pybmad.osc_write_rectpipe_grn
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_cartesian_map

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L6209)

::: pybmad.parse_cartesian_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_cylindrical_map

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L6408)

::: pybmad.parse_cylindrical_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_gen_grad_map

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7011)

::: pybmad.parse_gen_grad_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_grid_field

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L6631)

::: pybmad.parse_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_integer_list

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7262)

::: pybmad.parse_integer_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_integer_list2

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7328)

::: pybmad.parse_integer_list2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_real_list

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7439)

::: pybmad.parse_real_list
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parse_real_list2

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7605)

::: pybmad.parse_real_list2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_add_constant

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L1284)

::: pybmad.parser_add_constant
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_call_check

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L182)

::: pybmad.parser_call_check
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_fast_complex_read

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L8457)

::: pybmad.parser_fast_complex_read
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_fast_integer_read

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L8390)

::: pybmad.parser_fast_integer_read
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_fast_real_read

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L8559)

::: pybmad.parser_fast_real_read
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_file_stack

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L364)

::: pybmad.parser_file_stack
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_get_integer

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7717)

::: pybmad.parser_get_integer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_get_logical

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L7747)

::: pybmad.parser_get_logical
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_identify_fork_to_element

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L5524)

::: pybmad.parser_identify_fork_to_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_init_custom_elements

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L8059)

::: pybmad.parser_init_custom_elements
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_print_line

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L8010)

::: pybmad.parser_print_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_read_lr_wake

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L1688)

::: pybmad.parser_read_lr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_read_old_format_lr_wake

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L1826)

::: pybmad.parser_read_old_format_lr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_read_old_format_sr_wake

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L1943)

::: pybmad.parser_read_old_format_sr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_read_sr_wake

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L1450)

::: pybmad.parser_read_sr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### parser_transfer_control_struct

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L8286)

::: pybmad.parser_transfer_control_struct
    options:
      show_root_heading: false
      show_root_toc_entry: false

### particle_in_global_frame

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/time_tracker_mod.f90#L815)

::: pybmad.particle_in_global_frame
    options:
      show_root_heading: false
      show_root_toc_entry: false

### particle_is_moving_backwards

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2071)

::: pybmad.particle_is_moving_backwards
    options:
      show_root_heading: false
      show_root_toc_entry: false

### particle_is_moving_forward

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2078)

::: pybmad.particle_is_moving_forward
    options:
      show_root_heading: false
      show_root_toc_entry: false

### particle_rf_time

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2086)

::: pybmad.particle_rf_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### patch_flips_propagation_direction

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2096)

::: pybmad.patch_flips_propagation_direction
    options:
      show_root_heading: false
      show_root_toc_entry: false

### patch_length

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2103)

::: pybmad.patch_length
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_absorption_and_phase_shift

Fortran source: [`bmad/interface/xraylib_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/xraylib_interface.f90#L34)

::: pybmad.photon_absorption_and_phase_shift
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_add_to_detector_statistics

Fortran source: [`bmad/photon/photon_target_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_target_mod.f90#L239)

::: pybmad.photon_add_to_detector_statistics
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_reflection

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_reflection_mod.f90#L767)

::: pybmad.photon_reflection
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_reflection_std_surface_init

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_reflection_mod.f90#L57)

::: pybmad.photon_reflection_std_surface_init
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_reflectivity

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_reflection_mod.f90#L626)

::: pybmad.photon_reflectivity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_target_corner_calc

Fortran source: [`bmad/photon/photon_target_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_target_mod.f90#L173)

::: pybmad.photon_target_corner_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_target_setup

Fortran source: [`bmad/photon/photon_target_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_target_mod.f90#L29)

::: pybmad.photon_target_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### photon_type

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_utils_mod.f90#L57)

::: pybmad.photon_type
    options:
      show_root_heading: false
      show_root_toc_entry: false

### physical_ele_end

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2111)

::: pybmad.physical_ele_end
    options:
      show_root_heading: false
      show_root_toc_entry: false

### point_photon_emission

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L304)

::: pybmad.point_photon_emission
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_branch

Fortran sources (overloaded):

- `pointer_to_branch_given_name`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L75)
- `pointer_to_branch_given_ele`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L85)

::: pybmad.pointer_to_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_ele

Fortran sources (overloaded):

- `pointer_to_ele1`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4047)
- `pointer_to_ele2`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4095)
- `pointer_to_ele3`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4123)
- `pointer_to_ele4`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4156)

::: pybmad.pointer_to_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_element_at_s

Fortran source: [`bmad/modules/element_at_s_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/element_at_s_mod.f90#L267)

::: pybmad.pointer_to_element_at_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_fibre

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2139)

::: pybmad.pointer_to_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_field_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2146)

::: pybmad.pointer_to_field_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_girder

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2155)

::: pybmad.pointer_to_girder
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_lord

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2173)

::: pybmad.pointer_to_lord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_multipass_lord

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2183)

::: pybmad.pointer_to_multipass_lord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_next_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2192)

::: pybmad.pointer_to_next_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_slave

Fortran source: [`bmad/modules/bmad_struct.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_struct.f90#L2724)

::: pybmad.pointer_to_slave
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_super_lord

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2201)

::: pybmad.pointer_to_super_lord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_surface_displacement_pt

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_utils_mod.f90#L586)

::: pybmad.pointer_to_surface_displacement_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_surface_segmented_pt

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_utils_mod.f90#L484)

::: pybmad.pointer_to_surface_segmented_pt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_wake_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2210)

::: pybmad.pointer_to_wake_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pointer_to_wall3d

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/wall3d_mod.f90#L1089)

::: pybmad.pointer_to_wall3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### polar_to_spinor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2232)

::: pybmad.polar_to_spinor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### polar_to_vec

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2239)

::: pybmad.polar_to_vec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### project_emit_to_xyz

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L1191)

::: pybmad.project_emit_to_xyz
    options:
      show_root_heading: false
      show_root_toc_entry: false

### psi_prime_sca

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/longitudinal_profile_mod.f90#L86)

::: pybmad.psi_prime_sca
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_bookkeeper

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2246)

::: pybmad.ptc_bookkeeper
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_calculate_tracking_step_size

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L1431)

::: pybmad.ptc_calculate_tracking_step_size
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_check_for_lost_particle

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L1550)

::: pybmad.ptc_check_for_lost_particle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_closed_orbit_calc

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L657)

::: pybmad.ptc_closed_orbit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_emit_calc

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L347)

::: pybmad.ptc_emit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_layouts_resplit

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L1507)

::: pybmad.ptc_layouts_resplit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_one_turn_mat_and_closed_orbit_calc

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L305)

::: pybmad.ptc_one_turn_mat_and_closed_orbit_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_ran_seed_put

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2259)

::: pybmad.ptc_ran_seed_put
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_read_flat_file

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2264)

::: pybmad.ptc_read_flat_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_set_rf_state_for_c_normal

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2273)

::: pybmad.ptc_set_rf_state_for_c_normal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_set_taylor_order_if_needed

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L47)

::: pybmad.ptc_set_taylor_order_if_needed
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_spin_calc

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L457)

::: pybmad.ptc_spin_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_track_all

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L587)

::: pybmad.ptc_track_all
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ptc_transfer_map_with_spin

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2285)

::: pybmad.ptc_transfer_map_with_spin
    options:
      show_root_heading: false
      show_root_toc_entry: false

### pwd_mat

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/longitudinal_profile_mod.f90#L679)

::: pybmad.pwd_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rad1_damp_and_stoc_mats

Fortran source: [`bmad/modules/rad_6d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/rad_6d_mod.f90#L406)

::: pybmad.rad1_damp_and_stoc_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rad_damp_and_stoc_mats

Fortran source: [`bmad/modules/rad_6d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/rad_6d_mod.f90#L261)

::: pybmad.rad_damp_and_stoc_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rad_g_integrals

Fortran source: [`bmad/modules/rad_6d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/rad_6d_mod.f90#L744)

::: pybmad.rad_g_integrals
    options:
      show_root_heading: false
      show_root_toc_entry: false

### radiation_integrals

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2304)

::: pybmad.radiation_integrals
    options:
      show_root_heading: false
      show_root_toc_entry: false

### radiation_map_setup

Fortran source: [`bmad/modules/radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/radiation_mod.f90#L166)

::: pybmad.radiation_map_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ramper_slave_setup

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L614)

::: pybmad.ramper_slave_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ramper_value

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L621)

::: pybmad.ramper_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### randomize_lr_wake_frequencies

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/wake_mod.f90#L29)

::: pybmad.randomize_lr_wake_frequencies
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rchomp

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/output/write_lattice_file_mod.f90#L270)

::: pybmad.rchomp
    options:
      show_root_heading: false
      show_root_toc_entry: false

### re_allocate

Fortran sources (overloaded):

- `re_allocate_wall3d_vertex_array`: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/wall3d_mod.f90#L36)
- `re_allocate_wall3d_section_array`: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/wall3d_mod.f90#L84)

::: pybmad.re_allocate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### re_allocate_eles

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2314)

::: pybmad.re_allocate_eles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### re_associate_node_array

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L836)

::: pybmad.re_associate_node_array
    options:
      show_root_heading: false
      show_root_toc_entry: false

### re_str

Fortran sources (overloaded):

- `re_str_rp`: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/output/write_lattice_file_mod.f90#L184)
- `re_str_qp`: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/output/write_lattice_file_mod.f90#L205)

::: pybmad.re_str
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_beam_ascii

Fortran source: [`bmad/multiparticle/beam_file_io.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_file_io.f90#L834)

::: pybmad.read_beam_ascii
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_beam_file

Fortran source: [`bmad/multiparticle/beam_file_io.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_file_io.f90#L380)

::: pybmad.read_beam_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_binary_cartesian_map

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/binary_parser_mod.f90#L85)

::: pybmad.read_binary_cartesian_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_binary_cylindrical_map

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/binary_parser_mod.f90#L211)

::: pybmad.read_binary_cylindrical_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_binary_grid_field

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/binary_parser_mod.f90#L338)

::: pybmad.read_binary_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_digested_bmad_file

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2358)

::: pybmad.read_digested_bmad_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### read_surface_reflection_file

Fortran source: [`bmad/photon/photon_reflection_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_reflection_mod.f90#L403)

::: pybmad.read_surface_reflection_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_beam

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2368)

::: pybmad.reallocate_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_bp_com_const

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L3112)

::: pybmad.reallocate_bp_com_const
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_bunch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2377)

::: pybmad.reallocate_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_control

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2322)

::: pybmad.reallocate_control
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_coord

Fortran sources (overloaded):

- `reallocate_coord_n`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L168)
- `reallocate_coord_lat`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L175)
- `reallocate_coord_array`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L183)

::: pybmad.reallocate_coord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reallocate_expression_stack

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2329)

::: pybmad.reallocate_expression_stack
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rel_tracking_charge_to_mass

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2337)

::: pybmad.rel_tracking_charge_to_mass
    options:
      show_root_heading: false
      show_root_toc_entry: false

### relative_mode_flip

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2385)

::: pybmad.relative_mode_flip
    options:
      show_root_heading: false
      show_root_toc_entry: false

### release_rad_int_cache

Fortran source: [`bmad/modules/radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/radiation_mod.f90#L23)

::: pybmad.release_rad_int_cache
    options:
      show_root_heading: false
      show_root_toc_entry: false

### remove_constant_taylor

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L1915)

::: pybmad.remove_constant_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### remove_dead_from_bunch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L630)

::: pybmad.remove_dead_from_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### remove_eles_from_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2345)

::: pybmad.remove_eles_from_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### remove_lord_slave_link

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2352)

::: pybmad.remove_lord_slave_link
    options:
      show_root_heading: false
      show_root_toc_entry: false

### reverse_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2393)

::: pybmad.reverse_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rf_coupler_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2409)

::: pybmad.rf_coupler_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rf_is_on

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2421)

::: pybmad.rf_is_on
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rf_ref_time_offset

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2429)

::: pybmad.rf_ref_time_offset
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rfun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L723)

::: pybmad.rfun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rk_adaptive_time_step

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/time_tracker_mod.f90#L392)

::: pybmad.rk_adaptive_time_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rk_time_step1

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/time_tracker_mod.f90#L518)

::: pybmad.rk_time_step1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate3

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/astra_interface_mod.f90#L362)

::: pybmad.rotate3
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_em_field

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/em_field_mod.f90#L173)

::: pybmad.rotate_em_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_field_zx

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/gpt_interface_mod.f90#L1588)

::: pybmad.rotate_field_zx
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_for_curved_surface

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2437)

::: pybmad.rotate_for_curved_surface
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_spin

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2446)

::: pybmad.rotate_spin
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_spin_a_step

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2453)

::: pybmad.rotate_spin_a_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### rotate_spin_given_field

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2462)

::: pybmad.rotate_spin_given_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### s_body_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2470)

::: pybmad.s_body_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### s_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2478)

::: pybmad.s_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sad_mult_hard_bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/fringe_mod.f90#L447)

::: pybmad.sad_mult_hard_bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sad_soft_bend_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/fringe_mod.f90#L892)

::: pybmad.sad_soft_bend_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### save_a_beam_step

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2484)

::: pybmad.save_a_beam_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### save_a_bunch_step

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2494)

::: pybmad.save_a_bunch_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### save_a_step

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2504)

::: pybmad.save_a_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sbend_body_with_k1_map

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2518)

::: pybmad.sbend_body_with_k1_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sc_adaptive_step

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/space_charge_mod.f90#L290)

::: pybmad.sc_adaptive_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sc_step

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/space_charge_mod.f90#L220)

::: pybmad.sc_step
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_active_fixer

Fortran source: [`bmad/modules/fixer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/fixer_mod.f90#L32)

::: pybmad.set_active_fixer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_custom_attribute_name

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L2863)

::: pybmad.set_custom_attribute_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_attribute

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2530)

::: pybmad.set_ele_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_defaults

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2642)

::: pybmad.set_ele_defaults
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_name

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2540)

::: pybmad.set_ele_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_real_attribute

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2547)

::: pybmad.set_ele_real_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ele_status_stale

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2557)

::: pybmad.set_ele_status_stale
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_flags_for_changed_attribute

Fortran sources (overloaded):

- `set_flags_for_changed_integer_attribute`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L335)
- `set_flags_for_changed_logical_attribute`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L343)
- `set_flags_for_changed_lat_attribute`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L351)
- `set_flags_for_changed_real_attribute`: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L358)

::: pybmad.set_flags_for_changed_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_fringe_on_off

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2566)

::: pybmad.set_fringe_on_off
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_lords_status_stale

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2573)

::: pybmad.set_lords_status_stale
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_on_off

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2582)

::: pybmad.set_on_off
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_orbit_to_zero

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2594)

::: pybmad.set_orbit_to_zero
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ptc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2602)

::: pybmad.set_ptc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ptc_base_state

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2610)

::: pybmad.set_ptc_base_state
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ptc_com_pointers

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L1005)

::: pybmad.set_ptc_com_pointers
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ptc_quiet

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L3365)

::: pybmad.set_ptc_quiet
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_ptc_verbose

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L1176)

::: pybmad.set_ptc_verbose
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_pwd_ele

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/longitudinal_profile_mod.f90#L547)

::: pybmad.set_pwd_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_status_flags

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2617)

::: pybmad.set_status_flags
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_tune

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2649)

::: pybmad.set_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2624)

::: pybmad.set_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### set_z_tune

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2634)

::: pybmad.set_z_tune
    options:
      show_root_heading: false
      show_root_toc_entry: false

### settable_dep_var_bookkeeping

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L4979)

::: pybmad.settable_dep_var_bookkeeping
    options:
      show_root_heading: false
      show_root_toc_entry: false

### setup_high_energy_space_charge_calc

Fortran source: [`bmad/space_charge/high_energy_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/high_energy_space_charge_mod.f90#L32)

::: pybmad.setup_high_energy_space_charge_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sigma_mat_ptc_to_bmad

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L1327)

::: pybmad.sigma_mat_ptc_to_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### significant_difference

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2673)

::: pybmad.significant_difference
    options:
      show_root_heading: false
      show_root_toc_entry: false

### skip_ele_blender

Fortran source: [`bmad/interface/blender_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/blender_interface_mod.f90#L87)

::: pybmad.skip_ele_blender
    options:
      show_root_heading: false
      show_root_toc_entry: false

### slice_lattice

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2681)

::: pybmad.slice_lattice
    options:
      show_root_heading: false
      show_root_toc_entry: false

### soft_quadrupole_edge_kick

Fortran source: [`bmad/modules/fringe_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/fringe_mod.f90#L555)

::: pybmad.soft_quadrupole_edge_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sol_quad_mat6_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2690)

::: pybmad.sol_quad_mat6_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### solve_psi_adaptive

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/longitudinal_profile_mod.f90#L168)

::: pybmad.solve_psi_adaptive
    options:
      show_root_heading: false
      show_root_toc_entry: false

### solve_psi_fixed_steps

Fortran source: [`bmad/multiparticle/longitudinal_profile_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/longitudinal_profile_mod.f90#L227)

::: pybmad.solve_psi_fixed_steps
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sort_complex_taylor_terms

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L604)

::: pybmad.sort_complex_taylor_terms
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_dn_dpz_from_mat8

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2731)

::: pybmad.spin_dn_dpz_from_mat8
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_dn_dpz_from_qmap

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2739)

::: pybmad.spin_dn_dpz_from_qmap
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_map1_normalize

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2748)

::: pybmad.spin_map1_normalize
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_mat8_resonance_strengths

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2762)

::: pybmad.spin_mat8_resonance_strengths
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_mat_to_eigen

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2754)

::: pybmad.spin_mat_to_eigen
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_omega

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2769)

::: pybmad.spin_omega
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_quat_resonance_strengths

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2779)

::: pybmad.spin_quat_resonance_strengths
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spin_taylor_to_linear

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2786)

::: pybmad.spin_taylor_to_linear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spinor_to_polar

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2794)

::: pybmad.spinor_to_polar
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spinor_to_vec

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2801)

::: pybmad.spinor_to_vec
    options:
      show_root_heading: false
      show_root_toc_entry: false

### spline_fit_orbit

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2808)

::: pybmad.spline_fit_orbit
    options:
      show_root_heading: false
      show_root_toc_entry: false

### split_expression_string

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L1907)

::: pybmad.split_expression_string
    options:
      show_root_heading: false
      show_root_toc_entry: false

### split_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2815)

::: pybmad.split_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sprint_spin_taylor_map

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2828)

::: pybmad.sprint_spin_taylor_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sr_longitudinal_wake_particle

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/wake_mod.f90#L298)

::: pybmad.sr_longitudinal_wake_particle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sr_transverse_wake_particle

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/wake_mod.f90#L394)

::: pybmad.sr_transverse_wake_particle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### sr_z_long_wake

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/wake_mod.f90#L501)

::: pybmad.sr_z_long_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### srdt_calc

Fortran source: [`bmad/modules/srdt_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/srdt_mod.f90#L81)

::: pybmad.srdt_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### srdt_lsq_solution

Fortran source: [`bmad/modules/srdt_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/srdt_mod.f90#L656)

::: pybmad.srdt_lsq_solution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### start_branch_at

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2835)

::: pybmad.start_branch_at
    options:
      show_root_heading: false
      show_root_toc_entry: false

### stream_ele_end

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2843)

::: pybmad.stream_ele_end
    options:
      show_root_heading: false
      show_root_toc_entry: false

### string_attrib

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L2350)

::: pybmad.string_attrib
    options:
      show_root_heading: false
      show_root_toc_entry: false

### strong_beam_sigma_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2849)

::: pybmad.strong_beam_sigma_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### strong_beam_strength

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2856)

::: pybmad.strong_beam_strength
    options:
      show_root_heading: false
      show_root_toc_entry: false

### surface_grid_displacement

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_utils_mod.f90#L186)

::: pybmad.surface_grid_displacement
    options:
      show_root_heading: false
      show_root_toc_entry: false

### switch_attrib_value_name

Fortran source: [`bmad/modules/attribute_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/attribute_mod.f90#L2431)

::: pybmad.switch_attrib_value_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### symp_lie_bmad

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2863)

::: pybmad.symp_lie_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### t6_to_b123

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L49)

::: pybmad.t6_to_b123
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taper_mag_strengths

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2874)

::: pybmad.taper_mag_strengths
    options:
      show_root_heading: false
      show_root_toc_entry: false

### target_min_max_calc

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L1190)

::: pybmad.target_min_max_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### target_rot_mats

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L1141)

::: pybmad.target_rot_mats
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taylor_equal_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4919)

::: pybmad.taylor_equal_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taylor_inverse

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L1981)

::: pybmad.taylor_inverse
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taylor_propagate1

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L2490)

::: pybmad.taylor_propagate1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taylor_to_mad_map

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1785)

::: pybmad.taylor_to_mad_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### taylors_equal_taylors

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L4956)

::: pybmad.taylors_equal_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tilt_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2883)

::: pybmad.tilt_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tilt_coords_photon

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2892)

::: pybmad.tilt_coords_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tilt_mat6

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2899)

::: pybmad.tilt_mat6
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_eta_reading

Fortran source: [`bmad/modules/measurement_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/measurement_mod.f90#L209)

::: pybmad.to_eta_reading
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_fieldmap_coords

Fortran source: [`bmad/modules/em_field_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/em_field_mod.f90#L101)

::: pybmad.to_fieldmap_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_orbit_reading

Fortran source: [`bmad/modules/measurement_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/measurement_mod.f90#L127)

::: pybmad.to_orbit_reading
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_phase_and_coupling_reading

Fortran source: [`bmad/modules/measurement_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/measurement_mod.f90#L288)

::: pybmad.to_phase_and_coupling_reading
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_photon_angle_coords

Fortran source: [`bmad/photon/photon_target_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_target_mod.f90#L355)

::: pybmad.to_photon_angle_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### to_surface_coords

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2905)

::: pybmad.to_surface_coords
    options:
      show_root_heading: false
      show_root_toc_entry: false

### touschek_lifetime

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/touschek_mod.f90#L82)

::: pybmad.touschek_lifetime
    options:
      show_root_heading: false
      show_root_toc_entry: false

### touschek_rate1

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/touschek_mod.f90#L427)

::: pybmad.touschek_rate1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### touschek_rate1_zap

Fortran source: [`bmad/multiparticle/touschek_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/touschek_mod.f90#L269)

::: pybmad.touschek_rate1_zap
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3170)

::: pybmad.track1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_beam

Fortran source: [`bmad/multiparticle/beam_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_mod.f90#L193)

::: pybmad.track1_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bmad

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3183)

::: pybmad.track1_bmad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bmad_photon

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3195)

::: pybmad.track1_bmad_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bunch

Fortran source: [`bmad/multiparticle/beam_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_mod.f90#L238)

::: pybmad.track1_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bunch_csr

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/csr_and_space_charge_mod.f90#L129)

::: pybmad.track1_bunch_csr
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bunch_csr3d

Fortran source: [`bmad/space_charge/csr_and_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/csr_and_space_charge_mod.f90#L1862)

::: pybmad.track1_bunch_csr3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bunch_hom

Fortran source: [`bmad/multiparticle/beam_utils.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_utils.f90#L32)

::: pybmad.track1_bunch_hom
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_bunch_space_charge

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3204)

::: pybmad.track1_bunch_space_charge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_crystal

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L845)

::: pybmad.track1_crystal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_diffraction_plate_or_mask

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L145)

::: pybmad.track1_diffraction_plate_or_mask
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_high_energy_space_charge

Fortran source: [`bmad/space_charge/high_energy_space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/high_energy_space_charge_mod.f90#L169)

::: pybmad.track1_high_energy_space_charge
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_lens

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L27)

::: pybmad.track1_lens
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_linear

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3214)

::: pybmad.track1_linear
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_lr_wake

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/wake_mod.f90#L110)

::: pybmad.track1_lr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_mad

Fortran source: [`bmad/modules/mad_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mad_mod.f90#L1597)

::: pybmad.track1_mad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_mirror

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L435)

::: pybmad.track1_mirror
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_mosaic_crystal

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L618)

::: pybmad.track1_mosaic_crystal
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_multilayer_mirror

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L482)

::: pybmad.track1_multilayer_mirror
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_radiation

Fortran source: [`bmad/modules/radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/radiation_mod.f90#L61)

::: pybmad.track1_radiation
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_radiation_center

Fortran source: [`bmad/modules/radiation_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/radiation_mod.f90#L271)

::: pybmad.track1_radiation_center
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_runge_kutta

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3222)

::: pybmad.track1_runge_kutta
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_sample

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L235)

::: pybmad.track1_sample
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_spin

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3234)

::: pybmad.track1_spin
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_spin_integration

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3243)

::: pybmad.track1_spin_integration
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_spin_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3252)

::: pybmad.track1_spin_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_sr_wake

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/wake_mod.f90#L703)

::: pybmad.track1_sr_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_symp_lie_ptc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3260)

::: pybmad.track1_symp_lie_ptc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3269)

::: pybmad.track1_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track1_time_runge_kutta

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3279)

::: pybmad.track1_time_runge_kutta
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_beambeam

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2912)

::: pybmad.track_a_beambeam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_bend

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2923)

::: pybmad.track_a_bend
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_bend_photon

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L1246)

::: pybmad.track_a_bend_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_capillary

Fortran source: [`bmad/photon/capillary_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/capillary_mod.f90#L33)

::: pybmad.track_a_capillary
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_converter

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2933)

::: pybmad.track_a_converter
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_crab_cavity

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2943)

::: pybmad.track_a_crab_cavity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_drift

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2953)

::: pybmad.track_a_drift
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_drift_photon

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2963)

::: pybmad.track_a_drift_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_foil

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3080)

::: pybmad.track_a_foil
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_gkicker

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2971)

::: pybmad.track_a_gkicker
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_lcavity

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2981)

::: pybmad.track_a_lcavity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_lcavity_old

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L2991)

::: pybmad.track_a_lcavity_old
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_mask

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3001)

::: pybmad.track_a_mask
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_match

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3011)

::: pybmad.track_a_match
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_patch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3031)

::: pybmad.track_a_patch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_patch_photon

Fortran source: [`bmad/photon/track1_photon_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/track1_photon_mod.f90#L65)

::: pybmad.track_a_patch_photon
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_pickup

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3021)

::: pybmad.track_a_pickup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_quadrupole

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3040)

::: pybmad.track_a_quadrupole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_rfcavity

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3050)

::: pybmad.track_a_rfcavity
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_sad_mult

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3060)

::: pybmad.track_a_sad_mult
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_sol_quad

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3070)

::: pybmad.track_a_sol_quad
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_thick_multipole

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3090)

::: pybmad.track_a_thick_multipole
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_wiggler

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3100)

::: pybmad.track_a_wiggler
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_a_zero_length_element

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3110)

::: pybmad.track_a_zero_length_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_all

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3120)

::: pybmad.track_all
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_beam

Fortran source: [`bmad/multiparticle/beam_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_mod.f90#L40)

::: pybmad.track_beam
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_bunch

Fortran source: [`bmad/multiparticle/beam_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_mod.f90#L103)

::: pybmad.track_bunch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_bunch_time

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3130)

::: pybmad.track_bunch_time
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_bunch_to_s

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/space_charge_mod.f90#L429)

::: pybmad.track_bunch_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_bunch_to_t

Fortran source: [`bmad/space_charge/space_charge_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/space_charge_mod.f90#L497)

::: pybmad.track_bunch_to_t
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_complex_taylor

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L828)

::: pybmad.track_complex_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_from_s_to_s

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3140)

::: pybmad.track_from_s_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_many

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3150)

::: pybmad.track_many
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_to_surface

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3161)

::: pybmad.track_to_surface
    options:
      show_root_heading: false
      show_root_toc_entry: false

### track_until_dead

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/time_tracker_mod.f90#L1141)

::: pybmad.track_until_dead
    options:
      show_root_heading: false
      show_root_toc_entry: false

### tracking_rad_map_setup

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3290)

::: pybmad.tracking_rad_map_setup
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_ac_kick

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3300)

::: pybmad.transfer_ac_kick
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_branch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3306)

::: pybmad.transfer_branch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_branch_parameters

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3313)

::: pybmad.transfer_branch_parameters
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_branches

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3320)

::: pybmad.transfer_branches
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3327)

::: pybmad.transfer_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_ele_taylor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3335)

::: pybmad.transfer_ele_taylor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_eles

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3342)

::: pybmad.transfer_eles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_fieldmap

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3349)

::: pybmad.transfer_fieldmap
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_fixer_params

Fortran source: [`bmad/modules/fixer_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/fixer_mod.f90#L113)

::: pybmad.transfer_fixer_params
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_lat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3356)

::: pybmad.transfer_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_lat_parameters

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3363)

::: pybmad.transfer_lat_parameters
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_map_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3370)

::: pybmad.transfer_map_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_map_from_s_to_s

Fortran source: [`bmad/modules/transfer_map_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/transfer_map_mod.f90#L59)

::: pybmad.transfer_map_from_s_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_mat2_from_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3391)

::: pybmad.transfer_mat2_from_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_mat_from_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3383)

::: pybmad.transfer_mat_from_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_matrix_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3398)

::: pybmad.transfer_matrix_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3408)

::: pybmad.transfer_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### transfer_wake

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3415)

::: pybmad.transfer_wake
    options:
      show_root_heading: false
      show_root_toc_entry: false

### truncate_complex_taylor_to_order

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L900)

::: pybmad.truncate_complex_taylor_to_order
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss1_propagate

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3470)

::: pybmad.twiss1_propagate
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss3_at_start

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L1384)

::: pybmad.twiss3_at_start
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss3_from_twiss2

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L1332)

::: pybmad.twiss3_from_twiss2
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss3_propagate1

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L1265)

::: pybmad.twiss3_propagate1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss3_propagate_all

Fortran source: [`bmad/modules/mode3_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/mode3_mod.f90#L1236)

::: pybmad.twiss3_propagate_all
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_and_track

Fortran sources (overloaded):

- `twiss_and_track_branch`: [`bmad/modules/twiss_and_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/twiss_and_track_mod.f90#L88)
- `twiss_and_track_all`: [`bmad/modules/twiss_and_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/twiss_and_track_mod.f90#L133)

::: pybmad.twiss_and_track
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_and_track_at_s

Fortran source: [`bmad/modules/twiss_and_track_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/twiss_and_track_mod.f90#L347)

::: pybmad.twiss_and_track_at_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_and_track_from_s_to_s

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3427)

::: pybmad.twiss_and_track_from_s_to_s
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_and_track_intra_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3439)

::: pybmad.twiss_and_track_intra_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_at_element

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3452)

::: pybmad.twiss_at_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_at_start

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3461)

::: pybmad.twiss_at_start
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_from_tracking

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3498)

::: pybmad.twiss_from_tracking
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_propagate1

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3508)

::: pybmad.twiss_propagate1
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_propagate_all

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3515)

::: pybmad.twiss_propagate_all
    options:
      show_root_heading: false
      show_root_toc_entry: false

### twiss_to_1_turn_mat

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3523)

::: pybmad.twiss_to_1_turn_mat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_complex_taylors

Fortran source: [`bmad/modules/complex_taylor_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/complex_taylor_mod.f90#L281)

::: pybmad.type_complex_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_coord

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3530)

::: pybmad.type_coord
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3536)

::: pybmad.type_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_end_stuff

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L909)

::: pybmad.type_end_stuff
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_expression_tree

Fortran source: [`bmad/modules/expression_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/expression_mod.f90#L683)

::: pybmad.type_expression_tree
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_ptc_fibre

Fortran source: [`bmad/ptc/ptc_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_interface_mod.f90#L385)

::: pybmad.type_ptc_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_ptc_layout

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L28)

::: pybmad.type_ptc_layout
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_taylors

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3549)

::: pybmad.type_taylors
    options:
      show_root_heading: false
      show_root_toc_entry: false

### type_twiss

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3560)

::: pybmad.type_twiss
    options:
      show_root_heading: false
      show_root_toc_entry: false

### update_ele_from_fibre

Fortran source: [`bmad/ptc/ptc_layout_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/ptc/ptc_layout_mod.f90#L1204)

::: pybmad.update_ele_from_fibre
    options:
      show_root_heading: false
      show_root_toc_entry: false

### update_fibre_from_ele

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3592)

::: pybmad.update_fibre_from_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### update_floor_angles

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3585)

::: pybmad.update_floor_angles
    options:
      show_root_heading: false
      show_root_toc_entry: false

### valid_field_calc

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3599)

::: pybmad.valid_field_calc
    options:
      show_root_heading: false
      show_root_toc_entry: false

### valid_fringe_type

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3607)

::: pybmad.valid_fringe_type
    options:
      show_root_heading: false
      show_root_toc_entry: false

### valid_mat6_calc_method

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3615)

::: pybmad.valid_mat6_calc_method
    options:
      show_root_heading: false
      show_root_toc_entry: false

### valid_spin_tracking_method

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3623)

::: pybmad.valid_spin_tracking_method
    options:
      show_root_heading: false
      show_root_toc_entry: false

### valid_tracking_method

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3631)

::: pybmad.valid_tracking_method
    options:
      show_root_heading: false
      show_root_toc_entry: false

### value_of_attribute

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3639)

::: pybmad.value_of_attribute
    options:
      show_root_heading: false
      show_root_toc_entry: false

### value_to_line

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/output/write_lattice_file_mod.f90#L434)

::: pybmad.value_to_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### vec_to_polar

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3650)

::: pybmad.vec_to_polar
    options:
      show_root_heading: false
      show_root_toc_entry: false

### vec_to_spinor

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3658)

::: pybmad.vec_to_spinor
    options:
      show_root_heading: false
      show_root_toc_entry: false

### verify_valid_name

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L2582)

::: pybmad.verify_valid_name
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_for_bend_angle

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3666)

::: pybmad.w_mat_for_bend_angle
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_for_tilt

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3689)

::: pybmad.w_mat_for_tilt
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_for_x_pitch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3673)

::: pybmad.w_mat_for_x_pitch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### w_mat_for_y_pitch

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3681)

::: pybmad.w_mat_for_y_pitch
    options:
      show_root_heading: false
      show_root_toc_entry: false

### wall3d_d_radius

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/wall3d_mod.f90#L659)

::: pybmad.wall3d_d_radius
    options:
      show_root_heading: false
      show_root_toc_entry: false

### wall3d_initializer

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/wall3d_mod.f90#L135)

::: pybmad.wall3d_initializer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### wall3d_section_initializer

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/wall3d_mod.f90#L207)

::: pybmad.wall3d_section_initializer
    options:
      show_root_heading: false
      show_root_toc_entry: false

### wall3d_to_position

Fortran source: [`bmad/modules/wall3d_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/wall3d_mod.f90#L1153)

::: pybmad.wall3d_to_position
    options:
      show_root_heading: false
      show_root_toc_entry: false

### word_to_value

Fortran source: [`bmad/parsing/bmad_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/bmad_parser_mod.f90#L1076)

::: pybmad.word_to_value
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_ascii_beam_file

Fortran source: [`bmad/multiparticle/beam_file_io.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_file_io.f90#L123)

::: pybmad.write_ascii_beam_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_astra_bend

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/astra_interface_mod.f90#L344)

::: pybmad.write_astra_bend
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_astra_field_grid_file

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/astra_interface_mod.f90#L504)

::: pybmad.write_astra_field_grid_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_astra_field_grid_file_3d

Fortran source: [`bmad/interface/astra_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/astra_interface_mod.f90#L718)

::: pybmad.write_astra_field_grid_file_3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_beam_file

Fortran source: [`bmad/multiparticle/beam_file_io.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/beam_file_io.f90#L29)

::: pybmad.write_beam_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_beam_floor_positions

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3697)

::: pybmad.write_beam_floor_positions
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_binary_cartesian_map

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/binary_parser_mod.f90#L25)

::: pybmad.write_binary_cartesian_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_binary_cylindrical_map

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/binary_parser_mod.f90#L149)

::: pybmad.write_binary_cylindrical_map
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_binary_grid_field

Fortran source: [`bmad/parsing/binary_parser_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/parsing/binary_parser_mod.f90#L275)

::: pybmad.write_binary_grid_field
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_blender_ele

Fortran source: [`bmad/interface/blender_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/blender_interface_mod.f90#L112)

::: pybmad.write_blender_ele
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_blender_lat_layout

Fortran source: [`bmad/interface/blender_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/blender_interface_mod.f90#L15)

::: pybmad.write_blender_lat_layout
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_bmad_lattice_file

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3717)

::: pybmad.write_bmad_lattice_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_gpt_field_grid_file_1d

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/gpt_interface_mod.f90#L797)

::: pybmad.write_gpt_field_grid_file_1d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_gpt_field_grid_file_2d

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/gpt_interface_mod.f90#L1008)

::: pybmad.write_gpt_field_grid_file_2d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_gpt_field_grid_file_3d

Fortran source: [`bmad/interface/gpt_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/gpt_interface_mod.f90#L1264)

::: pybmad.write_gpt_field_grid_file_3d
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lat_line

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/output/write_lattice_file_mod.f90#L331)

::: pybmad.write_lat_line
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_elegant_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3727)

::: pybmad.write_lattice_elegant_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_foreign_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3739)

::: pybmad.write_lattice_foreign_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_mad_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3751)

::: pybmad.write_lattice_mad_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_pals_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3763)

::: pybmad.write_lattice_pals_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_sad_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3771)

::: pybmad.write_lattice_sad_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_lattice_scibmad_format

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3780)

::: pybmad.write_lattice_scibmad_format
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_line_element

Fortran source: [`bmad/output/write_lattice_file_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/output/write_lattice_file_mod.f90#L142)

::: pybmad.write_line_element
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_opal_field_grid_file

Fortran source: [`bmad/interface/opal_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/opal_interface_mod.f90#L419)

::: pybmad.write_opal_field_grid_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_opal_lattice_file

Fortran source: [`bmad/interface/opal_interface_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/opal_interface_mod.f90#L26)

::: pybmad.write_opal_lattice_file
    options:
      show_root_heading: false
      show_root_toc_entry: false

### write_time_particle_distribution

Fortran source: [`bmad/modules/time_tracker_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/time_tracker_mod.f90#L963)

::: pybmad.write_time_particle_distribution
    options:
      show_root_heading: false
      show_root_toc_entry: false

### xlafun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L377)

::: pybmad.xlafun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### xraylib_nist_compound

Fortran source: [`bmad/interface/xraylib_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/interface/xraylib_interface.f90#L457)

::: pybmad.xraylib_nist_compound
    options:
      show_root_heading: false
      show_root_toc_entry: false

### ylafun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L397)

::: pybmad.ylafun
    options:
      show_root_heading: false
      show_root_toc_entry: false

### z_at_surface

Fortran source: [`bmad/photon/photon_utils_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/photon/photon_utils_mod.f90#L96)

::: pybmad.z_at_surface
    options:
      show_root_heading: false
      show_root_toc_entry: false

### zero_ele_kicks

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3799)

::: pybmad.zero_ele_kicks
    options:
      show_root_heading: false
      show_root_toc_entry: false

### zero_ele_offsets

Fortran source: [`bmad/modules/bmad_routine_interface.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/modules/bmad_routine_interface.f90#L3805)

::: pybmad.zero_ele_offsets
    options:
      show_root_heading: false
      show_root_toc_entry: false

### zero_lr_wakes_in_lat

Fortran source: [`bmad/multiparticle/wake_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/multiparticle/wake_mod.f90#L75)

::: pybmad.zero_lr_wakes_in_lat
    options:
      show_root_heading: false
      show_root_toc_entry: false

### zlafun

Fortran source: [`bmad/space_charge/open_spacecharge_core_mod.f90`](https://github.com/bmad-sim/bmad-ecosystem/blob/4e3330c0a436e793938d2510a585f4f322235bac/bmad/space_charge/open_spacecharge_core_mod.f90#L412)

::: pybmad.zlafun
    options:
      show_root_heading: false
      show_root_toc_entry: false
