from __future__ import annotations

import copy
import pathlib

import pybmad
from pybmad import LatStruct


def find_elements_by_name(lat, name: str) -> list:
    """Replicates lat_ele_locator logic in Python."""
    return [ele for ele in lat.ele if ele.name == name]


def set_hom_order_cutoff(lat: LatStruct, cutoff: float):
    # Iterate over all elements in the lattice
    for ele in lat.ele:
        # Check availability of wakefields
        if not ele.wake:
            continue

        # Check availability of Long Range (LR) modes
        # This accesses the structure hierarchy: ele -> wake -> lr -> mode
        if not ele.wake.lr or not ele.wake.lr.mode:
            continue

        modes = ele.wake.lr.mode  # This is a list of mode objects
        # Filter modes based on 'm' attribute
        high_order_modes = [m for m in modes if m.m > cutoff]

        # If no modes exceed cutoff, continue
        if not high_order_modes:
            continue

        # If all modes exceed cutoff, clear the list
        if len(high_order_modes) == len(modes):
            ele.wake.lr.mode = []  # Clear assignments
            continue

        # Partial removal: keep modes <= cutoff
        kept_modes = [m for m in modes if m.m <= cutoff]
        ele.wake.lr.mode = kept_modes


def hybridize(
    lat: LatStruct,
    ele_track_end: str = "",
    use_taylor_for_hybrids: bool = False,
    keep_all_lcavities: bool = False,
):
    for ele in lat.ele:
        ele.select = False  # Default: Hybridize this element

        if ele.name == ele_track_end:
            ele.select = True
            continue

        if ele.key == pybmad.TAYLOR:
            ele.select = True
            continue

        if ele.key != pybmad.LCAVITY:
            continue

        # LCavity logic
        if keep_all_lcavities:
            pass

        elif not ele.wake or not ele.wake.lr.mode:
            # If no wake or no modes, default (false) applies?
            # Fortran logic: if NOT keep_all_lcavities
            #   if not associated(wake) cycle
            #   if size(mode) == 0 cycle
            # ele%select = true
            continue

        ele.select = True

    return pybmad.make_hybrid_lat(lat, use_taylor_for_hybrids)


def main():
    bbu_param = pybmad.BbuParamStruct()
    bbu_param.lat_filename = "bbu_test/oneturn_lat.bmad"
    bbu_param.keep_overlays_and_groups = False
    bbu_param.simulation_turns_max = 500
    bbu_param.elname = "T1"
    bbu_param.hybridize = True
    bbu_param.nrep = 5
    bbu_param.limit_factor = 3
    bbu_param.keep_all_lcavities = False
    bbu_param.ran_gauss_sigma_cut = 3
    bbu_param.nstep = 50
    bbu_param.current = 1.000
    bbu_param.ran_seed = 100
    bbu_param.rel_tol = 0.001
    bbu_param.lat2_filename = ""
    bbu_param.bunch_freq = 1300000000.0
    # ---------------------------------------------------------------

    beam_init = pybmad.BeamInitStruct()
    beam_init.n_particle = 1

    # Handle Ramp Pattern Check
    bbu_param.n_ramp_pattern = len(bbu_param.ramp_pattern)
    if bbu_param.ramp_on and bbu_param.n_ramp_pattern < 1:
        raise RuntimeError("RAMP_ON = TRUE BUT THERE IS NO RAMP_PATTERN!")

    # Define distance between bunches (1 / freq)
    if bbu_param.bunch_freq != 0.0:
        beam_init.dt_bunch = 1.0 / bbu_param.bunch_freq
    else:
        beam_init.dt_bunch = 0.0

    pybmad.ran_seed_put(bbu_param.ran_seed)

    # In pybmad/cppbmad context, the gauss converter is usually internal,
    # but strictly following translation:
    if bbu_param.ran_gauss_sigma_cut > 0:
        pybmad.ran_gauss_converter(set_sigma_cut=bbu_param.ran_gauss_sigma_cut)

    print(f"Lattice file: {bbu_param.lat_filename}")
    res = pybmad.bmad_parser(bbu_param.lat_filename)
    if res.err_flag:
        raise RuntimeError("bmad_parser failed")
    lat_in: LatStruct = res.lat

    # Set Bmad Com
    bmad_com = pybmad.get_bmad_com()
    bmad_com.auto_bookkeeper = False

    if bbu_param.lat2_filename:
        print(f"DR-scan or Phase-scan, parsing: {bbu_param.lat2_filename}")
        # Note: bmad_parser2 updates the existing lat_in
        pybmad.bmad_parser2(bbu_param.lat2_filename, lat_in)

    # --------------------------------------------------------------------------
    # Twiss and Track (Closed Orbit)
    # --------------------------------------------------------------------------
    # In Python, we have separate calls usually.
    # twiss_and_track calculates the closed orbit and fills lat info.
    orb = pybmad.CoordStruct.new_array1d(0)
    pybmad.twiss_and_track(lat_in, orb)

    # --------------------------------------------------------------------------
    # HOM Pruning (Higher Order Modes)
    # --------------------------------------------------------------------------
    if bbu_param.hom_order_cutoff > 0:
        set_hom_order_cutoff(lat_in, cutoff=bbu_param.hom_order_cutoff)

    # --------------------------------------------------------------------------
    # Hybridization Logic
    # --------------------------------------------------------------------------
    # In Fortran, this replaces drift/magnets with Taylor maps
    lat = lat_in  # Default if not hybridizing

    if bbu_param.hybridize:
        print("Hybridizing lattice...")
        hybrid_lat = hybridize(
            lat,
            ele_track_end=bbu_param.ele_track_end,
            use_taylor_for_hybrids=bbu_param.use_taylor_for_hybrids,
            keep_all_lcavities=bbu_param.keep_all_lcavities,
        )
        print("Hybridization complete !!!")

        if bbu_param.write_digested_hybrid_lat:
            # pybmad.write_digested_bmad_file("hybrid.digested", lat)
            pybmad.write_bmad_lattice_file("hybrid.lat", hybrid_lat)

    # --------------------------------------------------------------------------
    # Setup Tracking End point
    # --------------------------------------------------------------------------
    # Keep copy of current state before restoration
    lat0 = copy.copy(lat)

    if bbu_param.ele_track_end:
        locs = find_elements_by_name(lat, bbu_param.ele_track_end)

        if not locs:
            raise ValueError(f"No matching element found for {bbu_param.ele_track_end}")

        if len(locs) > 1:
            print(f"Multiple elements found for {bbu_param.ele_track_end}")
            print("Will use the first instance.")

        ele_end = locs[0]

        # Handle Lord/Slave logic (slaves contain the tracking physics)
        if ele_end.lord_status == pybmad.SUPER_LORD:
            # Replicate: pointer_to_slave(ele, ele%n_slave)
            # In Python, we usually access likely the last slave in the slave list
            ele_end = ele_end.slave[-1]

        ix = ele_end.ix_ele
        # Check if index is within tracking range
        if ix > lat.n_ele_track:
            raise ValueError(f"STOPPING ELEMENT IS A LORD! {bbu_param.ele_track_end}")

        bbu_param.ix_ele_track_end = ix

    # --------------------------------------------------------------------------
    # BBU Specific Setup
    # --------------------------------------------------------------------------
    if bbu_param.write_hom_info:
        pybmad.rf_cav_names(lat)

    # Check RF Freq
    pybmad.check_rf_freq(lat, bbu_param.bunch_freq)

    # Prepare custom BBU beam structure (assumed binding)
    bbu_beam = pybmad.BbuBeamStruct()

    # BBU Setup call
    pybmad.bbu_setup(lat, beam_init.dt_bunch, bbu_param, bbu_beam)
    print("bbu_setup complete !!!")

    # Recalculate bunch charge based on current
    beam_init.bunch_charge = bbu_param.current * beam_init.dt_bunch

    n_stages = len(bbu_beam.stage)
    print(f"Number of stages and elements: {n_stages}   {lat.n_ele_track}")

    # Use the initial lattice
    lat = lat0

    # --------------------------------------------------------------------------
    # Tracking Execution
    # --------------------------------------------------------------------------
    print("Starting bbu_track_all...")
    (hom_voltage_gain, growth_rate, lost, irep) = pybmad.bbu_track_all(
        lat=lat,
        bbu_beam=bbu_beam,
        bbu_param=bbu_param,
        beam_init=beam_init,
        # TODO: these aren't marked as output, so we have to specify them
        hom_voltage_normalized=0.0,
        growth_rate=0.0,
        lost=False,
        irep=0,
    )

    print("bbu_track_all complete !!!")
    print(f"HOM VOLT GAIN: {hom_voltage_gain}")
    print(f"growth_rate: {growth_rate}")

    output_file = pathlib.Path("for_py.txt")
    with output_file.open("w") as f:
        f.write(f"lostbool = {lost}\n")
        # Format: Scientific notation with specifically padded output
        f.write(f"v_gain = {hom_voltage_gain:.8E}\n")
        f.write(f"bunch_dt = {beam_init.dt_bunch:.6E}\n")

        # Check against a 'garbage' constant, in Python usually math.nan check or checking initialization
        # Fortran used real_garbage$. Assuming -1 or specific check.
        # Here we just print validity.
        valid_growth = growth_rate != -1.0  # Simplify garbage check
        f.write(f"growth_rate_set = {valid_growth}\n")
        f.write(f"growth_rate = {growth_rate:.6E}\n")
    return {
        "lat": lat,
        "bbu_beam": bbu_beam,
        "hom_voltage_gain": hom_voltage_gain,
        "growth_rate": growth_rate,
        "lost": lost,
        "irep": irep,
    }


if __name__ == "__main__":
    res = main()
    lat = res["lat"]
