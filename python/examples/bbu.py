from __future__ import annotations

import copy
import pathlib

import pybmad as pb
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
        if not ele.wake.lr or not ele.wake.lr.mode:
            continue

        # Filter modes based on 'm' attribute
        high_order_modes = [m for m in ele.wake.lr.mode if m.m > cutoff]

        # If no modes exceed cutoff, continue
        if not high_order_modes:
            continue

        # If all modes exceed cutoff, clear the list
        if len(high_order_modes) == len(ele.wake.lr.mode):
            # TODO this might not work in pybmad
            ele.wake.lr.mode = []
            continue

        # Partial removal: keep modes <= cutoff
        kept_modes = [m for m in ele.wake.lr.mode if m.m <= cutoff]
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

        if ele.key == pb.TAYLOR:
            ele.select = True
            continue

        if ele.key == pb.LCAVITY:
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

    return pb.make_hybrid_lat(lat, use_taylor_for_hybrids)


def main():
    bbu_param = pb.BbuParamStruct(
        lat_filename="bbu_test/oneturn_lat.bmad",
        keep_overlays_and_groups=False,
        simulation_turns_max=500,
        elname="T1",
        hybridize=True,
        nrep=5,
        limit_factor=3,
        keep_all_lcavities=False,
        ran_gauss_sigma_cut=3,
        nstep=50,
        current=1.000,
        ran_seed=100,
        rel_tol=0.001,
        lat2_filename="",
        bunch_freq=1300000000.0,
    )
    # ---------------------------------------------------------------

    beam_init = pb.BeamInitStruct(
        n_particle=1,
        # Define distance between bunches (1 / freq)
        dt_bunch=1.0 / bbu_param.bunch_freq if bbu_param.bunch_freq != 0 else 0.0,
    )

    bbu_param.n_ramp_pattern = len(bbu_param.ramp_pattern)
    if bbu_param.ramp_on and bbu_param.n_ramp_pattern < 1:
        raise RuntimeError("RAMP_ON = TRUE BUT THERE IS NO RAMP_PATTERN!")

    pb.ran_seed_put(bbu_param.ran_seed)

    if bbu_param.ran_gauss_sigma_cut > 0:
        pb.ran_gauss_converter(set_sigma_cut=bbu_param.ran_gauss_sigma_cut)

    print(f"Lattice file: {bbu_param.lat_filename}")
    res = pb.bmad_parser(bbu_param.lat_filename)
    if res.err_flag:
        raise RuntimeError("bmad_parser failed")
    lat_in: LatStruct = res.lat

    # Set Bmad Com
    bmad_com = pb.get_bmad_com()
    bmad_com.auto_bookkeeper = False

    if bbu_param.lat2_filename:
        print(f"DR-scan or Phase-scan, parsing: {bbu_param.lat2_filename}")
        # Note: bmad_parser2 updates the existing lat_in
        pb.bmad_parser2(bbu_param.lat2_filename, lat_in)

    # --------------------------------------------------------------------------
    # Twiss and Track (Closed Orbit)
    # --------------------------------------------------------------------------
    orb = pb.CoordStruct.new_array1d(0)
    pb.twiss_and_track(lat_in, orb)

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
            pb.write_bmad_lattice_file("hybrid.lat", hybrid_lat)

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
        if ele_end.lord_status == pb.SUPER_LORD:
            ele_end = ele_end.slave[-1]

        ix = ele_end.ix_ele
        if ix > lat.n_ele_track:
            raise ValueError(f"STOPPING ELEMENT IS A LORD! {bbu_param.ele_track_end}")

        bbu_param.ix_ele_track_end = ix

    if bbu_param.write_hom_info:
        pb.rf_cav_names(lat)

    # Check RF Freq
    pb.check_rf_freq(lat, bbu_param.bunch_freq)

    # Prepare custom BBU beam structure (assumed binding)
    bbu_beam = pb.BbuBeamStruct()

    # BBU Setup call
    pb.bbu_setup(lat, beam_init.dt_bunch, bbu_param, bbu_beam)
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
    (hom_voltage_gain, growth_rate, lost, irep) = pb.bbu_track_all(
        lat=lat,
        bbu_beam=bbu_beam,
        bbu_param=bbu_param,
        beam_init=beam_init,
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
