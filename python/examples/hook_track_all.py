"""
Particle tracking example including some hooks.
"""

from __future__ import annotations

import argparse
import logging

import pybmad as pb


def particle_track(
    lat_file: str,
    ix_branch: int = 0,
    start_vec: list[float] | None = None,
    start_spin: list[float] | None = None,
    ran_seed: int = 1234,
):
    if start_vec is None:
        start_vec = [1e-3, 0, 0, 0, 0, 0]
    if start_spin is None:
        start_spin = [0.0, 1.0, 0.0]

    pb.ran_seed_put(ran_seed)

    print(f'Parsing lattice: "{lat_file}"')
    parsed = pb.bmad_parser(lat_file)
    if parsed.err_flag:
        raise RuntimeError(f"Failed to parse lattice: {lat_file}")
    lat = parsed.lat

    # Set bmad_com parameters (after parsing, so these override lattice settings).
    bmad_com = pb.get_bmad_com()
    bmad_com.radiation_damping_on = True
    bmad_com.radiation_fluctuations_on = False
    # bmad_com.spin_tracking_on = True

    branch = lat.branch[ix_branch]

    orbit = pb.CoordStruct.new_array1d(0)
    pb.reallocate_coord(orbit, lat, ix_branch)

    pb.init_coord(orbit[0], start_vec, branch.ele[0], pb.DOWNSTREAM_END, spin=start_spin)

    pb.lattice_bookkeeper(lat)

    seen = []

    def preprocess(
        start_orb: pb.CoordStruct,
        ele: pb.EleStruct,
        param: pb.LatParamStruct,
        err_flag: bool,
        finished: bool,
        radiation_included: bool,
        track: pb.TrackStruct | None,
    ):
        # Observe each element just before it is tracked. Returning None leaves all
        # out-parameters unchanged, so Bmad's standard tracking proceeds.
        logging.debug(
            f"Ele: {ele.name} Starting orbit: {start_orb.vec} {param=} {err_flag=} {finished=} {radiation_included=} {track=}"
        )
        seen.append((ele.name, float(start_orb.vec[0])))
        # Return None (leave finished=False) so Bmad tracks the element normally and
        # track1_postprocess still fires. Returning finished=True makes track1 return
        # early -- before postprocess -- and skips the actual tracking.
        # return (err_flag, finished, radiation_included)

    def postprocess(
        start_orb: pb.CoordStruct, ele: pb.EleStruct, param: pb.LatParamStruct, end_orb: pb.CoordStruct
    ):
        # Observe the exit orbit after each element. end_orb is a live proxy.
        logging.debug(f"Ele: {ele.name} Starting orbit: {start_orb.vec} {param=} {end_orb=}")
        seen[-1] = (*seen[-1], float(end_orb.vec[0]))

    pb.hooks.track1_preprocess = preprocess
    pb.hooks.track1_postprocess = postprocess

    result = pb.track_all(lat, orbit, ix_branch)
    track_state = result.track_state

    print(f"\nHooks observed {len(seen)} element(s):\n")

    for item in seen:
        if len(item) == 3:
            name, start_orbit_x, end_orbit_x = item
            print(f"- {name}: {start_orbit_x=:g} {end_orbit_x=:g}")
        else:
            print("not just 1 pre/post hook?", len(item), item)

    if track_state != pb.MOVING_FORWARD:
        print(f"Particle lost at element: {track_state}")
        ix_end = track_state
    else:
        ix_end = branch.n_ele_track - 1

    for ie in range(ix_end + 1):
        ele = branch.ele[ie]
        vec = list(orbit[ie].vec)
        print(ele.name, vec)
        break

    # NOTE: if we're doing multiple passes, we might:
    # set initial orbit = final orbit for next turn.
    # orbit[0] = orbit[n_track]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Track a particle through a lattice for multiple turns.")
    parser.add_argument(
        "--lat",
        dest="lat_file",
        default="$ACC_ROOT_DIR/bmad-doc/lattices/small_ring/small_ring.bmad",
        help="Bmad lattice file",
    )
    parser.add_argument("--ix-branch", type=int, default=0, help="Branch index")
    parser.add_argument("--ran-seed", type=int, default=1234, help="Random seed (0 = system clock)")
    parser.add_argument(
        "--start-vec",
        type=float,
        nargs=6,
        default=[1e-3, 0, 0, 0, 0, 0],
        help="Starting (x, px, y, py, z, pz)",
    )
    parser.add_argument(
        "--start-spin",
        type=float,
        nargs=3,
        default=[0.0, 1.0, 0.0],
        help="Starting spin (Sx, Sy, Sz)",
    )

    # To see all the debug log messages:
    # logging.basicConfig(level="DEBUG")
    args = parser.parse_args()

    res = particle_track(
        lat_file=args.lat_file,
        ran_seed=args.ran_seed,
        ix_branch=args.ix_branch,
        start_vec=args.start_vec,
        start_spin=args.start_spin,
    )
