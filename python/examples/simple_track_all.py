"""
Simple particle tracking example.
"""

from __future__ import annotations

import argparse

import pybmad


def particle_track(
    lat_file: str = "small_ring.bmad",
    ix_branch: int = 0,
    start_vec: list[float] | None = None,
    start_spin: list[float] | None = None,
    ran_seed: int = 1234,
):
    if start_vec is None:
        start_vec = [1e-3, 0, 0, 0, 0, 0]
    if start_spin is None:
        start_spin = [0.0, 1.0, 0.0]

    pybmad.ran_seed_put(ran_seed)

    print(f'Parsing lattice: "{lat_file}"')
    parsed = pybmad.bmad_parser(lat_file)
    if parsed.err_flag:
        raise RuntimeError(f"Failed to parse lattice: {lat_file}")
    lat = parsed.lat

    # Set bmad_com parameters (after parsing, so these override lattice settings).
    bmad_com = pybmad.get_bmad_com()
    bmad_com.radiation_damping_on = True
    bmad_com.radiation_fluctuations_on = False
    # bmad_com.spin_tracking_on = True

    branch = lat.branch[ix_branch]

    orbit = pybmad.CoordStruct.new_array1d(0)
    pybmad.reallocate_coord(orbit, lat, ix_branch)

    pybmad.init_coord(orbit[0], start_vec, branch.ele[0], pybmad.DOWNSTREAM_END, spin=start_spin)

    result = pybmad.track_all(lat, orbit, ix_branch)
    track_state = result.track_state

    if track_state != pybmad.MOVING_FORWARD:
        print(f"Particle lost at element: {track_state}")
        ix_end = track_state
    else:
        ix_end = branch.n_ele_track - 1

    for ie in range(ix_end + 1):
        ele = branch.ele[ie]
        vec = list(orbit[ie].vec)
        print(ele.name, vec)

    # NOTE: if we're doing multiple passes, we might:
    # set initial orbit = final orbit for next turn.
    # orbit[0] = orbit[n_track]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Track a particle through a lattice for multiple turns.")
    parser.add_argument("--lat", dest="lat_file", default="small_ring.bmad", help="Bmad lattice file")
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
    args = parser.parse_args()

    res = particle_track(
        lat_file=args.lat_file,
        ran_seed=args.ran_seed,
        ix_branch=args.ix_branch,
        start_vec=args.start_vec,
        start_spin=args.start_spin,
    )
