"""
Program wall_generator

Parses lattice file and outputs wall.out file with points that can be linearly
interpolated to make wall geometry.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import sys
from typing import IO

import pybmad as pb
from pybmad import BranchStruct, CoordStruct, RealAlloc1D


def write_header(fp: IO[str]) -> None:
    """
    Write the header to the output file.

    Parameters
    ----------
    fp : IO[str]
        The open file object to write to.
    """
    # Mimicking Fortran (9a18) format
    headers = ["x", "normal_x", "y", "normal_y", "z", "normal_z", "ix_ele", "angle_index", "s"]
    units = ["m", "1", "m", "1", "m", "1", "1", "1", "m"]

    # Format string for 18-character width columns
    fmt = "{:>18}" * 9

    fp.write("#" + fmt.format(*headers) + "\n")
    fp.write("#" + fmt.format(*units) + "\n")


def write_line(f_out: IO[str], orb_out: CoordStruct, ix_ele: int, i_angle: int, s: float):
    vec_str = "".join([f"{val:18.10E}" for val in orb_out.vec])
    f_out.write(f"{vec_str}{ix_ele:>18}{i_angle:>18}{s:18.10E}\n")


def print_logo() -> None:
    print(r"")
    print(r"                         __   ___       ___  __       ___  __   __  ")
    print(r"|  |  /\  |    |        / _` |__  |\ | |__  |__)  /\   |  /  \ |__) ")
    print(r"|/\| /~~\ |___ |___ ___ \__> |___ | \| |___ |  \ /~~\  |  \__/ |  \ ")
    print(r"")


def get_wall(branch: BranchStruct, *, n_angles: int = 2):
    wall3d = branch.wall3d[0]
    for section in wall3d.section:
        ele = branch.ele[section.ix_ele]

        # Calculate relative s position within the element
        # section.s is absolute s, ele.s is total s at exit of element
        s_rel = section.s - (ele.s - ele.value[pb.EleAttribute.L])

        for i_angle in range(n_angles):
            theta = i_angle * 2 * math.pi / n_angles

            if ele.key == pb.SBEND:
                theta = theta - ele.value[pb.EleAttribute.REF_TILT_TOT]

            r_wall, _dr_dtheta, _ix_vertex = pb.calc_wall_radius(
                section.v.view(), math.cos(theta), math.sin(theta)
            )

            dummy_orb = pb.CoordStruct(
                ix_ele=ele.ix_ele,
                # WARNING: NORMAL IS NOT COMPUTED!!!
                # x, px, y, py, z, pz
                vec=[
                    r_wall * math.cos(theta),
                    0.0,
                    r_wall * math.sin(theta),
                    0.0,
                    s_rel,
                    0.0,
                ],
            )

            orb_out = pb.particle_in_global_frame(
                dummy_orb,
                branch,
                in_time_coordinates=True,
                in_body_frame=False,
            )

            yield section, i_angle, ele, orb_out


def get_wall_contour(
    branch: BranchStruct,
    n_angles: int = 2,
    ds: float = 0.0,
    r0=1e-3,
):
    point_s = 0.0

    point = pb.CoordStruct()
    point.s = 0.0
    # Initial dummy vector setup from Fortran: vec(6) = 1.0 (index 5)
    # x, px, y, py, z, pz
    X, _PX, Y, _PY, Z, PZ = range(6)
    point_vec = RealAlloc1D()
    point_vec.resize_bounds(0, 5)
    point_vec[PZ] = 1.0
    last_ele_s = branch.ele[branch.n_ele_track].s

    while True:
        if point_s > last_ele_s:
            break

        ix_ele = pb.element_at_s(branch, point_s, True).ix_ele

        ele = branch.ele[ix_ele]

        point_vec[Z] = point_s - (ele.s - ele.value[pb.L])

        for i_angle in range(n_angles):
            theta = i_angle * 2 * math.pi / n_angles
            if ele.key == pb.SBEND:
                theta = theta - ele.value[pb.EleAttribute.REF_TILT_TOT]

            point_vec[X] = r0 * math.cos(theta)
            point_vec[Y] = r0 * math.sin(theta)

            # Calculate d_radius
            w_res = pb.wall3d_d_radius(point_vec.view(), ele, 1)
            d_radius = w_res.d_radius
            norm_x, norm_y, norm_z = w_res.perp

            r_wall = r0 - d_radius

            dummy_orb = pb.CoordStruct(
                ix_ele=ele.ix_ele,
                # Indices in global frame: [x, normal_x, y, normal_y, z, normal_z]
                vec=[
                    # x
                    r_wall * math.cos(theta),
                    norm_x,
                    # y
                    r_wall * math.sin(theta),
                    norm_y,
                    # z
                    point_vec[Z],  # The longitudinal pos
                    norm_z,
                ],
            )

            orb_out = pb.particle_in_global_frame(
                dummy_orb, branch, in_time_coordinates=True, in_body_frame=False
            )
            yield point_s, i_angle, ele, orb_out

        point_s += ds


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Wall Generator: Create wall geometry from lattice file.",
        prog="wall_generator",
        usage="%(prog)s <lattice> [n_angles] [ds]",
    )
    parser.add_argument("lattice", help="The input Bmad lattice file")
    parser.add_argument(
        "n_angles", type=int, nargs="?", default=2, help="Number of angles in a section (default: 2)"
    )
    parser.add_argument("ds", type=float, nargs="?", default=0.0, help="Step size for contour (default: 0)")
    parser.add_argument("--terse", action="store_true")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    lat_name: str = args.lattice
    n_angles: int = args.n_angles
    ds: float | None = args.ds
    verbose: bool = not args.terse

    if verbose:
        print_logo()

    if n_angles < 1:
        print(f"Not enough angles: {n_angles}", file=sys.stderr)
        sys.exit(1)

    if verbose:
        print(f"Creating wall for lattice file: {lat_name}")
        print(f"Using number of angles: {n_angles}")
        if ds != 0:
            print(f"using ds: {ds}")

    lat = pb.bmad_parser(lat_name).lat

    with pathlib.Path("wall.out").open("w") as f_out:
        branch = lat.branch[0]

        if not branch.wall3d:
            print("No branch.wall3d, exiting...", file=sys.stderr)
            sys.exit(1)

        write_header(f_out)

        for section, i_angle, ele, orb_out in get_wall(branch, n_angles=n_angles):
            write_line(f_out, orb_out=orb_out, ix_ele=ele.ix_ele, i_angle=i_angle, s=section.s)

    if verbose:
        print("Wrote: wall.out")

    if not ds:
        return

    with pathlib.Path("wall_contour.out").open("w") as f_out:
        write_header(f_out)

        for point_s, i_angle, ele, orb_out in get_wall_contour(branch, n_angles=n_angles, ds=ds):
            write_line(f_out, orb_out=orb_out, ix_ele=ele.ix_ele, i_angle=i_angle, s=point_s)

    if verbose:
        print("Wrote: wall_contour.out")


if __name__ == "__main__":
    main()
