#!/usr/bin/env python3
from __future__ import annotations

import math
import os
from collections.abc import Iterable

import matplotlib.axes
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from pytao import Tao

import pybmad
from pybmad import (
    GIRDER,
    GROUP,
    MULTIPASS_LORD,
    MULTIPASS_SLAVE,
    OVERLAY,
    RBEND,
    # PATCH,
    SBEND,
    SUPER_SLAVE,
    EleAttribute,
    EleStruct,
    FloorPositionStruct,
    TaoEleShapeStruct,
    TaoGraphStruct,
    TaoPlotStruct,
)

TAO_INIT = os.environ.get("TAO_INIT", "$ACC_ROOT_DIR/bmad-doc/tao_examples/optics_matching/tao.init")

tao_color_to_mpl: dict[str, str | None] = {
    "not_set": None,
    "white": "white",
    "black": "black",
    "red": "red",
    "green": "green",
    "blue": "blue",
    "cyan": "cyan",
    "magenta": "magenta",
    "yellow": "yellow",
    "orange": "orange",
    "yellow_green": "yellowgreen",
    "light_green": "lightgreen",
    "navy_blue": "navy",
    "purple": "purple",
    "reddish_purple": "mediumpurple",  # ?
    "dark_grey": "darkgrey",
    "light_grey": "lightgrey",
    "transparent": "none",
}


def cross_product(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def rotate_vec(vec: Iterable[float], axis_idx: int, angle: float) -> list[float]:
    """Rotate vector around axis index (0=x, 1=y, 2=z) by angle."""
    c = math.cos(angle)
    s = math.sin(angle)
    v = list(vec)
    if axis_idx == 0:  # x-axis
        new_y = c * v[1] - s * v[2]
        new_z = s * v[1] + c * v[2]
        return [v[0], new_y, new_z]
    if axis_idx == 1:  # y-axis
        new_x = c * v[0] + s * v[2]
        new_z = -s * v[0] + c * v[2]
        return [new_x, v[1], new_z]
    if axis_idx == 2:  # z-axis
        new_x = c * v[0] - s * v[1]
        new_y = s * v[0] + c * v[1]
        return [new_x, new_y, v[2]]
    return v


def floor_to_screen(graph, r_floor):
    """Wrapper for tao_floor_to_screen_coords returning (x,y)."""
    if not isinstance(r_floor, FloorPositionStruct):
        f = FloorPositionStruct()
        f.r = r_floor
    else:
        f = r_floor

    scr = pybmad.tao_floor_to_screen_coords(graph, f)
    return list(scr.r)[:2]


def tao_draw_floor_plan_orbit(graph: TaoGraphStruct, tao_lat, ax: matplotlib.axes.Axes):
    """
    Draw the particle orbit.
    """
    lat = tao_lat.lat

    for branch in lat.branch:
        x_path = []
        y_path = []

        f0 = FloorPositionStruct()
        f0.r = [0.0, 0.0, 0.0]
        # In Bmad, ele 0 is beginning marker.
        # We trace through all elements.

        for idx, ele in enumerate(branch.ele):
            # Skip multipass slaves pass > 1 if configured
            if ele.slave_status == MULTIPASS_SLAVE and graph.floor_plan.draw_only_first_pass:
                # Logic for multipass visibility similar to elements
                # If we only draw first pass
                ix_pass = pybmad.multipass_chain(ele).ix_pass
                if ix_pass > 1:
                    continue

            # Step through element
            n_seg = 1
            if ele.key in {SBEND, RBEND}:
                n_seg = 10

            l_ele = ele.value[EleAttribute.L]
            ds = l_ele / max(1, n_seg)

            # Start point
            if idx == 0:
                # Element 0
                f_glob = pybmad.coords_local_curvilinear_to_floor(f0, ele, True).global_position
                xs, ys = floor_to_screen(graph, f_glob)
                x_path.append(xs)
                y_path.append(ys)

            for j in range(1, n_seg + 1):
                s_pos = j * ds
                f_local = FloorPositionStruct()
                f_local.r = [0.0, 0.0, s_pos]
                f_glob = pybmad.coords_local_curvilinear_to_floor(f_local, ele, True).global_position
                xs, ys = floor_to_screen(graph, f_glob)
                x_path.append(xs)
                y_path.append(ys)

        # Plot
        ax.plot(x_path, y_path, color="blue", linestyle="--", linewidth=0.5, alpha=0.5, label="Orbit")


def tao_draw_ele_for_floor_plan(
    graph: TaoGraphStruct,
    ele: EleStruct,
    ele_shape: TaoEleShapeStruct,
    label_name: str,
    offset1: float,
    offset2: float,
    ax: matplotlib.axes.Axes,
):
    """
    Draw a single element for the floor plan.
    Port of tao_draw_ele_for_floor_plan from tao_plot_mod.f90.
    """

    # Check element ends
    ele1, _ele2 = pybmad.find_element_ends(ele)
    if not ele1:
        return

    # Check if data or var (which don't have physical extent usually)
    is_data_or_var = False
    if ele_shape and (ele_shape.ele_id.startswith("data::") or ele_shape.ele_id.startswith("var::")):
        is_data_or_var = True

    is_bend = ele.key in {SBEND, RBEND} and not is_data_or_var

    # Calculate floor coordinates for element ends
    # We compute local curvilinear to floor for start (s=0) and end (s=L)

    f0 = FloorPositionStruct()
    f0.r = [0.0, 0.0, 0.0]

    f1 = FloorPositionStruct()
    l_val = ele.value[EleAttribute.L]
    f1.r = [0.0, 0.0, l_val]

    # Transform to global floor coords
    floor1 = pybmad.coords_local_curvilinear_to_floor(f0, ele, True).global_position
    floor2 = pybmad.coords_local_curvilinear_to_floor(f1, ele, True).global_position

    if is_data_or_var:
        # Data/Var often drawn at a single point
        floor1 = floor2

    end1_scr_x, end1_scr_y = floor_to_screen(graph, floor1)
    end2_scr_x, end2_scr_y = floor_to_screen(graph, floor2)

    # Drawing Logic

    # Bend Geometry
    # If it's a bend, we might want to draw an arc.
    x_bend_cen = []
    y_bend_cen = []
    x_bend_curr_minus = []  # -offset2
    y_bend_curr_minus = []
    x_bend_curr_plus = []  # +offset1
    y_bend_curr_plus = []

    plot_page = pybmad.get_super_universe().plot_page
    python_scale = 10.0
    scale = plot_page.floor_plan_shape_scale / python_scale

    if is_bend and ele_shape and ele_shape.draw:
        # Calculate intermediate points
        n_seg = 11  # Segments
        ds_step = l_val / (n_seg - 1)

        off1_s = offset1 * scale
        off2_s = offset2 * scale

        for j in range(n_seg):
            s_pos = j * ds_step

            # Center
            f_local = FloorPositionStruct()
            f_local.r = [0.0, 0.0, s_pos]  # x, y, s

            f_global = pybmad.coords_local_curvilinear_to_floor(f_local, ele, True).global_position
            xs, ys = floor_to_screen(graph, f_global)
            x_bend_cen.append(xs)
            y_bend_cen.append(ys)

            # +offset1 (Left usually)
            f_loc_p = FloorPositionStruct()
            f_loc_p.r = [off1_s, 0.0, s_pos]
            f_glob_p = pybmad.coords_local_curvilinear_to_floor(f_loc_p, ele, True).global_position
            xp, yp = floor_to_screen(graph, f_glob_p)
            x_bend_curr_plus.append(xp)
            y_bend_curr_plus.append(yp)

            # -offset2 (Right usually)
            f_loc_m = FloorPositionStruct()
            f_loc_m.r = [-off2_s, 0.0, s_pos]
            f_glob_m = pybmad.coords_local_curvilinear_to_floor(f_loc_m, ele, True).global_position
            xm, ym = floor_to_screen(graph, f_glob_m)
            x_bend_curr_minus.append(xm)
            y_bend_curr_minus.append(ym)

    if not ele_shape or not ele_shape.draw:
        # Default draw: straight line
        if is_bend and x_bend_cen:
            ax.plot(x_bend_cen, y_bend_cen, color="black", linewidth=1)
        else:
            ax.plot([end1_scr_x, end2_scr_x], [end1_scr_y, end2_scr_y], color="black", linewidth=1)
        return

    # Scale factors
    # Factor to match visual expectation if unit difference exists
    # If Fortran implementation scales differently, we might need adjustment.

    size = ele_shape.size * scale
    off1 = offset1 * scale
    off2 = offset2 * scale

    # Direction vector on screen
    dx = end2_scr_x - end1_scr_x
    dy = end2_scr_y - end1_scr_y
    length = math.hypot(dx, dy)

    if length == 0:
        nx, ny = 0.0, 0.0
    else:
        # Normal vector (rotated 90 deg)
        nx = -dy / length
        ny = dx / length

    # Offsets for width
    dx1 = off1 * nx
    dy1 = off1 * ny
    dx2 = off2 * nx
    dy2 = off2 * ny

    # Shape specific drawing
    shape_full = ele_shape.shape
    if ":" in shape_full:
        _prefix, shape = shape_full.split(":", 1)
    else:
        shape = shape_full

    color = tao_color_to_mpl[ele_shape.color]

    lw = ele_shape.line_width  # in points probably
    if lw == 0:
        lw = 1

    # Helpers for plotting
    def draw_line(xpts, ypts):
        ax.plot(xpts, ypts, color=color, linewidth=lw)

    def fill_poly(xpts, ypts):
        # Matplotlib uses (N, 2) array
        verts = list(zip(xpts, ypts, strict=False))
        poly = Polygon(verts, closed=True, facecolor=color, edgecolor=color, linewidth=lw)
        ax.add_patch(poly)

    # 4 corners relative to the centerline on screen
    # Top Left (Start, +offset1)
    xtl = end1_scr_x + dx1
    ytl = end1_scr_y + dy1
    # Bottom Left (Start, -offset2)
    xbl = end1_scr_x - dx2
    ybl = end1_scr_y - dy2
    # Top Right (End, +offset1)
    xtr = end2_scr_x + dx1
    ytr = end2_scr_y + dy1
    # Bottom Right (End, -offset2)
    xbr = end2_scr_x - dx2
    ybr = end2_scr_y - dy2

    if shape == "box":
        if is_bend and x_bend_curr_plus:
            # Draw curved box
            # Top edge: plus offsets
            # Bottom edge: minus offsets reversed
            # Ends: straight lines connecting them

            # X path: [plus_0...plus_N, minus_N...minus_0, plus_0]
            x_path = x_bend_curr_plus + x_bend_curr_minus[::-1] + [x_bend_curr_plus[0]]
            y_path = y_bend_curr_plus + y_bend_curr_minus[::-1] + [y_bend_curr_plus[0]]

            draw_line(x_path, y_path)
        else:
            # Draw box
            draw_line([xtl, xtr, xbr, xbl, xtl], [ytl, ytr, ybr, ybl, ytl])

    elif shape == "xbox":
        if is_bend and x_bend_curr_plus:
            x_path = x_bend_curr_plus + x_bend_curr_minus[::-1] + [x_bend_curr_plus[0]]
            y_path = y_bend_curr_plus + y_bend_curr_minus[::-1] + [y_bend_curr_plus[0]]
            draw_line(x_path, y_path)
            # X
            draw_line(
                [x_bend_curr_plus[0], x_bend_curr_minus[-1]], [y_bend_curr_plus[0], y_bend_curr_minus[-1]]
            )
            draw_line(
                [x_bend_curr_minus[0], x_bend_curr_plus[-1]], [y_bend_curr_minus[0], y_bend_curr_plus[-1]]
            )
        else:
            draw_line([xtl, xtr, xbr, xbl, xtl], [ytl, ytr, ybr, ybl, ytl])
        # Cross
        draw_line([xtl, xbr], [ytl, ybr])
        draw_line([xbl, xtr], [ybl, ytr])

    elif shape == "bow_tie":
        draw_line([xtl, xbl], [ytl, ybl])  # Left vertical
        draw_line([xtr, xbr], [ytr, ybr])  # Right vertical
        draw_line([xtl, xbr], [ytl, ybr])  # Diagonal
        draw_line([xbl, xtr], [ybl, ytr])  # Diagonal

    elif shape == "diamond":
        # Midpoints
        x_mid = (end1_scr_x + end2_scr_x) / 2
        y_mid = (end1_scr_y + end2_scr_y) / 2

        # Start point: (end1_scr_x, end1_scr_y)
        # End point: (end2_scr_x, end2_scr_y)
        # Top point: Midline shifted by offset1
        dtx = dx1
        dty = dy1
        dbx = -dx2
        dby = -dy2  # Bottom vector

        xm_top = x_mid + dtx
        ym_top = y_mid + dty
        xm_bot = x_mid + dbx
        ym_bot = y_mid + dby

        draw_line(
            [end1_scr_x, xm_top, end2_scr_x, xm_bot, end1_scr_x],
            [end1_scr_y, ym_top, end2_scr_y, ym_bot, end1_scr_y],
        )

    elif shape == "circle":
        # Circle at center
        x_mid = (end1_scr_x + end2_scr_x) / 2
        y_mid = (end1_scr_y + end2_scr_y) / 2
        radius = size / 2.0
        circle = Circle((x_mid, y_mid), radius, color=color, fill=False, linewidth=lw)
        ax.add_patch(circle)

    elif "triangle" in shape:
        # d_triangle, u_triangle, l, r etc.
        # Simplified handling as generic marker or oriented triangle
        # Center
        x_mid = (end1_scr_x + end2_scr_x) / 2
        y_mid = (end1_scr_y + end2_scr_y) / 2

    else:
        # Fallback to line
        draw_line([end1_scr_x, end2_scr_x], [end1_scr_y, end2_scr_y])

    # Label
    if label_name:
        # Center position
        x_mid = (end1_scr_x + end2_scr_x) / 2
        y_mid = (end1_scr_y + end2_scr_y) / 2

        # Shift perpendicular based on offset2 and text height approximation
        # Fortran:
        # call tao_floor_to_screen_coords (graph, floor_mid - offset2 * 1.5, label_cen)

        text_shift = 1.5 * scale * offset2  # Heuristic
        lx = x_mid - nx * text_shift
        ly = y_mid - ny * text_shift

        rotation = math.degrees(math.atan2(dy, dx))

        if graph.floor_plan.flip_label_side:
            lx = x_mid + nx * text_shift
            ly = y_mid + ny * text_shift

        ax.text(
            lx,
            ly,
            label_name,
            rotation=rotation,
            ha="center",
            va="center",
            color=color,
            clip_on=True,
            fontsize=8,
        )


def tao_draw_building_wall(graph: TaoGraphStruct, ax: matplotlib.axes.Axes):
    """
    Draw building walls.
    Port of tao_draw_building_wall from tao_plot_mod.f90
    """

    s = pybmad.get_super_universe()
    if not graph.floor_plan.draw_building_wall or not s.building_wall.section:
        return

    for sec in s.building_wall.section:
        # Wall shape
        wall_shape = pybmad.tao_pointer_to_building_wall_shape(sec.name)
        if not wall_shape:
            continue

        color = wall_shape.color
        lw = wall_shape.line_width
        if lw == 0:
            lw = 1

        n_pts = len(sec.point)
        # sec.point is array of TaoBuildingWallPointStruct

        for i in range(1, n_pts):
            pt0_raw = sec.point[i - 1]
            pt1_raw = sec.point[i]

            # Orient points (handle Global vs Local)
            w_pt0 = pybmad.tao_oreint_building_wall_pt(pt0_raw)
            w_pt1 = pybmad.tao_oreint_building_wall_pt(pt1_raw)

            f0 = FloorPositionStruct()
            f0.r = [w_pt0.x, 0.0, w_pt0.z]
            f1 = FloorPositionStruct()
            f1.r = [w_pt1.x, 0.0, w_pt1.z]

            x0, y0 = floor_to_screen(graph, f0)
            x1, y1 = floor_to_screen(graph, f1)

            radius = w_pt1.radius

            if radius == 0:
                # Straight line
                ax.plot([x0, x1], [y0, y1], color=color, linewidth=lw)
            else:
                # Arc
                # We need center
                f_cen = FloorPositionStruct()
                f_cen.r = [w_pt1.x_center, 0.0, w_pt1.z_center]
                _xc, _yc = floor_to_screen(graph, f_cen)

                # Angles
                # dx0 = x0 - xc; dy0 = y0 - yc
                # dx1 = x1 - xc; dy1 = y1 - yc
                # ang0 = math.degrees(math.atan2(dy0, dx0))
                # ang1 = math.degrees(math.atan2(dy1, dx1))

                # Simple segmentation
                n_seg = 20
                theta0 = math.atan2(w_pt0.z - w_pt1.z_center, w_pt0.x - w_pt1.x_center)
                theta1 = math.atan2(w_pt1.z - w_pt1.z_center, w_pt1.x - w_pt1.x_center)

                # Ensure correct direction for arc
                d_theta = theta1 - theta0
                if d_theta > math.pi:
                    d_theta -= 2 * math.pi
                if d_theta < -math.pi:
                    d_theta += 2 * math.pi

                x_arc = []
                y_arc = []

                for j in range(n_seg + 1):
                    frac = j / float(n_seg)
                    th = theta0 + frac * d_theta

                    x_w = w_pt1.x_center + abs(radius) * math.cos(th)
                    z_w = w_pt1.z_center + abs(radius) * math.sin(th)

                    f_pt = FloorPositionStruct()
                    f_pt.r = [x_w, 0.0, z_w]
                    xs, ys = floor_to_screen(graph, f_pt)
                    x_arc.append(xs)
                    y_arc.append(ys)

                ax.plot(x_arc, y_arc, color=color, linewidth=lw)


def get_ele_slaves(ele: pybmad.EleStruct) -> list[pybmad.EleStruct]:
    res = []
    for j in range(1, ele.n_slave + 1):
        slave_ptr = pybmad.pointer_to_slave(ele, j)
        slave = slave_ptr.slave_ptr
        if slave is not None:
            res.append(slave)
    return res


def draw_this_floor_plan(
    uni: pybmad.TaoUniverseStruct,
    plot: TaoPlotStruct,  # noqa: ARG001
    graph: TaoGraphStruct,
    ax: matplotlib.axes.Axes,
):
    """
    Draw floor plan for a specific universe
    """
    s = pybmad.get_super_universe()
    orbit_lat = graph.floor_plan.orbit_lattice

    if orbit_lat == "model":
        tao_lat = uni.model
    elif orbit_lat == "design":
        tao_lat = uni.design
    elif orbit_lat == "base":
        tao_lat = uni.base
    else:
        tao_lat = uni.model

    if tao_lat is None:
        return

    lat = tao_lat.lat

    for branch in lat.branch:
        for ele in branch.ele:
            # Skip super_slave as they are part of a lord which might be drawn
            if ele.slave_status == SUPER_SLAVE:
                continue

            ix_shape_min = 1

            while True:
                shape_info = pybmad.tao_ele_shape_info(
                    uni.ix_uni, ele, s.plot_page.floor_plan.ele_shape.view(), ix_shape_min
                )

                e_shape = None
                if shape_info and shape_info.e_shape:
                    e_shape = shape_info.e_shape
                    label_name = shape_info.label_name
                    offset1 = shape_info.y1
                    offset2 = shape_info.y2
                    ix_shape_min = shape_info.ix_shape_min
                else:
                    # No shape found
                    # Check special elements
                    if ele.key in (OVERLAY, GROUP, GIRDER):
                        # Logic to continue
                        pass
                    break

                if graph.floor_plan.draw_only_first_pass and ele.slave_status == MULTIPASS_SLAVE:
                    ix_pass = pybmad.multipass_chain(ele).ix_pass
                    if ix_pass > 1:
                        break

                if ele.lord_status == MULTIPASS_LORD:
                    for j, slave in enumerate(get_ele_slaves(ele)):
                        if graph.floor_plan.draw_only_first_pass and j > 0:
                            break

                        s_info2 = pybmad.tao_ele_shape_info(
                            uni.ix_uni, slave, s.plot_page.floor_plan.ele_shape
                        )
                        if s_info2 and s_info2.e_shape:
                            continue  # Slave handles it

                        tao_draw_ele_for_floor_plan(graph, slave, e_shape, label_name, offset1, offset2, ax)
                else:
                    tao_draw_ele_for_floor_plan(graph, ele, e_shape, label_name, offset1, offset2, ax)

                if not e_shape.multi:
                    break

    tao_draw_floor_plan_orbit(graph, tao_lat, ax)

    if s.building_wall.section and graph.floor_plan.draw_building_wall:
        tao_draw_building_wall(graph, ax)


def tao_draw_floor_plan():
    """
    Main entry to draw floor plan.
    Iterates over regions and graphs, finding 'floor_plan' type graphs.
    """

    s = pybmad.get_super_universe()

    regions = []
    for r in s.plot_page.region:
        if r.visible and r.plot and r.plot.name:  # Basic check
            regions.append(r)

    for region in regions:
        plot = region.plot
        for graph in plot.graph:
            if graph.type != "floor_plan":
                continue
            if not graph.is_valid:
                continue

            title = graph.title if graph.title else "Floor Plan"
            _fig, ax = plt.subplots()
            ax.set_title(title)
            ax.set_aspect("equal")

            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            # graph.floor_plan.view might be 'xz' (default) or 'xy' etc.
            # tao_floor_to_screen_coords handles the projection.

            ix_univ = graph.ix_universe
            if ix_univ == -2:  # All
                # Python 0-based, Fortran 1-based universe array idx
                universes = list(s.u)
            else:
                universes = [s.u[pybmad.tao_universe_index(ix_univ) - 1]]

            for uni in universes:
                draw_this_floor_plan(uni, plot, graph, ax)

            if graph.draw_grid:
                ax.grid(True)

            plt.show()


def main():
    tao = Tao(init_file=TAO_INIT, plot="mpl")
    print(tao)
    plt.ion()
    tao.plot("floor_plan")
    plt.xlim(0, 5.0)
    # plt.ylim(-1, 0.5)
    print("Initialized Tao.")

    tao_draw_floor_plan()
    plt.xlim(0, 5.0)
    # plt.ylim(-1, 0.5)
    plt.ioff()
    plt.show()


if __name__ == "__main__":
    main()
