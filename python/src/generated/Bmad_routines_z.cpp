#include "pybmad/generated/Bmad_routines_z.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_z(nb::module_ &m) {
  nb::class_<Bmad::ZAtSurface>(m, "ZAtSurface", "z_at_surface return type")
      .def_ro("err_flag", &Bmad::ZAtSurface::err_flag)
      .def_ro("dz_dxy", &Bmad::ZAtSurface::dz_dxy)
      .def_ro("z", &Bmad::ZAtSurface::z)
      .def("__len__", [](const Bmad::ZAtSurface &) { return 3; })
      .def("__getitem__", [](const Bmad::ZAtSurface &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.dz_dxy);
        if (i == 2)
          return nb::cast(s.z);
        throw nb::index_error();
      });
  m.def(
      "z_at_surface",
      &Bmad::z_at_surface,
      nb::arg("ele"),
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("extend_grid") = nb::none(),
      R"""(Routine return the height (z) of the surface for a particular (x,y) position.
Remember: +z points into the element.

Parameters
----------
ele : EleStruct
    Element

x : float
    Photon coordinates on surface.

y : float
    Photon coordinates on surface.

extend_grid : bool, optional
    If a grid is involved and (x, y) is outside of the grid, and extend_grid = True: Pretend (x, y) is at
    edge. Default is False.

Returns
-------
err_flag : bool
    Set True if cannot compute z due to, say, point being outside of ellipseoid or grid bounds.

z : float
    z coordinate.

dz_dxy : 1D array of float (shape: 2), optional
    Surface slope at (x, y).
)"""
  );
  m.def(
      "zero_ele_kicks",
      &Bmad::zero_ele_kicks,
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine zero_ele_kicks

Parameters
----------
ele : EleStruct
    Element with possible nonzero kicks.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with no kicks.
)"""
  );
  m.def(
      "zero_ele_offsets",
      &Bmad::zero_ele_offsets,
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine zero_ele_offsets

Parameters
----------
ele : EleStruct
    Element with possible nonzero offsets, etc.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with no (mis)orientation.
)"""
  );
  m.def(
      "zero_lr_wakes_in_lat",
      &Bmad::zero_lr_wakes_in_lat,
      nb::arg("lat"),
      R"""(Routine to zero the long range wake amplitudes for the elements that have
long range wakes in a lattice.

Parameters
----------
lat : LatStruct
    Lattice
)"""
  );
  m.def(
      "zlafun",
      &Bmad::zlafun,
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("z"),
      R"""(Wrapper for Fortran routine zlafun

Parameters
----------
x : float

y : float

z : float

Returns
-------
res : float
)"""
  );
}
