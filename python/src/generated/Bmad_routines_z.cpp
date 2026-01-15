#include "pybmad/generated/Bmad_routines_z.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_z(py::module &m) {
  py::class_<Bmad::ZAtSurface, std::unique_ptr<Bmad::ZAtSurface>>(
      m,
      "ZAtSurface",
      "z_at_surface return type"
  )
      .def_readonly("err_flag", &Bmad::ZAtSurface::err_flag)
      .def_readonly("dz_dxy", &Bmad::ZAtSurface::dz_dxy)
      .def_readonly("z", &Bmad::ZAtSurface::z)
      .def("__len__", [](const Bmad::ZAtSurface &) { return 3; })
      .def("__getitem__", [](const Bmad::ZAtSurface &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.dz_dxy);
        if (i == 2)
          return py::cast(s.z);
        throw py::index_error();
      });
  m.def(
      "z_at_surface",
      &Bmad::z_at_surface,
      py::arg("ele"),
      py::arg("x"),
      py::arg("y"),
      py::arg("extend_grid") = py::none(),
      R"""(Function z_at_surface (ele, x, y, err_flag, extend_grid, dz_dxy) result (z)

Routine return the height (z) of the surface for a particular (x,y) position.
Remember: +z points into the element.

Parameters
----------
ele : EleStruct
    Element x, y        -- real(rp): Photon coordinates on surface.
extend_grid : bool, optional
    If a grid is involved and (x, y) is outside of the grid, and extend_grid = True: Pretend (x, y) is at
    edge. Default is False.

Returns
-------
z : float
    z coordinate.
err_flag : bool
    Set True if cannot compute z due to, say, point being outside of ellipseoid or grid bounds.
dz_dxy : float
    Surface slope at (x, y).

Notes
-----
Remember: +z points into the element.
)"""
  );
  m.def(
      "zero_ele_kicks",
      &Bmad::zero_ele_kicks,
      R"""(Parameters
----------
ele : EleStruct
    Element with no kicks.
)"""
  );
  m.def(
      "zero_ele_offsets",
      &Bmad::zero_ele_offsets,
      R"""(Parameters
----------
ele : EleStruct
    Element with no (mis)orientation.
)"""
  );
  m.def(
      "zero_lr_wakes_in_lat",
      &Bmad::zero_lr_wakes_in_lat,
      py::arg("lat"),
      R"""(Subroutine zero_lr_wakes_in_lat (lat)

Routine to zero the long range wake amplitudes for the elements that have
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
      py::arg("x"),
      py::arg("y"),
      py::arg("z"),
      py::arg("res"),
      R"""(Parameters
----------
x : 
y : 
z : 
res : 
)"""
  );
}
