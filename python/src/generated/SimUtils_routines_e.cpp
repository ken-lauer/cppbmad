#include "pybmad/generated/SimUtils_routines_e.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_e(py::module &m) {
  py::class_<SimUtils::Elbd, std::unique_ptr<SimUtils::Elbd>>(m, "Elbd", "elbd return type")
      .def_readonly("b", &SimUtils::Elbd::b)
      .def_readonly("d", &SimUtils::Elbd::d)
      .def("__len__", [](const SimUtils::Elbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Elbd &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.b);
        if (i == 1)
          return py::cast(s.d);
        throw py::index_error();
      });
  m.def(
      "elbd",
      &SimUtils::elbd,
      py::arg("phi"),
      py::arg("phic"),
      py::arg("mc"),
      R"""(Wrapper for Fortran routine elbd

Parameters
----------
phi : float

phic : float

mc : float

Returns
-------
b : float

d : float
)"""
  );
  py::class_<SimUtils::Elcbd, std::unique_ptr<SimUtils::Elcbd>>(m, "Elcbd", "elcbd return type")
      .def_readonly("b", &SimUtils::Elcbd::b)
      .def_readonly("dx", &SimUtils::Elcbd::dx)
      .def("__len__", [](const SimUtils::Elcbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Elcbd &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.b);
        if (i == 1)
          return py::cast(s.dx);
        throw py::index_error();
      });
  m.def(
      "elcbd",
      &SimUtils::elcbd,
      py::arg("c0"),
      py::arg("mc"),
      R"""(Wrapper for Fortran routine elcbd

Parameters
----------
c0 : float

mc : float

Returns
-------
b : float

dx : float
)"""
  );
  py::class_<SimUtils::Ellipinc, std::unique_ptr<SimUtils::Ellipinc>>(
      m,
      "Ellipinc",
      "ellipinc return type"
  )
      .def_readonly("ellipkinc", &SimUtils::Ellipinc::ellipkinc)
      .def_readonly("ellipeinc", &SimUtils::Ellipinc::ellipeinc)
      .def("__len__", [](const SimUtils::Ellipinc &) { return 2; })
      .def("__getitem__", [](const SimUtils::Ellipinc &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.ellipkinc);
        if (i == 1)
          return py::cast(s.ellipeinc);
        throw py::index_error();
      });
  m.def(
      "ellipinc",
      &SimUtils::ellipinc,
      py::arg("phi"),
      py::arg("m"),
      R"""(subroutine ellipinc(phi, m, ellipkinc, ellipeinc)

Calculates the first and second incomplete elliptic integrals,
using methods from T. Fukushima, (2011, 2018)

Uses classical transformations to handle negative m.
This package needs a function for the third kind to use the new 2018 transformations.
)"""
  );
  py::class_<SimUtils::Elsbd, std::unique_ptr<SimUtils::Elsbd>>(m, "Elsbd", "elsbd return type")
      .def_readonly("b", &SimUtils::Elsbd::b)
      .def_readonly("d", &SimUtils::Elsbd::d)
      .def("__len__", [](const SimUtils::Elsbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Elsbd &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.b);
        if (i == 1)
          return py::cast(s.d);
        throw py::index_error();
      });
  m.def(
      "elsbd",
      &SimUtils::elsbd,
      py::arg("s0"),
      py::arg("mc"),
      R"""(Wrapper for Fortran routine elsbd

Parameters
----------
s0 : float

mc : float

Returns
-------
b : float

d : float
)"""
  );
  m.def(
      "end_akima_spline_calc",
      &SimUtils::end_akima_spline_calc,
      py::arg("spline"),
      py::arg("which_end"),
      R"""(Subroutine end_akima_spline_calc (spline, which_end)

Routine to calculate the slopes at the ends of a spline array

Parameters
----------
spline : 1D array of SplineStruct
    Array of splines.
    This parameter is an input/output and is modified in-place.
    As an output, spline: Array with slopes at end calculated.

which_end : int
    0 => calculate slopes for the start end of the array. 1 => calculate slopes for the end end of the array.
)"""
  );
  m.def(
      "err_exit",
      &SimUtils::err_exit,
      py::arg("err_str") = py::none(),
      R"""(Wrapper for Fortran routine err_exit

Parameters
----------
err_str : character, optional
)"""
  );
}
