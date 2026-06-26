#include "pybmad/generated/SimUtils_routines_e.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_e(nb::module_ &m) {
  nb::class_<SimUtils::Elbd>(m, "Elbd", "elbd return type")
      .def_ro("b", &SimUtils::Elbd::b)
      .def_ro("d", &SimUtils::Elbd::d)
      .def("__len__", [](const SimUtils::Elbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Elbd &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.b);
        if (i == 1)
          return nb::cast(s.d);
        throw nb::index_error();
      });
  m.def(
      "elbd",
      &SimUtils::elbd,
      nb::arg("phi"),
      nb::arg("phic"),
      nb::arg("mc"),
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
  nb::class_<SimUtils::Elcbd>(m, "Elcbd", "elcbd return type")
      .def_ro("b", &SimUtils::Elcbd::b)
      .def_ro("dx", &SimUtils::Elcbd::dx)
      .def("__len__", [](const SimUtils::Elcbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Elcbd &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.b);
        if (i == 1)
          return nb::cast(s.dx);
        throw nb::index_error();
      });
  m.def(
      "elcbd",
      &SimUtils::elcbd,
      nb::arg("c0"),
      nb::arg("mc"),
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
  nb::class_<SimUtils::Ellipinc>(m, "Ellipinc", "ellipinc return type")
      .def_ro("ellipkinc", &SimUtils::Ellipinc::ellipkinc)
      .def_ro("ellipeinc", &SimUtils::Ellipinc::ellipeinc)
      .def("__len__", [](const SimUtils::Ellipinc &) { return 2; })
      .def("__getitem__", [](const SimUtils::Ellipinc &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.ellipkinc);
        if (i == 1)
          return nb::cast(s.ellipeinc);
        throw nb::index_error();
      });
  m.def(
      "ellipinc",
      &SimUtils::ellipinc,
      nb::arg("phi"),
      nb::arg("m"),
      R"""(Calculates the first and second incomplete elliptic integrals,
using methods from T. Fukushima, (2011, 2018)

Uses classical transformations to handle negative m.
This package needs a function for the third kind to use the new 2018 transformations.
)"""
  );
  nb::class_<SimUtils::Elsbd>(m, "Elsbd", "elsbd return type")
      .def_ro("b", &SimUtils::Elsbd::b)
      .def_ro("d", &SimUtils::Elsbd::d)
      .def("__len__", [](const SimUtils::Elsbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Elsbd &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.b);
        if (i == 1)
          return nb::cast(s.d);
        throw nb::index_error();
      });
  m.def(
      "elsbd",
      &SimUtils::elsbd,
      nb::arg("s0"),
      nb::arg("mc"),
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
      nb::arg("spline"),
      nb::arg("which_end"),
      R"""(Routine to calculate the slopes at the ends of a spline array

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
      nb::arg("err_str") = nb::none(),
      R"""(Wrapper for Fortran routine err_exit

Parameters
----------
err_str : str, optional
)"""
  );
}
