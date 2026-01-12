#include "pybmad/generated/SimUtils_routines_t.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyToStr python_to_str(
    double num,
    std::optional<int> max_signif,
    std::string string) {
  SimUtils::to_str(num, make_opt_ref(max_signif), string);
  auto py_result{PyToStr{num, max_signif, string}};
  return py_result;
}
PyTypeThisFile python_type_this_file(std::string filename) {
  SimUtils::type_this_file(filename);
  auto py_result{PyTypeThisFile{filename}};
  return py_result;
}

void init_SimUtils_routines_t(py::module& m) {
  py::class_<PyToStr, std::unique_ptr<PyToStr>>(
      m, "ToStr", "Fortran routine to_str return value")
      .def_readonly("num", &PyToStr::num)
      .def_readonly("max_signif", &PyToStr::max_signif)
      .def_readonly("string", &PyToStr::string)
      .def("__len__", [](const PyToStr&) { return 3; })
      .def("__getitem__", [](const PyToStr& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.num);
        if (i == 1)
          return py::cast(s.max_signif);
        if (i == 2)
          return py::cast(s.string);
        throw py::index_error();
      });
  m.def(
      "to_str",
      &python_to_str,
      py::arg("num"),
      py::arg("max_signif") = py::none(),
      py::arg("string"),
      R"""(no longer exists
  subroutine test_tune_tracker_lock (tracker_locked)
    implicit none
    logical tracker_locked(2)
  end subroutine

  )""");
  py::class_<
      SimUtils::TricubicCmplxEval,
      std::unique_ptr<SimUtils::TricubicCmplxEval>>(
      m,
      "TricubicCmplxEval",
      "Fortran routine tricubic_cmplx_eval return value")
      .def_readonly("df_dx", &SimUtils::TricubicCmplxEval::df_dx)
      .def_readonly("df_dy", &SimUtils::TricubicCmplxEval::df_dy)
      .def_readonly("df_dz", &SimUtils::TricubicCmplxEval::df_dz)
      .def_readonly("f_val", &SimUtils::TricubicCmplxEval::f_val)
      .def("__len__", [](const SimUtils::TricubicCmplxEval&) { return 4; })
      .def(
          "__getitem__",
          [](const SimUtils::TricubicCmplxEval& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.df_dx);
            if (i == 1)
              return py::cast(s.df_dy);
            if (i == 2)
              return py::cast(s.df_dz);
            if (i == 3)
              return py::cast(s.f_val);
            throw py::index_error();
          });
  m.def(
      "tricubic_cmplx_eval",
      &SimUtils::tricubic_cmplx_eval,
      py::arg("x_norm"),
      py::arg("y_norm"),
      py::arg("z_norm"),
      py::arg("tri_coef"),
      R"""(Function tricubic_cmplx_eval (x_norm, y_norm, z_norm, tri_coef, df_dx, df_dy, df_dz) result (f_val)

  Routine to evaluate a tricubic interpolating complex function.

  Use the routine tricubic_interpolation_cmplx_coefs to generate tri_coef.

  Note: In the equations below, the eight points of the grid box being interpolated range
  from (x0, y0, z0) to (x0+dx, y0+dy, z0+dz).

  Parameters
  ----------
  x_norm : float
      x_norm = (x - x0) / dx
  y_norm : float
      y_norm = (y - y0) / dy
  z_norm : float
      z_norm = (z - z0) / dz
  tri_coef : TricubicCmplxCoefStruct
      Coefficients.

  Returns
  -------
  f_val : complex
      Value of f.
  df_dx : complex
      Normalized first derivative: True df/dx = df_dx * dx
  df_dy : complex
      Normalized first derivative: True df/dy = df_dy * dy
  df_dz : complex
      Normalized first derivative: True df/dz = df_dz * dz
  )""");
  py::class_<PyTypeThisFile, std::unique_ptr<PyTypeThisFile>>(
      m, "TypeThisFile", "Fortran routine type_this_file return value")
      .def_readonly("filename", &PyTypeThisFile::filename)
      .def("__len__", [](const PyTypeThisFile&) { return 1; })
      .def("__getitem__", [](const PyTypeThisFile& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.filename);
        throw py::index_error();
      });
  m.def(
      "type_this_file",
      &python_type_this_file,
      py::arg("filename"),
      R"""(Parameters
  ----------
  filename : 
  )""");
}
