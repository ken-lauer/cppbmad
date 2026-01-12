#include "pybmad/generated/SimUtils_routines_c.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyCalcFileNumber {
  std::string file_name;
  int num_in;
  int num_out;
  bool err_flag;
};
PyCalcFileNumber python_calc_file_number(
    std::string file_name,
    int num_in,
    int num_out,
    bool err_flag) {
  SimUtils::calc_file_number(file_name, num_in, num_out, err_flag);
  auto py_result{PyCalcFileNumber{file_name, num_in, num_out, err_flag}};
  return py_result;
}
struct PyChangeFileNumber {
  std::string file_name;
  int change;
};
PyChangeFileNumber python_change_file_number(
    std::string file_name,
    int change) {
  SimUtils::change_file_number(file_name, change);
  auto py_result{PyChangeFileNumber{file_name, change}};
  return py_result;
}

struct PyCoarseFrequencyEstimate {
  double frequency;
  std::optional<bool> error;
};
PyCoarseFrequencyEstimate python_coarse_frequency_estimate(
    RealAlloc1D& data,
    std::optional<bool> error = std::nullopt) {
  auto _result = SimUtils::coarse_frequency_estimate(data, make_opt_ref(error));
  auto py_result{PyCoarseFrequencyEstimate{_result, error}};
  return py_result;
}
struct PyComplexErrorFunction {
  double wr;
  double wi;
  double zr;
  double zi;
};
PyComplexErrorFunction python_complex_error_function(
    double wr,
    double wi,
    double zr,
    double zi) {
  SimUtils::complex_error_function(wr, wi, zr, zi);
  auto py_result{PyComplexErrorFunction{wr, wi, zr, zi}};
  return py_result;
}
struct PyCosOne {
  double cos1;
};
PyCosOne python_cos_one(double angle, double cos1) {
  SimUtils::cos_one(angle, cos1);
  auto py_result{PyCosOne{cos1}};
  return py_result;
}
struct PyCosc {
  double y;
};
PyCosc python_cosc(double x, std::optional<int> nd, double y) {
  SimUtils::cosc(x, nd, y);
  auto py_result{PyCosc{y}};
  return py_result;
}

void init_SimUtils_routines_c(py::module& m) {
  m.def(
      "calc_file_number",
      &python_calc_file_number,
      py::arg("file_name"),
      py::arg("num_in"),
      py::arg("num_out"),
      py::arg("err_flag"),
      R"""(Parameters
  ----------
  file_name : 
  num_in : 
  num_out : 
  err_flag : 
  )""");
  py::class_<PyCalcFileNumber, std::unique_ptr<PyCalcFileNumber>>(
      m, "CalcFileNumber", "Fortran routine calc_file_number return value")
      .def_readonly("file_name", &PyCalcFileNumber::file_name)
      .def_readonly("num_in", &PyCalcFileNumber::num_in)
      .def_readonly("num_out", &PyCalcFileNumber::num_out)
      .def_readonly("err_flag", &PyCalcFileNumber::err_flag)
      .def("__len__", [](const PyCalcFileNumber&) { return 4; })
      .def("__getitem__", [](const PyCalcFileNumber& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.file_name);
        if (i == 1)
          return py::cast(s.num_in);
        if (i == 2)
          return py::cast(s.num_out);
        if (i == 3)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "change_file_number",
      &python_change_file_number,
      py::arg("file_name"),
      py::arg("change"),
      R"""(Parameters
  ----------
  file_name : 
  change : 
  )""");
  py::class_<PyChangeFileNumber, std::unique_ptr<PyChangeFileNumber>>(
      m, "ChangeFileNumber", "Fortran routine change_file_number return value")
      .def_readonly("file_name", &PyChangeFileNumber::file_name)
      .def_readonly("change", &PyChangeFileNumber::change)
      .def("__len__", [](const PyChangeFileNumber&) { return 2; })
      .def("__getitem__", [](const PyChangeFileNumber& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.file_name);
        if (i == 1)
          return py::cast(s.change);
        throw py::index_error();
      });
  m.def(
      "charge_of",
      &SimUtils::charge_of,
      py::arg("species"),
      py::arg("default_") = py::none(),
      R"""(Function charge_of (species, default) result (charge)

  Routine to return the charge, in units of e+, of a particle.

  Parameters
  ----------
  species : int
      Species ID.
  default : int, optional
      If present then use default value if species = not_set$.

  Returns
  -------
  charge : int
      particle charge.
  )""");
  m.def(
      "charge_to_mass_of",
      &SimUtils::charge_to_mass_of,
      py::arg("species"),
      R"""(Function charge_to_mass_of (species) result (charge_mass_ratio)

  Routine to return the charge (in units of e+) to mass (in units of eV) ratio of a particle.

  Parameters
  ----------
  species : int
      Species ID.

  Returns
  -------
  charge_mass_ratio : float
      particle charge to mass ratio. (1/eV)
  )""");
  m.def(
      "coarse_frequency_estimate",
      &python_coarse_frequency_estimate,
      py::arg("data"),
      py::arg("error") = py::none(),
      R"""(Function coarse_frequency_estimate(data, error) result(frequency)

  Simple function to take periodic data and estimate
  the most dominant frequency by FFT.

  Parameters
  ----------
  data : float
      data to analyze. Preferably size(data) is a power of 2 Otherwise the data is padded with zeros.

  Returns
  -------
  frequency : float
      Frequency corresponding to the largest FFT amplitude
  err : bool
      Error: not enough data. Frequency is near 0 or 0.5
  )""");
  py::class_<
      PyCoarseFrequencyEstimate,
      std::unique_ptr<PyCoarseFrequencyEstimate>>(
      m,
      "CoarseFrequencyEstimate",
      "Fortran routine coarse_frequency_estimate return value")
      .def_readonly("frequency", &PyCoarseFrequencyEstimate::frequency)
      .def_readonly("error", &PyCoarseFrequencyEstimate::error)
      .def("__len__", [](const PyCoarseFrequencyEstimate&) { return 2; })
      .def(
          "__getitem__",
          [](const PyCoarseFrequencyEstimate& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.frequency);
            if (i == 1)
              return py::cast(s.error);
            throw py::index_error();
          });
  m.def(
      "complex_error_function",
      &python_complex_error_function,
      py::arg("wr"),
      py::arg("wi"),
      py::arg("zr"),
      py::arg("zi"),
      R"""(Parameters
  ----------
  wr : 
  wi : 
  zr : 
  zi : 
  )""");
  py::class_<PyComplexErrorFunction, std::unique_ptr<PyComplexErrorFunction>>(
      m,
      "ComplexErrorFunction",
      "Fortran routine complex_error_function return value")
      .def_readonly("wr", &PyComplexErrorFunction::wr)
      .def_readonly("wi", &PyComplexErrorFunction::wi)
      .def_readonly("zr", &PyComplexErrorFunction::zr)
      .def_readonly("zi", &PyComplexErrorFunction::zi)
      .def("__len__", [](const PyComplexErrorFunction&) { return 4; })
      .def(
          "__getitem__",
          [](const PyComplexErrorFunction& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.wr);
            if (i == 1)
              return py::cast(s.wi);
            if (i == 2)
              return py::cast(s.zr);
            if (i == 3)
              return py::cast(s.zi);
            throw py::index_error();
          });
  m.def(
      "cos_one",
      &python_cos_one,
      py::arg("angle"),
      py::arg("cos1"),
      R"""(Parameters
  ----------
  angle : 
  cos1 : 
  )""");
  py::class_<PyCosOne, std::unique_ptr<PyCosOne>>(
      m, "CosOne", "Fortran routine cos_one return value")
      .def_readonly("cos1", &PyCosOne::cos1)
      .def("__len__", [](const PyCosOne&) { return 1; })
      .def("__getitem__", [](const PyCosOne& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.cos1);
        throw py::index_error();
      });
  m.def(
      "cosc",
      &python_cosc,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::arg("y"),
      R"""(Parameters
  ----------
  x : 
  nd : 
  y : 
  )""");
  py::class_<PyCosc, std::unique_ptr<PyCosc>>(
      m, "Cosc", "Fortran routine cosc return value")
      .def_readonly("y", &PyCosc::y)
      .def("__len__", [](const PyCosc&) { return 1; })
      .def("__getitem__", [](const PyCosc& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.y);
        throw py::index_error();
      });
  m.def(
      "create_a_spline",
      &SimUtils::create_a_spline,
      py::arg("r0"),
      py::arg("r1"),
      py::arg("slope0"),
      py::arg("slope1"),
      R"""(Function create_a_spline (r0, r1, slope0, slope1) result (spline)

  Routine to create a single spline given end point positions and slopes.
  The spline will pass through the data points and have the given slopes
  at these points.

  Modules used:
    use spline_mod

  Parameters
  ----------
  r0 : float
      Start (x, y) point.
  r1 : float
      End (x, y) point.
  slope0 : float
      Starting slope.
  slope1 : float
      End slope.

  Returns
  -------
  spline : SplineStruct
      Spline.
  )""");
  m.def(
      "cross_product",
      &SimUtils::cross_product,
      py::arg("a"),
      py::arg("b"),
      py::arg("c"),
      R"""(Parameters
  ----------
  a : float
      Input vectors.
  b : 
  c : 
  )""");
}
