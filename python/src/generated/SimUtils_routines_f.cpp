#include "pybmad/generated/SimUtils_routines_f.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyFactorial python_factorial(int n, double fact) {
  SimUtils::factorial(n, fact);
  auto py_result{PyFactorial{n, fact}};
  return py_result;
}
PyFileDirectorizer python_file_directorizer(
    std::string in_file,
    std::string out_file,
    std::string directory,
    bool add_switch) {
  SimUtils::file_directorizer(in_file, out_file, directory, add_switch);
  auto py_result{PyFileDirectorizer{in_file, out_file, directory, add_switch}};
  return py_result;
}
PyFileGet python_file_get(
    std::string string,
    std::string dflt_file_name,
    std::string file_name) {
  SimUtils::file_get(string, dflt_file_name, file_name);
  auto py_result{PyFileGet{string, dflt_file_name, file_name}};
  return py_result;
}
PyFileGetOpen python_file_get_open(
    std::string string,
    std::string dflt_file_name,
    std::string file_name,
    int file_unit,
    bool readonly) {
  SimUtils::file_get_open(
      string, dflt_file_name, file_name, file_unit, readonly);
  auto py_result{
      PyFileGetOpen{string, dflt_file_name, file_name, file_unit, readonly}};
  return py_result;
}
PyFileSuffixer python_file_suffixer(
    std::string in_file_name,
    std::string out_file_name,
    std::string suffix,
    bool add_switch) {
  SimUtils::file_suffixer(in_file_name, out_file_name, suffix, add_switch);
  auto py_result{
      PyFileSuffixer{in_file_name, out_file_name, suffix, add_switch}};
  return py_result;
}
PyFindLocationInt python_find_location_int(
    IntAlloc1D& arr,
    int value,
    int ix_match) {
  SimUtils::find_location(arr, value, ix_match);
  auto py_result{PyFindLocationInt{value, ix_match}};
  return py_result;
}
PyFindLocationLogic python_find_location_logic(
    BoolAlloc1D& arr,
    bool value,
    int ix_match) {
  SimUtils::find_location(arr, value, ix_match);
  auto py_result{PyFindLocationLogic{value, ix_match}};
  return py_result;
}
PyFindLocationReal python_find_location_real(
    RealAlloc1D& arr,
    double value,
    int ix_match) {
  SimUtils::find_location(arr, value, ix_match);
  auto py_result{PyFindLocationReal{ix_match}};
  return py_result;
}
PyFixedwindowls python_fixedwindowls(double ynew, int id, double z) {
  SimUtils::fixedwindowls(ynew, id, z);
  auto py_result{PyFixedwindowls{z}};
  return py_result;
}

void init_SimUtils_routines_f(py::module& m) {
  py::class_<PyFactorial, std::unique_ptr<PyFactorial>>(
      m, "Factorial", "factorial return type")
      .def_readonly("n", &PyFactorial::n)
      .def_readonly("fact", &PyFactorial::fact)
      .def("__len__", [](const PyFactorial&) { return 2; })
      .def("__getitem__", [](const PyFactorial& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.n);
        if (i == 1)
          return py::cast(s.fact);
        throw py::index_error();
      });
  m.def(
      "factorial",
      &python_factorial,
      py::arg("n"),
      py::arg("fact"),
      R"""(Parameters
  ----------
  n : 
  fact : 
  )""");
  m.def(
      "faddeeva_function",
      &SimUtils::faddeeva_function,
      py::arg("z"),
      py::arg("w"),
      py::arg("dw"),
      R"""(Parameters
  ----------
  z : 
  w : 
  dw : 
  )""");
  m.def(
      "fft_1d",
      &SimUtils::fft_1d,
      py::arg("arr"),
      py::arg("isign"),
      R"""(no longer exists
  subroutine fff_sub(line, error)
    implicit none
    character(*) line
    logical error
  end subroutine

  Parameters
  ----------
  arr : complex
      Input array.
      This parameter is an input/output and is modified in-place. As an output: FFT of array.
  isign : int
      -1 => "Forward" transform, +1 => "Backwards" transform.
  )""");
  py::class_<PyFileDirectorizer, std::unique_ptr<PyFileDirectorizer>>(
      m, "FileDirectorizer", "file_directorizer return type")
      .def_readonly("in_file", &PyFileDirectorizer::in_file)
      .def_readonly("out_file", &PyFileDirectorizer::out_file)
      .def_readonly("directory", &PyFileDirectorizer::directory)
      .def_readonly("add_switch", &PyFileDirectorizer::add_switch)
      .def("__len__", [](const PyFileDirectorizer&) { return 4; })
      .def("__getitem__", [](const PyFileDirectorizer& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.in_file);
        if (i == 1)
          return py::cast(s.out_file);
        if (i == 2)
          return py::cast(s.directory);
        if (i == 3)
          return py::cast(s.add_switch);
        throw py::index_error();
      });
  m.def(
      "file_directorizer",
      &python_file_directorizer,
      py::arg("in_file"),
      py::arg("out_file"),
      py::arg("directory"),
      py::arg("add_switch"),
      R"""(Parameters
  ----------
  in_file : 
  out_file : 
  directory : 
  add_switch : 
  )""");
  py::class_<PyFileGet, std::unique_ptr<PyFileGet>>(
      m, "FileGet", "file_get return type")
      .def_readonly("string", &PyFileGet::string)
      .def_readonly("dflt_file_name", &PyFileGet::dflt_file_name)
      .def_readonly("file_name", &PyFileGet::file_name)
      .def("__len__", [](const PyFileGet&) { return 3; })
      .def("__getitem__", [](const PyFileGet& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.string);
        if (i == 1)
          return py::cast(s.dflt_file_name);
        if (i == 2)
          return py::cast(s.file_name);
        throw py::index_error();
      });
  m.def(
      "file_get",
      &python_file_get,
      py::arg("string"),
      py::arg("dflt_file_name"),
      py::arg("file_name"),
      R"""(Parameters
  ----------
  string : 
  dflt_file_name : 
  file_name : 
  )""");
  py::class_<PyFileGetOpen, std::unique_ptr<PyFileGetOpen>>(
      m, "FileGetOpen", "file_get_open return type")
      .def_readonly("string", &PyFileGetOpen::string)
      .def_readonly("dflt_file_name", &PyFileGetOpen::dflt_file_name)
      .def_readonly("file_name", &PyFileGetOpen::file_name)
      .def_readonly("file_unit", &PyFileGetOpen::file_unit)
      .def_readonly("readonly", &PyFileGetOpen::readonly)
      .def("__len__", [](const PyFileGetOpen&) { return 5; })
      .def("__getitem__", [](const PyFileGetOpen& s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.string);
        if (i == 1)
          return py::cast(s.dflt_file_name);
        if (i == 2)
          return py::cast(s.file_name);
        if (i == 3)
          return py::cast(s.file_unit);
        if (i == 4)
          return py::cast(s.readonly);
        throw py::index_error();
      });
  m.def(
      "file_get_open",
      &python_file_get_open,
      py::arg("string"),
      py::arg("dflt_file_name"),
      py::arg("file_name"),
      py::arg("file_unit"),
      py::arg("readonly"),
      R"""(Parameters
  ----------
  string : 
  dflt_file_name : 
  file_name : 
  file_unit : 
  readonly : 
  )""");
  py::class_<PyFileSuffixer, std::unique_ptr<PyFileSuffixer>>(
      m, "FileSuffixer", "file_suffixer return type")
      .def_readonly("in_file_name", &PyFileSuffixer::in_file_name)
      .def_readonly("out_file_name", &PyFileSuffixer::out_file_name)
      .def_readonly("suffix", &PyFileSuffixer::suffix)
      .def_readonly("add_switch", &PyFileSuffixer::add_switch)
      .def("__len__", [](const PyFileSuffixer&) { return 4; })
      .def("__getitem__", [](const PyFileSuffixer& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.in_file_name);
        if (i == 1)
          return py::cast(s.out_file_name);
        if (i == 2)
          return py::cast(s.suffix);
        if (i == 3)
          return py::cast(s.add_switch);
        throw py::index_error();
      });
  m.def(
      "file_suffixer",
      &python_file_suffixer,
      py::arg("in_file_name"),
      py::arg("out_file_name"),
      py::arg("suffix"),
      py::arg("add_switch"),
      R"""(Parameters
  ----------
  in_file_name : 
  out_file_name : 
  suffix : 
  add_switch : 
  )""");
  py::class_<PyFindLocationInt, std::unique_ptr<PyFindLocationInt>>(
      m, "FindLocationInt", "find_location_int return type")
      .def_readonly("value", &PyFindLocationInt::value)
      .def_readonly("ix_match", &PyFindLocationInt::ix_match)
      .def("__len__", [](const PyFindLocationInt&) { return 2; })
      .def("__getitem__", [](const PyFindLocationInt& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.value);
        if (i == 1)
          return py::cast(s.ix_match);
        throw py::index_error();
      });
  m.def(
      "find_location",
      py::overload_cast<IntAlloc1D&, int, int>(&python_find_location_int),
      py::arg("arr"),
      py::arg("value"),
      py::arg("ix_match"),
      R"""(Parameters
  ----------
  arr : 
  value : 
  ix_match : 
  )""");
  py::class_<PyFindLocationLogic, std::unique_ptr<PyFindLocationLogic>>(
      m, "FindLocationLogic", "find_location_logic return type")
      .def_readonly("value", &PyFindLocationLogic::value)
      .def_readonly("ix_match", &PyFindLocationLogic::ix_match)
      .def("__len__", [](const PyFindLocationLogic&) { return 2; })
      .def(
          "__getitem__", [](const PyFindLocationLogic& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.value);
            if (i == 1)
              return py::cast(s.ix_match);
            throw py::index_error();
          });
  m.def(
      "find_location",
      py::overload_cast<BoolAlloc1D&, bool, int>(&python_find_location_logic),
      py::arg("arr"),
      py::arg("value"),
      py::arg("ix_match"),
      R"""(Parameters
  ----------
  arr : 
  value : 
  ix_match : 
  )""");
  py::class_<PyFindLocationReal, std::unique_ptr<PyFindLocationReal>>(
      m, "FindLocationReal", "find_location_real return type")
      .def_readonly("ix_match", &PyFindLocationReal::ix_match)
      .def("__len__", [](const PyFindLocationReal&) { return 1; })
      .def("__getitem__", [](const PyFindLocationReal& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.ix_match);
        throw py::index_error();
      });
  m.def(
      "find_location",
      py::overload_cast<RealAlloc1D&, double, int>(&python_find_location_real),
      py::arg("arr"),
      py::arg("value"),
      py::arg("ix_match"),
      R"""(Parameters
  ----------
  arr : 
      real(rp), logical, or integer
  value : unknown
      :).
  ix_match : 
  )""");
  m.def(
      "fine_frequency_estimate",
      &SimUtils::fine_frequency_estimate,
      py::arg("data"),
      R"""(Function fine_frequency_estimate(data) result(frequency)

  Uses Laskar's method to accurately find the most dominant frequency
  A coarse estimate is first made by FFT.

  Parameters
  ----------
  data : float
      data to analyze

  Returns
  -------
  frequency : float
      Frequency corresponding to the largest FFT amplitude
  )""");
  py::class_<PyFixedwindowls, std::unique_ptr<PyFixedwindowls>>(
      m, "Fixedwindowls", "fixedwindowls return type")
      .def_readonly("z", &PyFixedwindowls::z)
      .def("__len__", [](const PyFixedwindowls&) { return 1; })
      .def("__getitem__", [](const PyFixedwindowls& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.z);
        throw py::index_error();
      });
  m.def(
      "fixedwindowls",
      &python_fixedwindowls,
      py::arg("ynew"),
      py::arg("id"),
      py::arg("z"),
      R"""(Function fixedWindowLS

  Main function of the windowLS modult.  Each call to this function adds a data point to the fit
  and returns the derivative evaluated at the end of the window.  It is assumed that all data points
  are separeted by the same interval.
  This module is initialized with zeros for all data points, and so the results are unreliable until
  a number of data points equal to N has been entered.

  initFixedWindowLS must be called prior to calling this function.  destFixedWindowLS should be
  called when the instance is no longer needed.

  )""");
  py::class_<
      SimUtils::FourierAmplitude,
      std::unique_ptr<SimUtils::FourierAmplitude>>(
      m, "FourierAmplitude", "fourier_amplitude return type")
      .def_readonly("cos_amp", &SimUtils::FourierAmplitude::cos_amp)
      .def_readonly("sin_amp", &SimUtils::FourierAmplitude::sin_amp)
      .def_readonly("dcos_amp", &SimUtils::FourierAmplitude::dcos_amp)
      .def_readonly("dsin_amp", &SimUtils::FourierAmplitude::dsin_amp)
      .def("__len__", [](const SimUtils::FourierAmplitude&) { return 4; })
      .def(
          "__getitem__",
          [](const SimUtils::FourierAmplitude& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.cos_amp);
            if (i == 1)
              return py::cast(s.sin_amp);
            if (i == 2)
              return py::cast(s.dcos_amp);
            if (i == 3)
              return py::cast(s.dsin_amp);
            throw py::index_error();
          });
  m.def(
      "fourier_amplitude",
      &SimUtils::fourier_amplitude,
      py::arg("data"),
      py::arg("frequency"),
      R"""(Subroutine fourier_amplitude(data, frequency, cos_amp, sin_amp, dcos_amp, dsin_amp)

  Computes cos_amp = (1/N) * sum_n=0^{N-1} data(n-1) cos(twopi*frequency*n)
      and  sin_amp = (1/N) * sum_n=0^{N-1} data(n-1) sin(twopi*frequency*n)
      and optionally dcos_amp = d/dfrequency cos_amp
                     dsin_amp = d/dfrequency sin_amp

  Parameters
  ----------
  data : float
      data to analyze
  frequency : float
      frequency

  Returns
  -------
  cos_amp : float
      cosine amplitude
  sin_amp : float
      sine amplitude
  dcos_amp : float
      cosine amplitude derivative
  dsin_amp : float
      sine amplitude derivative
  )""");
}
