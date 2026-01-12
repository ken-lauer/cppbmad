#include "pybmad/generated/SimUtils_routines_i.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyIBessel {
  int m;
  double arg;
  double i_bes;
};
PyIBessel python_i_bessel(int m, double arg, double i_bes) {
  SimUtils::i_bessel(m, arg, i_bes);
  auto py_result{PyIBessel{m, arg, i_bes}};
  return py_result;
}
struct PyIBesselExtended {
  int m;
  double arg;
  std::complex<double> i_bes;
};
PyIBesselExtended python_i_bessel_extended(
    int m,
    double arg,
    std::complex<double> i_bes) {
  SimUtils::i_bessel_extended(m, arg, i_bes);
  auto py_result{PyIBesselExtended{m, arg, i_bes}};
  return py_result;
}
struct PyIncrementFileNumber {
  std::string file_name;
  int digits;
  int number;
  std::string cnumber;
};
PyIncrementFileNumber python_increment_file_number(
    std::string file_name,
    int digits,
    int number,
    std::string cnumber) {
  SimUtils::increment_file_number(file_name, digits, number, cnumber);
  auto py_result{PyIncrementFileNumber{file_name, digits, number, cnumber}};
  return py_result;
}
struct PyIndexNocase {
  std::string string1;
  std::string string2;
  int indx;
};
PyIndexNocase python_index_nocase(
    std::string string1,
    std::string string2,
    int indx) {
  SimUtils::index_nocase(string1, string2, indx);
  auto py_result{PyIndexNocase{string1, string2, indx}};
  return py_result;
}
struct PyIntStr {
  int int_;
  std::optional<int> width;
  std::string str;
};
PyIntStr python_int_str(int int_, std::optional<int> width, std::string str) {
  SimUtils::int_str(int_, make_opt_ref(width), str);
  auto py_result{PyIntStr{int_, width, str}};
  return py_result;
}
struct PyInterpolatedFft {
  bool calc_ok;
  std::optional<int> opt_dump_spectrum;
  std::optional<int> opt_dump_index;
  double this_fft;
};
PyInterpolatedFft python_interpolated_fft(
    ComplexAlloc1D& cdata,
    bool calc_ok,
    std::optional<int> opt_dump_spectrum,
    std::optional<int> opt_dump_index,
    double this_fft) {
  SimUtils::interpolated_fft(
      cdata,
      calc_ok,
      make_opt_ref(opt_dump_spectrum),
      make_opt_ref(opt_dump_index),
      this_fft);
  auto py_result{
      PyInterpolatedFft{calc_ok, opt_dump_spectrum, opt_dump_index, this_fft}};
  return py_result;
}
struct PyInterpolatedFftGsl {
  bool calc_ok;
  std::optional<int> opt_dump_spectrum;
  std::optional<int> opt_dump_index;
  double this_fft;
};
PyInterpolatedFftGsl python_interpolated_fft_gsl(
    ComplexAlloc1D& cdata,
    bool calc_ok,
    std::optional<int> opt_dump_spectrum,
    std::optional<int> opt_dump_index,
    double this_fft) {
  SimUtils::interpolated_fft_gsl(
      cdata,
      calc_ok,
      make_opt_ref(opt_dump_spectrum),
      make_opt_ref(opt_dump_index),
      this_fft);
  auto py_result{PyInterpolatedFftGsl{
      calc_ok, opt_dump_spectrum, opt_dump_index, this_fft}};
  return py_result;
}
struct PyIsAlphabetic {
  std::string string;
  std::optional<std::string> valid_chars;
  bool is_alpha;
};
PyIsAlphabetic python_is_alphabetic(
    std::string string,
    std::optional<std::string> valid_chars,
    bool is_alpha) {
  SimUtils::is_alphabetic(string, make_opt_ref(valid_chars), is_alpha);
  auto py_result{PyIsAlphabetic{string, valid_chars, is_alpha}};
  return py_result;
}
struct PyIsDecreasingSequence {
  bool is_decreasing;
};
PyIsDecreasingSequence python_is_decreasing_sequence(
    RealAlloc1D& array,
    std::optional<bool> strict,
    bool is_decreasing) {
  SimUtils::is_decreasing_sequence(array, strict, is_decreasing);
  auto py_result{PyIsDecreasingSequence{is_decreasing}};
  return py_result;
}
struct PyIsIncreasingSequence {
  bool is_increasing;
};
PyIsIncreasingSequence python_is_increasing_sequence(
    RealAlloc1D& array,
    std::optional<bool> strict,
    bool is_increasing) {
  SimUtils::is_increasing_sequence(array, strict, is_increasing);
  auto py_result{PyIsIncreasingSequence{is_increasing}};
  return py_result;
}
struct PyIsInteger {
  std::string string;
  std::optional<int> int_;
  std::optional<std::string> delims;
  std::optional<int> ix_word;
  bool valid;
};
PyIsInteger python_is_integer(
    std::string string,
    std::optional<int> int_,
    std::optional<std::string> delims,
    std::optional<int> ix_word,
    bool valid) {
  SimUtils::is_integer(
      string,
      make_opt_ref(int_),
      make_opt_ref(delims),
      make_opt_ref(ix_word),
      valid);
  auto py_result{PyIsInteger{string, int_, delims, ix_word, valid}};
  return py_result;
}
struct PyIsLogical {
  std::string string;
  std::optional<bool> ignore;
  bool valid;
};
PyIsLogical python_is_logical(
    std::string string,
    std::optional<bool> ignore,
    bool valid) {
  SimUtils::is_logical(string, make_opt_ref(ignore), valid);
  auto py_result{PyIsLogical{string, ignore, valid}};
  return py_result;
}
struct PyIsReal {
  std::string string;
  std::optional<bool> ignore;
  std::optional<double> real_num;
  bool valid;
};
PyIsReal python_is_real(
    std::string string,
    std::optional<bool> ignore,
    std::optional<double> real_num,
    bool valid) {
  SimUtils::is_real(
      string, make_opt_ref(ignore), make_opt_ref(real_num), valid);
  auto py_result{PyIsReal{string, ignore, real_num, valid}};
  return py_result;
}

void init_SimUtils_routines_i(py::module& m) {
  m.def(
      "i_bessel",
      &python_i_bessel,
      py::arg("m"),
      py::arg("arg"),
      py::arg("i_bes"),
      R"""(Parameters
  ----------
  m : 
  arg : 
  i_bes : 
  )""");
  py::class_<PyIBessel, std::unique_ptr<PyIBessel>>(
      m, "IBessel", "Fortran routine i_bessel return value")
      .def_readonly("m", &PyIBessel::m)
      .def_readonly("arg", &PyIBessel::arg)
      .def_readonly("i_bes", &PyIBessel::i_bes)
      .def("__len__", [](const PyIBessel&) { return 3; })
      .def("__getitem__", [](const PyIBessel& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.m);
        if (i == 1)
          return py::cast(s.arg);
        if (i == 2)
          return py::cast(s.i_bes);
        throw py::index_error();
      });
  m.def(
      "i_bessel_extended",
      &python_i_bessel_extended,
      py::arg("m"),
      py::arg("arg"),
      py::arg("i_bes"),
      R"""(Parameters
  ----------
  m : 
  arg : 
  i_bes : 
  )""");
  py::class_<PyIBesselExtended, std::unique_ptr<PyIBesselExtended>>(
      m, "IBesselExtended", "Fortran routine i_bessel_extended return value")
      .def_readonly("m", &PyIBesselExtended::m)
      .def_readonly("arg", &PyIBesselExtended::arg)
      .def_readonly("i_bes", &PyIBesselExtended::i_bes)
      .def("__len__", [](const PyIBesselExtended&) { return 3; })
      .def("__getitem__", [](const PyIBesselExtended& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.m);
        if (i == 1)
          return py::cast(s.arg);
        if (i == 2)
          return py::cast(s.i_bes);
        throw py::index_error();
      });
  m.def(
      "increment_file_number",
      &python_increment_file_number,
      py::arg("file_name"),
      py::arg("digits"),
      py::arg("number"),
      py::arg("cnumber"),
      R"""(Parameters
  ----------
  file_name : 
  digits : 
  number : 
  cnumber : 
  )""");
  py::class_<PyIncrementFileNumber, std::unique_ptr<PyIncrementFileNumber>>(
      m,
      "IncrementFileNumber",
      "Fortran routine increment_file_number return value")
      .def_readonly("file_name", &PyIncrementFileNumber::file_name)
      .def_readonly("digits", &PyIncrementFileNumber::digits)
      .def_readonly("number", &PyIncrementFileNumber::number)
      .def_readonly("cnumber", &PyIncrementFileNumber::cnumber)
      .def("__len__", [](const PyIncrementFileNumber&) { return 4; })
      .def(
          "__getitem__",
          [](const PyIncrementFileNumber& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.file_name);
            if (i == 1)
              return py::cast(s.digits);
            if (i == 2)
              return py::cast(s.number);
            if (i == 3)
              return py::cast(s.cnumber);
            throw py::index_error();
          });
  m.def(
      "index_nocase",
      &python_index_nocase,
      py::arg("string1"),
      py::arg("string2"),
      py::arg("indx"),
      R"""(Parameters
  ----------
  string1 : 
  string2 : 
  indx : 
  )""");
  py::class_<PyIndexNocase, std::unique_ptr<PyIndexNocase>>(
      m, "IndexNocase", "Fortran routine index_nocase return value")
      .def_readonly("string1", &PyIndexNocase::string1)
      .def_readonly("string2", &PyIndexNocase::string2)
      .def_readonly("indx", &PyIndexNocase::indx)
      .def("__len__", [](const PyIndexNocase&) { return 3; })
      .def("__getitem__", [](const PyIndexNocase& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.string1);
        if (i == 1)
          return py::cast(s.string2);
        if (i == 2)
          return py::cast(s.indx);
        throw py::index_error();
      });
  m.def(
      "initfixedwindowls",
      &SimUtils::initfixedwindowls,
      py::arg("N"),
      py::arg("dt"),
      py::arg("order"),
      py::arg("der"),
      R"""(Function initFixedWindowLS

  Initializes an instance of the fixed window least squares module.
  See module documentation (getf windowLS_mod) for use details.
  Any instance of windowLS created with this module should be destroyed with destFixedWindowLS.

  Parameters
  ----------
  N : int
      Number of data points to fit over. aka window size.
  dt : float
      Time interval between data points. It is assumed that the data is separated by fixed time intervals.
  order : int
      Order of fit polynomial.  Must be greater than or equal to der.
  der : int
      Order of derivative to be returned. Set der=0 to obtain the fit. <return value>  -- INTEGER: id of
      windowLS instance created.
  )""");
  m.def(
      "int_str",
      &python_int_str,
      py::arg("int_"),
      py::arg("width") = py::none(),
      py::arg("str"),
      R"""(Parameters
  ----------
  int : 
  width : 
  str : 
  )""");
  py::class_<PyIntStr, std::unique_ptr<PyIntStr>>(
      m, "IntStr", "Fortran routine int_str return value")
      .def_readonly("int_", &PyIntStr::int_)
      .def_readonly("width", &PyIntStr::width)
      .def_readonly("str", &PyIntStr::str)
      .def("__len__", [](const PyIntStr&) { return 3; })
      .def("__getitem__", [](const PyIntStr& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.int_);
        if (i == 1)
          return py::cast(s.width);
        if (i == 2)
          return py::cast(s.str);
        throw py::index_error();
      });
  m.def(
      "interpolated_fft",
      &python_interpolated_fft,
      py::arg("cdata"),
      py::arg("calc_ok"),
      py::arg("opt_dump_spectrum") = py::none(),
      py::arg("opt_dump_index") = py::none(),
      py::arg("this_fft"),
      R"""(Function interpolated_fft (cdata, calc_ok, opt_dump_spectrum, opt_dump_index) result (this_fft)

  Windows the complex data and used Numerical Recipes four1 to find the peak in the spectrum.
  The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
  available.


  Returns
  -------
  this_fft
  )""");
  py::class_<PyInterpolatedFft, std::unique_ptr<PyInterpolatedFft>>(
      m, "InterpolatedFft", "Fortran routine interpolated_fft return value")
      .def_readonly("calc_ok", &PyInterpolatedFft::calc_ok)
      .def_readonly("opt_dump_spectrum", &PyInterpolatedFft::opt_dump_spectrum)
      .def_readonly("opt_dump_index", &PyInterpolatedFft::opt_dump_index)
      .def_readonly("this_fft", &PyInterpolatedFft::this_fft)
      .def("__len__", [](const PyInterpolatedFft&) { return 4; })
      .def("__getitem__", [](const PyInterpolatedFft& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.calc_ok);
        if (i == 1)
          return py::cast(s.opt_dump_spectrum);
        if (i == 2)
          return py::cast(s.opt_dump_index);
        if (i == 3)
          return py::cast(s.this_fft);
        throw py::index_error();
      });
  m.def(
      "interpolated_fft_gsl",
      &python_interpolated_fft_gsl,
      py::arg("cdata"),
      py::arg("calc_ok"),
      py::arg("opt_dump_spectrum") = py::none(),
      py::arg("opt_dump_index") = py::none(),
      py::arg("this_fft"),
      R"""(function interpolated_fft_gsl

  Windows the complex data and uses a mixed-radix GSL routine to find the peak in the spectrum.
  The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
  available.

  )""");
  py::class_<PyInterpolatedFftGsl, std::unique_ptr<PyInterpolatedFftGsl>>(
      m,
      "InterpolatedFftGsl",
      "Fortran routine interpolated_fft_gsl return value")
      .def_readonly("calc_ok", &PyInterpolatedFftGsl::calc_ok)
      .def_readonly(
          "opt_dump_spectrum", &PyInterpolatedFftGsl::opt_dump_spectrum)
      .def_readonly("opt_dump_index", &PyInterpolatedFftGsl::opt_dump_index)
      .def_readonly("this_fft", &PyInterpolatedFftGsl::this_fft)
      .def("__len__", [](const PyInterpolatedFftGsl&) { return 4; })
      .def(
          "__getitem__",
          [](const PyInterpolatedFftGsl& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.calc_ok);
            if (i == 1)
              return py::cast(s.opt_dump_spectrum);
            if (i == 2)
              return py::cast(s.opt_dump_index);
            if (i == 3)
              return py::cast(s.this_fft);
            throw py::index_error();
          });
  m.def(
      "is_alphabetic",
      &python_is_alphabetic,
      py::arg("string"),
      py::arg("valid_chars") = py::none(),
      py::arg("is_alpha"),
      R"""(no longer exists
  function inverse_prob (val) result (prob)
    import
    implicit none
    real(rp) prob
    real(rp) val
  end function


  Returns
  -------
  prob
  )""");
  py::class_<PyIsAlphabetic, std::unique_ptr<PyIsAlphabetic>>(
      m, "IsAlphabetic", "Fortran routine is_alphabetic return value")
      .def_readonly("string", &PyIsAlphabetic::string)
      .def_readonly("valid_chars", &PyIsAlphabetic::valid_chars)
      .def_readonly("is_alpha", &PyIsAlphabetic::is_alpha)
      .def("__len__", [](const PyIsAlphabetic&) { return 3; })
      .def("__getitem__", [](const PyIsAlphabetic& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.string);
        if (i == 1)
          return py::cast(s.valid_chars);
        if (i == 2)
          return py::cast(s.is_alpha);
        throw py::index_error();
      });
  m.def(
      "is_decreasing_sequence",
      &python_is_decreasing_sequence,
      py::arg("array"),
      py::arg("strict") = py::none(),
      py::arg("is_decreasing"),
      R"""(Parameters
  ----------
  array : float
      Sequence.
  strict : bool, optional
      If True (default) sequence must be strictly decreasing.
  is_decreasing : 
  )""");
  py::class_<PyIsDecreasingSequence, std::unique_ptr<PyIsDecreasingSequence>>(
      m,
      "IsDecreasingSequence",
      "Fortran routine is_decreasing_sequence return value")
      .def_readonly("is_decreasing", &PyIsDecreasingSequence::is_decreasing)
      .def("__len__", [](const PyIsDecreasingSequence&) { return 1; })
      .def(
          "__getitem__",
          [](const PyIsDecreasingSequence& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_decreasing);
            throw py::index_error();
          });
  m.def(
      "is_false",
      &SimUtils::is_false,
      py::arg("param"),
      R"""(Function is_false (param) result (this_false)

  Routine to translate from a real number to a boolian True or False.
  Translation: 0 = False, nonzero = True

  Also see: is_true and int_logic

  The typical use of this routine is for parameters in ele_struct%value(:) which
  is a real array. Some of the elements in the %value array are used to specify
  boolian attributes. For example, quadrupoles use ele%value(scale_multipoles$).

  Parameters
  ----------
  param : float
      Real number to be translated

  Returns
  -------
  this_false : bool
      Set True if param is zero. False otherwise.

  Notes
  -----
  Related routines:
  is_true int_logic ) which is a real array. Some of the elements in the %value array are used to specify
  boolian attributes. For example quadrupoles use ele%value(scale_multipoles$).
  )""");
  m.def(
      "is_increasing_sequence",
      &python_is_increasing_sequence,
      py::arg("array"),
      py::arg("strict") = py::none(),
      py::arg("is_increasing"),
      R"""(Parameters
  ----------
  array : float
      Sequence.
  strict : bool, optional
      If True (default) sequence must be strictly increasing.
  is_increasing : 
  )""");
  py::class_<PyIsIncreasingSequence, std::unique_ptr<PyIsIncreasingSequence>>(
      m,
      "IsIncreasingSequence",
      "Fortran routine is_increasing_sequence return value")
      .def_readonly("is_increasing", &PyIsIncreasingSequence::is_increasing)
      .def("__len__", [](const PyIsIncreasingSequence&) { return 1; })
      .def(
          "__getitem__",
          [](const PyIsIncreasingSequence& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_increasing);
            throw py::index_error();
          });
  m.def(
      "is_integer",
      &python_is_integer,
      py::arg("string"),
      py::arg("int_") = py::none(),
      py::arg("delims") = py::none(),
      py::arg("ix_word") = py::none(),
      py::arg("valid"),
      R"""(Parameters
  ----------
  string : 
  int : 
  delims : 
  ix_word : 
  valid : 
  )""");
  py::class_<PyIsInteger, std::unique_ptr<PyIsInteger>>(
      m, "IsInteger", "Fortran routine is_integer return value")
      .def_readonly("string", &PyIsInteger::string)
      .def_readonly("int_", &PyIsInteger::int_)
      .def_readonly("delims", &PyIsInteger::delims)
      .def_readonly("ix_word", &PyIsInteger::ix_word)
      .def_readonly("valid", &PyIsInteger::valid)
      .def("__len__", [](const PyIsInteger&) { return 5; })
      .def("__getitem__", [](const PyIsInteger& s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.string);
        if (i == 1)
          return py::cast(s.int_);
        if (i == 2)
          return py::cast(s.delims);
        if (i == 3)
          return py::cast(s.ix_word);
        if (i == 4)
          return py::cast(s.valid);
        throw py::index_error();
      });
  m.def(
      "is_logical",
      &python_is_logical,
      py::arg("string"),
      py::arg("ignore") = py::none(),
      py::arg("valid"),
      R"""(Parameters
  ----------
  string : 
  ignore : 
  valid : 
  )""");
  py::class_<PyIsLogical, std::unique_ptr<PyIsLogical>>(
      m, "IsLogical", "Fortran routine is_logical return value")
      .def_readonly("string", &PyIsLogical::string)
      .def_readonly("ignore", &PyIsLogical::ignore)
      .def_readonly("valid", &PyIsLogical::valid)
      .def("__len__", [](const PyIsLogical&) { return 3; })
      .def("__getitem__", [](const PyIsLogical& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.string);
        if (i == 1)
          return py::cast(s.ignore);
        if (i == 2)
          return py::cast(s.valid);
        throw py::index_error();
      });
  m.def(
      "is_real",
      &python_is_real,
      py::arg("string"),
      py::arg("ignore") = py::none(),
      py::arg("real_num") = py::none(),
      py::arg("valid"),
      R"""(Parameters
  ----------
  string : 
  ignore : 
  real_num : 
  valid : 
  )""");
  py::class_<PyIsReal, std::unique_ptr<PyIsReal>>(
      m, "IsReal", "Fortran routine is_real return value")
      .def_readonly("string", &PyIsReal::string)
      .def_readonly("ignore", &PyIsReal::ignore)
      .def_readonly("real_num", &PyIsReal::real_num)
      .def_readonly("valid", &PyIsReal::valid)
      .def("__len__", [](const PyIsReal&) { return 4; })
      .def("__getitem__", [](const PyIsReal& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.string);
        if (i == 1)
          return py::cast(s.ignore);
        if (i == 2)
          return py::cast(s.real_num);
        if (i == 3)
          return py::cast(s.valid);
        throw py::index_error();
      });
  m.def(
      "is_subatomic_species",
      &SimUtils::is_subatomic_species,
      py::arg("species"),
      R"""(Function is_subatomic_species(species) result (is_subatomic)

  Routine to return True if species argument corresponds to a subatomic particle.

  Parameters
  ----------
  species : int
      Spicies ID.

  Returns
  -------
  is_subatomic : bool
      Set True if species corresponds to a subatomic particle.
  )""");
  m.def(
      "is_true",
      &SimUtils::is_true,
      py::arg("param"),
      R"""(Function is_true (param) result (this_true)

  Routine to translate from a real number to a boolian True or False.
  Translation: 0 = False, nonzero = True

  Also see: is_false and int_logic

  The typical use of this routine is for parameters in ele_struct%value(:) which
  is a real array. Some of the elements in the %value array are used to specify
  boolian attributes. For example, quadrupoles use ele%value(scale_multipoles$).

  Parameters
  ----------
  param : float
      Real number to be translated

  Returns
  -------
  this_true : bool
      Set False if param is zero. True otherwise.

  Notes
  -----
  Related routines:
  is_false int_logic ) which is a real array. Some of the elements in the %value array are used to specify
  boolian attributes. For example quadrupoles use ele%value(scale_multipoles$).
  )""");
}
