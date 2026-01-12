#include "pybmad/generated/SimUtils_routines_s.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PySetParameterInt {
  int param_val;
  int set_val;
  int save_val;
};
PySetParameterInt python_set_parameter_int(
    int param_val,
    int set_val,
    int save_val) {
  SimUtils::set_parameter(param_val, set_val, save_val);
  auto py_result{PySetParameterInt{param_val, set_val, save_val}};
  return py_result;
}
struct PySetParameterLogic {
  bool param_val;
  bool set_val;
  bool save_val;
};
PySetParameterLogic python_set_parameter_logic(
    bool param_val,
    bool set_val,
    bool save_val) {
  SimUtils::set_parameter(param_val, set_val, save_val);
  auto py_result{PySetParameterLogic{param_val, set_val, save_val}};
  return py_result;
}
struct PySetParameterReal {
  double param_val;
  double set_val;
  double save_val;
};
PySetParameterReal python_set_parameter_real(
    double param_val,
    double set_val,
    double save_val) {
  SimUtils::set_parameter(param_val, set_val, save_val);
  auto py_result{PySetParameterReal{param_val, set_val, save_val}};
  return py_result;
}
struct PySinc {
  double y;
};
PySinc python_sinc(double x, std::optional<int> nd, double y) {
  SimUtils::sinc(x, nd, y);
  auto py_result{PySinc{y}};
  return py_result;
}
struct PySincc {
  double y;
};
PySincc python_sincc(double x, std::optional<int> nd, double y) {
  SimUtils::sincc(x, nd, y);
  auto py_result{PySincc{y}};
  return py_result;
}
struct PySinhxX {
  double y;
};
PySinhxX python_sinhx_x(double x, std::optional<int> nd, double y) {
  SimUtils::sinhx_x(x, nd, y);
  auto py_result{PySinhxX{y}};
  return py_result;
}
struct PySkipHeader {
  int ix_unit;
  bool error_flag;
};
PySkipHeader python_skip_header(int ix_unit, bool error_flag) {
  SimUtils::skip_header(ix_unit, error_flag);
  auto py_result{PySkipHeader{ix_unit, error_flag}};
  return py_result;
}
struct PySqrtAlpha {
  double y;
};
PySqrtAlpha python_sqrt_alpha(double alpha, double x, double y) {
  SimUtils::sqrt_alpha(alpha, x, y);
  auto py_result{PySqrtAlpha{y}};
  return py_result;
}
struct PySqrtOne {
  double ds1;
};
PySqrtOne python_sqrt_one(double x, std::optional<int> nd, double ds1) {
  SimUtils::sqrt_one(x, nd, ds1);
  auto py_result{PySqrtOne{ds1}};
  return py_result;
}
struct PyStrCount {
  std::string str;
  std::string match;
  int num;
};
PyStrCount python_str_count(std::string str, std::string match, int num) {
  SimUtils::str_count(str, match, num);
  auto py_result{PyStrCount{str, match, num}};
  return py_result;
}
struct PyStrFirstInSet {
  std::string line;
  std::string set;
  std::optional<bool> ignore_clauses;
  int ix_match;
};
PyStrFirstInSet python_str_first_in_set(
    std::string line,
    std::string set,
    std::optional<bool> ignore_clauses,
    int ix_match) {
  SimUtils::str_first_in_set(line, set, make_opt_ref(ignore_clauses), ix_match);
  auto py_result{PyStrFirstInSet{line, set, ignore_clauses, ix_match}};
  return py_result;
}
struct PyStrFirstNotInSet {
  std::string line;
  std::string set;
  int ix_match;
};
PyStrFirstNotInSet python_str_first_not_in_set(
    std::string line,
    std::string set,
    int ix_match) {
  SimUtils::str_first_not_in_set(line, set, ix_match);
  auto py_result{PyStrFirstNotInSet{line, set, ix_match}};
  return py_result;
}
struct PyStrLastInSet {
  std::string line;
  std::string set;
  int ix_match;
};
PyStrLastInSet python_str_last_in_set(
    std::string line,
    std::string set,
    int ix_match) {
  SimUtils::str_last_in_set(line, set, ix_match);
  auto py_result{PyStrLastInSet{line, set, ix_match}};
  return py_result;
}
struct PyStrLastNotInSet {
  std::string line;
  std::string set;
  int ix_match;
};
PyStrLastNotInSet python_str_last_not_in_set(
    std::string line,
    std::string set,
    int ix_match) {
  SimUtils::str_last_not_in_set(line, set, ix_match);
  auto py_result{PyStrLastNotInSet{line, set, ix_match}};
  return py_result;
}
struct PyStrMatchWild {
  std::string str;
  std::string pat;
  bool a_match;
};
PyStrMatchWild python_str_match_wild(
    std::string str,
    std::string pat,
    bool a_match) {
  SimUtils::str_match_wild(str, pat, a_match);
  auto py_result{PyStrMatchWild{str, pat, a_match}};
  return py_result;
}
struct PyStrSubstitute {
  std::string string;
  std::optional<std::string> str_match;
  std::optional<std::string> str_replace;
  std::optional<bool> do_trim;
  std::optional<bool> ignore_escaped;
};
PyStrSubstitute python_str_substitute(
    std::string string,
    std::optional<std::string> str_match = std::nullopt,
    std::optional<std::string> str_replace = std::nullopt,
    std::optional<bool> do_trim = std::nullopt,
    std::optional<bool> ignore_escaped = std::nullopt) {
  SimUtils::str_substitute(
      string,
      make_opt_ref(str_match),
      make_opt_ref(str_replace),
      make_opt_ref(do_trim),
      make_opt_ref(ignore_escaped));
  auto py_result{
      PyStrSubstitute{string, str_match, str_replace, do_trim, ignore_escaped}};
  return py_result;
}
struct PyStringToInt {
  std::string line;
  int default_;
  bool err_flag;
  std::optional<bool> err_print_flag;
  int value;
};
PyStringToInt python_string_to_int(
    std::string line,
    int default_,
    bool err_flag,
    std::optional<bool> err_print_flag,
    int value) {
  SimUtils::string_to_int(
      line, default_, err_flag, make_opt_ref(err_print_flag), value);
  auto py_result{
      PyStringToInt{line, default_, err_flag, err_print_flag, value}};
  return py_result;
}
struct PyStringToReal {
  std::string line;
  double default_;
  bool err_flag;
  std::optional<bool> err_print_flag;
  double value;
};
PyStringToReal python_string_to_real(
    std::string line,
    double default_,
    bool err_flag,
    std::optional<bool> err_print_flag,
    double value) {
  SimUtils::string_to_real(
      line, default_, err_flag, make_opt_ref(err_print_flag), value);
  auto py_result{
      PyStringToReal{line, default_, err_flag, err_print_flag, value}};
  return py_result;
}
struct PyStringTrim {
  std::string in_string;
  std::string out_string;
  int word_len;
};
PyStringTrim python_string_trim(
    std::string in_string,
    std::string out_string,
    int word_len) {
  SimUtils::string_trim(in_string, out_string, word_len);
  auto py_result{PyStringTrim{in_string, out_string, word_len}};
  return py_result;
}
struct PyStringTrim2 {
  std::string in_str;
  std::string delimitors;
  std::string out_str;
  int ix_word;
  std::string delim;
  int ix_next;
};
PyStringTrim2 python_string_trim2(
    std::string in_str,
    std::string delimitors,
    std::string out_str,
    int ix_word,
    std::string delim,
    int ix_next) {
  SimUtils::string_trim2(in_str, delimitors, out_str, ix_word, delim, ix_next);
  auto py_result{
      PyStringTrim2{in_str, delimitors, out_str, ix_word, delim, ix_next}};
  return py_result;
}
struct PySystemCommand {
  std::string line;
  std::optional<bool> err_flag;
};
PySystemCommand python_system_command(
    std::string line,
    std::optional<bool> err_flag = std::nullopt) {
  SimUtils::system_command(line, make_opt_ref(err_flag));
  auto py_result{PySystemCommand{line, err_flag}};
  return py_result;
}

void init_SimUtils_routines_s(py::module& m) {
  m.def(
      "set_parameter",
      py::overload_cast<int, int, int>(&python_set_parameter_int),
      py::arg("param_val"),
      py::arg("set_val"),
      py::arg("save_val"),
      R"""(Parameters
  ----------
  param_val : 
  set_val : 
  save_val : 
  )""");
  py::class_<PySetParameterInt, std::unique_ptr<PySetParameterInt>>(
      m, "SetParameterInt", "Fortran routine set_parameter_int return value")
      .def_readonly("param_val", &PySetParameterInt::param_val)
      .def_readonly("set_val", &PySetParameterInt::set_val)
      .def_readonly("save_val", &PySetParameterInt::save_val)
      .def("__len__", [](const PySetParameterInt&) { return 3; })
      .def("__getitem__", [](const PySetParameterInt& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.param_val);
        if (i == 1)
          return py::cast(s.set_val);
        if (i == 2)
          return py::cast(s.save_val);
        throw py::index_error();
      });
  m.def(
      "set_parameter",
      py::overload_cast<bool, bool, bool>(&python_set_parameter_logic),
      py::arg("param_val"),
      py::arg("set_val"),
      py::arg("save_val"),
      R"""(Parameters
  ----------
  param_val : 
  set_val : 
  save_val : 
  )""");
  py::class_<PySetParameterLogic, std::unique_ptr<PySetParameterLogic>>(
      m,
      "SetParameterLogic",
      "Fortran routine set_parameter_logic return value")
      .def_readonly("param_val", &PySetParameterLogic::param_val)
      .def_readonly("set_val", &PySetParameterLogic::set_val)
      .def_readonly("save_val", &PySetParameterLogic::save_val)
      .def("__len__", [](const PySetParameterLogic&) { return 3; })
      .def(
          "__getitem__", [](const PySetParameterLogic& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.param_val);
            if (i == 1)
              return py::cast(s.set_val);
            if (i == 2)
              return py::cast(s.save_val);
            throw py::index_error();
          });
  m.def(
      "set_parameter",
      py::overload_cast<double, double, double>(&python_set_parameter_real),
      py::arg("param_val"),
      py::arg("set_val"),
      py::arg("save_val"),
      R"""(Parameters
  ----------
  param_val : 
  set_val : 
  save_val : 
  )""");
  py::class_<PySetParameterReal, std::unique_ptr<PySetParameterReal>>(
      m, "SetParameterReal", "Fortran routine set_parameter_real return value")
      .def_readonly("param_val", &PySetParameterReal::param_val)
      .def_readonly("set_val", &PySetParameterReal::set_val)
      .def_readonly("save_val", &PySetParameterReal::save_val)
      .def("__len__", [](const PySetParameterReal&) { return 3; })
      .def("__getitem__", [](const PySetParameterReal& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.param_val);
        if (i == 1)
          return py::cast(s.set_val);
        if (i == 2)
          return py::cast(s.save_val);
        throw py::index_error();
      });
  m.def(
      "set_species_charge",
      &SimUtils::set_species_charge,
      py::arg("species_in"),
      py::arg("charge"),
      R"""(Function set_species_charge(species_in, charge) result(species_charged)

  Routine to return the ID for a particle of the same type as species_in but with a different charge.
  Exception: If species_in corresponds to a subatomic particle, the charge argument is ignored and
  species_charged will be set equal to species_in.

  Parameters
  ----------
  species_in : int
      Input species.
  charge : int
      Charge to set species_charged to.

  Returns
  -------
  species_charged : int
      Species of the same type as species_in but with different charge.
  )""");
  m.def(
      "sinc",
      &python_sinc,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::arg("y"),
      R"""(Parameters
  ----------
  x : 
  nd : 
  y : 
  )""");
  py::class_<PySinc, std::unique_ptr<PySinc>>(
      m, "Sinc", "Fortran routine sinc return value")
      .def_readonly("y", &PySinc::y)
      .def("__len__", [](const PySinc&) { return 1; })
      .def("__getitem__", [](const PySinc& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.y);
        throw py::index_error();
      });
  m.def(
      "sincc",
      &python_sincc,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::arg("y"),
      R"""(Parameters
  ----------
  x : 
  nd : 
  y : 
  )""");
  py::class_<PySincc, std::unique_ptr<PySincc>>(
      m, "Sincc", "Fortran routine sincc return value")
      .def_readonly("y", &PySincc::y)
      .def("__len__", [](const PySincc&) { return 1; })
      .def("__getitem__", [](const PySincc& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.y);
        throw py::index_error();
      });
  m.def(
      "sinhx_x",
      &python_sinhx_x,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::arg("y"),
      R"""(Parameters
  ----------
  x : 
  nd : 
  y : 
  )""");
  py::class_<PySinhxX, std::unique_ptr<PySinhxX>>(
      m, "SinhxX", "Fortran routine sinhx_x return value")
      .def_readonly("y", &PySinhxX::y)
      .def("__len__", [](const PySinhxX&) { return 1; })
      .def("__getitem__", [](const PySinhxX& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.y);
        throw py::index_error();
      });
  m.def(
      "skip_header",
      &python_skip_header,
      py::arg("ix_unit"),
      py::arg("error_flag"),
      R"""(Parameters
  ----------
  ix_unit : 
  error_flag : 
  )""");
  py::class_<PySkipHeader, std::unique_ptr<PySkipHeader>>(
      m, "SkipHeader", "Fortran routine skip_header return value")
      .def_readonly("ix_unit", &PySkipHeader::ix_unit)
      .def_readonly("error_flag", &PySkipHeader::error_flag)
      .def("__len__", [](const PySkipHeader&) { return 2; })
      .def("__getitem__", [](const PySkipHeader& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.ix_unit);
        if (i == 1)
          return py::cast(s.error_flag);
        throw py::index_error();
      });
  m.def(
      "species_id",
      &SimUtils::species_id,
      py::arg("name"),
      py::arg("default_") = py::none(),
      py::arg("print_err") = py::none(),
      R"""(Function species_id (name, default, print_err) result(species)

  Routine to return the integer ID index of a particle species given the name.

  For subatomic particles, the case does not matter.
  For all other types of particles, the case does matter.

  Parameters
  ----------
  name : unknown
      Name of the species.
  default : int, optional
      Default species to use if name is blank or 'ref_species'. If not present, a blank name is an error.
  print_err : bool, optional
      Print error message? Default is True. If False, return species = invalid$,

  Returns
  -------
  species : int
      Species ID. Will return invalid$ if name is not valid. Will return not_set$ if name is blank
  )""");
  m.def(
      "species_id_from_openpmd",
      &SimUtils::species_id_from_openpmd,
      py::arg("pmd_name"),
      py::arg("charge"),
      R"""(Function species_id_from_openpmd (pmd_name, charge) result(species)

  Routine to return the Bmad species ID given the openPMD species name and given particle charge.
  Note: If pmd_name corresponds to a subatomic particle, the charge argument is ignored.

  Parameters
  ----------
  pmd_name : unknown
      OpenPMD species name.
  charge : int
      Species charge. Ignored for subatomic particles.

  Returns
  -------
  species : int
      Bmad spicies ID number.
  )""");
  m.def(
      "species_name",
      &SimUtils::species_name,
      py::arg("species"),
      R"""(Function species_name (species) result(name)

  Routine to return the name of a particle species given the integer index.

  Parameters
  ----------
  species : int
      Species ID.

  Returns
  -------
  name : unknown
      Name of the species. Will return 'INVALID!' (= invalid_name) if index is not valid.
  )""");
  m.def(
      "species_of",
      &SimUtils::species_of,
      py::arg("mass"),
      py::arg("charge"),
      R"""(Function species_of (mass, charge) result (species)

  Routine to return the integer ID index of a particle species given the mass and charge.
  Note: Currently this routine only works for subatomic particles and is used for decoding PTC flat files.

  Parameters
  ----------
  mass : float
      Mass of the particle
  charge : int
      Charge of the particle.

  Returns
  -------
  species : int
      Species ID. Will return invalid$ if name is not valid.
  )""");
  m.def(
      "spin_of",
      &SimUtils::spin_of,
      py::arg("species"),
      py::arg("non_subatomic_default") = py::none(),
      R"""(Function spin_of (species, non_subatomic_default) result (spin)

  Routine to return the spin, in units of hbar, of a particle.
  This routine is only valid for subatomic particles.
  For all other particles, the returned spin value will be the value of non_subatomic_default.

  Parameters
  ----------
  species : int
      Species ID.
  non_subatomic_default : float, optional
      Default value to be used for non-subatomic species. Default value of this argument is zero.

  Returns
  -------
  spin : float
      Particle spin.
  )""");
  m.def(
      "spline1",
      &SimUtils::spline1,
      py::arg("a_spline"),
      py::arg("x"),
      py::arg("n") = py::none(),
      R"""(Function spline1 (a_spline, x, n) result (y)

  Function for spline evaluation using a single spline (instead of a spline array).

  Parameters
  ----------
  a_spline : SplineStruct
      Single spline structure.
  x : float
      Point for evaluation.
  n : int, optional
      Output derivative order. May be -1, 0, 1, 2, or 3. Default is 0. n = -1 => output is integral of y from
      a_spline.x0 to x. n = 1 => output is dy/dx, n = 2 => output is d^2y/dx^2, etc.

  Returns
  -------
  y : float
      Interpolated spline value or derivative.

  Notes
  -----
  Related routines:
  spline_evaluate spline_akima_interpolate use spline_mod
  )""");
  m.def(
      "spline_akima",
      &SimUtils::spline_akima,
      py::arg("spline"),
      R"""(Subroutine spline_akima (spline, ok)

  Given a set of (x,y) points we want to interpolate between the points.
  This subroutine computes the semi-hermite cubic spline developed by
  Hiroshi Akima. The spline goes thorugh all the points (that is, it is
  not a smoothing spline). For interpolation use:
    spline_evaluate
    spline_akima_interpolate ! You do not need to call spline_akima if you use this routine.

  Reference:
    H Akima, "A New Method of Interpolation and Smooth Curve Fitting Based
    on Local Procedures", J. Assoc. Comp. Mach., Vol 17(4), 589-602 (1970).

  Modules used:
    use spline_mod

  Parameters
  ----------
  spline : SplineStruct
      .x0  -- X-component of a point. Note: points must be in assending order. .y0  -- Y-component of a point.

  Returns
  -------
  ok : bool
      Set .false. if something is wrong (like less than 2 points used).
  )""");
  m.def(
      "spline_akima_interpolate",
      &SimUtils::spline_akima_interpolate,
      py::arg("x_knot"),
      py::arg("y_knot"),
      py::arg("x"),
      R"""(Subroutine spline_akima_interpolate (x_knot, y_knot, x, ok, y, dy)

  Routine to interpolate using an akima spline.

  When evaluating at enough points, this routine is slower than calling spline_akima to
  first evaluate the spline coefficients and then repeatedly calling spline_evaluate.

  The advantage of this routine is that only the (x, y) knot points need to be stored
  and it will be faster if the number of evaluations is small.

  This routine will extrapolate past the range of x_knot(:) up to a distance equal to the
  length between an end point and the point just inside the end point.

  Parameters
  ----------
  x_knot : float
      Array of x values for the knot points. Must have more than 2 points and be in asending order.
  y_knot : float
      Array of y values for the knot points. Must be same size as x_knot(:).
  x : float
      Point to evaluate at.

  Returns
  -------
  ok : bool
      Set .true. if everything ok, That is, x is within the spline range.
  y : float
      Spline interpolation.
  dy : float
      Spline derivative interpolation.
  )""");
  py::class_<
      SimUtils::SplineAkimaInterpolate,
      std::unique_ptr<SimUtils::SplineAkimaInterpolate>>(
      m,
      "SplineAkimaInterpolate",
      "Fortran routine spline_akima_interpolate return value")
      .def_readonly("ok", &SimUtils::SplineAkimaInterpolate::ok)
      .def_readonly("y", &SimUtils::SplineAkimaInterpolate::y)
      .def_readonly("dy", &SimUtils::SplineAkimaInterpolate::dy)
      .def("__len__", [](const SimUtils::SplineAkimaInterpolate&) { return 3; })
      .def(
          "__getitem__",
          [](const SimUtils::SplineAkimaInterpolate& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.ok);
            if (i == 1)
              return py::cast(s.y);
            if (i == 2)
              return py::cast(s.dy);
            throw py::index_error();
          });
  m.def(
      "spline_evaluate",
      &SimUtils::spline_evaluate,
      py::arg("spline"),
      py::arg("x"),
      R"""(Subroutine spline_evaluate (spline, x, ok, y, dy)

  Subroutine to evalueate a spline at a set of points.

  A point outside of the range of knot points is an error.

  Parameters
  ----------
  spline : SplineStruct
      Spline structure.
  x : float
      point for evaluation.

  Returns
  -------
  ok : bool
      Set .true. if everything ok. That is, x is within the spline range.
  y : float
      Spline interpolation.
  dy : float
      Spline derivative interpolation.

  Notes
  -----
  Related routines:
  spline1 spline_akima_interpolate A spline may be generated using for example the spline_akima routine. use
  spline_mod
  )""");
  py::class_<
      SimUtils::SplineEvaluate,
      std::unique_ptr<SimUtils::SplineEvaluate>>(
      m, "SplineEvaluate", "Fortran routine spline_evaluate return value")
      .def_readonly("ok", &SimUtils::SplineEvaluate::ok)
      .def_readonly("y", &SimUtils::SplineEvaluate::y)
      .def_readonly("dy", &SimUtils::SplineEvaluate::dy)
      .def("__len__", [](const SimUtils::SplineEvaluate&) { return 3; })
      .def(
          "__getitem__",
          [](const SimUtils::SplineEvaluate& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.ok);
            if (i == 1)
              return py::cast(s.y);
            if (i == 2)
              return py::cast(s.dy);
            throw py::index_error();
          });
  m.def(
      "sqrt_alpha",
      &python_sqrt_alpha,
      py::arg("alpha"),
      py::arg("x"),
      py::arg("y"),
      R"""(Parameters
  ----------
  alpha : 
  x : 
  y : 
  )""");
  py::class_<PySqrtAlpha, std::unique_ptr<PySqrtAlpha>>(
      m, "SqrtAlpha", "Fortran routine sqrt_alpha return value")
      .def_readonly("y", &PySqrtAlpha::y)
      .def("__len__", [](const PySqrtAlpha&) { return 1; })
      .def("__getitem__", [](const PySqrtAlpha& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.y);
        throw py::index_error();
      });
  m.def(
      "sqrt_one",
      &python_sqrt_one,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::arg("ds1"),
      R"""(Parameters
  ----------
  x : 
  nd : 
  ds1 : 
  )""");
  py::class_<PySqrtOne, std::unique_ptr<PySqrtOne>>(
      m, "SqrtOne", "Fortran routine sqrt_one return value")
      .def_readonly("ds1", &PySqrtOne::ds1)
      .def("__len__", [](const PySqrtOne&) { return 1; })
      .def("__getitem__", [](const PySqrtOne& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.ds1);
        throw py::index_error();
      });
  m.def(
      "str_count",
      &python_str_count,
      py::arg("str"),
      py::arg("match"),
      py::arg("num"),
      R"""(Parameters
  ----------
  str : 
  match : 
  num : 
  )""");
  py::class_<PyStrCount, std::unique_ptr<PyStrCount>>(
      m, "StrCount", "Fortran routine str_count return value")
      .def_readonly("str", &PyStrCount::str)
      .def_readonly("match", &PyStrCount::match)
      .def_readonly("num", &PyStrCount::num)
      .def("__len__", [](const PyStrCount&) { return 3; })
      .def("__getitem__", [](const PyStrCount& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.str);
        if (i == 1)
          return py::cast(s.match);
        if (i == 2)
          return py::cast(s.num);
        throw py::index_error();
      });
  m.def(
      "str_downcase",
      &SimUtils::str_downcase,
      py::arg("src"),
      R"""(Parameters
  ----------
  dst : 
  src : 
  )""");
  m.def(
      "str_first_in_set",
      &python_str_first_in_set,
      py::arg("line"),
      py::arg("set"),
      py::arg("ignore_clauses") = py::none(),
      py::arg("ix_match"),
      R"""(Parameters
  ----------
  line : 
  set : 
  ignore_clauses : 
  ix_match : 
  )""");
  py::class_<PyStrFirstInSet, std::unique_ptr<PyStrFirstInSet>>(
      m, "StrFirstInSet", "Fortran routine str_first_in_set return value")
      .def_readonly("line", &PyStrFirstInSet::line)
      .def_readonly("set", &PyStrFirstInSet::set)
      .def_readonly("ignore_clauses", &PyStrFirstInSet::ignore_clauses)
      .def_readonly("ix_match", &PyStrFirstInSet::ix_match)
      .def("__len__", [](const PyStrFirstInSet&) { return 4; })
      .def("__getitem__", [](const PyStrFirstInSet& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.line);
        if (i == 1)
          return py::cast(s.set);
        if (i == 2)
          return py::cast(s.ignore_clauses);
        if (i == 3)
          return py::cast(s.ix_match);
        throw py::index_error();
      });
  m.def(
      "str_first_not_in_set",
      &python_str_first_not_in_set,
      py::arg("line"),
      py::arg("set"),
      py::arg("ix_match"),
      R"""(Parameters
  ----------
  line : 
  set : 
  ix_match : 
  )""");
  py::class_<PyStrFirstNotInSet, std::unique_ptr<PyStrFirstNotInSet>>(
      m,
      "StrFirstNotInSet",
      "Fortran routine str_first_not_in_set return value")
      .def_readonly("line", &PyStrFirstNotInSet::line)
      .def_readonly("set", &PyStrFirstNotInSet::set)
      .def_readonly("ix_match", &PyStrFirstNotInSet::ix_match)
      .def("__len__", [](const PyStrFirstNotInSet&) { return 3; })
      .def("__getitem__", [](const PyStrFirstNotInSet& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.line);
        if (i == 1)
          return py::cast(s.set);
        if (i == 2)
          return py::cast(s.ix_match);
        throw py::index_error();
      });
  m.def(
      "str_last_in_set",
      &python_str_last_in_set,
      py::arg("line"),
      py::arg("set"),
      py::arg("ix_match"),
      R"""(Parameters
  ----------
  line : 
  set : 
  ix_match : 
  )""");
  py::class_<PyStrLastInSet, std::unique_ptr<PyStrLastInSet>>(
      m, "StrLastInSet", "Fortran routine str_last_in_set return value")
      .def_readonly("line", &PyStrLastInSet::line)
      .def_readonly("set", &PyStrLastInSet::set)
      .def_readonly("ix_match", &PyStrLastInSet::ix_match)
      .def("__len__", [](const PyStrLastInSet&) { return 3; })
      .def("__getitem__", [](const PyStrLastInSet& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.line);
        if (i == 1)
          return py::cast(s.set);
        if (i == 2)
          return py::cast(s.ix_match);
        throw py::index_error();
      });
  m.def(
      "str_last_not_in_set",
      &python_str_last_not_in_set,
      py::arg("line"),
      py::arg("set"),
      py::arg("ix_match"),
      R"""(Parameters
  ----------
  line : 
  set : 
  ix_match : 
  )""");
  py::class_<PyStrLastNotInSet, std::unique_ptr<PyStrLastNotInSet>>(
      m, "StrLastNotInSet", "Fortran routine str_last_not_in_set return value")
      .def_readonly("line", &PyStrLastNotInSet::line)
      .def_readonly("set", &PyStrLastNotInSet::set)
      .def_readonly("ix_match", &PyStrLastNotInSet::ix_match)
      .def("__len__", [](const PyStrLastNotInSet&) { return 3; })
      .def("__getitem__", [](const PyStrLastNotInSet& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.line);
        if (i == 1)
          return py::cast(s.set);
        if (i == 2)
          return py::cast(s.ix_match);
        throw py::index_error();
      });
  m.def(
      "str_match_wild",
      &python_str_match_wild,
      py::arg("str"),
      py::arg("pat"),
      py::arg("a_match"),
      R"""(Parameters
  ----------
  str : 
  pat : 
  a_match : 
  )""");
  py::class_<PyStrMatchWild, std::unique_ptr<PyStrMatchWild>>(
      m, "StrMatchWild", "Fortran routine str_match_wild return value")
      .def_readonly("str", &PyStrMatchWild::str)
      .def_readonly("pat", &PyStrMatchWild::pat)
      .def_readonly("a_match", &PyStrMatchWild::a_match)
      .def("__len__", [](const PyStrMatchWild&) { return 3; })
      .def("__getitem__", [](const PyStrMatchWild& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.str);
        if (i == 1)
          return py::cast(s.pat);
        if (i == 2)
          return py::cast(s.a_match);
        throw py::index_error();
      });
  m.def(
      "str_substitute",
      &python_str_substitute,
      py::arg("string"),
      py::arg("str_match") = py::none(),
      py::arg("str_replace") = py::none(),
      py::arg("do_trim") = py::none(),
      py::arg("ignore_escaped") = py::none(),
      R"""(Parameters
  ----------
  string : 
  str_match : 
  str_replace : 
  do_trim : 
  ignore_escaped : 
  )""");
  py::class_<PyStrSubstitute, std::unique_ptr<PyStrSubstitute>>(
      m, "StrSubstitute", "Fortran routine str_substitute return value")
      .def_readonly("string", &PyStrSubstitute::string)
      .def_readonly("str_match", &PyStrSubstitute::str_match)
      .def_readonly("str_replace", &PyStrSubstitute::str_replace)
      .def_readonly("do_trim", &PyStrSubstitute::do_trim)
      .def_readonly("ignore_escaped", &PyStrSubstitute::ignore_escaped)
      .def("__len__", [](const PyStrSubstitute&) { return 5; })
      .def("__getitem__", [](const PyStrSubstitute& s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.string);
        if (i == 1)
          return py::cast(s.str_match);
        if (i == 2)
          return py::cast(s.str_replace);
        if (i == 3)
          return py::cast(s.do_trim);
        if (i == 4)
          return py::cast(s.ignore_escaped);
        throw py::index_error();
      });
  m.def(
      "str_upcase",
      &SimUtils::str_upcase,
      py::arg("src"),
      R"""(Parameters
  ----------
  dst : 
  src : 
  )""");
  m.def(
      "string_to_int",
      &python_string_to_int,
      py::arg("line"),
      py::arg("default_"),
      py::arg("err_flag"),
      py::arg("err_print_flag") = py::none(),
      py::arg("value"),
      R"""(Parameters
  ----------
  line : 
  default : 
  err_flag : 
  err_print_flag : 
  value : 
  )""");
  py::class_<PyStringToInt, std::unique_ptr<PyStringToInt>>(
      m, "StringToInt", "Fortran routine string_to_int return value")
      .def_readonly("line", &PyStringToInt::line)
      .def_readonly("default_", &PyStringToInt::default_)
      .def_readonly("err_flag", &PyStringToInt::err_flag)
      .def_readonly("err_print_flag", &PyStringToInt::err_print_flag)
      .def_readonly("value", &PyStringToInt::value)
      .def("__len__", [](const PyStringToInt&) { return 5; })
      .def("__getitem__", [](const PyStringToInt& s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.line);
        if (i == 1)
          return py::cast(s.default_);
        if (i == 2)
          return py::cast(s.err_flag);
        if (i == 3)
          return py::cast(s.err_print_flag);
        if (i == 4)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "string_to_real",
      &python_string_to_real,
      py::arg("line"),
      py::arg("default_"),
      py::arg("err_flag"),
      py::arg("err_print_flag") = py::none(),
      py::arg("value"),
      R"""(Parameters
  ----------
  line : 
  default : 
  err_flag : 
  err_print_flag : 
  value : 
  )""");
  py::class_<PyStringToReal, std::unique_ptr<PyStringToReal>>(
      m, "StringToReal", "Fortran routine string_to_real return value")
      .def_readonly("line", &PyStringToReal::line)
      .def_readonly("default_", &PyStringToReal::default_)
      .def_readonly("err_flag", &PyStringToReal::err_flag)
      .def_readonly("err_print_flag", &PyStringToReal::err_print_flag)
      .def_readonly("value", &PyStringToReal::value)
      .def("__len__", [](const PyStringToReal&) { return 5; })
      .def("__getitem__", [](const PyStringToReal& s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.line);
        if (i == 1)
          return py::cast(s.default_);
        if (i == 2)
          return py::cast(s.err_flag);
        if (i == 3)
          return py::cast(s.err_print_flag);
        if (i == 4)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "string_trim",
      &python_string_trim,
      py::arg("in_string"),
      py::arg("out_string"),
      py::arg("word_len"),
      R"""(Parameters
  ----------
  in_string : 
  out_string : 
  word_len : 
  )""");
  py::class_<PyStringTrim, std::unique_ptr<PyStringTrim>>(
      m, "StringTrim", "Fortran routine string_trim return value")
      .def_readonly("in_string", &PyStringTrim::in_string)
      .def_readonly("out_string", &PyStringTrim::out_string)
      .def_readonly("word_len", &PyStringTrim::word_len)
      .def("__len__", [](const PyStringTrim&) { return 3; })
      .def("__getitem__", [](const PyStringTrim& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.in_string);
        if (i == 1)
          return py::cast(s.out_string);
        if (i == 2)
          return py::cast(s.word_len);
        throw py::index_error();
      });
  m.def(
      "string_trim2",
      &python_string_trim2,
      py::arg("in_str"),
      py::arg("delimitors"),
      py::arg("out_str"),
      py::arg("ix_word"),
      py::arg("delim"),
      py::arg("ix_next"),
      R"""(Parameters
  ----------
  in_str : 
  delimitors : 
  out_str : 
  ix_word : 
  delim : 
  ix_next : 
  )""");
  py::class_<PyStringTrim2, std::unique_ptr<PyStringTrim2>>(
      m, "StringTrim2", "Fortran routine string_trim2 return value")
      .def_readonly("in_str", &PyStringTrim2::in_str)
      .def_readonly("delimitors", &PyStringTrim2::delimitors)
      .def_readonly("out_str", &PyStringTrim2::out_str)
      .def_readonly("ix_word", &PyStringTrim2::ix_word)
      .def_readonly("delim", &PyStringTrim2::delim)
      .def_readonly("ix_next", &PyStringTrim2::ix_next)
      .def("__len__", [](const PyStringTrim2&) { return 6; })
      .def("__getitem__", [](const PyStringTrim2& s, int i) -> py::object {
        if (i < 0)
          i += 6;
        if (i == 0)
          return py::cast(s.in_str);
        if (i == 1)
          return py::cast(s.delimitors);
        if (i == 2)
          return py::cast(s.out_str);
        if (i == 3)
          return py::cast(s.ix_word);
        if (i == 4)
          return py::cast(s.delim);
        if (i == 5)
          return py::cast(s.ix_next);
        throw py::index_error();
      });
  m.def(
      "super_bicubic_coef",
      &SimUtils::super_bicubic_coef,
      py::arg("y"),
      py::arg("y1"),
      py::arg("y2"),
      py::arg("y12"),
      py::arg("d1"),
      py::arg("d2"),
      R"""(Subroutine super_bicubic_coef(y, y1, y2, y12, d1, d2, c)

  Routine to compute coefficients for bicubic interpolation.
  This is from NR bcucof.

  Parameters
  ----------
  y : float
      Function values at grid points.
  y1 : float
      dy/dx1 derivatives.
  y2 : float
      dy/dx2 derivatives.
  y12 : float
      d2y/dx1*dx2 second derivatives.
  d1 : float
      Grid width in 1-direction.
  d2 : float
      Grid width in 2-direction.

  Returns
  -------
  c : float
      Coefficients.
  )""");
  m.def(
      "super_bicubic_interpolation",
      &SimUtils::super_bicubic_interpolation,
      py::arg("y"),
      py::arg("y1"),
      py::arg("y2"),
      py::arg("y12"),
      py::arg("x1l"),
      py::arg("x1u"),
      py::arg("x2l"),
      py::arg("x2u"),
      py::arg("x1"),
      py::arg("x2"),
      R"""(Subroutine super_bicubic_interpolation(y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, ansy, ansy1, ansy2)

  Routine to do bicubic interpolation.
  This is from NR bcuint.

  Note! The four grid points are arrayed in counter-clockwise order beginning from the lower left.
  So, for example, y = [y_ll, y_lu, y_uu, y_ul] where "l" = lower, "u" = upper index.

  Parameters
  ----------
  y : float
      Function values at grid points.
  y1 : float
      dy/dx1 derivatives.
  y2 : float
      dy/dx2 derivatives.
  y12 : float
      d2y/dx1*dx2 second derivatives.
  x1l : float
      1-direction coordinate at lower points.
  x1u : float
      1-direction coordinate at upper points
  x2l : float
      2-direction coordinate at lower points.
  x2u : float
      2-direction coordinate at upper points
  x1 : float
      1-direction coordinate at point to evaluate.
  x2 : float
      2-direction coordinate at point to evaluate.

  Returns
  -------
  ansy : float
      Interpolation value.
  ansy1 : float
      1-direction derivative at interpolation point.
  ansy2 : float
      2-direction derivative at interpolation point.
  )""");
  py::class_<
      SimUtils::SuperBicubicInterpolation,
      std::unique_ptr<SimUtils::SuperBicubicInterpolation>>(
      m,
      "SuperBicubicInterpolation",
      "Fortran routine super_bicubic_interpolation return value")
      .def_readonly("ansy", &SimUtils::SuperBicubicInterpolation::ansy)
      .def_readonly("ansy1", &SimUtils::SuperBicubicInterpolation::ansy1)
      .def_readonly("ansy2", &SimUtils::SuperBicubicInterpolation::ansy2)
      .def(
          "__len__",
          [](const SimUtils::SuperBicubicInterpolation&) { return 3; })
      .def(
          "__getitem__",
          [](const SimUtils::SuperBicubicInterpolation& s,
             int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.ansy);
            if (i == 1)
              return py::cast(s.ansy1);
            if (i == 2)
              return py::cast(s.ansy2);
            throw py::index_error();
          });
  m.def(
      "super_polint",
      &SimUtils::super_polint,
      py::arg("xa"),
      py::arg("ya"),
      py::arg("x"),
      R"""(Function super_polint (xa, ya, x, y, dy)

  This is essentially polint from Numerical Recipes.

  Parameters
  ----------
  xa : float
  x : float

  Returns
  -------
  y : float
  dy : float
  )""");
  py::class_<SimUtils::SuperPolint, std::unique_ptr<SimUtils::SuperPolint>>(
      m, "SuperPolint", "Fortran routine super_polint return value")
      .def_readonly("y", &SimUtils::SuperPolint::y)
      .def_readonly("dy", &SimUtils::SuperPolint::dy)
      .def("__len__", [](const SimUtils::SuperPolint&) { return 2; })
      .def(
          "__getitem__",
          [](const SimUtils::SuperPolint& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.y);
            if (i == 1)
              return py::cast(s.dy);
            throw py::index_error();
          });
  m.def(
      "super_poly",
      &SimUtils::super_poly,
      py::arg("x"),
      py::arg("coeffs"),
      R"""(Function super_poly (x, coef) result (value)

  Routine to compute Sum: coef(i)*x^i

  Parameters
  ----------
  x : float
      Variable.
  coef : float
      Coefficients.

  Returns
  -------
  value : float
      Polynomial value.
  )""");
  m.def(
      "super_sobseq",
      &SimUtils::super_sobseq,
      py::arg("ran_state") = py::none(),
      R"""(Subroutine super_sobseq (x, ran_state)

  Routine patterened after sobseq in Numerical Recipes.
  Difference is that this version has an argument for the internal state.

  Parameters
  ----------
  ran_state : RandomStateStruct, optional
      Generator state. See the ran_seed_put documentation for more details.

  Returns
  -------
  x : float
      Random vector.
  )""");
  m.def(
      "super_sort",
      &SimUtils::super_sort,
      py::arg("arr"),
      R"""(Subroutine super_sort(arr)

  Routine to sort an integer array in place.
  This is the NR routine sort modified to sort integers.

  Parameters
  ----------
  arr : int
      Array of integers.
      This parameter is an input/output and is modified in-place. As an output: Sorted array.
  )""");
  m.def(
      "system_command",
      &python_system_command,
      py::arg("line"),
      py::arg("err_flag") = py::none(),
      R"""(Parameters
  ----------
  line : 
  err_flag : 
  )""");
  py::class_<PySystemCommand, std::unique_ptr<PySystemCommand>>(
      m, "SystemCommand", "Fortran routine system_command return value")
      .def_readonly("line", &PySystemCommand::line)
      .def_readonly("err_flag", &PySystemCommand::err_flag)
      .def("__len__", [](const PySystemCommand&) { return 2; })
      .def("__getitem__", [](const PySystemCommand& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.line);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
}
