#include "pybmad/generated/Bmad_routines_v.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyValidFieldCalc python_valid_field_calc(
    EleProxy& ele,
    int field_calc,
    bool is_valid) {
  Bmad::valid_field_calc(ele, field_calc, is_valid);
  auto py_result{PyValidFieldCalc{is_valid}};
  return py_result;
}
PyValidFringeType python_valid_fringe_type(
    EleProxy& ele,
    int fringe_type,
    bool is_valid) {
  Bmad::valid_fringe_type(ele, fringe_type, is_valid);
  auto py_result{PyValidFringeType{is_valid}};
  return py_result;
}
PyValidMat6CalcMethod python_valid_mat6_calc_method(
    EleProxy& ele,
    int species,
    int mat6_calc_method,
    bool is_valid) {
  Bmad::valid_mat6_calc_method(ele, species, mat6_calc_method, is_valid);
  auto py_result{PyValidMat6CalcMethod{is_valid}};
  return py_result;
}
PyValidSpinTrackingMethod python_valid_spin_tracking_method(
    EleProxy& ele,
    int spin_tracking_method,
    bool is_valid) {
  Bmad::valid_spin_tracking_method(ele, spin_tracking_method, is_valid);
  auto py_result{PyValidSpinTrackingMethod{is_valid}};
  return py_result;
}
PyValidTrackingMethod python_valid_tracking_method(
    EleProxy& ele,
    int species,
    int tracking_method,
    bool is_valid) {
  Bmad::valid_tracking_method(ele, species, tracking_method, is_valid);
  auto py_result{PyValidTrackingMethod{is_valid}};
  return py_result;
}
PyValueOfAttribute python_value_of_attribute(
    EleProxy& ele,
    std::string attrib_name,
    std::optional<bool> err_print_flag,
    std::optional<double> err_value,
    double value) {
  auto _result = Bmad::value_of_attribute(
      ele, attrib_name, err_print_flag, err_value, value);
  auto py_result{PyValueOfAttribute{_result, value}};
  return py_result;
}
PyValueToLine python_value_to_line(
    std::string line,
    double value,
    std::string str,
    std::string typ,
    std::optional<bool> ignore_if_zero = std::nullopt,
    std::optional<bool> use_comma = std::nullopt) {
  Bmad::value_to_line(
      line,
      value,
      str,
      typ,
      make_opt_ref(ignore_if_zero),
      make_opt_ref(use_comma));
  auto py_result{
      PyValueToLine{line, value, str, typ, ignore_if_zero, use_comma}};
  return py_result;
}

void init_Bmad_routines_v(py::module& m) {
  py::class_<PyValidFieldCalc, std::unique_ptr<PyValidFieldCalc>>(
      m, "ValidFieldCalc", "valid_field_calc return type")
      .def_readonly("is_valid", &PyValidFieldCalc::is_valid)
      .def("__len__", [](const PyValidFieldCalc&) { return 1; })
      .def("__getitem__", [](const PyValidFieldCalc& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_valid);
        throw py::index_error();
      });
  m.def(
      "valid_field_calc",
      &python_valid_field_calc,
      py::arg("ele"),
      py::arg("field_calc"),
      py::arg("is_valid"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element.
  field_calc : int
      bmad_standard$, etc.
  is_valid : 
  )""");
  py::class_<PyValidFringeType, std::unique_ptr<PyValidFringeType>>(
      m, "ValidFringeType", "valid_fringe_type return type")
      .def_readonly("is_valid", &PyValidFringeType::is_valid)
      .def("__len__", [](const PyValidFringeType&) { return 1; })
      .def("__getitem__", [](const PyValidFringeType& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_valid);
        throw py::index_error();
      });
  m.def(
      "valid_fringe_type",
      &python_valid_fringe_type,
      py::arg("ele"),
      py::arg("fringe_type"),
      py::arg("is_valid"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element.
  fringe_type : int
      bmad_standard$, etc.
  is_valid : 
  )""");
  py::class_<PyValidMat6CalcMethod, std::unique_ptr<PyValidMat6CalcMethod>>(
      m, "ValidMat6CalcMethod", "valid_mat6_calc_method return type")
      .def_readonly("is_valid", &PyValidMat6CalcMethod::is_valid)
      .def("__len__", [](const PyValidMat6CalcMethod&) { return 1; })
      .def(
          "__getitem__",
          [](const PyValidMat6CalcMethod& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_valid);
            throw py::index_error();
          });
  m.def(
      "valid_mat6_calc_method",
      &python_valid_mat6_calc_method,
      py::arg("ele"),
      py::arg("species"),
      py::arg("mat6_calc_method"),
      py::arg("is_valid"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element.
  species : 
      Type of particle being tracked. electron$, etc. or not_set$
  mat6_calc_method : int
      bmad_standard$, etc.
  is_valid : 
  )""");
  py::class_<
      PyValidSpinTrackingMethod,
      std::unique_ptr<PyValidSpinTrackingMethod>>(
      m, "ValidSpinTrackingMethod", "valid_spin_tracking_method return type")
      .def_readonly("is_valid", &PyValidSpinTrackingMethod::is_valid)
      .def("__len__", [](const PyValidSpinTrackingMethod&) { return 1; })
      .def(
          "__getitem__",
          [](const PyValidSpinTrackingMethod& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_valid);
            throw py::index_error();
          });
  m.def(
      "valid_spin_tracking_method",
      &python_valid_spin_tracking_method,
      py::arg("ele"),
      py::arg("spin_tracking_method"),
      py::arg("is_valid"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element.
  spin_tracking_method : int
      bmad_standard$, etc.
  is_valid : 
  )""");
  py::class_<PyValidTrackingMethod, std::unique_ptr<PyValidTrackingMethod>>(
      m, "ValidTrackingMethod", "valid_tracking_method return type")
      .def_readonly("is_valid", &PyValidTrackingMethod::is_valid)
      .def("__len__", [](const PyValidTrackingMethod&) { return 1; })
      .def(
          "__getitem__",
          [](const PyValidTrackingMethod& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_valid);
            throw py::index_error();
          });
  m.def(
      "valid_tracking_method",
      &python_valid_tracking_method,
      py::arg("ele"),
      py::arg("species"),
      py::arg("tracking_method"),
      py::arg("is_valid"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element.
  species : 
      Type of particle being tracked. electron$, etc. or not_set$
  tracking_method : int
      bmad_standard$, etc.
  is_valid : 
  )""");
  py::class_<PyValueOfAttribute, std::unique_ptr<PyValueOfAttribute>>(
      m, "ValueOfAttribute", "value_of_attribute return type")
      .def_readonly("err_flag", &PyValueOfAttribute::err_flag)
      .def_readonly("value", &PyValueOfAttribute::value)
      .def("__len__", [](const PyValueOfAttribute&) { return 2; })
      .def("__getitem__", [](const PyValueOfAttribute& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "value_of_attribute",
      &python_value_of_attribute,
      py::arg("ele"),
      py::arg("attrib_name"),
      py::arg("err_print_flag") = py::none(),
      py::arg("err_value") = py::none(),
      py::arg("value"),
      R"""(Parameters
  ----------
  ele : EleStruct
      After this routine finishes Ptr_attrib will point to a variable within this element.
  attrib_name : unknown
      Name of attribute. Must be uppercase. For example: "HKICK".
  err_flag : bool
      Set True if attribtute not found. False otherwise.
  err_print_flag : bool, optional
      If present and True then print an error message if there is an  error.
  err_value : float, optional
      Value to set value argument if there is an error. Default is 0.
  value : 
  )""");
  py::class_<PyValueToLine, std::unique_ptr<PyValueToLine>>(
      m, "ValueToLine", "value_to_line return type")
      .def_readonly("line", &PyValueToLine::line)
      .def_readonly("value", &PyValueToLine::value)
      .def_readonly("str", &PyValueToLine::str)
      .def_readonly("typ", &PyValueToLine::typ)
      .def_readonly("ignore_if_zero", &PyValueToLine::ignore_if_zero)
      .def_readonly("use_comma", &PyValueToLine::use_comma)
      .def("__len__", [](const PyValueToLine&) { return 6; })
      .def("__getitem__", [](const PyValueToLine& s, int i) -> py::object {
        if (i < 0)
          i += 6;
        if (i == 0)
          return py::cast(s.line);
        if (i == 1)
          return py::cast(s.value);
        if (i == 2)
          return py::cast(s.str);
        if (i == 3)
          return py::cast(s.typ);
        if (i == 4)
          return py::cast(s.ignore_if_zero);
        if (i == 5)
          return py::cast(s.use_comma);
        throw py::index_error();
      });
  m.def(
      "value_to_line",
      &python_value_to_line,
      py::arg("line"),
      py::arg("value"),
      py::arg("str"),
      py::arg("typ"),
      py::arg("ignore_if_zero") = py::none(),
      py::arg("use_comma") = py::none(),
      R"""(Parameters
  ----------
  line : 
  value : 
  str : 
  typ : 
  ignore_if_zero : 
  use_comma : 
  )""");
  m.def(
      "vec_to_polar",
      &Bmad::vec_to_polar,
      py::arg("vec"),
      py::arg("phase") = py::none(),
      py::arg("polar"),
      R"""(Parameters
  ----------
  vec : float
      unitary spin vector
  phase : float, optional
      Phase of the spinor, if not given then set to zero
  polar : 
  )""");
  m.def(
      "vec_to_spinor",
      &Bmad::vec_to_spinor,
      py::arg("vec"),
      py::arg("phase") = py::none(),
      py::arg("spinor"),
      R"""(Parameters
  ----------
  vec : float
      Spin vector in cartesian coordinates
  phase : float
      Phase of the spinor, if not given then set to zero
  spinor : 
  )""");
  m.def(
      "verify_valid_name",
      &Bmad::verify_valid_name,
      py::arg("name"),
      py::arg("ix_name"),
      py::arg("pure_name") = py::none(),
      py::arg("include_wild") = py::none(),
      R"""(Function verify_valid_name (name, ix_name, pure_name, include_wild) result (is_valid)

  Routine to check if a name is well formed. Examples:
    "0>>Q0"                           -- Invalid (will only be valid after lattice expansion).
    "Q1##1"                           -- Invalid (double hash not accepted).
    "Q2A_C.\7#"                       -- Pure name (no "[", "]", "(", ")", "%" characters present).
    "Q3[GRID_FIELD(1)%FIELD_SCALE]"   -- Valid but not a pure name.
    "RFCAVITY::*"                     -- Valid if include_wild = True.

  This subroutine is used by bmad_parser and bmad_parser2.
  This subroutine is not intended for general use.


  Parameters
  ----------
  name : unknown
      Name(1:ix_name) is the string to check.
  ix_name : int
      Number of characters in the name.
  pure_name : bool, optional
      If True, reject names that contain "[", "]", "(", ")", "." characters. Default is False.
  include_wild : bool, optional
      Name can include wild card characters and additionally type prefixes like "QUAD::". Default is False.

  Returns
  -------
  is_valid : bool
      True if name is well formed. False otherwise.
  )""");
}
