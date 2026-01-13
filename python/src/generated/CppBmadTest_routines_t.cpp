#include "pybmad/generated/CppBmadTest_routines_t.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyTestCharacterScalar python_test_character_scalar(
    std::string val_in,
    std::string val_inout,
    std::optional<std::string> val_in_opt = std::nullopt,
    std::optional<std::string> val_inout_opt = std::nullopt) {
  auto _result = CppBmadTest::test_character_scalar(
      val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestCharacterScalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestComplexScalar python_test_complex_scalar(
    std::complex<double> val_in,
    std::complex<double> val_inout,
    std::optional<std::complex<double>> val_in_opt = std::nullopt,
    std::optional<std::complex<double>> val_inout_opt = std::nullopt) {
  auto _result = CppBmadTest::test_complex_scalar(
      val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestComplexScalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestInteger8Scalar python_test_integer8_scalar(
    int64_t val_in,
    int64_t val_inout,
    std::optional<int64_t> val_in_opt = std::nullopt,
    std::optional<int64_t> val_inout_opt = std::nullopt) {
  auto _result = CppBmadTest::test_integer8_scalar(
      val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestInteger8Scalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestIntegerScalar python_test_integer_scalar(
    int val_in,
    int val_inout,
    std::optional<int> val_in_opt = std::nullopt,
    std::optional<int> val_inout_opt = std::nullopt) {
  auto _result = CppBmadTest::test_integer_scalar(
      val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestIntegerScalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestLogicalScalar python_test_logical_scalar(
    bool val_in,
    bool val_inout,
    std::optional<bool> val_in_opt = std::nullopt,
    std::optional<bool> val_inout_opt = std::nullopt) {
  auto _result = CppBmadTest::test_logical_scalar(
      val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestLogicalScalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestReal16Scalar python_test_real16_scalar(
    long double val_in,
    long double val_inout,
    std::optional<long double> val_in_opt = std::nullopt,
    std::optional<long double> val_inout_opt = std::nullopt) {
  auto _result = CppBmadTest::test_real16_scalar(
      val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestReal16Scalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestRealScalar python_test_real_scalar(
    double val_in,
    double val_inout,
    std::optional<double> val_in_opt = std::nullopt,
    std::optional<double> val_inout_opt = std::nullopt) {
  auto _result = CppBmadTest::test_real_scalar(
      val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestRealScalar{_result, val_inout, val_inout_opt}};
  return py_result;
}

void init_CppBmadTest_routines_t(py::module& m) {
  py::class_<
      CppBmadTest::TestBunchStructArray,
      std::unique_ptr<CppBmadTest::TestBunchStructArray>>(
      m, "TestBunchStructArray", "test_bunch_struct_array return type")
      .def_readonly("arr_out", &CppBmadTest::TestBunchStructArray::arr_out)
      .def_readonly(
          "opt_status", &CppBmadTest::TestBunchStructArray::opt_status)
      .def(
          "__len__", [](const CppBmadTest::TestBunchStructArray&) { return 2; })
      .def(
          "__getitem__",
          [](const CppBmadTest::TestBunchStructArray& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.arr_out);
            if (i == 1)
              return py::cast(s.opt_status);
            throw py::index_error();
          });
  m.def(
      "test_bunch_struct_array",
      &CppBmadTest::test_bunch_struct_array,
      py::arg("arr_in"),
      py::arg("arr_inout"),
      py::arg("arr_in_opt") = py::none(),
      py::arg("arr_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  arr_in : 
  arr_inout : 
  arr_out : 
  opt_status : 
  arr_in_opt : 
  arr_inout_opt : 
  )""");
  py::class_<
      CppBmadTest::TestBunchStructScalar,
      std::unique_ptr<CppBmadTest::TestBunchStructScalar>>(
      m, "TestBunchStructScalar", "test_bunch_struct_scalar return type")
      .def_readonly("val_out", &CppBmadTest::TestBunchStructScalar::val_out)
      .def_readonly(
          "opt_status", &CppBmadTest::TestBunchStructScalar::opt_status)
      .def(
          "__len__",
          [](const CppBmadTest::TestBunchStructScalar&) { return 2; })
      .def(
          "__getitem__",
          [](const CppBmadTest::TestBunchStructScalar& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.val_out);
            if (i == 1)
              return py::cast(s.opt_status);
            throw py::index_error();
          });
  m.def(
      "test_bunch_struct_scalar",
      &CppBmadTest::test_bunch_struct_scalar,
      py::arg("val_in"),
      py::arg("val_inout"),
      py::arg("val_in_opt") = py::none(),
      py::arg("val_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  val_in : 
  val_inout : 
  val_out : 
  opt_status : 
  val_in_opt : 
  val_inout_opt : 
  )""");
  py::class_<PyTestCharacterScalar, std::unique_ptr<PyTestCharacterScalar>>(
      m, "TestCharacterScalar", "test_character_scalar return type")
      .def_readonly("val_out", &PyTestCharacterScalar::val_out)
      .def_readonly("opt_status", &PyTestCharacterScalar::opt_status)
      .def_readonly("val_inout", &PyTestCharacterScalar::val_inout)
      .def_readonly("val_inout_opt", &PyTestCharacterScalar::val_inout_opt)
      .def("__len__", [](const PyTestCharacterScalar&) { return 4; })
      .def(
          "__getitem__",
          [](const PyTestCharacterScalar& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.val_out);
            if (i == 1)
              return py::cast(s.opt_status);
            if (i == 2)
              return py::cast(s.val_inout);
            if (i == 3)
              return py::cast(s.val_inout_opt);
            throw py::index_error();
          });
  m.def(
      "test_character_scalar",
      &python_test_character_scalar,
      py::arg("val_in"),
      py::arg("val_inout"),
      py::arg("val_in_opt") = py::none(),
      py::arg("val_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  val_in : 
  val_inout : 
  val_out : 
  opt_status : 
  val_in_opt : 
  val_inout_opt : 
  )""");
  py::class_<
      CppBmadTest::TestComplexArray,
      std::unique_ptr<CppBmadTest::TestComplexArray>>(
      m, "TestComplexArray", "test_complex_array return type")
      .def_readonly("arr_out", &CppBmadTest::TestComplexArray::arr_out)
      .def_readonly("opt_status", &CppBmadTest::TestComplexArray::opt_status)
      .def("__len__", [](const CppBmadTest::TestComplexArray&) { return 2; })
      .def(
          "__getitem__",
          [](const CppBmadTest::TestComplexArray& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.arr_out);
            if (i == 1)
              return py::cast(s.opt_status);
            throw py::index_error();
          });
  m.def(
      "test_complex_array",
      &CppBmadTest::test_complex_array,
      py::arg("arr_in"),
      py::arg("arr_inout"),
      py::arg("arr_in_opt") = py::none(),
      py::arg("arr_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  arr_in : 
  arr_inout : 
  arr_out : 
  opt_status : 
  arr_in_opt : 
  arr_inout_opt : 
  )""");
  py::class_<PyTestComplexScalar, std::unique_ptr<PyTestComplexScalar>>(
      m, "TestComplexScalar", "test_complex_scalar return type")
      .def_readonly("val_out", &PyTestComplexScalar::val_out)
      .def_readonly("opt_status", &PyTestComplexScalar::opt_status)
      .def_readonly("val_inout", &PyTestComplexScalar::val_inout)
      .def_readonly("val_inout_opt", &PyTestComplexScalar::val_inout_opt)
      .def("__len__", [](const PyTestComplexScalar&) { return 4; })
      .def(
          "__getitem__", [](const PyTestComplexScalar& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.val_out);
            if (i == 1)
              return py::cast(s.opt_status);
            if (i == 2)
              return py::cast(s.val_inout);
            if (i == 3)
              return py::cast(s.val_inout_opt);
            throw py::index_error();
          });
  m.def(
      "test_complex_scalar",
      &python_test_complex_scalar,
      py::arg("val_in"),
      py::arg("val_inout"),
      py::arg("val_in_opt") = py::none(),
      py::arg("val_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  val_in : 
  val_inout : 
  val_out : 
  opt_status : 
  val_in_opt : 
  val_inout_opt : 
  )""");
  py::class_<
      CppBmadTest::TestInteger8Array,
      std::unique_ptr<CppBmadTest::TestInteger8Array>>(
      m, "TestInteger8Array", "test_integer8_array return type")
      .def_readonly("arr_out", &CppBmadTest::TestInteger8Array::arr_out)
      .def_readonly("opt_status", &CppBmadTest::TestInteger8Array::opt_status)
      .def("__len__", [](const CppBmadTest::TestInteger8Array&) { return 2; })
      .def(
          "__getitem__",
          [](const CppBmadTest::TestInteger8Array& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.arr_out);
            if (i == 1)
              return py::cast(s.opt_status);
            throw py::index_error();
          });
  m.def(
      "test_integer8_array",
      &CppBmadTest::test_integer8_array,
      py::arg("arr_in"),
      py::arg("arr_inout"),
      py::arg("arr_in_opt") = py::none(),
      py::arg("arr_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  arr_in : 
  arr_inout : 
  arr_out : 
  opt_status : 
  arr_in_opt : 
  arr_inout_opt : 
  )""");
  py::class_<PyTestInteger8Scalar, std::unique_ptr<PyTestInteger8Scalar>>(
      m, "TestInteger8Scalar", "test_integer8_scalar return type")
      .def_readonly("val_out", &PyTestInteger8Scalar::val_out)
      .def_readonly("opt_status", &PyTestInteger8Scalar::opt_status)
      .def_readonly("val_inout", &PyTestInteger8Scalar::val_inout)
      .def_readonly("val_inout_opt", &PyTestInteger8Scalar::val_inout_opt)
      .def("__len__", [](const PyTestInteger8Scalar&) { return 4; })
      .def(
          "__getitem__",
          [](const PyTestInteger8Scalar& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.val_out);
            if (i == 1)
              return py::cast(s.opt_status);
            if (i == 2)
              return py::cast(s.val_inout);
            if (i == 3)
              return py::cast(s.val_inout_opt);
            throw py::index_error();
          });
  m.def(
      "test_integer8_scalar",
      &python_test_integer8_scalar,
      py::arg("val_in"),
      py::arg("val_inout"),
      py::arg("val_in_opt") = py::none(),
      py::arg("val_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  val_in : 
  val_inout : 
  val_out : 
  opt_status : 
  val_in_opt : 
  val_inout_opt : 
  )""");
  py::class_<
      CppBmadTest::TestIntegerArray,
      std::unique_ptr<CppBmadTest::TestIntegerArray>>(
      m, "TestIntegerArray", "test_integer_array return type")
      .def_readonly("arr_out", &CppBmadTest::TestIntegerArray::arr_out)
      .def_readonly("opt_status", &CppBmadTest::TestIntegerArray::opt_status)
      .def("__len__", [](const CppBmadTest::TestIntegerArray&) { return 2; })
      .def(
          "__getitem__",
          [](const CppBmadTest::TestIntegerArray& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.arr_out);
            if (i == 1)
              return py::cast(s.opt_status);
            throw py::index_error();
          });
  m.def(
      "test_integer_array",
      &CppBmadTest::test_integer_array,
      py::arg("arr_in"),
      py::arg("arr_inout"),
      py::arg("arr_in_opt") = py::none(),
      py::arg("arr_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  arr_in : 
  arr_inout : 
  arr_out : 
  opt_status : 
  arr_in_opt : 
  arr_inout_opt : 
  )""");
  py::class_<PyTestIntegerScalar, std::unique_ptr<PyTestIntegerScalar>>(
      m, "TestIntegerScalar", "test_integer_scalar return type")
      .def_readonly("val_out", &PyTestIntegerScalar::val_out)
      .def_readonly("opt_status", &PyTestIntegerScalar::opt_status)
      .def_readonly("val_inout", &PyTestIntegerScalar::val_inout)
      .def_readonly("val_inout_opt", &PyTestIntegerScalar::val_inout_opt)
      .def("__len__", [](const PyTestIntegerScalar&) { return 4; })
      .def(
          "__getitem__", [](const PyTestIntegerScalar& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.val_out);
            if (i == 1)
              return py::cast(s.opt_status);
            if (i == 2)
              return py::cast(s.val_inout);
            if (i == 3)
              return py::cast(s.val_inout_opt);
            throw py::index_error();
          });
  m.def(
      "test_integer_scalar",
      &python_test_integer_scalar,
      py::arg("val_in"),
      py::arg("val_inout"),
      py::arg("val_in_opt") = py::none(),
      py::arg("val_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  val_in : 
  val_inout : 
  val_out : 
  opt_status : 
  val_in_opt : 
  val_inout_opt : 
  )""");
  py::class_<
      CppBmadTest::TestLogicalArray,
      std::unique_ptr<CppBmadTest::TestLogicalArray>>(
      m, "TestLogicalArray", "test_logical_array return type")
      .def_readonly("arr_out", &CppBmadTest::TestLogicalArray::arr_out)
      .def_readonly("opt_status", &CppBmadTest::TestLogicalArray::opt_status)
      .def("__len__", [](const CppBmadTest::TestLogicalArray&) { return 2; })
      .def(
          "__getitem__",
          [](const CppBmadTest::TestLogicalArray& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.arr_out);
            if (i == 1)
              return py::cast(s.opt_status);
            throw py::index_error();
          });
  m.def(
      "test_logical_array",
      &CppBmadTest::test_logical_array,
      py::arg("arr_in"),
      py::arg("arr_inout"),
      py::arg("arr_in_opt") = py::none(),
      py::arg("arr_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  arr_in : 
  arr_inout : 
  arr_out : 
  opt_status : 
  arr_in_opt : 
  arr_inout_opt : 
  )""");
  py::class_<PyTestLogicalScalar, std::unique_ptr<PyTestLogicalScalar>>(
      m, "TestLogicalScalar", "test_logical_scalar return type")
      .def_readonly("val_out", &PyTestLogicalScalar::val_out)
      .def_readonly("opt_status", &PyTestLogicalScalar::opt_status)
      .def_readonly("val_inout", &PyTestLogicalScalar::val_inout)
      .def_readonly("val_inout_opt", &PyTestLogicalScalar::val_inout_opt)
      .def("__len__", [](const PyTestLogicalScalar&) { return 4; })
      .def(
          "__getitem__", [](const PyTestLogicalScalar& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.val_out);
            if (i == 1)
              return py::cast(s.opt_status);
            if (i == 2)
              return py::cast(s.val_inout);
            if (i == 3)
              return py::cast(s.val_inout_opt);
            throw py::index_error();
          });
  m.def(
      "test_logical_scalar",
      &python_test_logical_scalar,
      py::arg("val_in"),
      py::arg("val_inout"),
      py::arg("val_in_opt") = py::none(),
      py::arg("val_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  val_in : 
  val_inout : 
  val_out : 
  opt_status : 
  val_in_opt : 
  val_inout_opt : 
  )""");
  py::class_<
      CppBmadTest::TestReal16Array,
      std::unique_ptr<CppBmadTest::TestReal16Array>>(
      m, "TestReal16Array", "test_real16_array return type")
      .def_readonly("arr_out", &CppBmadTest::TestReal16Array::arr_out)
      .def_readonly("opt_status", &CppBmadTest::TestReal16Array::opt_status)
      .def("__len__", [](const CppBmadTest::TestReal16Array&) { return 2; })
      .def(
          "__getitem__",
          [](const CppBmadTest::TestReal16Array& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.arr_out);
            if (i == 1)
              return py::cast(s.opt_status);
            throw py::index_error();
          });
  m.def(
      "test_real16_array",
      &CppBmadTest::test_real16_array,
      py::arg("arr_in"),
      py::arg("arr_inout"),
      py::arg("arr_in_opt") = py::none(),
      py::arg("arr_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  arr_in : 
  arr_inout : 
  arr_out : 
  opt_status : 
  arr_in_opt : 
  arr_inout_opt : 
  )""");
  py::class_<PyTestReal16Scalar, std::unique_ptr<PyTestReal16Scalar>>(
      m, "TestReal16Scalar", "test_real16_scalar return type")
      .def_readonly("val_out", &PyTestReal16Scalar::val_out)
      .def_readonly("opt_status", &PyTestReal16Scalar::opt_status)
      .def_readonly("val_inout", &PyTestReal16Scalar::val_inout)
      .def_readonly("val_inout_opt", &PyTestReal16Scalar::val_inout_opt)
      .def("__len__", [](const PyTestReal16Scalar&) { return 4; })
      .def("__getitem__", [](const PyTestReal16Scalar& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.val_out);
        if (i == 1)
          return py::cast(s.opt_status);
        if (i == 2)
          return py::cast(s.val_inout);
        if (i == 3)
          return py::cast(s.val_inout_opt);
        throw py::index_error();
      });
  m.def(
      "test_real16_scalar",
      &python_test_real16_scalar,
      py::arg("val_in"),
      py::arg("val_inout"),
      py::arg("val_in_opt") = py::none(),
      py::arg("val_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  val_in : 
  val_inout : 
  val_out : 
  opt_status : 
  val_in_opt : 
  val_inout_opt : 
  )""");
  py::class_<
      CppBmadTest::TestRealArray,
      std::unique_ptr<CppBmadTest::TestRealArray>>(
      m, "TestRealArray", "test_real_array return type")
      .def_readonly("arr_out", &CppBmadTest::TestRealArray::arr_out)
      .def_readonly("opt_status", &CppBmadTest::TestRealArray::opt_status)
      .def("__len__", [](const CppBmadTest::TestRealArray&) { return 2; })
      .def(
          "__getitem__",
          [](const CppBmadTest::TestRealArray& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.arr_out);
            if (i == 1)
              return py::cast(s.opt_status);
            throw py::index_error();
          });
  m.def(
      "test_real_array",
      &CppBmadTest::test_real_array,
      py::arg("arr_in"),
      py::arg("arr_inout"),
      py::arg("arr_in_opt") = py::none(),
      py::arg("arr_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  arr_in : 
  arr_inout : 
  arr_out : 
  opt_status : 
  arr_in_opt : 
  arr_inout_opt : 
  )""");
  py::class_<PyTestRealScalar, std::unique_ptr<PyTestRealScalar>>(
      m, "TestRealScalar", "test_real_scalar return type")
      .def_readonly("val_out", &PyTestRealScalar::val_out)
      .def_readonly("opt_status", &PyTestRealScalar::opt_status)
      .def_readonly("val_inout", &PyTestRealScalar::val_inout)
      .def_readonly("val_inout_opt", &PyTestRealScalar::val_inout_opt)
      .def("__len__", [](const PyTestRealScalar&) { return 4; })
      .def("__getitem__", [](const PyTestRealScalar& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.val_out);
        if (i == 1)
          return py::cast(s.opt_status);
        if (i == 2)
          return py::cast(s.val_inout);
        if (i == 3)
          return py::cast(s.val_inout_opt);
        throw py::index_error();
      });
  m.def(
      "test_real_scalar",
      &python_test_real_scalar,
      py::arg("val_in"),
      py::arg("val_inout"),
      py::arg("val_in_opt") = py::none(),
      py::arg("val_inout_opt") = py::none(),
      R"""(Parameters
  ----------
  val_in : 
  val_inout : 
  val_out : 
  opt_status : 
  val_in_opt : 
  val_inout_opt : 
  )""");
}
