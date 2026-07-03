#include "pybmad/generated/CppBmadTest_routines_t.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

PyTestCharacterScalar python_test_character_scalar(
    std::string val_in,
    std::string val_inout,
    std::optional<std::string> val_in_opt = std::nullopt,
    std::optional<std::string> val_inout_opt = std::nullopt
) {
  auto _result = CppBmadTest::test_character_scalar(
      val_in,
      val_inout,
      val_in_opt,
      make_opt_ref(val_inout_opt)
  );
  auto py_result{PyTestCharacterScalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestComplexScalar python_test_complex_scalar(
    std::complex<double> val_in,
    std::complex<double> val_inout,
    std::optional<std::complex<double>> val_in_opt = std::nullopt,
    std::optional<std::complex<double>> val_inout_opt = std::nullopt
) {
  auto _result =
      CppBmadTest::test_complex_scalar(val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestComplexScalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestInteger8Scalar python_test_integer8_scalar(
    int64_t val_in,
    int64_t val_inout,
    std::optional<int64_t> val_in_opt = std::nullopt,
    std::optional<int64_t> val_inout_opt = std::nullopt
) {
  auto _result =
      CppBmadTest::test_integer8_scalar(val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestInteger8Scalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestIntegerScalar python_test_integer_scalar(
    int val_in,
    int val_inout,
    std::optional<int> val_in_opt = std::nullopt,
    std::optional<int> val_inout_opt = std::nullopt
) {
  auto _result =
      CppBmadTest::test_integer_scalar(val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestIntegerScalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestLogicalScalar python_test_logical_scalar(
    bool val_in,
    bool val_inout,
    std::optional<bool> val_in_opt = std::nullopt,
    std::optional<bool> val_inout_opt = std::nullopt
) {
  auto _result =
      CppBmadTest::test_logical_scalar(val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestLogicalScalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestReal16Scalar python_test_real16_scalar(
    long double val_in,
    long double val_inout,
    std::optional<long double> val_in_opt = std::nullopt,
    std::optional<long double> val_inout_opt = std::nullopt
) {
  auto _result =
      CppBmadTest::test_real16_scalar(val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestReal16Scalar{_result, val_inout, val_inout_opt}};
  return py_result;
}
PyTestRealScalar python_test_real_scalar(
    double val_in,
    double val_inout,
    std::optional<double> val_in_opt = std::nullopt,
    std::optional<double> val_inout_opt = std::nullopt
) {
  auto _result =
      CppBmadTest::test_real_scalar(val_in, val_inout, val_in_opt, make_opt_ref(val_inout_opt));
  auto py_result{PyTestRealScalar{_result, val_inout, val_inout_opt}};
  return py_result;
}

void init_CppBmadTest_routines_t(nb::module_ &m) {
  nb::class_<CppBmadTest::TestBunchStructArray>(
      m,
      "TestBunchStructArray",
      "test_bunch_struct_array return type"
  )
      .def_ro("arr_out", &CppBmadTest::TestBunchStructArray::arr_out)
      .def_ro("opt_status", &CppBmadTest::TestBunchStructArray::opt_status)
      .def("__len__", [](const CppBmadTest::TestBunchStructArray &) { return 2; })
      .def("__getitem__", [](const CppBmadTest::TestBunchStructArray &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.arr_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        throw nb::index_error();
      });
  m.def(
      "test_bunch_struct_array",
      &CppBmadTest::test_bunch_struct_array,
      nb::arg("arr_in"),
      nb::arg("arr_inout"),
      nb::arg("arr_in_opt") = nb::none(),
      nb::arg("arr_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_bunch_struct_array

Parameters
----------
arr_in : 1D array of BunchStruct

arr_inout : 1D array of BunchStruct

arr_in_opt : 1D array of BunchStruct, optional

arr_inout_opt : 1D array of BunchStruct, optional

Returns
-------
arr_out : 1D array of BunchStruct

opt_status : 1D array of int (shape: 2)
)"""
  );
  nb::class_<CppBmadTest::TestBunchStructScalar>(
      m,
      "TestBunchStructScalar",
      "test_bunch_struct_scalar return type"
  )
      .def_ro("val_out", &CppBmadTest::TestBunchStructScalar::val_out)
      .def_ro("opt_status", &CppBmadTest::TestBunchStructScalar::opt_status)
      .def("__len__", [](const CppBmadTest::TestBunchStructScalar &) { return 2; })
      .def("__getitem__", [](const CppBmadTest::TestBunchStructScalar &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.val_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        throw nb::index_error();
      });
  m.def(
      "test_bunch_struct_scalar",
      [](BunchStruct &val_in,
         BunchStruct &val_inout,
         BunchStruct *val_in_opt,
         BunchStruct *val_inout_opt) {
        auto fn = static_cast<
            CppBmadTest::
                TestBunchStructScalar (*)(BunchStruct &, BunchStruct &, optional_ref<BunchStruct>, optional_ref<BunchStruct>)>(
            &CppBmadTest::test_bunch_struct_scalar
        );
        return fn(val_in, val_inout, ptr_to_opt_ref(val_in_opt), ptr_to_opt_ref(val_inout_opt));
      },
      nb::arg("val_in"),
      nb::arg("val_inout"),
      nb::arg("val_in_opt") = nb::none(),
      nb::arg("val_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_bunch_struct_scalar

Parameters
----------
val_in : BunchStruct

val_inout : BunchStruct

val_in_opt : BunchStruct, optional

val_inout_opt : BunchStruct, optional

Returns
-------
val_out : BunchStruct

opt_status : 1D array of int (shape: 2)
)"""
  );
  nb::class_<CppBmadTest::TestCharacterArray>(
      m,
      "TestCharacterArray",
      "test_character_array return type"
  )
      .def_ro("arr_out", &CppBmadTest::TestCharacterArray::arr_out)
      .def_ro("opt_status", &CppBmadTest::TestCharacterArray::opt_status)
      .def("__len__", [](const CppBmadTest::TestCharacterArray &) { return 2; })
      .def("__getitem__", [](const CppBmadTest::TestCharacterArray &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.arr_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        throw nb::index_error();
      });
  m.def(
      "test_character_array",
      [](CharacterAlloc1D &arr_in,
         CharacterAlloc1D &arr_inout,
         CharacterAlloc1D *arr_in_opt,
         CharacterAlloc1D *arr_inout_opt) {
        auto fn = static_cast<
            CppBmadTest::
                TestCharacterArray (*)(CharacterAlloc1D &, CharacterAlloc1D &, optional_ref<CharacterAlloc1D>, optional_ref<CharacterAlloc1D>)>(
            &CppBmadTest::test_character_array
        );
        return fn(arr_in, arr_inout, ptr_to_opt_ref(arr_in_opt), ptr_to_opt_ref(arr_inout_opt));
      },
      nb::arg("arr_in"),
      nb::arg("arr_inout"),
      nb::arg("arr_in_opt") = nb::none(),
      nb::arg("arr_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_character_array

Parameters
----------
arr_in : 1D array of str

arr_inout : 1D array of str

arr_in_opt : 1D array of str, optional

arr_inout_opt : 1D array of str, optional

Returns
-------
arr_out : 1D array of str

opt_status : 1D array of int (shape: 2)
)"""
  );
  nb::class_<PyTestCharacterScalar>(m, "TestCharacterScalar", "test_character_scalar return type")
      .def_ro("val_out", &PyTestCharacterScalar::val_out)
      .def_ro("opt_status", &PyTestCharacterScalar::opt_status)
      .def_ro("val_inout", &PyTestCharacterScalar::val_inout)
      .def_ro("val_inout_opt", &PyTestCharacterScalar::val_inout_opt)
      .def("__len__", [](const PyTestCharacterScalar &) { return 4; })
      .def("__getitem__", [](const PyTestCharacterScalar &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.val_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        if (i == 2)
          return nb::cast(s.val_inout);
        if (i == 3)
          return nb::cast(s.val_inout_opt);
        throw nb::index_error();
      });
  m.def(
      "test_character_scalar",
      &python_test_character_scalar,
      nb::arg("val_in"),
      nb::arg("val_inout"),
      nb::arg("val_in_opt") = nb::none(),
      nb::arg("val_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_character_scalar

Parameters
----------
val_in : str

val_inout : str

val_in_opt : str, optional

val_inout_opt : str, optional

Returns
-------
val_inout : str

val_out : str

opt_status : 1D array of int (shape: 2)

val_inout_opt : str, optional
)"""
  );
  nb::class_<CppBmadTest::TestComplexArray>(m, "TestComplexArray", "test_complex_array return type")
      .def_ro("arr_out", &CppBmadTest::TestComplexArray::arr_out)
      .def_ro("opt_status", &CppBmadTest::TestComplexArray::opt_status)
      .def("__len__", [](const CppBmadTest::TestComplexArray &) { return 2; })
      .def("__getitem__", [](const CppBmadTest::TestComplexArray &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.arr_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        throw nb::index_error();
      });
  m.def(
      "test_complex_array",
      &CppBmadTest::test_complex_array,
      nb::arg("arr_in"),
      nb::arg("arr_inout"),
      nb::arg("arr_in_opt") = nb::none(),
      nb::arg("arr_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_complex_array

Parameters
----------
arr_in : 1D array of complex

arr_inout : 1D array of complex

arr_in_opt : 1D array of complex, optional

arr_inout_opt : 1D array of complex, optional

Returns
-------
arr_out : 1D array of complex

opt_status : 1D array of int (shape: 2)
)"""
  );
  nb::class_<PyTestComplexScalar>(m, "TestComplexScalar", "test_complex_scalar return type")
      .def_ro("val_out", &PyTestComplexScalar::val_out)
      .def_ro("opt_status", &PyTestComplexScalar::opt_status)
      .def_ro("val_inout", &PyTestComplexScalar::val_inout)
      .def_ro("val_inout_opt", &PyTestComplexScalar::val_inout_opt)
      .def("__len__", [](const PyTestComplexScalar &) { return 4; })
      .def("__getitem__", [](const PyTestComplexScalar &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.val_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        if (i == 2)
          return nb::cast(s.val_inout);
        if (i == 3)
          return nb::cast(s.val_inout_opt);
        throw nb::index_error();
      });
  m.def(
      "test_complex_scalar",
      &python_test_complex_scalar,
      nb::arg("val_in"),
      nb::arg("val_inout"),
      nb::arg("val_in_opt") = nb::none(),
      nb::arg("val_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_complex_scalar

Parameters
----------
val_in : complex

val_inout : complex

val_in_opt : complex, optional

val_inout_opt : complex, optional

Returns
-------
val_inout : complex

val_out : complex

opt_status : 1D array of int (shape: 2)

val_inout_opt : complex, optional
)"""
  );
  nb::class_<CppBmadTest::TestInteger8Array>(
      m,
      "TestInteger8Array",
      "test_integer8_array return type"
  )
      .def_ro("arr_out", &CppBmadTest::TestInteger8Array::arr_out)
      .def_ro("opt_status", &CppBmadTest::TestInteger8Array::opt_status)
      .def("__len__", [](const CppBmadTest::TestInteger8Array &) { return 2; })
      .def("__getitem__", [](const CppBmadTest::TestInteger8Array &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.arr_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        throw nb::index_error();
      });
  m.def(
      "test_integer8_array",
      &CppBmadTest::test_integer8_array,
      nb::arg("arr_in"),
      nb::arg("arr_inout"),
      nb::arg("arr_in_opt") = nb::none(),
      nb::arg("arr_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_integer8_array

Parameters
----------
arr_in : 1D array of int

arr_inout : 1D array of int

arr_in_opt : 1D array of int, optional

arr_inout_opt : 1D array of int, optional

Returns
-------
arr_out : 1D array of int

opt_status : 1D array of int (shape: 2)
)"""
  );
  nb::class_<PyTestInteger8Scalar>(m, "TestInteger8Scalar", "test_integer8_scalar return type")
      .def_ro("val_out", &PyTestInteger8Scalar::val_out)
      .def_ro("opt_status", &PyTestInteger8Scalar::opt_status)
      .def_ro("val_inout", &PyTestInteger8Scalar::val_inout)
      .def_ro("val_inout_opt", &PyTestInteger8Scalar::val_inout_opt)
      .def("__len__", [](const PyTestInteger8Scalar &) { return 4; })
      .def("__getitem__", [](const PyTestInteger8Scalar &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.val_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        if (i == 2)
          return nb::cast(s.val_inout);
        if (i == 3)
          return nb::cast(s.val_inout_opt);
        throw nb::index_error();
      });
  m.def(
      "test_integer8_scalar",
      &python_test_integer8_scalar,
      nb::arg("val_in"),
      nb::arg("val_inout"),
      nb::arg("val_in_opt") = nb::none(),
      nb::arg("val_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_integer8_scalar

Parameters
----------
val_in : int

val_inout : int

val_in_opt : int, optional

val_inout_opt : int, optional

Returns
-------
val_inout : int

val_out : int

opt_status : 1D array of int (shape: 2)

val_inout_opt : int, optional
)"""
  );
  nb::class_<CppBmadTest::TestIntegerArray>(m, "TestIntegerArray", "test_integer_array return type")
      .def_ro("arr_out", &CppBmadTest::TestIntegerArray::arr_out)
      .def_ro("opt_status", &CppBmadTest::TestIntegerArray::opt_status)
      .def("__len__", [](const CppBmadTest::TestIntegerArray &) { return 2; })
      .def("__getitem__", [](const CppBmadTest::TestIntegerArray &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.arr_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        throw nb::index_error();
      });
  m.def(
      "test_integer_array",
      &CppBmadTest::test_integer_array,
      nb::arg("arr_in"),
      nb::arg("arr_inout"),
      nb::arg("arr_in_opt") = nb::none(),
      nb::arg("arr_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_integer_array

Parameters
----------
arr_in : 1D array of int

arr_inout : 1D array of int

arr_in_opt : 1D array of int, optional

arr_inout_opt : 1D array of int, optional

Returns
-------
arr_out : 1D array of int

opt_status : 1D array of int (shape: 2)
)"""
  );
  nb::class_<PyTestIntegerScalar>(m, "TestIntegerScalar", "test_integer_scalar return type")
      .def_ro("val_out", &PyTestIntegerScalar::val_out)
      .def_ro("opt_status", &PyTestIntegerScalar::opt_status)
      .def_ro("val_inout", &PyTestIntegerScalar::val_inout)
      .def_ro("val_inout_opt", &PyTestIntegerScalar::val_inout_opt)
      .def("__len__", [](const PyTestIntegerScalar &) { return 4; })
      .def("__getitem__", [](const PyTestIntegerScalar &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.val_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        if (i == 2)
          return nb::cast(s.val_inout);
        if (i == 3)
          return nb::cast(s.val_inout_opt);
        throw nb::index_error();
      });
  m.def(
      "test_integer_scalar",
      &python_test_integer_scalar,
      nb::arg("val_in"),
      nb::arg("val_inout"),
      nb::arg("val_in_opt") = nb::none(),
      nb::arg("val_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_integer_scalar

Parameters
----------
val_in : int

val_inout : int

val_in_opt : int, optional

val_inout_opt : int, optional

Returns
-------
val_inout : int

val_out : int

opt_status : 1D array of int (shape: 2)

val_inout_opt : int, optional
)"""
  );
  nb::class_<CppBmadTest::TestLogicalArray>(m, "TestLogicalArray", "test_logical_array return type")
      .def_ro("arr_out", &CppBmadTest::TestLogicalArray::arr_out)
      .def_ro("opt_status", &CppBmadTest::TestLogicalArray::opt_status)
      .def("__len__", [](const CppBmadTest::TestLogicalArray &) { return 2; })
      .def("__getitem__", [](const CppBmadTest::TestLogicalArray &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.arr_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        throw nb::index_error();
      });
  m.def(
      "test_logical_array",
      [](BoolAlloc1D &arr_in,
         BoolAlloc1D &arr_inout,
         BoolAlloc1D *arr_in_opt,
         BoolAlloc1D *arr_inout_opt) {
        auto fn = static_cast<
            CppBmadTest::
                TestLogicalArray (*)(BoolAlloc1D &, BoolAlloc1D &, optional_ref<BoolAlloc1D>, optional_ref<BoolAlloc1D>)>(
            &CppBmadTest::test_logical_array
        );
        return fn(arr_in, arr_inout, ptr_to_opt_ref(arr_in_opt), ptr_to_opt_ref(arr_inout_opt));
      },
      nb::arg("arr_in"),
      nb::arg("arr_inout"),
      nb::arg("arr_in_opt") = nb::none(),
      nb::arg("arr_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_logical_array

Parameters
----------
arr_in : 1D array of bool

arr_inout : 1D array of bool

arr_in_opt : 1D array of bool, optional

arr_inout_opt : 1D array of bool, optional

Returns
-------
arr_out : 1D array of bool

opt_status : 1D array of int (shape: 2)
)"""
  );
  nb::class_<PyTestLogicalScalar>(m, "TestLogicalScalar", "test_logical_scalar return type")
      .def_ro("val_out", &PyTestLogicalScalar::val_out)
      .def_ro("opt_status", &PyTestLogicalScalar::opt_status)
      .def_ro("val_inout", &PyTestLogicalScalar::val_inout)
      .def_ro("val_inout_opt", &PyTestLogicalScalar::val_inout_opt)
      .def("__len__", [](const PyTestLogicalScalar &) { return 4; })
      .def("__getitem__", [](const PyTestLogicalScalar &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.val_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        if (i == 2)
          return nb::cast(s.val_inout);
        if (i == 3)
          return nb::cast(s.val_inout_opt);
        throw nb::index_error();
      });
  m.def(
      "test_logical_scalar",
      &python_test_logical_scalar,
      nb::arg("val_in"),
      nb::arg("val_inout"),
      nb::arg("val_in_opt") = nb::none(),
      nb::arg("val_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_logical_scalar

Parameters
----------
val_in : bool

val_inout : bool

val_in_opt : bool, optional

val_inout_opt : bool, optional

Returns
-------
val_inout : bool

val_out : bool

opt_status : 1D array of int (shape: 2)

val_inout_opt : bool, optional
)"""
  );
  nb::class_<PyTestReal16Scalar>(m, "TestReal16Scalar", "test_real16_scalar return type")
      .def_ro("val_out", &PyTestReal16Scalar::val_out)
      .def_ro("opt_status", &PyTestReal16Scalar::opt_status)
      .def_ro("val_inout", &PyTestReal16Scalar::val_inout)
      .def_ro("val_inout_opt", &PyTestReal16Scalar::val_inout_opt)
      .def("__len__", [](const PyTestReal16Scalar &) { return 4; })
      .def("__getitem__", [](const PyTestReal16Scalar &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.val_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        if (i == 2)
          return nb::cast(s.val_inout);
        if (i == 3)
          return nb::cast(s.val_inout_opt);
        throw nb::index_error();
      });
  m.def(
      "test_real16_scalar",
      &python_test_real16_scalar,
      nb::arg("val_in"),
      nb::arg("val_inout"),
      nb::arg("val_in_opt") = nb::none(),
      nb::arg("val_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_real16_scalar

Parameters
----------
val_in : float

val_inout : float

val_in_opt : float, optional

val_inout_opt : float, optional

Returns
-------
val_inout : float

val_out : float

opt_status : 1D array of int (shape: 2)

val_inout_opt : float, optional
)"""
  );
  nb::class_<CppBmadTest::TestRealArray>(m, "TestRealArray", "test_real_array return type")
      .def_ro("arr_out", &CppBmadTest::TestRealArray::arr_out)
      .def_ro("opt_status", &CppBmadTest::TestRealArray::opt_status)
      .def("__len__", [](const CppBmadTest::TestRealArray &) { return 2; })
      .def("__getitem__", [](const CppBmadTest::TestRealArray &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.arr_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        throw nb::index_error();
      });
  m.def(
      "test_real_array",
      &CppBmadTest::test_real_array,
      nb::arg("arr_in"),
      nb::arg("arr_inout"),
      nb::arg("arr_in_opt") = nb::none(),
      nb::arg("arr_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_real_array

Parameters
----------
arr_in : 1D array of float

arr_inout : 1D array of float

arr_in_opt : 1D array of float, optional

arr_inout_opt : 1D array of float, optional

Returns
-------
arr_out : 1D array of float

opt_status : 1D array of int (shape: 2)
)"""
  );
  nb::class_<PyTestRealScalar>(m, "TestRealScalar", "test_real_scalar return type")
      .def_ro("val_out", &PyTestRealScalar::val_out)
      .def_ro("opt_status", &PyTestRealScalar::opt_status)
      .def_ro("val_inout", &PyTestRealScalar::val_inout)
      .def_ro("val_inout_opt", &PyTestRealScalar::val_inout_opt)
      .def("__len__", [](const PyTestRealScalar &) { return 4; })
      .def("__getitem__", [](const PyTestRealScalar &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.val_out);
        if (i == 1)
          return nb::cast(s.opt_status);
        if (i == 2)
          return nb::cast(s.val_inout);
        if (i == 3)
          return nb::cast(s.val_inout_opt);
        throw nb::index_error();
      });
  m.def(
      "test_real_scalar",
      &python_test_real_scalar,
      nb::arg("val_in"),
      nb::arg("val_inout"),
      nb::arg("val_in_opt") = nb::none(),
      nb::arg("val_inout_opt") = nb::none(),
      R"""(Wrapper for Fortran routine test_real_scalar

Parameters
----------
val_in : float

val_inout : float

val_in_opt : float, optional

val_inout_opt : float, optional

Returns
-------
val_inout : float

val_out : float

opt_status : 1D array of int (shape: 2)

val_inout_opt : float, optional
)"""
  );
}
