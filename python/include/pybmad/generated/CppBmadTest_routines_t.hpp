#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_CppBmadTest_routines_t(py::module &m);

struct PyTestCharacterScalar : public CppBmadTest::TestCharacterScalar {
  std::string val_inout;
  std::optional<std::string> val_inout_opt;
  PyTestCharacterScalar(
      CppBmadTest::TestCharacterScalar _base,
      std::string val_inout,
      std::optional<std::string> val_inout_opt
  )
      : CppBmadTest::TestCharacterScalar(std::move(_base))
      , val_inout(val_inout)
      , val_inout_opt(val_inout_opt) {}
};
struct PyTestComplexScalar : public CppBmadTest::TestComplexScalar {
  std::complex<double> val_inout;
  std::optional<std::complex<double>> val_inout_opt;
  PyTestComplexScalar(
      CppBmadTest::TestComplexScalar _base,
      std::complex<double> val_inout,
      std::optional<std::complex<double>> val_inout_opt
  )
      : CppBmadTest::TestComplexScalar(std::move(_base))
      , val_inout(val_inout)
      , val_inout_opt(val_inout_opt) {}
};
struct PyTestInteger8Scalar : public CppBmadTest::TestInteger8Scalar {
  int64_t val_inout;
  std::optional<int64_t> val_inout_opt;
  PyTestInteger8Scalar(
      CppBmadTest::TestInteger8Scalar _base,
      int64_t val_inout,
      std::optional<int64_t> val_inout_opt
  )
      : CppBmadTest::TestInteger8Scalar(std::move(_base))
      , val_inout(val_inout)
      , val_inout_opt(val_inout_opt) {}
};
struct PyTestIntegerScalar : public CppBmadTest::TestIntegerScalar {
  int val_inout;
  std::optional<int> val_inout_opt;
  PyTestIntegerScalar(
      CppBmadTest::TestIntegerScalar _base,
      int val_inout,
      std::optional<int> val_inout_opt
  )
      : CppBmadTest::TestIntegerScalar(std::move(_base))
      , val_inout(val_inout)
      , val_inout_opt(val_inout_opt) {}
};
struct PyTestLogicalScalar : public CppBmadTest::TestLogicalScalar {
  bool val_inout;
  std::optional<bool> val_inout_opt;
  PyTestLogicalScalar(
      CppBmadTest::TestLogicalScalar _base,
      bool val_inout,
      std::optional<bool> val_inout_opt
  )
      : CppBmadTest::TestLogicalScalar(std::move(_base))
      , val_inout(val_inout)
      , val_inout_opt(val_inout_opt) {}
};
struct PyTestReal16Scalar : public CppBmadTest::TestReal16Scalar {
  long double val_inout;
  std::optional<long double> val_inout_opt;
  PyTestReal16Scalar(
      CppBmadTest::TestReal16Scalar _base,
      long double val_inout,
      std::optional<long double> val_inout_opt
  )
      : CppBmadTest::TestReal16Scalar(std::move(_base))
      , val_inout(val_inout)
      , val_inout_opt(val_inout_opt) {}
};
struct PyTestRealScalar : public CppBmadTest::TestRealScalar {
  double val_inout;
  std::optional<double> val_inout_opt;
  PyTestRealScalar(
      CppBmadTest::TestRealScalar _base,
      double val_inout,
      std::optional<double> val_inout_opt
  )
      : CppBmadTest::TestRealScalar(std::move(_base))
      , val_inout(val_inout)
      , val_inout_opt(val_inout_opt) {}
};
