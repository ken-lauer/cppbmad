#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "bmad/generated/cppbmad_test_routines.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/json.hpp"
#include "bmad/types.h"

using namespace Bmad;

using json = nlohmann::json;
CppBmadTest::TestBunchStructArray CppBmadTest::test_bunch_struct_array(
    BunchProxyAlloc1D& arr_in,
    BunchProxyAlloc1D& arr_inout,
    optional_ref<BunchProxyAlloc1D> arr_in_opt,
    optional_ref<BunchProxyAlloc1D> arr_inout_opt) {
  // intent=in allocatable type array
  // intent=inout allocatable type array
  // intent=out allocatable type array
  auto arr_out{BunchProxyAlloc1D()};
  FixedArray1D<Int, 2> _opt_status;
  // intent=in allocatable type array
  auto* _arr_in_opt = arr_in_opt.has_value()
      ? arr_in_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  // intent=inout allocatable type array
  auto* _arr_inout_opt = arr_inout_opt.has_value()
      ? arr_inout_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_test_bunch_struct_array(
      /* void* */ arr_in.get_fortran_ptr(),
      /* void* */ arr_inout.get_fortran_ptr(),
      /* void* */ arr_out.get_fortran_ptr(),
      /* int* */ _opt_status.data(),
      /* void* */ _arr_in_opt,
      /* void* */ _arr_inout_opt);
  return TestBunchStructArray{std::move(arr_out), _opt_status};
}
CppBmadTest::TestBunchStructScalar CppBmadTest::test_bunch_struct_scalar(
    BunchProxy& val_in,
    BunchProxy& val_inout,
    optional_ref<BunchProxy> val_in_opt,
    optional_ref<BunchProxy> val_inout_opt) {
  BunchProxy _val_out;
  FixedArray1D<Int, 2> _opt_status;
  auto* _val_in_opt = val_in_opt.has_value()
      ? val_in_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  auto* _val_inout_opt = val_inout_opt.has_value()
      ? val_inout_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_test_bunch_struct_scalar(
      /* void* */ val_in.get_fortran_ptr(),
      /* void* */ val_inout.get_fortran_ptr(),
      /* void* */ _val_out.get_fortran_ptr(),
      /* int* */ _opt_status.data(),
      /* void* */ _val_in_opt,
      /* void* */ _val_inout_opt);
  return TestBunchStructScalar{std::move(_val_out), _opt_status};
}
CppBmadTest::TestCharacterScalar CppBmadTest::test_character_scalar(
    std::string val_in,
    std::string& val_inout,
    std::optional<std::string> val_in_opt,
    optional_ref<std::string> val_inout_opt) {
  auto _val_in = val_in.c_str();
  auto _val_inout = val_inout.c_str(); // ptr, inout, required
  char _val_out[4096];
  FixedArray1D<Int, 2> _opt_status;
  const char* _val_in_opt =
      val_in_opt.has_value() ? val_in_opt->c_str() : nullptr;
  const char* _val_inout_opt =
      val_inout_opt.has_value() ? val_inout_opt->get().c_str() : nullptr;
  fortran_test_character_scalar(
      /* const char* */ _val_in,
      /* const char* */ _val_inout,
      /* const char* */ _val_out,
      /* int* */ _opt_status.data(),
      /* const char* */ _val_in_opt,
      /* const char* */ _val_inout_opt);
  return TestCharacterScalar{_val_out, _opt_status};
}
CppBmadTest::TestComplexArray CppBmadTest::test_complex_array(
    ComplexAlloc1D& arr_in,
    ComplexAlloc1D& arr_inout,
    optional_ref<ComplexAlloc1D> arr_in_opt,
    optional_ref<ComplexAlloc1D> arr_inout_opt) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  // intent=out allocatable general array
  auto arr_out{ComplexAlloc1D()};
  FixedArray1D<Int, 2> _opt_status;
  // intent=in allocatable general array
  auto* _arr_in_opt = arr_in_opt.has_value()
      ? arr_in_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  // intent=inout allocatable general array
  auto* _arr_inout_opt = arr_inout_opt.has_value()
      ? arr_inout_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_test_complex_array(
      /* void* */ arr_in.get_fortran_ptr(),
      /* void* */ arr_inout.get_fortran_ptr(),
      /* void* */ arr_out.get_fortran_ptr(),
      /* int* */ _opt_status.data(),
      /* void* */ _arr_in_opt,
      /* void* */ _arr_inout_opt);
  return TestComplexArray{std::move(arr_out), _opt_status};
}
CppBmadTest::TestComplexScalar CppBmadTest::test_complex_scalar(
    std::complex<double> val_in,
    std::complex<double>& val_inout,
    std::optional<std::complex<double>> val_in_opt,
    optional_ref<std::complex<double>> val_inout_opt) {
  std::complex<double> _val_out{};
  FixedArray1D<Int, 2> _opt_status;
  std::complex<double> val_in_opt_lvalue;
  auto* _val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto* _val_inout_opt = val_inout_opt.has_value() ? &val_inout_opt->get()
                                                   : nullptr; // inout, optional
  fortran_test_complex_scalar(
      /* std::complex<double>& */ val_in,
      /* std::complex<double>& */ val_inout,
      /* std::complex<double>& */ _val_out,
      /* int* */ _opt_status.data(),
      /* std::complex<double>* */ _val_in_opt,
      /* std::complex<double>* */ _val_inout_opt);
  return TestComplexScalar{_val_out, _opt_status};
}
CppBmadTest::TestInteger8Array CppBmadTest::test_integer8_array(
    Int8Alloc1D& arr_in,
    Int8Alloc1D& arr_inout,
    optional_ref<Int8Alloc1D> arr_in_opt,
    optional_ref<Int8Alloc1D> arr_inout_opt) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  // intent=out allocatable general array
  auto arr_out{Int8Alloc1D()};
  FixedArray1D<Int, 2> _opt_status;
  // intent=in allocatable general array
  auto* _arr_in_opt = arr_in_opt.has_value()
      ? arr_in_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  // intent=inout allocatable general array
  auto* _arr_inout_opt = arr_inout_opt.has_value()
      ? arr_inout_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_test_integer8_array(
      /* void* */ arr_in.get_fortran_ptr(),
      /* void* */ arr_inout.get_fortran_ptr(),
      /* void* */ arr_out.get_fortran_ptr(),
      /* int* */ _opt_status.data(),
      /* void* */ _arr_in_opt,
      /* void* */ _arr_inout_opt);
  return TestInteger8Array{std::move(arr_out), _opt_status};
}
CppBmadTest::TestInteger8Scalar CppBmadTest::test_integer8_scalar(
    int64_t val_in,
    int64_t& val_inout,
    std::optional<int64_t> val_in_opt,
    optional_ref<int64_t> val_inout_opt) {
  int64_t _val_out{};
  FixedArray1D<Int, 2> _opt_status;
  int64_t val_in_opt_lvalue;
  auto* _val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto* _val_inout_opt = val_inout_opt.has_value() ? &val_inout_opt->get()
                                                   : nullptr; // inout, optional
  fortran_test_integer8_scalar(
      /* int64_t& */ val_in,
      /* int64_t& */ val_inout,
      /* int64_t& */ _val_out,
      /* int* */ _opt_status.data(),
      /* int64_t* */ _val_in_opt,
      /* int64_t* */ _val_inout_opt);
  return TestInteger8Scalar{_val_out, _opt_status};
}
CppBmadTest::TestIntegerArray CppBmadTest::test_integer_array(
    IntAlloc1D& arr_in,
    IntAlloc1D& arr_inout,
    optional_ref<IntAlloc1D> arr_in_opt,
    optional_ref<IntAlloc1D> arr_inout_opt) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  // intent=out allocatable general array
  auto arr_out{IntAlloc1D()};
  FixedArray1D<Int, 2> _opt_status;
  // intent=in allocatable general array
  auto* _arr_in_opt = arr_in_opt.has_value()
      ? arr_in_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  // intent=inout allocatable general array
  auto* _arr_inout_opt = arr_inout_opt.has_value()
      ? arr_inout_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_test_integer_array(
      /* void* */ arr_in.get_fortran_ptr(),
      /* void* */ arr_inout.get_fortran_ptr(),
      /* void* */ arr_out.get_fortran_ptr(),
      /* int* */ _opt_status.data(),
      /* void* */ _arr_in_opt,
      /* void* */ _arr_inout_opt);
  return TestIntegerArray{std::move(arr_out), _opt_status};
}
CppBmadTest::TestIntegerScalar CppBmadTest::test_integer_scalar(
    int val_in,
    int& val_inout,
    std::optional<int> val_in_opt,
    optional_ref<int> val_inout_opt) {
  int _val_out{};
  FixedArray1D<Int, 2> _opt_status;
  int val_in_opt_lvalue;
  auto* _val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto* _val_inout_opt = val_inout_opt.has_value() ? &val_inout_opt->get()
                                                   : nullptr; // inout, optional
  fortran_test_integer_scalar(
      /* int& */ val_in,
      /* int& */ val_inout,
      /* int& */ _val_out,
      /* int* */ _opt_status.data(),
      /* int* */ _val_in_opt,
      /* int* */ _val_inout_opt);
  return TestIntegerScalar{_val_out, _opt_status};
}
CppBmadTest::TestLogicalArray CppBmadTest::test_logical_array(
    BoolAlloc1D& arr_in,
    BoolAlloc1D& arr_inout,
    optional_ref<BoolAlloc1D> arr_in_opt,
    optional_ref<BoolAlloc1D> arr_inout_opt) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  // intent=out allocatable general array
  auto arr_out{BoolAlloc1D()};
  FixedArray1D<Int, 2> _opt_status;
  // intent=in allocatable general array
  auto* _arr_in_opt = arr_in_opt.has_value()
      ? arr_in_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  // intent=inout allocatable general array
  auto* _arr_inout_opt = arr_inout_opt.has_value()
      ? arr_inout_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_test_logical_array(
      /* void* */ arr_in.get_fortran_ptr(),
      /* void* */ arr_inout.get_fortran_ptr(),
      /* void* */ arr_out.get_fortran_ptr(),
      /* int* */ _opt_status.data(),
      /* void* */ _arr_in_opt,
      /* void* */ _arr_inout_opt);
  return TestLogicalArray{std::move(arr_out), _opt_status};
}
CppBmadTest::TestLogicalScalar CppBmadTest::test_logical_scalar(
    bool val_in,
    bool& val_inout,
    std::optional<bool> val_in_opt,
    optional_ref<bool> val_inout_opt) {
  bool _val_out{};
  FixedArray1D<Int, 2> _opt_status;
  bool val_in_opt_lvalue;
  auto* _val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto* _val_inout_opt = val_inout_opt.has_value() ? &val_inout_opt->get()
                                                   : nullptr; // inout, optional
  fortran_test_logical_scalar(
      /* bool& */ val_in,
      /* bool& */ val_inout,
      /* bool& */ _val_out,
      /* int* */ _opt_status.data(),
      /* bool* */ _val_in_opt,
      /* bool* */ _val_inout_opt);
  return TestLogicalScalar{_val_out, _opt_status};
}
CppBmadTest::TestReal16Array CppBmadTest::test_real16_array(
    Real16Alloc1D& arr_in,
    Real16Alloc1D& arr_inout,
    optional_ref<Real16Alloc1D> arr_in_opt,
    optional_ref<Real16Alloc1D> arr_inout_opt) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  // intent=out allocatable general array
  auto arr_out{Real16Alloc1D()};
  FixedArray1D<Int, 2> _opt_status;
  // intent=in allocatable general array
  auto* _arr_in_opt = arr_in_opt.has_value()
      ? arr_in_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  // intent=inout allocatable general array
  auto* _arr_inout_opt = arr_inout_opt.has_value()
      ? arr_inout_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_test_real16_array(
      /* void* */ arr_in.get_fortran_ptr(),
      /* void* */ arr_inout.get_fortran_ptr(),
      /* void* */ arr_out.get_fortran_ptr(),
      /* int* */ _opt_status.data(),
      /* void* */ _arr_in_opt,
      /* void* */ _arr_inout_opt);
  return TestReal16Array{std::move(arr_out), _opt_status};
}
CppBmadTest::TestReal16Scalar CppBmadTest::test_real16_scalar(
    long double val_in,
    long double& val_inout,
    std::optional<long double> val_in_opt,
    optional_ref<long double> val_inout_opt) {
  long double _val_out{};
  FixedArray1D<Int, 2> _opt_status;
  long double val_in_opt_lvalue;
  auto* _val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto* _val_inout_opt = val_inout_opt.has_value() ? &val_inout_opt->get()
                                                   : nullptr; // inout, optional
  fortran_test_real16_scalar(
      /* long double& */ val_in,
      /* long double& */ val_inout,
      /* long double& */ _val_out,
      /* int* */ _opt_status.data(),
      /* long double* */ _val_in_opt,
      /* long double* */ _val_inout_opt);
  return TestReal16Scalar{_val_out, _opt_status};
}
CppBmadTest::TestRealArray CppBmadTest::test_real_array(
    RealAlloc1D& arr_in,
    RealAlloc1D& arr_inout,
    optional_ref<RealAlloc1D> arr_in_opt,
    optional_ref<RealAlloc1D> arr_inout_opt) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  // intent=out allocatable general array
  auto arr_out{RealAlloc1D()};
  FixedArray1D<Int, 2> _opt_status;
  // intent=in allocatable general array
  auto* _arr_in_opt = arr_in_opt.has_value()
      ? arr_in_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  // intent=inout allocatable general array
  auto* _arr_inout_opt = arr_inout_opt.has_value()
      ? arr_inout_opt->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_test_real_array(
      /* void* */ arr_in.get_fortran_ptr(),
      /* void* */ arr_inout.get_fortran_ptr(),
      /* void* */ arr_out.get_fortran_ptr(),
      /* int* */ _opt_status.data(),
      /* void* */ _arr_in_opt,
      /* void* */ _arr_inout_opt);
  return TestRealArray{std::move(arr_out), _opt_status};
}
CppBmadTest::TestRealScalar CppBmadTest::test_real_scalar(
    double val_in,
    double& val_inout,
    std::optional<double> val_in_opt,
    optional_ref<double> val_inout_opt) {
  double _val_out{};
  FixedArray1D<Int, 2> _opt_status;
  double val_in_opt_lvalue;
  auto* _val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto* _val_inout_opt = val_inout_opt.has_value() ? &val_inout_opt->get()
                                                   : nullptr; // inout, optional
  fortran_test_real_scalar(
      /* double& */ val_in,
      /* double& */ val_inout,
      /* double& */ _val_out,
      /* int* */ _opt_status.data(),
      /* double* */ _val_in_opt,
      /* double* */ _val_inout_opt);
  return TestRealScalar{_val_out, _opt_status};
}