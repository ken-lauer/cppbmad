#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.hpp"

using namespace Bmad;

namespace CppBmadTest {
extern "C" void fortran_test_bunch_struct_array(
    Bmad::array_descriptor_t& arr_in /* 1D_NOT_type in */,
    Bmad::array_descriptor_t& arr_inout /* 1D_NOT_type inout */,
    void* arr_out /* 1D_ALLOC_type inout */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    Bmad::array_descriptor_t& arr_in_opt /* 1D_NOT_type in */,
    Bmad::array_descriptor_t& arr_inout_opt /* 1D_NOT_type inout */);
FixedArray1D<Int, 2> test_bunch_struct_array(
    BunchProxyArray1D& arr_in,
    BunchProxyArray1D& arr_inout,
    BunchProxyAlloc1D& arr_out,
    optional_ref<BunchProxyArray1D> arr_in_opt = std::nullopt,
    optional_ref<BunchProxyArray1D> arr_inout_opt = std::nullopt);
extern "C" void fortran_test_bunch_struct_scalar(
    void* val_in /* 0D_NOT_type in */,
    void* val_inout /* 0D_NOT_type inout */,
    void* val_out /* 0D_NOT_type out */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    void* val_in_opt /* 0D_NOT_type in */,
    void* val_inout_opt /* 0D_NOT_type inout */);
struct TestBunchStructScalar {
  BunchProxy val_out;
  FixedArray1D<Int, 2> opt_status;
};
CppBmadTest::TestBunchStructScalar test_bunch_struct_scalar(
    BunchProxy& val_in,
    BunchProxy& val_inout,
    optional_ref<BunchProxy> val_in_opt = std::nullopt,
    optional_ref<BunchProxy> val_inout_opt = std::nullopt);

// Skipped unusable routine test_character_array:
// - Variable-sized in character array: 1D_NOT_character
// - Variable-sized inout character array: 1D_NOT_character
// - Variable-sized inout character array: 1D_ALLOC_character
// - Variable-sized in character array: 1D_NOT_character
// - Variable-sized inout character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_test_character_scalar(
    const char* val_in /* 0D_NOT_character in */,
    const char* val_inout /* 0D_NOT_character inout */,
    const char* val_out /* 0D_NOT_character out */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    const char* val_in_opt /* 0D_NOT_character in */,
    const char* val_inout_opt /* 0D_NOT_character inout */);
struct TestCharacterScalar {
  std::string val_out;
  FixedArray1D<Int, 2> opt_status;
};
CppBmadTest::TestCharacterScalar test_character_scalar(
    std::string val_in,
    std::string& val_inout,
    std::optional<std::string> val_in_opt = std::nullopt,
    optional_ref<std::string> val_inout_opt = std::nullopt);
extern "C" void fortran_test_complex_array(
    Bmad::array_descriptor_t& arr_in /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t& arr_inout /* 1D_NOT_complex inout */,
    void* arr_out /* 1D_ALLOC_complex inout */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    Bmad::array_descriptor_t& arr_in_opt /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t& arr_inout_opt /* 1D_NOT_complex inout */);
FixedArray1D<Int, 2> test_complex_array(
    FArray1D<Complex>& arr_in,
    FArray1D<Complex>& arr_inout,
    ComplexAlloc1D& arr_out,
    std::optional<FArray1D<Complex>> arr_in_opt = std::nullopt,
    std::optional<FArray1D<Complex>> arr_inout_opt = std::nullopt);
extern "C" void fortran_test_complex_scalar(
    std::complex<double>& val_in /* 0D_NOT_complex in */,
    std::complex<double>& val_inout /* 0D_NOT_complex inout */,
    std::complex<double>& val_out /* 0D_NOT_complex out */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    std::complex<double>* val_in_opt /* 0D_NOT_complex in */,
    std::complex<double>* val_inout_opt /* 0D_NOT_complex inout */);
struct TestComplexScalar {
  std::complex<double> val_out;
  FixedArray1D<Int, 2> opt_status;
};
CppBmadTest::TestComplexScalar test_complex_scalar(
    std::complex<double> val_in,
    std::complex<double>& val_inout,
    std::optional<std::complex<double>> val_in_opt = std::nullopt,
    optional_ref<std::complex<double>> val_inout_opt = std::nullopt);
extern "C" void fortran_test_integer8_array(
    Bmad::array_descriptor_t& arr_in /* 1D_NOT_integer8 in */,
    Bmad::array_descriptor_t& arr_inout /* 1D_NOT_integer8 inout */,
    void* arr_out /* 1D_ALLOC_integer8 inout */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    Bmad::array_descriptor_t& arr_in_opt /* 1D_NOT_integer8 in */,
    Bmad::array_descriptor_t& arr_inout_opt /* 1D_NOT_integer8 inout */);
FixedArray1D<Int, 2> test_integer8_array(
    FArray1D<Int8>& arr_in,
    FArray1D<Int8>& arr_inout,
    Int8Alloc1D& arr_out,
    std::optional<FArray1D<Int8>> arr_in_opt = std::nullopt,
    std::optional<FArray1D<Int8>> arr_inout_opt = std::nullopt);
extern "C" void fortran_test_integer8_scalar(
    int64_t& val_in /* 0D_NOT_integer8 in */,
    int64_t& val_inout /* 0D_NOT_integer8 inout */,
    int64_t& val_out /* 0D_NOT_integer8 out */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    int64_t* val_in_opt /* 0D_NOT_integer8 in */,
    int64_t* val_inout_opt /* 0D_NOT_integer8 inout */);
struct TestInteger8Scalar {
  int64_t val_out;
  FixedArray1D<Int, 2> opt_status;
};
CppBmadTest::TestInteger8Scalar test_integer8_scalar(
    int64_t val_in,
    int64_t& val_inout,
    std::optional<int64_t> val_in_opt = std::nullopt,
    optional_ref<int64_t> val_inout_opt = std::nullopt);
extern "C" void fortran_test_integer_array(
    Bmad::array_descriptor_t& arr_in /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t& arr_inout /* 1D_NOT_integer inout */,
    void* arr_out /* 1D_ALLOC_integer inout */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    Bmad::array_descriptor_t& arr_in_opt /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t& arr_inout_opt /* 1D_NOT_integer inout */);
FixedArray1D<Int, 2> test_integer_array(
    FArray1D<Int>& arr_in,
    FArray1D<Int>& arr_inout,
    IntAlloc1D& arr_out,
    std::optional<FArray1D<Int>> arr_in_opt = std::nullopt,
    std::optional<FArray1D<Int>> arr_inout_opt = std::nullopt);
extern "C" void fortran_test_integer_scalar(
    int& val_in /* 0D_NOT_integer in */,
    int& val_inout /* 0D_NOT_integer inout */,
    int& val_out /* 0D_NOT_integer out */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    int* val_in_opt /* 0D_NOT_integer in */,
    int* val_inout_opt /* 0D_NOT_integer inout */);
struct TestIntegerScalar {
  int val_out;
  FixedArray1D<Int, 2> opt_status;
};
CppBmadTest::TestIntegerScalar test_integer_scalar(
    int val_in,
    int& val_inout,
    std::optional<int> val_in_opt = std::nullopt,
    optional_ref<int> val_inout_opt = std::nullopt);
extern "C" void fortran_test_logical_array(
    void* arr_in /* 1D_ALLOC_logical in */,
    void* arr_inout /* 1D_ALLOC_logical inout */,
    void* arr_out /* 1D_ALLOC_logical inout */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    void* arr_in_opt /* 1D_ALLOC_logical in */,
    void* arr_inout_opt /* 1D_ALLOC_logical inout */);
FixedArray1D<Int, 2> test_logical_array(
    BoolAlloc1D& arr_in,
    BoolAlloc1D& arr_inout,
    BoolAlloc1D& arr_out,
    optional_ref<BoolAlloc1D> arr_in_opt = std::nullopt,
    optional_ref<BoolAlloc1D> arr_inout_opt = std::nullopt);
extern "C" void fortran_test_logical_scalar(
    bool& val_in /* 0D_NOT_logical in */,
    bool& val_inout /* 0D_NOT_logical inout */,
    bool& val_out /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    bool* val_in_opt /* 0D_NOT_logical in */,
    bool* val_inout_opt /* 0D_NOT_logical inout */);
struct TestLogicalScalar {
  bool val_out;
  FixedArray1D<Int, 2> opt_status;
};
CppBmadTest::TestLogicalScalar test_logical_scalar(
    bool val_in,
    bool& val_inout,
    std::optional<bool> val_in_opt = std::nullopt,
    optional_ref<bool> val_inout_opt = std::nullopt);
extern "C" void fortran_test_real16_array(
    void* arr_in /* 1D_ALLOC_real16 in */,
    void* arr_inout /* 1D_ALLOC_real16 inout */,
    void* arr_out /* 1D_ALLOC_real16 inout */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    void* arr_in_opt /* 1D_ALLOC_real16 in */,
    void* arr_inout_opt /* 1D_ALLOC_real16 inout */);
FixedArray1D<Int, 2> test_real16_array(
    Real16Alloc1D& arr_in,
    Real16Alloc1D& arr_inout,
    Real16Alloc1D& arr_out,
    optional_ref<Real16Alloc1D> arr_in_opt = std::nullopt,
    optional_ref<Real16Alloc1D> arr_inout_opt = std::nullopt);
extern "C" void fortran_test_real16_scalar(
    long double& val_in /* 0D_NOT_real16 in */,
    long double& val_inout /* 0D_NOT_real16 inout */,
    long double& val_out /* 0D_NOT_real16 out */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    long double* val_in_opt /* 0D_NOT_real16 in */,
    long double* val_inout_opt /* 0D_NOT_real16 inout */);
struct TestReal16Scalar {
  long double val_out;
  FixedArray1D<Int, 2> opt_status;
};
CppBmadTest::TestReal16Scalar test_real16_scalar(
    long double val_in,
    long double& val_inout,
    std::optional<long double> val_in_opt = std::nullopt,
    optional_ref<long double> val_inout_opt = std::nullopt);
extern "C" void fortran_test_real_array(
    Bmad::array_descriptor_t& arr_in /* 1D_NOT_real in */,
    Bmad::array_descriptor_t& arr_inout /* 1D_NOT_real inout */,
    void* arr_out /* 1D_ALLOC_real inout */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    Bmad::array_descriptor_t& arr_in_opt /* 1D_NOT_real in */,
    Bmad::array_descriptor_t& arr_inout_opt /* 1D_NOT_real inout */);
FixedArray1D<Int, 2> test_real_array(
    FArray1D<Real>& arr_in,
    FArray1D<Real>& arr_inout,
    RealAlloc1D& arr_out,
    std::optional<FArray1D<Real>> arr_in_opt = std::nullopt,
    std::optional<FArray1D<Real>> arr_inout_opt = std::nullopt);
extern "C" void fortran_test_real_scalar(
    double& val_in /* 0D_NOT_real in */,
    double& val_inout /* 0D_NOT_real inout */,
    double& val_out /* 0D_NOT_real out */,
    Bmad::array_descriptor_t& opt_status /* 1D_NOT_integer out */,
    double* val_in_opt /* 0D_NOT_real in */,
    double* val_inout_opt /* 0D_NOT_real inout */);
struct TestRealScalar {
  double val_out;
  FixedArray1D<Int, 2> opt_status;
};
CppBmadTest::TestRealScalar test_real_scalar(
    double val_in,
    double& val_inout,
    std::optional<double> val_in_opt = std::nullopt,
    optional_ref<double> val_inout_opt = std::nullopt);
} // namespace CppBmadTest
