#include "bmad/generated/cppbmad_test_routines.hpp"

#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "bmad/generated/proxy.hpp"
#include "bmad/json.hpp"
#include "bmad/types.hpp"

using namespace Bmad;

using json = nlohmann::json;
CppBmadTest::TestBunchStructArray CppBmadTest::test_bunch_struct_array(
    BunchStructArray1D arr_in,
    BunchStructArray1D arr_inout,
    std::optional<BunchStructArray1D> arr_in_opt,
    std::optional<BunchStructArray1D> arr_inout_opt
) {
  // arr_in: BunchStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _arr_in_desc;
  _arr_in_desc.rank = 1;
  _arr_in_desc.data_ptr = arr_in.data();
  _arr_in_desc.dims[0] = arr_in.size();
  _arr_in_desc.strides[0] = 1;
  // arr_inout: BunchStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _arr_inout_desc;
  _arr_inout_desc.rank = 1;
  _arr_inout_desc.data_ptr = arr_inout.data();
  _arr_inout_desc.dims[0] = arr_inout.size();
  _arr_inout_desc.strides[0] = 1;
  // intent=out allocatable type array
  auto arr_out{BunchStructAlloc1D()};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  // arr_in_opt: BunchStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _arr_in_opt_desc;
  _arr_in_opt_desc.rank = 1;
  if (arr_in_opt) {
    _arr_in_opt_desc.data_ptr = arr_in_opt->data();
    _arr_in_opt_desc.dims[0] = arr_in_opt->size();
  } else {
    _arr_in_opt_desc.data_ptr = nullptr;
    _arr_in_opt_desc.dims[0] = 0;
  }
  _arr_in_opt_desc.strides[0] = 1;
  // arr_inout_opt: BunchStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _arr_inout_opt_desc;
  _arr_inout_opt_desc.rank = 1;
  if (arr_inout_opt) {
    _arr_inout_opt_desc.data_ptr = arr_inout_opt->data();
    _arr_inout_opt_desc.dims[0] = arr_inout_opt->size();
  } else {
    _arr_inout_opt_desc.data_ptr = nullptr;
    _arr_inout_opt_desc.dims[0] = 0;
  }
  _arr_inout_opt_desc.strides[0] = 1;
  fortran_test_bunch_struct_array(
      /* Bmad::array_descriptor_t& */ _arr_in_desc,
      /* Bmad::array_descriptor_t& */ _arr_inout_desc,
      /* void* */ arr_out.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* Bmad::array_descriptor_t& */ _arr_in_opt_desc,
      /* Bmad::array_descriptor_t& */ _arr_inout_opt_desc
  );
  return TestBunchStructArray{std::move(arr_out), _opt_status};
}
CppBmadTest::TestBunchStructScalar CppBmadTest::test_bunch_struct_scalar(
    BunchStruct &val_in,
    BunchStruct &val_inout,
    optional_ref<BunchStruct> val_in_opt,
    optional_ref<BunchStruct> val_inout_opt
) {
  BunchStruct _val_out;
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  auto *_val_in_opt =
      val_in_opt.has_value() ? val_in_opt->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_val_inout_opt = val_inout_opt.has_value() ? val_inout_opt->get().get_fortran_ptr()
                                                   : nullptr; // input, optional
  fortran_test_bunch_struct_scalar(
      /* void* */ val_in.get_fortran_ptr(),
      /* void* */ val_inout.get_fortran_ptr(),
      /* void* */ _val_out.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* void* */ _val_in_opt,
      /* void* */ _val_inout_opt
  );
  return TestBunchStructScalar{std::move(_val_out), _opt_status};
}
CppBmadTest::TestCharacterScalar CppBmadTest::test_character_scalar(
    std::string val_in,
    std::string &val_inout,
    std::optional<std::string> val_in_opt,
    optional_ref<std::string> val_inout_opt
) {
  auto _val_in = val_in.c_str();
  auto _val_inout = val_inout.c_str(); // ptr, inout, required
  char _val_out[4096];
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  const char *_val_in_opt = val_in_opt.has_value() ? val_in_opt->c_str() : nullptr;
  const char *_val_inout_opt = val_inout_opt.has_value() ? val_inout_opt->get().c_str() : nullptr;
  fortran_test_character_scalar(
      /* const char* */ _val_in,
      /* const char* */ _val_inout,
      /* const char* */ _val_out,
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* const char* */ _val_in_opt,
      /* const char* */ _val_inout_opt
  );
  return TestCharacterScalar{_val_out, _opt_status};
}
CppBmadTest::TestComplexArray CppBmadTest::test_complex_array(
    FArray1D<Complex> &arr_in,
    FArray1D<Complex> &arr_inout,
    std::optional<FArray1D<Complex>> arr_in_opt,
    std::optional<FArray1D<Complex>> arr_inout_opt
) {
  // arr_in: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_in_desc;
  _arr_in_desc.rank = 1;
  _arr_in_desc.data_ptr = arr_in.data();
  _arr_in_desc.dims[0] = arr_in.size();
  // arr_inout: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_inout_desc;
  _arr_inout_desc.rank = 1;
  _arr_inout_desc.data_ptr = arr_inout.data();
  _arr_inout_desc.dims[0] = arr_inout.size();
  // intent=out allocatable general array
  auto arr_out{ComplexAlloc1D()};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  // arr_in_opt: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_in_opt_desc;
  _arr_in_opt_desc.rank = 1;
  if (arr_in_opt.has_value()) {
    _arr_in_opt_desc.data_ptr = arr_in_opt->data();
    _arr_in_opt_desc.dims[0] = arr_in_opt->size();
  } else {
    _arr_in_opt_desc.data_ptr = nullptr;
    _arr_in_opt_desc.dims[0] = 0;
  }
  // arr_inout_opt: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_inout_opt_desc;
  _arr_inout_opt_desc.rank = 1;
  if (arr_inout_opt.has_value()) {
    _arr_inout_opt_desc.data_ptr = arr_inout_opt->data();
    _arr_inout_opt_desc.dims[0] = arr_inout_opt->size();
  } else {
    _arr_inout_opt_desc.data_ptr = nullptr;
    _arr_inout_opt_desc.dims[0] = 0;
  }
  fortran_test_complex_array(
      /* Bmad::array_descriptor_t& */ _arr_in_desc,
      /* Bmad::array_descriptor_t& */ _arr_inout_desc,
      /* void* */ arr_out.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* Bmad::array_descriptor_t& */ _arr_in_opt_desc,
      /* Bmad::array_descriptor_t& */ _arr_inout_opt_desc
  );
  return TestComplexArray{std::move(arr_out), _opt_status};
}
CppBmadTest::TestComplexScalar CppBmadTest::test_complex_scalar(
    std::complex<double> val_in,
    std::complex<double> &val_inout,
    std::optional<std::complex<double>> val_in_opt,
    optional_ref<std::complex<double>> val_inout_opt
) {
  std::complex<double> _val_out{};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  std::complex<double> val_in_opt_lvalue;
  auto *_val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto *_val_inout_opt =
      val_inout_opt.has_value() ? &val_inout_opt->get() : nullptr; // inout, optional
  fortran_test_complex_scalar(
      /* std::complex<double>& */ val_in,
      /* std::complex<double>& */ val_inout,
      /* std::complex<double>& */ _val_out,
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* std::complex<double>* */ _val_in_opt,
      /* std::complex<double>* */ _val_inout_opt
  );
  return TestComplexScalar{_val_out, _opt_status};
}
CppBmadTest::TestInteger8Array CppBmadTest::test_integer8_array(
    FArray1D<Int8> &arr_in,
    FArray1D<Int8> &arr_inout,
    std::optional<FArray1D<Int8>> arr_in_opt,
    std::optional<FArray1D<Int8>> arr_inout_opt
) {
  // arr_in: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_in_desc;
  _arr_in_desc.rank = 1;
  _arr_in_desc.data_ptr = arr_in.data();
  _arr_in_desc.dims[0] = arr_in.size();
  // arr_inout: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_inout_desc;
  _arr_inout_desc.rank = 1;
  _arr_inout_desc.data_ptr = arr_inout.data();
  _arr_inout_desc.dims[0] = arr_inout.size();
  // intent=out allocatable general array
  auto arr_out{Int8Alloc1D()};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  // arr_in_opt: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_in_opt_desc;
  _arr_in_opt_desc.rank = 1;
  if (arr_in_opt.has_value()) {
    _arr_in_opt_desc.data_ptr = arr_in_opt->data();
    _arr_in_opt_desc.dims[0] = arr_in_opt->size();
  } else {
    _arr_in_opt_desc.data_ptr = nullptr;
    _arr_in_opt_desc.dims[0] = 0;
  }
  // arr_inout_opt: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_inout_opt_desc;
  _arr_inout_opt_desc.rank = 1;
  if (arr_inout_opt.has_value()) {
    _arr_inout_opt_desc.data_ptr = arr_inout_opt->data();
    _arr_inout_opt_desc.dims[0] = arr_inout_opt->size();
  } else {
    _arr_inout_opt_desc.data_ptr = nullptr;
    _arr_inout_opt_desc.dims[0] = 0;
  }
  fortran_test_integer8_array(
      /* Bmad::array_descriptor_t& */ _arr_in_desc,
      /* Bmad::array_descriptor_t& */ _arr_inout_desc,
      /* void* */ arr_out.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* Bmad::array_descriptor_t& */ _arr_in_opt_desc,
      /* Bmad::array_descriptor_t& */ _arr_inout_opt_desc
  );
  return TestInteger8Array{std::move(arr_out), _opt_status};
}
CppBmadTest::TestInteger8Scalar CppBmadTest::test_integer8_scalar(
    int64_t val_in,
    int64_t &val_inout,
    std::optional<int64_t> val_in_opt,
    optional_ref<int64_t> val_inout_opt
) {
  int64_t _val_out{};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  int64_t val_in_opt_lvalue;
  auto *_val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto *_val_inout_opt =
      val_inout_opt.has_value() ? &val_inout_opt->get() : nullptr; // inout, optional
  fortran_test_integer8_scalar(
      /* int64_t& */ val_in,
      /* int64_t& */ val_inout,
      /* int64_t& */ _val_out,
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* int64_t* */ _val_in_opt,
      /* int64_t* */ _val_inout_opt
  );
  return TestInteger8Scalar{_val_out, _opt_status};
}
CppBmadTest::TestIntegerArray CppBmadTest::test_integer_array(
    FArray1D<Int> &arr_in,
    FArray1D<Int> &arr_inout,
    std::optional<FArray1D<Int>> arr_in_opt,
    std::optional<FArray1D<Int>> arr_inout_opt
) {
  // arr_in: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_in_desc;
  _arr_in_desc.rank = 1;
  _arr_in_desc.data_ptr = arr_in.data();
  _arr_in_desc.dims[0] = arr_in.size();
  // arr_inout: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_inout_desc;
  _arr_inout_desc.rank = 1;
  _arr_inout_desc.data_ptr = arr_inout.data();
  _arr_inout_desc.dims[0] = arr_inout.size();
  // intent=out allocatable general array
  auto arr_out{IntAlloc1D()};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  // arr_in_opt: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_in_opt_desc;
  _arr_in_opt_desc.rank = 1;
  if (arr_in_opt.has_value()) {
    _arr_in_opt_desc.data_ptr = arr_in_opt->data();
    _arr_in_opt_desc.dims[0] = arr_in_opt->size();
  } else {
    _arr_in_opt_desc.data_ptr = nullptr;
    _arr_in_opt_desc.dims[0] = 0;
  }
  // arr_inout_opt: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_inout_opt_desc;
  _arr_inout_opt_desc.rank = 1;
  if (arr_inout_opt.has_value()) {
    _arr_inout_opt_desc.data_ptr = arr_inout_opt->data();
    _arr_inout_opt_desc.dims[0] = arr_inout_opt->size();
  } else {
    _arr_inout_opt_desc.data_ptr = nullptr;
    _arr_inout_opt_desc.dims[0] = 0;
  }
  fortran_test_integer_array(
      /* Bmad::array_descriptor_t& */ _arr_in_desc,
      /* Bmad::array_descriptor_t& */ _arr_inout_desc,
      /* void* */ arr_out.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* Bmad::array_descriptor_t& */ _arr_in_opt_desc,
      /* Bmad::array_descriptor_t& */ _arr_inout_opt_desc
  );
  return TestIntegerArray{std::move(arr_out), _opt_status};
}
CppBmadTest::TestIntegerScalar CppBmadTest::test_integer_scalar(
    int val_in,
    int &val_inout,
    std::optional<int> val_in_opt,
    optional_ref<int> val_inout_opt
) {
  int _val_out{};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  int val_in_opt_lvalue;
  auto *_val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto *_val_inout_opt =
      val_inout_opt.has_value() ? &val_inout_opt->get() : nullptr; // inout, optional
  fortran_test_integer_scalar(
      /* int& */ val_in,
      /* int& */ val_inout,
      /* int& */ _val_out,
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* int* */ _val_in_opt,
      /* int* */ _val_inout_opt
  );
  return TestIntegerScalar{_val_out, _opt_status};
}
CppBmadTest::TestLogicalArray CppBmadTest::test_logical_array(
    BoolAlloc1D &arr_in,
    BoolAlloc1D &arr_inout,
    optional_ref<BoolAlloc1D> arr_in_opt,
    optional_ref<BoolAlloc1D> arr_inout_opt
) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  // intent=out allocatable general array
  auto arr_out{BoolAlloc1D()};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  // intent=in allocatable general array
  auto *_arr_in_opt =
      arr_in_opt.has_value() ? arr_in_opt->get().get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable general array
  auto *_arr_inout_opt = arr_inout_opt.has_value() ? arr_inout_opt->get().get_fortran_ptr()
                                                   : nullptr; // input, optional
  fortran_test_logical_array(
      /* void* */ arr_in.get_fortran_ptr(),
      /* void* */ arr_inout.get_fortran_ptr(),
      /* void* */ arr_out.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* void* */ _arr_in_opt,
      /* void* */ _arr_inout_opt
  );
  return TestLogicalArray{std::move(arr_out), _opt_status};
}
CppBmadTest::TestLogicalScalar CppBmadTest::test_logical_scalar(
    bool val_in,
    bool &val_inout,
    std::optional<bool> val_in_opt,
    optional_ref<bool> val_inout_opt
) {
  bool _val_out{};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  bool val_in_opt_lvalue;
  auto *_val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto *_val_inout_opt =
      val_inout_opt.has_value() ? &val_inout_opt->get() : nullptr; // inout, optional
  fortran_test_logical_scalar(
      /* bool& */ val_in,
      /* bool& */ val_inout,
      /* bool& */ _val_out,
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* bool* */ _val_in_opt,
      /* bool* */ _val_inout_opt
  );
  return TestLogicalScalar{_val_out, _opt_status};
}
CppBmadTest::TestReal16Array CppBmadTest::test_real16_array(
    Real16Alloc1D &arr_in,
    Real16Alloc1D &arr_inout,
    optional_ref<Real16Alloc1D> arr_in_opt,
    optional_ref<Real16Alloc1D> arr_inout_opt
) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  // intent=out allocatable general array
  auto arr_out{Real16Alloc1D()};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  // intent=in allocatable general array
  auto *_arr_in_opt =
      arr_in_opt.has_value() ? arr_in_opt->get().get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable general array
  auto *_arr_inout_opt = arr_inout_opt.has_value() ? arr_inout_opt->get().get_fortran_ptr()
                                                   : nullptr; // input, optional
  fortran_test_real16_array(
      /* void* */ arr_in.get_fortran_ptr(),
      /* void* */ arr_inout.get_fortran_ptr(),
      /* void* */ arr_out.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* void* */ _arr_in_opt,
      /* void* */ _arr_inout_opt
  );
  return TestReal16Array{std::move(arr_out), _opt_status};
}
CppBmadTest::TestReal16Scalar CppBmadTest::test_real16_scalar(
    long double val_in,
    long double &val_inout,
    std::optional<long double> val_in_opt,
    optional_ref<long double> val_inout_opt
) {
  long double _val_out{};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  long double val_in_opt_lvalue;
  auto *_val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto *_val_inout_opt =
      val_inout_opt.has_value() ? &val_inout_opt->get() : nullptr; // inout, optional
  fortran_test_real16_scalar(
      /* long double& */ val_in,
      /* long double& */ val_inout,
      /* long double& */ _val_out,
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* long double* */ _val_in_opt,
      /* long double* */ _val_inout_opt
  );
  return TestReal16Scalar{_val_out, _opt_status};
}
CppBmadTest::TestRealArray CppBmadTest::test_real_array(
    FArray1D<Real> &arr_in,
    FArray1D<Real> &arr_inout,
    std::optional<FArray1D<Real>> arr_in_opt,
    std::optional<FArray1D<Real>> arr_inout_opt
) {
  // arr_in: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_in_desc;
  _arr_in_desc.rank = 1;
  _arr_in_desc.data_ptr = arr_in.data();
  _arr_in_desc.dims[0] = arr_in.size();
  // arr_inout: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_inout_desc;
  _arr_inout_desc.rank = 1;
  _arr_inout_desc.data_ptr = arr_inout.data();
  _arr_inout_desc.dims[0] = arr_inout.size();
  // intent=out allocatable general array
  auto arr_out{RealAlloc1D()};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  // arr_in_opt: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_in_opt_desc;
  _arr_in_opt_desc.rank = 1;
  if (arr_in_opt.has_value()) {
    _arr_in_opt_desc.data_ptr = arr_in_opt->data();
    _arr_in_opt_desc.dims[0] = arr_in_opt->size();
  } else {
    _arr_in_opt_desc.data_ptr = nullptr;
    _arr_in_opt_desc.dims[0] = 0;
  }
  // arr_inout_opt: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_inout_opt_desc;
  _arr_inout_opt_desc.rank = 1;
  if (arr_inout_opt.has_value()) {
    _arr_inout_opt_desc.data_ptr = arr_inout_opt->data();
    _arr_inout_opt_desc.dims[0] = arr_inout_opt->size();
  } else {
    _arr_inout_opt_desc.data_ptr = nullptr;
    _arr_inout_opt_desc.dims[0] = 0;
  }
  fortran_test_real_array(
      /* Bmad::array_descriptor_t& */ _arr_in_desc,
      /* Bmad::array_descriptor_t& */ _arr_inout_desc,
      /* void* */ arr_out.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* Bmad::array_descriptor_t& */ _arr_in_opt_desc,
      /* Bmad::array_descriptor_t& */ _arr_inout_opt_desc
  );
  return TestRealArray{std::move(arr_out), _opt_status};
}
CppBmadTest::TestRealScalar CppBmadTest::test_real_scalar(
    double val_in,
    double &val_inout,
    std::optional<double> val_in_opt,
    optional_ref<double> val_inout_opt
) {
  double _val_out{};
  // opt_status: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _opt_status_desc;
  _opt_status_desc.rank = 1;
  FixedArray1D<Int, 2> _opt_status;
  _opt_status_desc.data_ptr = _opt_status.data();
  _opt_status_desc.dims[0] = _opt_status.size();
  double val_in_opt_lvalue;
  auto *_val_in_opt{&val_in_opt_lvalue};
  if (val_in_opt.has_value()) {
    val_in_opt_lvalue = val_in_opt.value();
  } else {
    _val_in_opt = nullptr;
  }
  auto *_val_inout_opt =
      val_inout_opt.has_value() ? &val_inout_opt->get() : nullptr; // inout, optional
  fortran_test_real_scalar(
      /* double& */ val_in,
      /* double& */ val_inout,
      /* double& */ _val_out,
      /* Bmad::array_descriptor_t& */ _opt_status_desc,
      /* double* */ _val_in_opt,
      /* double* */ _val_inout_opt
  );
  return TestRealScalar{_val_out, _opt_status};
}