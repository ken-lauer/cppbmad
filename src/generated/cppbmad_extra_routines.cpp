#include "bmad/generated/cppbmad_extra_routines.hpp"

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
bool CppBmadExtra::set_ele_misalignments(
    EleStruct &ele,
    double x_offset,
    double y_offset,
    double z_offset,
    double x_pitch,
    double y_pitch,
    double tilt,
    std::optional<bool> check_free
) {
  bool check_free_lvalue;
  auto *_check_free{&check_free_lvalue};
  if (check_free.has_value()) {
    check_free_lvalue = check_free.value();
  } else {
    _check_free = nullptr;
  }
  bool _ok{};
  fortran_set_ele_misalignments(
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ x_offset,
      /* double& */ y_offset,
      /* double& */ z_offset,
      /* double& */ x_pitch,
      /* double& */ y_pitch,
      /* double& */ tilt,
      /* bool* */ _check_free,
      /* bool& */ _ok
  );
  return _ok;
}