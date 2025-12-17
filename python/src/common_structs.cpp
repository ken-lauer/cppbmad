#include "pybmad/common_structs.hpp"

void Pybmad::init_common_structs(py::module& m) {
  m.def(
      "get_bmad_com",
      &Bmad::get_bmad_com,
      "Get the shared BmadCommon structure");
  m.def(
      "get_space_charge_com",
      &Bmad::get_space_charge_com,
      "Get the shared SpaceChargeCommon structure");
  m.def(
      "get_super_universe",
      &Tao::get_super_universe,
      "Get the shared TaoSuperUniverse structure");
}
