#include "bmad/common_structs.hpp"

Bmad::BmadCommonProxy get_bmad_com() {
  return Bmad::BmadCommonProxy(Bmad::bmad_get_bmad_com());
}
Bmad::SpaceChargeCommonProxy get_space_charge_com() {
  return Bmad::SpaceChargeCommonProxy(Bmad::bmad_get_space_charge_com());
}
