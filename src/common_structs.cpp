#include "bmad/common_structs.hpp"

namespace Bmad {
Bmad::BmadCommonProxy get_bmad_com() {
  return Bmad::BmadCommonProxy(Bmad::bmad_get_bmad_com());
}
Bmad::SpaceChargeCommonProxy get_space_charge_com() {
  return Bmad::SpaceChargeCommonProxy(Bmad::bmad_get_space_charge_com());
}
} // namespace Bmad

namespace Tao {
Bmad::TaoSuperUniverseProxy get_super_universe() {
  return Bmad::TaoSuperUniverseProxy(Tao::tao_get_super_universe_ptr());
}
} // namespace Tao
