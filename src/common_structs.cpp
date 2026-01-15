#include "bmad/common_structs.hpp"

namespace Bmad {
Bmad::BmadCommonStruct get_bmad_com() { return Bmad::BmadCommonStruct(Bmad::bmad_get_bmad_com()); }
Bmad::SpaceChargeCommonStruct get_space_charge_com() {
  return Bmad::SpaceChargeCommonStruct(Bmad::bmad_get_space_charge_com());
}
} // namespace Bmad

namespace Tao {
Bmad::TaoSuperUniverseStruct get_super_universe() {
  return Bmad::TaoSuperUniverseStruct(Tao::tao_get_super_universe_ptr());
}
} // namespace Tao
