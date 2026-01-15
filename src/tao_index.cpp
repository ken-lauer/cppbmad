#include "bmad/tao_index.hpp"

using namespace Bmad;
using namespace Tao;

TaoUniverseStruct TaoUniverseIndexProxy::operator*() const {
  int n_universes = tao_get_n_universes();
  if (ix_uni_ < 1 || ix_uni_ > n_universes) {
    throw InvalidIndexException("universe", ix_uni_, n_universes);
  }
  return TaoUniverseStruct(tao_c_get_universe_ptr(ix_uni_));
}
