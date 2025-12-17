from __future__ import annotations

import pybmad
from pytao import Tao

tao = Tao(init_file="$ACC_ROOT_DIR/bmad-doc/tao_examples/optics_matching/tao.init", noplot=True)

print(tao)

print("\n".join(tao.cmd("show uni")))

s = pybmad.get_super_universe()

print("Lattice:", s.u[0].model.lat.use_name)
print("Model element #0 name", s.u[0].model.lat.ele[0].name)
print("Some orbit value", s.u[0].design.tao_branch[0].orbit[0].vec)
