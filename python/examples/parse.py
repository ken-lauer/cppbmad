from __future__ import annotations

import logging
import time

import pybmad

logger = logging.getLogger("pybmad-test")
logger.setLevel("DEBUG")
logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s")

t0 = time.monotonic()
res = pybmad.bmad_parser("${ACC_ROOT_DIR}/bmad-doc/tao_examples/fodo/fodo.bmad")
t1 = time.monotonic()
lat = res.lat
print("Parsed lattice {lat.use_name} in {t1 - t0} s (as seen from Python)", lat.use_name, t1 - t0)

print("Branch 0 elements:", lat.branch[0].ele)

for ele in lat.branch[0].ele:
    print(ele, ele.name, ele.lord, ele.n_lord)


res = pybmad.ele_to_taylor(lat.ele[1])
print(res.orbital_taylor)
print(res.spin_taylor)
# orb, spin = res
# print(orb, spin)
