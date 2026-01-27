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
ele = pybmad.pointer_to_ele(lat, 0)
print(ele)

ele_ptrs = pybmad.ElePointerStruct.new_array1d(0)

res = pybmad.lat_ele_locator(loc_str="BEGINNING", lat=lat, eles=ele_ptrs, n_loc=0)
assert len(ele_ptrs) == 1

res = pybmad.lat_ele_locator(loc_str="*", lat=lat, eles=ele_ptrs, n_loc=0)
assert len(ele_ptrs) == lat.n_ele_track + 2  # + beginning/end
