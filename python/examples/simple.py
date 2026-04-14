#!/usr/bin/env python3
"""
Port of the Fortran 'simple_bmad_program' program to pybmad.

Read in a lattice, modify CLEO_SOL solenoid KS strength, recalculate
transfer matrices, propagate Twiss, and print element info.
"""

from __future__ import annotations

import pybmad
from pybmad import EleKey

# Read in a lattice
res = pybmad.bmad_parser("$ACC_ROOT_DIR/code_examples/simple_bmad_program/lat.bmad")
if res.err_flag:
    raise RuntimeError("Failed to parse lattice")
lat = res.lat

# Find the CLEO_SOL element
(cleo,) = [ele for ele in lat.ele if ele.name == "CLEO_SOL"]
# See also: lat_ele_locator.py

# Modify the KS solenoid strength
ks_val = cleo.value[pybmad.KS]
cleo.value[pybmad.KS] = ks_val + 0.001

# Set flags and re-bookkeep
pybmad.set_flags_for_changed_attribute(cleo, cleo.value[pybmad.KS])
pybmad.lattice_bookkeeper(lat)

# Remake transfer matrix for the modified element
pybmad.lat_make_mat6(lat, ix_ele=cleo.ix_ele)

# Calculate starting Twiss params if the lattice is closed,
# and then propagate the Twiss parameters through the lattice.
assert lat.param is not None
if lat.param.geometry == pybmad.CLOSED:
    pybmad.twiss_at_start(lat)
pybmad.twiss_propagate_all(lat)

# Print info on the first 11 elements
print(f" {'Ix':>3}  {'Name':<16}  {'Ele_type':<16}  {'S':>12}  {'Beta_a':>12}")
for i in range(11):
    ele = lat.ele[i]
    key_name = EleKey(ele.key).name
    print(f"{i:4d}  {ele.name:<16}  {key_name:<16}  {ele.s:12.4f}  {ele.a.beta:12.4f}")


info = pybmad.type_ele(cleo)

print()
print("!---------------------------------------------------------")
print("! Information on element: CLEO_SOL")
print()
print("\n".join(info.lines))

# Alternatively, all of the information is accessible directly through the element:
key_name = EleKey(cleo.key).name
print(f"  Name:     {cleo.name}")
print(f"  Key:      {key_name}")
print(f"  ix_ele:   {cleo.ix_ele}")
print(f"  s:        {cleo.s:.6f}")
print(f"  KS:       {cleo.value[pybmad.KS]:.6f}")
print(f"  L:        {cleo.value[pybmad.L]:.6f}")
print(f"  beta_a:   {cleo.a.beta:.6f}")
print(f"  beta_b:   {cleo.b.beta:.6f}")
print(f"  alpha_a:  {cleo.a.alpha:.6f}")
print(f"  alpha_b:  {cleo.b.alpha:.6f}")
