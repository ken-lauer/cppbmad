#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import time

import pybmad

t0 = time.monotonic()

parsed = pybmad.bmad_parser("data/csr_example/lat.bmad")
if parsed.err_flag:
    raise RuntimeError("bmad_parser failed")

lat = parsed.lat
pybmad.ran_seed_put(123456)

beam_init = pybmad.BeamInitStruct(
    a_norm_emit=4e-12,
    b_norm_emit=4e-12,
    dPz_dz=0.0,
    sig_z=0.3e-3,
    sig_pz=0e-20,
    bunch_charge=0.01e-10,
    n_particle=1000,
    n_bunch=1,
)
bmad_com = pybmad.get_bmad_com()
bmad_com.csr_and_space_charge_on = True

space_charge_com = pybmad.get_space_charge_com()
space_charge_com.ds_track_step = 0.1
space_charge_com.n_bin = 400
space_charge_com.beam_chamber_height = 0.02
space_charge_com.n_shield_images = 0
space_charge_com.particle_bin_span = 8


assert lat.param is not None
init = pybmad.init_beam_distribution(lat.ele[0], lat.param, beam_init)
if init.err_flag:
    raise RuntimeError("init_beam_distribution failed (1)")

# First bunch and its particles
beam = init.beam
bunch = beam.bunch[0]
particles = bunch.particle
n_particles = len(particles)

# Calculate the average (centroid)
ave = [0.0] * 6
for i in range(6):
    total = sum(p.vec[i] for p in particles)
    ave[i] = (total / n_particles) if n_particles > 0 else 0.0

centroid = pybmad.CoordStruct.new_array1d(0)
pybmad.reallocate_coord(centroid, lat, 0)
pybmad.init_coord(centroid[0], ave, lat.ele[0], pybmad.DOWNSTREAM_END)

pybmad.track_all(lat, centroid)

init = pybmad.init_beam_distribution(lat.ele[0], lat.param, beam_init)
if init.err_flag:
    raise RuntimeError("init_beam_distribution failed (2)")

beam1 = init.beam
pybmad.track_beam(lat, beam1, ele1=None, ele2=None, centroid=centroid.view())

print("First particle coords at end of lattice:")
print(beam1.bunch[0].particle[0].vec.to_list())

with pathlib.Path("csr.dat").open("w") as out:
    idx = 1
    for part in beam1.bunch[0].particle:
        vec_str = " ".join([f"{val:.8f}" for val in part.vec])
        out.write(f"{idx:8d} ({vec_str})\n")
        idx += 1
t1 = time.monotonic()

print(f"Elapsed: {t1 - t0:0.3} s")
