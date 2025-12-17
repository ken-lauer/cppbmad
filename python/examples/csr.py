from __future__ import annotations

import pathlib
import time

import pybmad

t0 = time.monotonic()

res = pybmad.bmad_parser("data/csr_example/lat.bmad")
if res.err_flag:
    raise RuntimeError("bmad_parser failed")

lat = res.lat
pybmad.ran_seed_put(123456)

beam_init = pybmad.BeamInitStruct()
beam_init.a_norm_emit = 4e-12
beam_init.b_norm_emit = 4e-12
beam_init.dPz_dz = 0.0
beam_init.sig_z = 0.3e-3
beam_init.sig_pz = 0e-20
beam_init.bunch_charge = 0.01e-10
beam_init.n_particle = 1000
beam_init.n_bunch = 1

bmad_com = pybmad.get_bmad_com()
bmad_com.csr_and_space_charge_on = True

space_charge_com = pybmad.get_space_charge_com()
space_charge_com.ds_track_step = 0.1
space_charge_com.n_bin = 400
space_charge_com.beam_chamber_height = 0.02
space_charge_com.n_shield_images = 0
space_charge_com.particle_bin_span = 8

ele0 = lat.ele[0]
lat_param = lat.param
beam, err_flag, beam_init_set, _conserve_momentum = pybmad.init_beam_distribution(ele0, lat_param, beam_init)
if err_flag:
    raise RuntimeError("init_beam_distribution failed (1)")


# First bunch and its particles
bunch = beam.bunch[0]
particles = bunch.particle
n_particles = len(particles)

# Calculate the average (centroid)
ave = [0.0] * 6
for i in range(6):
    total = sum(p.vec[i] for p in particles)
    ave[i] = (total / n_particles) if n_particles > 0 else 0.0

centroid = pybmad.CoordProxyAlloc1D()
pybmad.reallocate_coord_lat(centroid, lat, 0)
pybmad.init_coord1(centroid[0], ave, ele0, pybmad.DOWNSTREAM_END)

pybmad.track_all(lat, centroid)

beam1, err_flag, beam_init_set, _conserve_momentum = pybmad.init_beam_distribution(ele0, lat_param, beam_init)
if err_flag:
    raise RuntimeError("init_beam_distribution failed (2)")

pybmad.track_beam(lat, beam1, ele1=None, ele2=None, centroid=centroid)

print("First particle coords at end of lattice:")
print(list(beam1.bunch[0].particle[0].vec))

with pathlib.Path("csr.dat").open("w") as out:
    idx = 1
    for part in beam1.bunch[0].particle:
        vec_str = " ".join([f"{val:.8f}" for val in part.vec])
        out.write(f"{idx:8d} ({vec_str})\n")
        idx += 1
t1 = time.monotonic()

print(f"Elapsed: {t1 - t0:0.3} s")
