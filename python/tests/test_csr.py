from __future__ import annotations

import sys

import numpy as np
import numpy.testing
import pytest
from conftest import CPPBMAD_REPO_ROOT, TESTS_ROOT

import pybmad as pb


def test_csr():
    lat = pb.bmad_parser(f"{CPPBMAD_REPO_ROOT}/data/csr_example/lat.bmad").lat

    pb.ran_seed_put(123456)

    beam_init = pb.BeamInitStruct(
        a_norm_emit=4e-12,
        b_norm_emit=4e-12,
        dPz_dz=0,
        sig_z=0.3e-3,
        sig_pz=0e-20,
        bunch_charge=0.01e-10,
        n_particle=1000,
        n_bunch=1,
    )

    bmad_com = pb.get_bmad_com()
    bmad_com.csr_and_space_charge_on = True

    space_charge_com = pb.get_space_charge_com()
    space_charge_com.ds_track_step = 0.1
    space_charge_com.n_bin = 400
    space_charge_com.beam_chamber_height = 0.02
    space_charge_com.n_shield_images = 0
    space_charge_com.particle_bin_span = 8

    ele0 = lat.ele[0]
    lat_param = lat.param

    if lat_param is None:
        raise RuntimeError("lat_param is None?")
    beam = pb.init_beam_distribution(ele0, lat_param, beam_init).beam

    # First bunch and its particles
    bunch = beam.bunch[0]
    particles = bunch.particle
    n_particles = len(particles)

    # Calculate the average (centroid)
    pmatrix = np.array([p.vec for p in particles])
    ave = np.zeros(6)
    if n_particles > 0:
        ave = np.mean(pmatrix, axis=0)

    centroid = pb.CoordStructAlloc1D()
    pb.reallocate_coord(centroid, lat, 0)
    centroid0 = centroid[0]
    pb.init_coord(centroid0, ave.tolist(), ele0, pb.DOWNSTREAM_END)

    pb.track_all(lat, centroid)

    beam1 = pb.init_beam_distribution(ele0, lat_param, beam_init).beam

    pb.track_beam(lat, beam1, centroid=centroid.view())

    first_particle_vec = np.asarray(beam1.bunch[0].particle[0].vec)

    expected = (
        1.44484e-07,
        8.51128e-09,
        -1.21621e-07,
        -4.29278e-09,
        2.06790e-04,
        2.52534e-08,
    )

    numpy.testing.assert_allclose(desired=expected, actual=first_particle_vec, atol=1e-6)

    with (TESTS_ROOT / "csr_python.dat").open("w") as out:
        for idx, part in enumerate(beam1.bunch[0].particle, start=1):
            formatted_vec = " ".join(f"{val:.8f}" for val in part.vec)
            out.write(f"{idx:8d} {formatted_vec}\n")


if __name__ == "__main__":
    sys.exit(pytest.main(["-v", __file__]))
