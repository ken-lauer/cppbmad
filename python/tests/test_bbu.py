from __future__ import annotations

import sys
from typing import NamedTuple

import pytest
from bbu import hybridize

import pybmad
from pybmad import LatStruct


class BbuExpectedResult(NamedTuple):
    """
    Holds the regression test configuration and expected results for a specific test case.
    """

    current: float
    expected_hom_voltage_gain: float
    expected_growth_rate: float
    expected_lost: bool
    # Tolerances
    hom_voltage_tol: float = 1e-8
    growth_rate_tol: float = 1e-8


# Define test cases individually
results_1A = BbuExpectedResult(
    current=1.0,
    expected_hom_voltage_gain=3.328330438634e-001,
    expected_growth_rate=-1.100114284509e000,
    expected_lost=False,
)

results_1mA = BbuExpectedResult(
    current=0.001,
    expected_hom_voltage_gain=3.589574300813e-001,
    expected_growth_rate=-1.024551476686e000,
    expected_lost=False,
)

# Note: The original code used a tolerance of 4e-8 for the 100A case's growth rate.
results_100A = BbuExpectedResult(
    current=100.0,
    expected_hom_voltage_gain=2.916793553899e-001,
    expected_growth_rate=-1.232100178045e000,
    expected_lost=False,
    growth_rate_tol=4e-8,
)


@pytest.fixture(scope="module")
def bbu_setup_data():
    """
    Performs the expensive parsings and setup once.
    Returns the (lat, bbu_beam, bbu_param, beam_init) tuple ready for tracking.
    """
    # 1. Initialize Parameters (mimicking bbu.init)
    bbu_param = pybmad.BbuParamStruct()
    bbu_param.lat_filename = "$ACC_ROOT_DIR/regression_tests/bbu_test/oneturn_lat.bmad"
    bbu_param.keep_overlays_and_groups = False
    bbu_param.simulation_turns_max = 500
    bbu_param.elname = "T1"
    bbu_param.hybridize = True
    bbu_param.nrep = 5
    bbu_param.limit_factor = 3
    bbu_param.keep_all_lcavities = False
    bbu_param.ran_gauss_sigma_cut = 3
    bbu_param.nstep = 50
    # bbu_param.current is set in the specific test function
    bbu_param.ran_seed = 100
    bbu_param.rel_tol = 0.001
    bbu_param.lat2_filename = ""
    bbu_param.bunch_freq = 1300000000.0

    # 2. Initialize Beam Init
    beam_init = pybmad.BeamInitStruct()
    beam_init.n_particle = 1

    # Calculate dt based on freq
    if bbu_param.bunch_freq != 0.0:
        beam_init.dt_bunch = 1.0 / bbu_param.bunch_freq
    else:
        beam_init.dt_bunch = 0.0

    # 3. Random Seed Setup
    pybmad.ran_seed_put(bbu_param.ran_seed)
    if bbu_param.ran_gauss_sigma_cut > 0:
        pybmad.ran_gauss_converter(set_sigma_cut=bbu_param.ran_gauss_sigma_cut)

    # 4. Parse Lattice
    res = pybmad.bmad_parser(bbu_param.lat_filename)

    lat_in: LatStruct = res.lat

    # 5. Twiss and Track (Closed Orbit)
    orb = pybmad.CoordStruct.new_array1d(0)
    pybmad.twiss_and_track(lat_in, orb)

    # 6. Hybridization
    lat = lat_in
    if bbu_param.hybridize:
        lat = hybridize(
            lat,
            ele_track_end=bbu_param.ele_track_end,
            use_taylor_for_hybrids=bbu_param.use_taylor_for_hybrids,
            keep_all_lcavities=bbu_param.keep_all_lcavities,
        )

    # 7. BBU Setup
    # Create the bbu_beam structure which holds the wakefield state
    bbu_beam = pybmad.BbuBeamStruct()
    pybmad.bbu_setup(lat, beam_init.dt_bunch, bbu_param, bbu_beam)
    return lat, bbu_beam, bbu_param, beam_init


@pytest.mark.parametrize(
    "case_data",
    [results_1A, results_1mA, results_100A],
    ids=["Current_1A", "Current_1mA", "Current_100A"],
)
def test_bbu_tracking(bbu_setup_data, case_data: BbuExpectedResult):
    """
    Runs the BBU tracking for different currents on the same setup.
    """
    lat, bbu_beam, bbu_param, beam_init = bbu_setup_data

    # Update state for specific test case
    bbu_param.current = case_data.current
    beam_init.bunch_charge = bbu_param.current * beam_init.dt_bunch

    (hom_voltage_gain, growth_rate, lost, _irep) = pybmad.bbu_track_all(
        lat=lat,
        bbu_beam=bbu_beam,
        bbu_param=bbu_param,
        beam_init=beam_init,
        # TODO: these aren't marked as output, so we have to specify them
        hom_voltage_normalized=0.0,
        growth_rate=0.0,
        lost=False,
        irep=0,
    )

    assert lost == case_data.expected_lost, f"Lost status mismatch for {case_data.current}A"

    assert hom_voltage_gain == pytest.approx(
        case_data.expected_hom_voltage_gain, abs=case_data.hom_voltage_tol
    ), f"HOM Voltage Gain mismatch for {case_data.current}A"

    assert growth_rate == pytest.approx(case_data.expected_growth_rate, abs=case_data.growth_rate_tol), (
        f"Growth Rate mismatch for {case_data.current}A"
    )


if __name__ == "__main__":
    sys.exit(pytest.main(["-v", __file__]))
