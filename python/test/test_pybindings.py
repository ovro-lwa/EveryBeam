# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

from everybeam import load_telescope, LOFAR, CoordinateSystem, Options
import pytest
import os
import numpy as np
from subprocess import check_call, call

DATADIR = os.environ["DATA_DIR"]


@pytest.fixture(scope="session", autouse=True)
def download_msets():
    """
    Download LBA and HBA measurement sets if they are not yet found.
    Fixture is run only once at start of session.
    """
    download_scripts = ["download_lofar_lba_ms.sh", "download_lofar_hba_ms.sh"]
    for script in download_scripts:
        call(f"sh {os.path.join(os.environ['SCRIPTS_DIR'], script)}", shell=True)


@pytest.fixture
def lba_setup():
    return {
        "filename": "LOFAR_LBA_MOCK.ms",
        "time": 4.92183348e09,
        # Frequency of channel 4
        "freq": 57884216.30859375,
        # Reference cpp solution for station 31, see cpp/test/tlofar_lba.test_hamaker
        "station_id": 31,
        "freq0": 5.78125e07,
        "direction": np.array([0.667806, -0.0770635, 0.740335]),
        "preapplied_beam_dir": np.array([0.655743, -0.0670973, 0.751996]),
        "station0": np.array([0.655743, -0.0670973, 0.751996]),
        "tile0": np.array([0.655743, -0.0670973, 0.751996]),
        "cpp_response": np.array(
            [
                [-0.71383788 + 0.00612506j, -0.4903527 + 0.00171652j],
                [-0.502122 + 0.00821683j, 0.7184408 - 0.00821723j],
            ]
        ),
    }


@pytest.fixture
def hba_setup():
    return {
        "filename": "LOFAR_HBA_MOCK.ms",
        "time": 4929192878.008341,
        "freq": 138476562.5,
        # Reference cpp solution for station 31, see cpp/test/tlofar_lba.test_hamaker
        "station_id": 63,
        "freq0": 138476562.5,
        "direction": np.array([0.42458804, 0.46299569, 0.77804112]),
        "preapplied_beam_dir": np.array([0.408326, 0.527345, 0.745102]),
        "station0": np.array([0.4083262, 0.52734471, 0.74510219]),
        "tile0": np.array([0.40832685, 0.52734421, 0.74510219]),
        "cpp_response": np.array(
            [
                [0.03259424 - 0.00023045994j, 0.12204104097 - 0.00091857865j],
                [0.13063535 - 0.0010039175j, -0.02934845 + 0.00023882818j],
            ]
        ),
    }


def check_consistency_station_response(telescope, time):
    """
    Check the internal consistency of various routes for computing the station
    response

    Parameters
    ----------
    telescope : everybeam.Telescope
    time : float
        Time MJD in s
    """
    # Checking internal consistency of the different methods
    stations_all = telescope.station_response(time)
    for station_idx in range(stations_all.shape[0]):
        channel_all = telescope.station_response(time, station_idx)
        np.testing.assert_allclose(stations_all[station_idx, :, :, :], channel_all)
        for channel_idx in range(channel_all.shape[0]):
            channel_single = telescope.station_response(time, station_idx, channel_idx)
            np.testing.assert_allclose(channel_all[channel_idx, :, :], channel_single)
            # Check against implementation where frequency is given as input
            frequency = telescope.channel_frequency(channel_idx)
            freq_single = telescope.station_response(time, station_idx, frequency)
            np.testing.assert_allclose(freq_single, channel_single)


def check_reference_solution(
    telescope, time, station_id, freq, direction, station0, tile0, reference_solution
):
    """
    Check computed response against provided reference solution.
    """
    response = telescope.station_response(
        time, station_id, freq, direction, station0, tile0,
    )
    np.testing.assert_allclose(response, reference_solution, atol=1e-6)

    # Check that the same response is obtained by multiplying the array factor with the
    # element response
    array_factor = telescope.array_factor(
        time, station_id, freq, direction, station0, tile0
    )
    element_response = telescope.element_response(time, station_id, freq, direction)
    np.testing.assert_allclose(
        np.matmul(array_factor, element_response), response, atol=1e-6
    )


@pytest.mark.parametrize(
    "ref, differential_beam",
    [
        (pytest.lazy_fixture("lba_setup"), False),
        (pytest.lazy_fixture("lba_setup"), True),
        (pytest.lazy_fixture("hba_setup"), False),
        (pytest.lazy_fixture("lba_setup"), True),
    ],
)
def test_lofar(ref, differential_beam):
    ms_path = os.path.join(DATADIR, ref["filename"])
    telescope = load_telescope(ms_path, use_differential_beam=differential_beam)

    assert isinstance(telescope, LOFAR)

    time = ref["time"]
    freq = ref["freq"]
    direction = ref["preapplied_beam_dir"] if differential_beam else ref["direction"]
    ref_solution = (
        np.eye(2, dtype=np.complex64) if differential_beam else ref["cpp_response"]
    )

    check_consistency_station_response(telescope, time)
    check_reference_solution(
        telescope,
        time,
        ref["station_id"],
        freq,
        direction,
        ref["station0"],
        ref["tile0"],
        ref_solution,
    )

    # Array factor in pointing direction should be unity
    array_factor_I = telescope.array_factor(
        time, ref["station_id"], freq, ref["station0"], ref["station0"], ref["station0"]
    )
    np.testing.assert_allclose(array_factor_I, np.eye(2, dtype=np.complex64), rtol=1e-6)

    # For equal station and tile direction, check that the same response is
    # obtained via two routes:
    # Response by explicitly specifying the direction
    response_1 = telescope.station_response(
        time,
        ref["station_id"],
        freq,
        ref["direction"],
        ref["station0"],
        ref["station0"],
    )

    # tile0 direction is implicit in station0 direction
    response_2 = telescope.station_response(
        time, ref["station_id"], freq, ref["direction"], ref["station0"],
    )
    np.testing.assert_allclose(response_1, response_2)

    # Check that we can reproduce this with array_factor * beam_response
    # Check that the same response is obtained by multiplying the array factor with the
    # element response
    array_factor = telescope.array_factor(
        time, ref["station_id"], freq, ref["direction"], ref["station0"]
    )
    element_response = telescope.element_response(
        time, ref["station_id"], freq, ref["direction"]
    )
    np.testing.assert_allclose(
        np.matmul(array_factor, element_response), response_1, rtol=1e-6
    )
