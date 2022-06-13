# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import everybeam as eb
import pytest
import os
from subprocess import check_call
import numpy as np
import warnings


EVERYBEAM_BASE_URL = "https://support.astron.nl/software/ci_data/EveryBeam/"
DATADIR = os.environ["DATA_DIR"]
RASCIL_FITS = "test_primary_beam_RADEC_MID_512.fits"


@pytest.fixture(autouse=True)
def download_rascil_fits():
    if not os.path.isfile(os.path.join(DATADIR, RASCIL_FITS)):
        print(os.getcwd())
        wget = f"wget -q {os.path.join(EVERYBEAM_BASE_URL, RASCIL_FITS)} -O {os.path.join(DATADIR, RASCIL_FITS)}"
        check_call(wget.split())


@pytest.fixture
def skamid_setup():
    gs = eb.GridSettings()
    gs.width = gs.height = 512
    gs.ra = 0.0
    gs.dec = -0.785398
    gs.dl = gs.dm = 8.0 / gs.width * np.pi / 180.0
    gs.l_shift = gs.m_shift = 0.0

    return {
        "filename": "SKA_MID_MOCK.ms",
        # Time not relevant for the analytical SKA-MID response model
        "time": 0.0,
        "freq": 1.3e9,
        "coordinate_system": gs,
        "station_id": 0,
        "response_model": "skamid_analytical",
    }


@pytest.mark.parametrize("settings", [pytest.lazy_fixture("skamid_setup")])
def test_skamid_analytical(settings):
    """
    Compare SKA-mid (analytical) gridded response against results obtained with
    RASCIL
    """

    try:
        from astropy.io import fits
    except:
        warnings.warn(
            UserWarning("Could not import astropy, so fits image checks are skipped.")
        )
        return

    ms_path = os.path.join(DATADIR, settings["filename"])
    rascil_fits = os.path.join(DATADIR, RASCIL_FITS)

    telescope = eb.load_telescope(
        ms_path, element_response_model=settings["response_model"]
    )
    assert type(telescope) is eb.SkaMid

    grid_response = telescope.gridded_response(
        settings["coordinate_system"],
        settings["time"],
        settings["freq"],
        settings["station_id"],
    )

    rascil_response = fits.open(rascil_fits)[0].data
    # In rascil, the voltage pattern is multiplied with its complex conjugate - even though the voltage
    # pattern is real-valued.
    dimage = (
        grid_response[..., 0, 0].flatten() * np.conj(grid_response[..., 0, 0].flatten())
        - rascil_response.flatten()
    )
    rms = np.sqrt(dimage.dot(dimage) / dimage.size)
    assert rms < 1e-5
