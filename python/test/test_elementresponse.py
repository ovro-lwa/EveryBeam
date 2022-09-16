# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

from everybeam import ElementResponse
from everybeam import ElementResponseModel
import pytest
import numpy as np


@pytest.mark.parametrize("use_string_names", [False, True])
def test_model_types(use_string_names):
    """
    Create differente ElementResponse types and check that the model is correct.
    """
    if use_string_names:
        hamaker_hba = ElementResponse.create("hamaker", "HBA")
        hamaker_lba = ElementResponse.create("hamaker", "LBA")
        lobes = ElementResponse.create("lobes", "CS002LBA")
        lobes_not_found = ElementResponse.create("lobes", "NOT_FOUND_LBA")
        oskar_dipole = ElementResponse.create("oskardipole")
        oskar_wave = ElementResponse.create("oskarsphericalwave")
    else:
        hamaker_hba = ElementResponse.create(
            ElementResponseModel.hamaker, "HBA"
        )
        hamaker_lba = ElementResponse.create(
            ElementResponseModel.hamaker_lba, "LBA"
        )
        lobes = ElementResponse.create(ElementResponseModel.lobes, "CS002LBA")
        lobes_not_found = ElementResponse.create(
            ElementResponseModel.lobes, "NOT_FOUND_LBA"
        )
        oskar_dipole = ElementResponse.create(
            ElementResponseModel.oskar_dipole
        )
        oskar_wave = ElementResponse.create(
            ElementResponseModel.skala40_spherical
        )

    assert hamaker_hba.model == ElementResponseModel.hamaker
    assert hamaker_lba.model == ElementResponseModel.hamaker
    assert lobes.model == ElementResponseModel.lobes
    assert lobes_not_found.model == ElementResponseModel.hamaker
    assert oskar_dipole.model == ElementResponseModel.oskar_dipole
    assert oskar_wave.model == ElementResponseModel.skala40_spherical


def test_response():
    """
    Test that response() supports both single values and arrays, and that the
    returned values are equal. Also test that the values with and without an
    element_id argument are equal.
    """
    element_response = ElementResponse.create("hamaker", "LBA")
    response_10 = element_response.response(10e6, 0.0, 0.0)
    response_11 = element_response.response(11e6, 0.1, 0.1)
    response_12 = element_response.response(12e6, 0.2, 0.1)
    response_10_array = element_response.response([10e6], [0.0], [0.0])
    response_all = element_response.response(
        [10e6, 11e6, 12e6], [0.0, 0.1, 0.2], [0.0, 0.1, 0.1]
    )

    response_10_id = element_response.response(0, 10e6, 0.0, 0.0)
    response_11_id = element_response.response(0, 11e6, 0.1, 0.1)
    response_12_id = element_response.response(0, 12e6, 0.2, 0.1)
    response_10_array_id = element_response.response([0], [10e6], [0.0], [0.0])
    response_all_id = element_response.response(
        [0, 0, 0], [10e6, 11e6, 12e6], [0.0, 0.1, 0.2], [0.0, 0.1, 0.1]
    )

    assert response_10.shape == (2, 2)
    assert response_11.shape == (2, 2)
    assert response_12.shape == (2, 2)
    assert response_10_array.shape == (1, 2, 2)
    assert response_all.shape == (3, 2, 2)

    assert np.array_equal(response_10_array[0], response_10)
    assert np.array_equal(response_all[0], response_10)
    assert np.array_equal(response_all[1], response_11)
    assert np.array_equal(response_all[2], response_12)

    assert np.array_equal(response_10, response_10_id)
    assert np.array_equal(response_11, response_11_id)
    assert np.array_equal(response_12, response_12_id)
    assert np.array_equal(response_10_array, response_10_array_id)
    assert np.array_equal(response_all, response_all_id)
