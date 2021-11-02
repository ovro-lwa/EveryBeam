.. _oskardemopointresponse:

OSKAR telescope: point response
===============================

This file explains in a step-by-step manner how to compute the beam response
for an OSKAR station. The complete code is shown at the bottom of this page.

It all starts with importing the ``everybeam`` python module, along with any other
libraries that are needed for your project.

::

    import everybeam as eb
    import numpy as np

In order to evaluate the beam response, we need to specify which Measurement Set will
be used, whether or not to apply the differential beam and which element response model
should be used. The element response should be eihter one of ``["hamaker", "lobes", "oskar_dipole", "skala40_wave"]``.
Even more important, the element response model needs to be chosen such that it is consistent with the provided telescope.
So for an OSKAR Measurement Set, the element response model should be either ``"oskar_dipole"`` or ``"skala40_wave"``.
In this demo we opt for ``"skala40_wave"`` - which is the element response model that is currently used in OSKAR simulations:

::

    # Set path to (OSKAR) MS
    ms_path = "../../test_data/OSKAR_MOCK.ms"

    # Use differential beam?
    use_differential_beam = False

    # Set element response model
    element_response_model = "skala40_wave"

Given the Measurement Set and the options, we can now load the telescope.
Given your Measurement Set, EveryBeam sorts out which telescope should be returned.
In this example we obviously expect to get an OSKAR telescope back.

::

    # Load the telescope
    telescope = eb.load_telescope(
        ms_path,
        use_differential_beam=use_differential_beam,
        element_response_model=element_response_model,
    )

    # Is this an OSKAR telescope?
    assert type(telescope) == eb.OSKAR

In order to evaluate the beam, we need to set some additional (observation related)
variables, such as the time, the station index, the frequency, and the point of interest (in ITRF-coordinates)

::
    # Set properties for the beam evaluation
    time = 4.45353e09
    station_id = 0
    freq = 5.0e07

    # Point of interest (given in ITRF coordinates)
    dir_itrf = np.array([-0.203137, 0.841818, -0.500078])

Now we are in principle all set to compute the beam response. For instance, computing
the array factor can be done as:

::

    array_factor_phase_centre = telescope.array_factor(
        time, station_id, freq, dir_itrf, dir_itrf
    )

In this specific example, the point of interest and the delay direction are
chosen identical. As a result, we expect that a unity (2x2) matrix is returned:

::

    np.testing.assert_allclose(
        array_factor_phase_centre, np.eye(2, dtype=np.complex64), rtol=1e-8
    )

To make the example a bit less trivial, we adapt the delay direction so that it
doesn't coincide with the point of interest, e.g.:

::

    delay_dir_itrf = dir_itrf + np.array([-0.1, 0.5, 0.025])

We now can calculate the array factor, element response and the full response
for a non trivial delay direction and check that the full response matches the
matrix-product of the array factor and the element response:

::

    array_factor = telescope.array_factor(time, station_id, freq, dir_itrf, delay_dir_itrf)

    # Element beam for station 0
    element_response = telescope.element_response(time, station_id, freq, dir_itrf)

    # Full beam for station 0
    response = telescope.station_response(time, station_id, freq, dir_itrf, delay_dir_itrf)

    # Full beam should match product of array_factor and element_response
    np.testing.assert_allclose(
        np.matmul(array_factor, element_response), response, rtol=1e-6
    )

A complete overview of the code is shown below:

.. literalinclude:: ../../../demo/python/oskarpointresponse.py
