.. _lofardemogriddedresponse:

LOFAR telescope: computing the gridded response
===============================================

The EveryBeam Python bindings can be used to evalute the response on a "grid" (e.g. an
image). This demo illustrates in a step-by-step manner a minimal example of doing this. A slightly more involved
(and graphically appealing) demo can be found in  :ref:`lofardemolobes`. The complete code is shown at the bottom of this page.

It all starts with importing the ``everybeam`` python module, along with any other
libraries that are needed for your project and specifying a path to the Measurement Set
that will be used:

::

    import everybeam as eb
    import numpy as np

    # Set path to LOFAR LBA MS and load telescope
    ms_path = "../../test_data/LOFAR_LBA_MOCK.ms"

The telescope can now be loaded with ``eb.load_telescope``. This should return an instance of a  ``LOFAR`` telescope.

::

    telescope = eb.load_telescope(ms_path)
    assert type(telescope) == eb.LOFAR

Please note that since no optional arguments were passed to ``load_telescope``, no differential beam will applied and the
default element response model (``"Hamaker"``) will be used.

In order to evaluate the beam response on a grid, we need to pass some user input, this includes a time, frequency and station index,

::

    time = 4.92183348e09
    freq = 57884216.30859375
    station_id = 23

as well as the settings for a grid. This can be done via an ``eb.GridSettings`` object:

::

    # Specify the settings of the grid (/image)
    gs = eb.GridSettings()
    gs.width = gs.height = 4
    gs.ra = -1.44194878
    gs.dec = 0.85078091
    gs.dl = gs.dm = 0.5 * np.pi / 180.0
    gs.l_shift = gs.m_shift = 0.0

So we defined a grid with 4x4 pixels, with its center located at ``(-1.44194878, 0.85078091)``.

We're now all set to compute the gridded response. This either can be done for all stations at once:

::

    # Get the gridded response for all stations at once
    grid_response_all = telescope.gridded_response(gs, time, freq)

Returning a 5-dimensional numpy array of shape ``(nstations, width, height, 2, 2)``, i.e. a Jones matrix
for each pixel in the grid, for each station. This can be easily confirmed with:

::

    # Check whether the returned numpy array has the expected shape
    assert grid_response_all.shape == (telescope.nr_stations, gs.width, gs.height, 2, 2)

We also could request the gridded beam response for a specific station ``station_id``:

::

    grid_response = telescope.gridded_response(gs, time, freq, station_id)

Of course, the returned response should match the corresponding entry in the array for all
gridded station responses:

::

    np.testing.assert_allclose(grid_response, grid_response_all[station_id, ...], rtol=1e-6)

In addition, we can check whether the gridded response at the center pixel - pixel (2, 2) in our example -
matches the station response at the facet center:

::

    station_response = telescope.station_response(time, station_id, freq, gs.ra, gs.dec)
    np.testing.assert_allclose(grid_response[2, 2, ...], station_response, rtol=1e-6)


A complete overview of the code is shown below:

.. literalinclude:: ../../../demo/python/lofargridresponse.py