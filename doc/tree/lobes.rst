.. _lobes:

LOBEs
=====

EveryBeam provides experimental support for the LOBEs beam model as an alternative to the Hamaker model for LOFAR observations.
Currently, LOBEs is supported for the following LOFAR LBA stations:

- The superterp LBA core stations (CS002LBA - CS007LBA)
- Core station CS302LBA (at low frequency resolution)
- International LBA station SE607LBA

.. note::
    To download and install the LOBEs coefficient files when building EveryBeam, please set the flag ``-DDOWNLOAD_LOBES=On``,
    see also the :ref:`buildinstructions`.
    To avoid redundant downloads of the large coefficient files, the value of this flag is ``Off`` by default.
    The downloaded coefficient files are stored in ``CMAKE_INSTALL_PREFIX/share/everybeam/lobes`` upon installation, with each
    available LOFAR station having its unique coefficient file.

.. note::
    The coefficients of the LOBEs model are quite large, up to 1GB in total. As a consequence,
    the Docker image that is pushed to the image registry and used for creating the Debian package
    does not include the LOBEs coefficients.



Fitting LOBEs coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~
EveryBeam contains a script to fit (Legendre polynomial) coefficients to simulated results, see ``scripts/coeff_scripts/convert_lobes.py``.
This script either takes one HDF5 file containing station coordinates and simulated result for all available frequencies, or a list of MATLAB files with simulated results per frequency.
In the latter case, station coordinates can be read from a separate MATLAB file or extracted using `lofarantpos <https://pypi.org/project/lofarantpos/>`_).
The order of the Legendre polynomial fit can be controlled by the user.

Using the LOBEs model
~~~~~~~~~~~~~~~~~~~~~
For a python demo on how to use the LOBEs model, please consult :ref:`lofardemolobes`.

.. note::
    EveryBeam falls back to the Hamaker model after spawning a warning for stations that have no
    corresponding LOBEs coefficient file. The search path for the LOBEs coefficient files
    is ``CMAKE_INSTALL_PREFIX/share/everybeam/lobes`` by default, but can be overridden via the
    ``Options.coeff_path`` field.