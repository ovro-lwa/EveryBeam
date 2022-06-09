.. _buildinstructions:

Build instructions
==================

Dependencies
~~~~~~~~~~~~
EveryBeam needs a number of dependencies in order to successfully compile. On a clean (ubuntu 20.04) system,
the dependencies can be installed with (see also the ``docker`` directory):

General packages:

::

    apt-get -y install wget git make cmake g++ doxygen \
    libboost-all-dev libhdf5-dev libfftw3-dev \
    libblas-dev liblapack-dev libgsl-dev libxml2-dev \
    libgtkmm-3.0-dev libpython3-dev python3-distutils

Astronomy-specific packages:

::

    apt-get -y install casacore-dev libcfitsio-dev wcslib-dev

In order to be able to build the documentation with ``make doc``, ``sphinx`` and some other documentation tools need to be installed:

::

    pip3 install sphinx sphinx_rtd_theme breathe myst-parser




Quick installation guide
~~~~~~~~~~~~~~~~~~~~~~~~

::

    git clone --recursive -j4 https://git.astron.nl/RD/EveryBeam.git
    cd EveryBeam
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=<everybeam_install_path> ..
    make install


Installation options
~~~~~~~~~~~~~~~~~~~~

(Use :code:`ccmake` or :code:`cmake -i` to configure all options.)

* :code:`BUILD_WITH_PYTHON`: build Python module 'everybeam' to use everybeam from Python
* :code:`BUILD_TESTING`: compile tests when building EveryBeam
* :code:`DOWNLOAD_LOBES`: download and install available LOBEs coefficients files (``OFF`` by default)

All other build options serve development purposes only, and can/should be left at the default values by a regular user.

All libraries are installed in :code:`<installpath>/lib`. The header files in
:code:`<installpath>/include`. The Python module in
:code:`<installpath>/lib/python{VERSION_MAJOR}.{VERSION_MINOR}/site-packages`. Make sure that your
:code:`LD_LIBRARY_PATH` and :code:`PYTHONPATH` are set as appropiate.
Data files, such as coefficient files for the Hamaker model, the SKALA4.0 model (OSKAR), and the LOBEs model in case ``DOWNLOAD_LOBES=On`` are
installed in ``<installpath>/share/everybeam``.
