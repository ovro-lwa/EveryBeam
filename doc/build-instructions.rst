Build instructions
==================

Dependencies:
~~~~~~~~~~~~~

To compile EveryBeam, the following dependencies are needed:

* Boost
* Casacore
* FFTW3 (float precision)
* CFITSIO

Quick installation guide:
~~~~~~~~~~~~~~~~~~~~~~~~~

::

    git clone --recursive -j4 https://git.astron.nl/RD/EveryBeam.git
    cd EveryBeam
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=<everybeam_install_path> ..
    make install


Installation options:
~~~~~~~~~~~~~~~~~~~~~

(Use :code:`ccmake` or :code:`cmake -i` to configure all options.)

* :code:`BUILD_LOBES`: include LOBEs support and download LOBEs coefficients
* :code:`BUILD_WITH_PYTHON`: build Python module 'everybeam' to use everybeam from Python
* :code:`BUILD_TESTING`: compile tests when building EveryBeam

All other build options are for development purposes only, and should be left at the default values by a regular user.

All libraries are installed in :code:`<installpath>/lib`. The header files in
:code:`<installpath>/include`. The Python module in
:code:`<installpath>/lib/python{VERSION_MAJOR}.{VERSION_MINOR}/site-packages`. Make sure that your
:code:`LD_LIBRARY_PATH` and :code:`PYTHONPATH` are set as appropiate.