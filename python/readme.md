Python wrappers
===============

This directory contains the python bindings for everybeam, along with a couple of utility functions for the LOBES model response. In order to compile the python bindings make sure that

    -DBUILD_WITH_PYTHON=On

Once installed, make sure that the installed `everybeam.cpython` shared object is on your `PYTHONPATH`. Check if everything works by running

    python3 -c "import everybeam"

from the commandline.