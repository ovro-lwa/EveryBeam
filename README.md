# EveryBeam library

This package can be used to compute the beam response for a variety of
radio telescopes, i.e.:

* LOFAR
* OSKAR
* MWA
* VLA
* ATCA

This package also provides an abstract interface to a selection of beam responses for apperture arrays (LOFAR/OSKAR), and beamformed versions thereof. Currently implemented are:

 * Hamaker LOFAR model
 * OSKAR spherical wave model
 * OSKAR-dipole: work in progress
 * LOBES: work in progress. A coefficient file is currently only available for LOFAR station CS302LBA. Selecting the LOBES model defaults back to Hamaker, except for the aforementioned station.

EveryBeam replaces the stand alone version of the LOFAR station response library (LOFARBeam).

EveryBeam is licensed under the terms of the GNU GPL3 license.

## Installation

On a clean ubuntu 18.04 computer, the following packages are required (see also the `docker` directory):

General packages

    apt-get -y install wget git make cmake g++ doxygen \
    libboost-all-dev libhdf5-dev libfftw3-dev \
    libblas-dev liblapack-dev libgsl-dev libxml2-dev \
    libgtkmm-3.0-dev libpython3-dev python3-distutils

Astronomy packages:

    apt-get -y install casacore-dev libcfitsio-dev wcslib-dev

Installation is then typically done as:

    mkdir build
    cd build
    cmake -DBUILD_WITH_PYTHON=On/Off ..
    make
    make install

The `BUILD_WITH_PYTHON` option indicates whether the python bindings for EveryBeam should be generated. If switched `On` (non-default), a shared `everybeam.cpython` should appear in the `[INSTALL_DIR]/lib/python[MAJOR].[MINOR]/site-packages` directory.

## Usage with DP3

To use Everybeam within DP3 - the streaming visibility framework, https://git.astron.nl/RD/DP3 - DP3 needs to be compiled against EveryBeam. To do so, make sure DP3 can find EveryBeam by adding the EveryBeam install dir to the `CMAKE_PREFIX_PATH`.

A test measurement set is included in DP3 (`tNDP3-generic.in_MS.tgz`).

To simulate visibilities with a certain element model, use `DP3 DP3.parset` with `DP3.parset` a parset file with the following contents:

    msin=tNDP3-generic.MS
    msout=.
    steps=[predict]
    predict.usebeammodel=True
    predict.elementmodel=oskardipole
    predict.sourcedb=tNDP3-generic.MS/sky  # sourcedb file

## Usage with WSClean

To use EveryBeam with WSClean (for A-term or primary beam corrections), WSClean needs to be compiled against EveryBeam. In order to do so, make sure WSClean can find EveryBeam by adding the EveryBeam install dir to the `CMAKE_PREFIX_PATH`.
