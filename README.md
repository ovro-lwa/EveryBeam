# EveryBeam library

This package can be used to compute the beam response for a variety of
radio telescopes, i.e.:

* LOFAR
* OSKAR
* MWA
* VLA
* ATCA

This package also provides an abstract interface to a selection of beam responses for apperture arrays (LOFAR/OSKAR), and beamformed versions thereof. Currently implemented are:

 * Hamaker
 * OSKAR spherical wave model
 * OSKAR-dipole: work in progress
 * LOBES: work in progress. A coefficient file is currently only available for LOFAR station CS302LBA. Selecting the LOBES model defaults back to Hamaker, except for the aforementioned station.

EveryBeam aims to gradually replace the stand alone version of the LOFAR station response library (LOFARBeam).

EveryBeam is licensed under the terms of the [GNU GPL3 license](LICENSE). 
 
## Installation

On a clean ubuntu 18.04 computer, the following packages are required, see also the [docker file](docker/Dockerfile-base):

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

A [`Dockerfile`](docker/Dockerfile-everybeam) to compile a working version of DPPP and EveryBeam is included in the `test` directory.

## Usage with DPPP

To use with DPPP, the streaming visibility framework used for LOFAR, currently the `development` branch from that repository is required, 
see https://git.astron.nl/RD/DP3/.

A test measurementset is included in DPPP (`tNDPPP-generic.in_MS.tgz`).

To simulate visibilities with a certain element model, use `DPPP DPPP.parset` where `DPPP.parset` contains:

    msin=tNDPPP-generic.MS
    msout=.
    steps=[predict]
    predict.usebeammodel=True
    predict.elementmodel=oskardipole
    predict.sourcedb=tNDPPP-generic.MS/sky  # sourcedb file

## Usage with WSClean

To use EveryBeam with WSClean (for A-term or primary beam corrections), WSClean needs to be compiled against EveryBeam. In order to do so, make sure WSClean can find EveryBeam by adding the EveryBeam install dir to the `CMAKE_PREFIX_PATH`.

## Design

See [docs/design.md](@ref designpage) for design considerations.

## Documentation
The `doxygen` documentation for EveryBeam can be found at https://www.astron.nl/citt/EveryBeam/

## Compatibility note
This package has undergone the same branching process as [the DP3 package](https://github.com/lofar-astron/DP3): it is a continuation of the LOFAR beam library at svn.astron.nl/LOFAR, and was extracted from and branched off LOFAR Release 3.2. The LOFAR beam package will likely not be maintained in the ASTRON repository.
