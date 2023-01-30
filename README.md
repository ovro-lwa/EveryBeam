# EveryBeam library

This package can be used to compute the beam response for a variety of
radio telescopes, i.e.:

* LOFAR
* SKA/OSKAR
* ATCA
* GMRT
* VLA
* MWA

This package also provides an abstract interface to a selection of beam responses for apperture arrays (LOFAR/OSKAR), and beamformed versions thereof. Currently implemented are:

 * Hamaker LOFAR model
 * OSKAR spherical wave model
 * OSKAR-dipole: work in progress
 * LOBEs: work in progress. A coefficient file is currently only available for a limited number of LOFAR stations. Selecting the LOBEs model defaults back to Hamaker, in case no coefficient file is available.

EveryBeam replaces the stand alone version of the LOFAR station response library (LOFARBeam).

EveryBeam is licensed under the terms of the GNU GPL3 license.

## Documentation and Installation Instructions

[Documentation](https://everybeam.readthedocs.io) along with [installation instructions](https://everybeam.readthedocs.io/en/latest/build-instructions.html) can be found at the provided links.

## Usage with DP3

To use Everybeam within [DP3](https://git.astron.nl/RD/DP3) - the streaming visibility framework - DP3 needs to be compiled against EveryBeam. To do so, make sure DP3 can find EveryBeam by adding the EveryBeam install dir to the `CMAKE_PREFIX_PATH`.

A test measurement set is included in DP3 (`tNDP3-generic.in_MS.tgz`).

To simulate visibilities with a certain element model, use `DP3 DP3.parset` with `DP3.parset` a parset file with the following contents:

    msin=tNDP3-generic.MS
    msout=.
    steps=[predict]
    predict.usebeammodel=True
    predict.elementmodel=oskardipole
    predict.sourcedb=tNDP3-generic.MS/sky  # sourcedb file

## Usage with WSClean

To use EveryBeam with [WSClean](https://gitlab.com/aroffringa/wsclean) (for A-term or primary beam corrections), WSClean needs to be compiled against EveryBeam. In order to do so, make sure WSClean can find EveryBeam by adding the EveryBeam install dir to the `CMAKE_PREFIX_PATH`.
