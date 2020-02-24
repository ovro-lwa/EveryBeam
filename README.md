# LOFAR beam library

Stand alone version of the LOFAR station response library.

This package provides an abstract interface to a selection of beam responses, and beamformed versions thereof. Currently implemented are:

 * Hamaker
 * OSKAR-dipole
 * LOBES (work in progress)
 * OSKAR spherical wave model

## Installation

On a clean ubuntu 18.04 computer, the following packages are required:

```
apt install build-essential casacore-dev cmake python3-distutils libblas-dev liblapack-dev libhdf5-dev
```

Installation is then typically done as:

```
mkdir build
cd build
cmake ..
make
make install
```

A `Dockerfile` to compile a working version of DPPP and LOFARBeam is included in the `test` directory.

## Use with DPPP

To use with DPPP, the streaming visibility framework used for LOFAR, currently the branch `oskar` from that repository is required.

A test measurementset is included in DPPP (`tNDPPP-generic.in_MS.tgz`).

To simulate visibilities with a certain element model, use `DPPP DPPP.parset` where `DPPP.parset` contains:

```
msin=tNDPPP-generic.MS
msout=.
steps=[predict]
predict.usebeammodel=True
predict.elementmodel=oskardipole
predict.sourcedb=tNDPPP-generic.MS/sky  # sourcedb file
```

## Compatibility note
This package has undergone the same branching process as [the DP3 package](https://github.com/lofar-astron/DP3): it is a continuation of the LOFAR beam library at svn.astron.nl/LOFAR, and was extracted from and branched off LOFAR Release 3.2. The LOFAR beam package will likely not be maintained in the ASTRON repository.
