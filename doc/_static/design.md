Design {#designpage}
===================

## Objectives:

 * Make Beam library instrument agnostic
 * Allow extensions to implement their own instrument specific antenna models
 * Allow extensions to implement their own instrument specific beamformers:
   discrete time-delay beamformer, spatial nulling beamformer



## Strategy to meet Objectives

 * Split off all LOFAR specific things from core library.
 * The core library implements `Stations`, `BeamFormers` and `Elements`,
   independent of the instrument.
 * The core library can use instrument specific implementations of
   `BeamFormers` and `Elements`.

LOFAR specific is [lofarreadutils](@ref lofarreadutils.h) which reads
metadata from a LOFAR MeasurementSet. Based on that data it composes Station
consisting of BeamFormers and Elements. OSKAR specific is [msv3readutils](@ref msv3readutils.h) for reading OSKAR simulated measurement sets.

The inputs of a BeamFormer are elements at a given positions.
A BeamFormer is also an Element. This allows composition of multilevel
beamformers. For example a HBA tile is a BeamFormer with HBA elements as inputs.
An HBA field is a beamformer with multiple tiles as inputs. An HBA station is
either a single HBA field (HBA0 or HBA1) or a beamformer with fields HBA0 and
HBA1 as inputs.

## Package outline

<img src="https://git.astron.nl/RD/EveryBeam/-/raw/6f8f7c3e9c2e50ad14e3b1afe538b2602460f97a/docs/everybeam_uml.png" alt="everybeam_uml" width="1200"/>

In addition, EveryBeam offers a suite of aterm calculation methods in the `aterms` namespace. This namespace rougly has the following layout:

<img src="https://git.astron.nl/RD/EveryBeam/-/raw/92f9aa123560fa2841468b2201d259b2c4b1ddd7/docs/everybeam_aterms.png" alt="everybeam_aterms" width="1200"/>

## Classes:

* [Station](@ref everybeam::Station) is a thin (but central!) wrapper around an
  Element object, with some convenience functions. 
  Constructor consumes:
    - `name`: name of station 
    - `position`: position of station
    - `model`: antenna model/element response model (Hamaker/OSKARDipole/OSKARSphericalWave/LOBES) 
  
  The all important attribute is the `antenna_`, which (somewhat confusingly) can denote either an `Element` or a `BeamFormer` pointer.
  `Station` class contains a couple of convenience specialisations for the `response` and `ArrayFactor` methods to facilitate implementation in dependencies such as `DP3` and `WSClean`.
* [Element](@ref everybeam::Element) - something for which a
  response can be computed, can be a single element (antenna) or a beamformer.
  Constructor consumes:
    - `CoordinateSystem`: coordinate system
    - `ElementResponse`: element response 
* [BeamFormer](@ref everybeam::BeamFormer) -
  Special kind of Element beamformer, has inputs that are Element objects.
  Inherits from `Antenna`, and includes `Element` (the latter for a not-immediately-obvious reason...)
  Constructor (optionally) consumes:
    - `CoordinateSystem`
    - `phase_reference_position`
* [ElementResponse](@ref everybeam::ElementResponse) - purely virtual class describing response of a single element (antenna), different models can be implemented.
  (Hamaker/OSKARDipole/OSKARSphericalWave/LOBES)
 See `hamaker`/`lobes`/`oskar` dirs for definitions of the inheriting classes
* [Antenna](@ref everybeam::ElementResponse) - Class describing a single antenna. `BeamFormer` and `Element` inherit `public` scope from `Antenna` (in particular the `response` method). 
Constructor (optionally) consumes:
    - `CoordinateSystem`
    - `phase_reference_position`
* [AtermConfig](@ref everybeam::aterms::ATermConfig) - Class to configure the aterms calculation. This class is typically interfaced from DP3 and WSClean in order to request aterm calculations from EveryBeam


