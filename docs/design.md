# Design


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

LOFAR specific is [LofarMetaDataUtil](@ref LofarMetaDataUtil.h) which reads
metadata from a LOFAR MeasurementSet. Based on that data it composes Station
consisting of BeamFormers and Elements.

The inputs of a BeamFormer are elements at a given positions.
A BeamFormer is also an Element. This allows composition of multilevel
beamformers. For example a HBA tile is a BeamFormer with HBA elements as inputs.
An HBA field is a beamformer with multiple tiles as inputs. An HBA station is
either a single HBA field (HBA0 or HBA1) or a beamformer with fields HBA0 and
HBA1 as inputs.





Classes:

* [Station](@ref everybeam::Station) is a thin wrapper around an
  Element object, with some convenience functions.
* [Element](@ref everybeam::Element) - something for which a
  response can be computed, can be a single element (antenna) or a beamformer.
* [BeamFormer](@ref everybeam::BeamFormer) -
  Special kind of Element beamformer, has inputs that are Element objects.
* [ElementResponse](@ref everybeam::ElementResponse) - response of
  a single element (antenna), different models can be implemented.


<!--Markdown | Less | Pretty
--- | --- | ---
*Still* | `renders` | **nicely**
1 | 2 | 3
-->







<!--```python
s = "Python syntax highlighting"
print s
```-->

