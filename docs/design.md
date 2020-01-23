# Design


## Objectives:

* Make Beam library instrument agnostic

* Allow extensions to implement their own instrument specific antenna models

* Allow extensions to implement their own instrument specific beamformers: discrete time-delay beamformer, spatial nulling beamformer



## Strategy to meet Objectives

* Split off all LOFAR specific things from core library

* The core library implements `Stations`, `Beamformers` and `Elements`, independent of the instrument

* The core library can use instrument specific implementations of `Beamformers` and `Elements`



LOFAR specific is LofarMetaDataUtil which reads metadata from a LOFAR MeasurementSet.
Based on that thata it composes Station consisting of Beamformers and Elements.
The inputs of a Beamformer are elements at a given positions.
A Beamformer is also an Element. This allows composition of multilevel beamformers

For example a HBA tile is a Beamformer with HBA elements as inputs.
A HBA field if a beamformer with multiple tiles as inputs.
An HBA station is either a single HBA field (HBA0 or HBA1) or a beamformer with fields HBA0 and HBA1 as inputs.


Markdown | Less | Pretty
--- | --- | ---
*Still* | `renders` | **nicely**
1 | 2 | 3








<!--```python
s = "Python syntax highlighting"
print s
```-->

