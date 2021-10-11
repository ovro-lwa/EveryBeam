.. _beamformers:

Beamformers and elements
========================

EveryBeam is designed such that it can handle an arbitrary number of nested beamformers within each station. The smallest entity in the
EveryBeam nomenclature is the :code:`Element`. Beamformers as well as the element inherit from the :code:`Antenna` base class.

For efficiency reasons, certain telescopes (in conjunction with their response models) have dedicated implementations for the
beamformer. This includes for example the :code:`BeamFormerIdenticalAntennas`, :code:`BeamFormerLofarLBA`, :code:`BeamFormerLofarHBA` classes.

.. doxygenclass:: everybeam::Antenna
    :members:
    :undoc-members:

.. doxygenclass:: everybeam::BeamFormer
    :members:
    :undoc-members:

.. doxygenclass:: everybeam::BeamFormerIdenticalAntennas
    :members:
    :undoc-members:

.. doxygenclass:: everybeam::BeamFormerLofarLBA
    :members:
    :undoc-members:

.. doxygenclass:: everybeam::BeamFormerLofarHBA
    :members:
    :undoc-members: