.. _cppinterface:

C++ interface
=============

Below the most important parts of the EveryBeam C++ API are documented. The :ref:`utilities` section
documents the tools to load a telescope from an input MeasurementSet and specify telescope or element
response model options. The section :ref:`everybeam::telescope` documents the actual telescope classes, with the :code:`Telescope` class being the mother of all telescopes.
A :code:`GriddedResponse` or :code:`PointResponse` object can be requested from the telescopes, to compute the beam response on a prescribed grid (i.e. image) or for
a prescribed direction, respectively. Documentation for the :code:`GriddedResponse` classes can be found in :ref:`everybeam::griddedresponse`, documentation for the
:code:`PointResponse` classes can be found in :ref:`everybeam::pointresponse`.

For advanced users that want to know more about the beam forming procedure, reference is made to :ref:`beamformers`

.. toctree::
   :maxdepth: 1

   cpp/utilities
   cpp/telescope
   cpp/pointresponse
   cpp/griddedresponse
   cpp/aterms
   cpp/beamformer
