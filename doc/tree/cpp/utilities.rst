.. _utilities:

utilities
=========

Utilities to read a measurement set into a :code:`Telescope` object,
set telescope related options or specify settings for the A-Term calculations.


.. doxygenenum:: everybeam::TelescopeType

.. doxygenfunction:: everybeam::GetTelescopeType

.. doxygenfunction:: everybeam::Load(const casacore::MeasurementSet &ms, const Options &options)

.. doxygenfunction:: everybeam::Load(const std::string &ms_name, const Options &options)

.. doxygenfunction:: everybeam::GetElementResponseEnum

.. doxygenstruct:: everybeam::Options

.. doxygenstruct:: everybeam::ATermSettings