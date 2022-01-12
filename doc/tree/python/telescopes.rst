Telescopes
==========

The EveryBeam python bindings support a number of Telescopes. Internally,
these inherit from the generic :code:`Telescope` class. This generic structure
should allow to define your own telescope by inheriting from the :code:`Telescope`
base class.

Classes
~~~~~~~
.. autoclass:: everybeam.Telescope
   :members:
   :inherited-members:

.. autoclass:: everybeam::PhasedArray
   :members:

.. autoclass:: everybeam.LOFAR
   :members:

.. autoclass:: everybeam.OSKAR
   :members:

.. autoclass:: everybeam.SkaMid
   :members: