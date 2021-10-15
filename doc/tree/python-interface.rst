Python interface
================

EveryBeam provides an interface via the :code:`everybeam` module. In order to build
this module, make sure the compile option :code:`BUILD_WITH_PYTHON` is turned :code:`ON`.

After successfully compiling and installing the python bindings, you should update your :code:`LD_LIBRARY_PATH`
and your :code:`PYTHONPATH` as follows:

::

    export LD_LIBRARY_PATH=<installpath>/lib/:$LD_LIBRARY_PATH
    export PYTHONPATH=<installpath>/lib/python<VERSION_MAJOR>.<VERSION_MINOR>/site-packages:$PYTHONPATH

where :code:`VERSION_MAJOR` and :code:`VERSION_MINOR` depend upon the specific version of python on your system.
The :code:`everybeam` module can now be imported in python with:

::

    import everybeam

For further examples on how the :code:`everybeam` module can be used, reference is made to the :ref:`pythondemos`.

.. toctree::
   :maxdepth: 1

   python/utils
   python/telescopes