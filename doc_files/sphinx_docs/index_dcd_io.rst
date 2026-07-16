.. sasmol documentation master file, created by
   sphinx-quickstart on Sun Aug 14 14:13:00 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DCD IO
======

Compatibility notes
-------------------

The Python DCD API preserves its existing compatibility contract:
``open_dcd_read(filename)`` returns the positional descriptor list used by
SASSIE/ZAZZIE callers, ``read_dcd_step(dcdfile, frame)`` consumes that list,
and ``close_dcd_read(dcdfile)`` continues to accept either the descriptor list
or the raw capsule.

The legacy C++ ``sasio`` wrapper also exposes a typed
``sasio::DcdReadHandle`` for streaming readers. The handle keeps the file
pointer and all header-derived state required by the low-level
``read_dcdstep()`` routine together for the lifetime of an open DCD:
``natoms``, ``nframes``, ``reverseEndian``, ``charmm``, ``istart``,
``nsavc``, ``delta``, and ``namnf``. C++ callers that need to stream frames
should prefer ``open_dcd_read_handle()``, the handle-based
``sasio::read_dcd_step(handle, frame, x, y, z)`` or
``Files::read_dcd_step(handle, frame)``, and ``close_dcd_read(handle)`` so
CHARMM/X-PLOR header state is not reconstructed or reset between open and
read.

:mod:`DCD IO`
-------------

.. automodule:: sasmol.dcd_io
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:
