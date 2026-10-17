.. sasmol documentation master file, created by
   sphinx-quickstart on Sun Aug 14 14:13:00 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Subset
======

Prepared selections
-------------------

``Mask.get_subset_mask_prepared(basis_filter)`` is an explicit alternative to
``Mask.get_subset_mask(basis_filter)``. The original method and its callers are
unchanged. The prepared method translates the expression once per call, then
uses the same molecular descriptors and atom traversal. It requires no native
compiler and caches neither expressions nor molecular results across calls.

.. code-block:: python

    import sasmol.system as system

    molecule = system.Molecule_Maker(3, name=['N', 'CA', 'C'])
    error, mask = molecule.get_subset_mask_prepared('name[i] == "CA"')
    if not error:
        indices = molecule.get_indices_from_mask(mask)

Developers own expression suitability and result handling, as with the existing
method. Existing selection messages, no-match/partial lists and successful
integer masks remain the library convention; there is no new error channel or
automatic fallback. The method prepares lazily at the first atom, so zero-atom
selections retain the existing behavior. Fresh descriptors are read on every
call. Subclasses overriding the old method must review adoption of the separate
method, which does not automatically call their override. Python frame/local
introspection and translation audit events are outside the parity claim.

API reference
-------------

:mod:`Subset`
----------------------

.. automodule:: sasmol.subset
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:
