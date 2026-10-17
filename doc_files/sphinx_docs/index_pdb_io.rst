.. sasmol documentation master file, created by
   sphinx-quickstart on Sun Aug 14 14:13:00 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PDB IO
======

:mod:`PDB IO`
-------------

.. automodule:: sasmol.pdb_io
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:

Unavailable psfgen atom serials
-------------------------------

An exact ``*****`` ATOM/HETATM serial is accepted as unavailable. ``index``
remains sequential and integer-valued. For that entry only, ``original_index``
uses the one-based atom position, not a recovered source serial. Readable
source serials remain unchanged. A substituted position colliding with a
preserved serial raises ``ValueError``. In connectivity-reading mode, a CONECT
endpoint depending on a substituted serial also raises ``ValueError``; no
bond identity is inferred. Other malformed serials, residue IDs and explicit
VMD-hex parsing retain their existing rules. Copy/subset operations preserve
the integer descriptor. The PDB writer is unchanged.
