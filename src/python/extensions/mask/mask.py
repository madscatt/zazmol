"""Compatibility wrapper for the compiled mask extension."""

try:
    from . import _mask
except ImportError:
    import _mask


def get_mask_array(farray, nname, resid, flexible_residues, nresidues, mtype):
    return _mask.get_mask_array(
        farray, nname, resid, flexible_residues, nresidues, mtype)
