"""
Jobflow-decorated point defects workflow functions.

These functions import core logic from the parent workflows module and
apply jobflow decorators.
"""

from jobflow import job
from .utils import ase_pmg_bridge
from .. import pointdefects as _core


@job
@ase_pmg_bridge
def create_interstitial(atoms, element, void_type, **kwargs):
    """Create interstitial defect (jobflow-decorated)."""
    # Remove KG parameter (not supported in jobflow)
    kwargs.pop("kg", None)
    return _core.create_interstitial(atoms, element, void_type, **kwargs)


@job
@ase_pmg_bridge
def create_substitutional(atoms, element, **kwargs):
    """Create substitutional defect (jobflow-decorated)."""
    # Remove KG parameter (not supported in jobflow)
    kwargs.pop("kg", None)
    return _core.create_substitutional(atoms, element, **kwargs)
