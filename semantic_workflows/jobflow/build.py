"""
Jobflow-decorated structure building workflow functions.

These functions import core logic from the parent workflows module and
apply jobflow decorators.
"""

from jobflow import job
from .utils import ase_pmg_bridge
from .. import build as _core


@job
@ase_pmg_bridge
def bulk(element, **kwargs):
    """Create bulk crystal structure (jobflow-decorated)."""
    # Remove cdict parameter if present (not supported in jobflow)
    kwargs.pop("cdict", None)
    return _core.bulk(element, **kwargs)


@job
@ase_pmg_bridge
def repeat(structure, rep, **kwargs):
    """Repeat structure (jobflow-decorated)."""
    kwargs.pop("cdict", None)
    return _core.repeat(structure, rep, **kwargs)


@job
@ase_pmg_bridge
def polycrystal(element, **kwargs):
    """Create polycrystal structure (jobflow-decorated)."""
    kwargs.pop("cdict", None)
    return _core.polycrystal(element, **kwargs)
