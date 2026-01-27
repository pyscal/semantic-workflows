"""
Jobflow-decorated elastic properties workflow functions.

These functions import core logic from the parent workflows module and
apply jobflow decorators.
"""

from jobflow import job
from .utils import ase_pmg_bridge
from .. import elastic as _core


@job
@ase_pmg_bridge
def calculate_elastic_constants(structure, pair_style, pair_coeff, **kwargs):
    """Calculate elastic constants (jobflow-decorated)."""
    # Remove cdict-related parameters (not supported in jobflow)
    kwargs.pop('cdict', None)
    kwargs.pop('potential_type', None)
    kwargs.pop('potential_doi', None)
    return _core.calculate_elastic_constants(structure, pair_style, pair_coeff, **kwargs)


@job
@ase_pmg_bridge
def mechanical_response_test(structure, pair_style, pair_coeff, **kwargs):
    """Perform mechanical response test (jobflow-decorated)."""
    kwargs.pop('cdict', None)
    kwargs.pop('potential_type', None)
    kwargs.pop('potential_doi', None)
    return _core.mechanical_response_test(structure, pair_style, pair_coeff, **kwargs)
