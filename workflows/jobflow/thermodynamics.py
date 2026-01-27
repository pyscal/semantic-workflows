"""
Jobflow-decorated thermodynamics workflow functions.

These functions import core logic from the parent workflows module and
apply jobflow decorators.
"""

from jobflow import job
from .utils import ase_pmg_bridge
from .. import thermodynamics as _core


@job
@ase_pmg_bridge
def calculate_thermal_properties(structure, pair_style, pair_coeff, **kwargs):
    """Calculate thermal properties (jobflow-decorated)."""
    # Remove cdict-related parameters (not supported in jobflow)
    kwargs.pop("cdict", None)
    kwargs.pop("potential_type", None)
    kwargs.pop("potential_doi", None)
    return _core.calculate_thermal_properties(
        structure, pair_style, pair_coeff, **kwargs
    )


@job
@ase_pmg_bridge
def calculate_free_energy(structure, pair_style, pair_coeff, **kwargs):
    """Calculate free energy (jobflow-decorated)."""
    # Remove cdict-related parameters (not supported in jobflow)
    kwargs.pop("cdict", None)
    kwargs.pop("potential_type", None)
    kwargs.pop("potential_doi", None)
    return _core.calculate_free_energy(structure, pair_style, pair_coeff, **kwargs)


# Re-export utility functions directly from core
calphy_default_input = _core.calphy_default_input
