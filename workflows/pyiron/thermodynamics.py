"""
Pyiron-decorated thermodynamics workflow functions.

These functions import core logic from the parent workflows module and
apply pyiron_workflow decorators.
"""

from pyiron_workflow import as_function_node
from .. import thermodynamics as _core


@as_function_node
def calculate_thermal_properties(structure, pair_style, pair_coeff, **kwargs):
    """Calculate thermal properties (pyiron-decorated)."""
    return _core.calculate_thermal_properties(
        structure, pair_style, pair_coeff, **kwargs
    )


@as_function_node
def calculate_free_energy(structure, pair_style, pair_coeff, **kwargs):
    """Calculate free energy (pyiron-decorated)."""
    return _core.calculate_free_energy(structure, pair_style, pair_coeff, **kwargs)


# Re-export utility functions directly from core
calphy_default_input = _core.calphy_default_input
