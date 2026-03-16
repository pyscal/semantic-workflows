"""
Pyiron-decorated elastic properties workflow functions.

These functions import core logic from the parent workflows module and
apply pyiron_workflow decorators.
"""

from pyiron_workflow import as_function_node
from .. import elastic as _core


@as_function_node
def calculate_elastic_constants(structure, pair_style, pair_coeff, **kwargs):
    """Calculate elastic constants (pyiron-decorated)."""
    return _core.calculate_elastic_constants(structure, pair_style, pair_coeff, **kwargs)


@as_function_node
def mechanical_response_test(structure, pair_style, pair_coeff, **kwargs):
    """Perform mechanical response test (pyiron-decorated)."""
    return _core.mechanical_response_test(structure, pair_style, pair_coeff, **kwargs)
