"""
Jobflow-decorated indentation workflow functions.

These functions import core logic from the parent workflows module and
apply jobflow decorators.
"""

from jobflow import job
from .utils import ase_pmg_bridge
from .. import indentation as _core


@job
@ase_pmg_bridge
def indentation_test(structure, pair_style, pair_coeff, **kwargs):
    """Perform nanoindentation simulation (jobflow-decorated)."""
    return _core.indentation_test(structure, pair_style, pair_coeff, **kwargs)


@job
@ase_pmg_bridge
def read_final_structure(structure):
    """Read final structure from indentation (jobflow-decorated)."""
    return _core.read_final_structure(structure)


@job
def plot_force_depth(results):
    """Plot force-depth curve (jobflow-decorated)."""
    return _core.plot_force_depth(results)


@job
def plot_temperature(results):
    """Plot temperature evolution (jobflow-decorated)."""
    return _core.plot_temperature(results)


@job
@ase_pmg_bridge
def plot_centrosymmetry(final_struct, n_neighbors=8):
    """Plot centrosymmetry parameter (jobflow-decorated)."""
    return _core.plot_centrosymmetry(final_struct, n_neighbors)
