"""
Jobflow-decorated EV curves workflow functions.

These functions import core logic from the parent workflows module and
apply jobflow decorators.
"""

from jobflow import job
from .utils import ase_pmg_bridge
from .. import evcurves as _core


@job
@ase_pmg_bridge
def calculate_ev_curves(structure, pair_style, pair_coeff, **kwargs):
    """Calculate energy-volume curves (jobflow-decorated)."""
    # Remove cdict-related parameters (not supported in jobflow)
    kwargs.pop("cdict", None)
    kwargs.pop("potential_type", None)
    kwargs.pop("potential_doi", None)
    return _core.calculate_ev_curves(structure, pair_style, pair_coeff, **kwargs)


@job
@ase_pmg_bridge
def relax_structure(structure, pair_style, pair_coeff, **kwargs):
    """Relax structure (jobflow-decorated)."""
    kwargs.pop("cdict", None)
    return _core.relax_structure(structure, pair_style, pair_coeff, **kwargs)


# Re-export utility functions (no decoration needed)
scale_atoms = _core.scale_atoms
fit_bm = _core.fit_bm
birch_murnaghan_eval = _core.birch_murnaghan_eval
calculate_energy = _core.calculate_energy
