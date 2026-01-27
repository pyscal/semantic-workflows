"""
Jobflow-decorated EV curves workflow functions.

These functions import core logic from the parent workflows module and
apply jobflow decorators.
"""

from jobflow import job
from .utils import ase_pmg_bridge
from .. import evcurves as _core

# Decorate core functions with jobflow decorators
# Note: jobflow functions don't use kg, potential_type, potential_doi parameters
calculate_ev_curves = ase_pmg_bridge(job(_core.calculate_ev_curves))
relax_structure = ase_pmg_bridge(job(_core.relax_structure))

# Re-export utility functions (no decoration needed)
scale_atoms = _core.scale_atoms
fit_bm = _core.fit_bm
birch_murnaghan_eval = _core.birch_murnaghan_eval
calculate_energy = _core.calculate_energy
