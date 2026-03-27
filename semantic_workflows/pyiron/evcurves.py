"""
Pyiron-decorated EV curves workflow functions.

These functions import core logic from the parent workflows module and
apply pyiron_workflow decorators.
"""

from pyiron_workflow import as_function_node
from .. import evcurves as _core

# Decorate core functions with pyiron_workflow decorator
calculate_ev_curves = as_function_node(_core.calculate_ev_curves)
relax_structure = as_function_node(_core.relax_structure)

# Re-export utility functions (no decoration needed)
scale_atoms = _core.scale_atoms
fit_bm = _core.fit_bm
birch_murnaghan_eval = _core.birch_murnaghan_eval
