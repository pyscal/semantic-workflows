"""
Jobflow-decorated structure building workflow functions.

These functions import core logic from the parent workflows module and
apply jobflow decorators.
"""

from jobflow import job
from .utils import ase_pmg_bridge
from .. import build as _core

# Decorate core functions with jobflow decorators
# Note: jobflow functions don't use kg parameter
bulk = ase_pmg_bridge(job(_core.bulk))
repeat = ase_pmg_bridge(job(_core.repeat))
polycrystal = ase_pmg_bridge(job(_core.polycrystal))
