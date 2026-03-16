"""
Pyiron-decorated structure building workflow functions.

These functions import core logic from the parent workflows module and
apply pyiron_workflow decorators.
"""

from pyiron_workflow import as_function_node
from .. import build as _core
from .templates import sample_template as template_dict

# Decorate core functions with pyiron_workflow decorator
bulk = as_function_node(_core.bulk)
repeat = as_function_node(_core.repeat)
polycrystal = as_function_node(_core.polycrystal)

# Re-export utility functions and constants
generate_id = _core.generate_id
bravais_lattice_dict = _core.bravais_lattice_dict
get_chemical_composition = _core.get_chemical_composition
get_cell_volume = _core.get_cell_volume
get_number_of_atoms = _core.get_number_of_atoms
get_simulation_cell_length = _core.get_simulation_cell_length
get_simulation_cell_vector = _core.get_simulation_cell_vector
get_simulation_cell_angle = _core.get_simulation_cell_angle
get_bravais_lattice = _core.get_bravais_lattice
get_spacegroup_symbol = _core.get_spacegroup_symbol
get_spacegroup_number = _core.get_spacegroup_number
update_attributes = _core.update_attributes
