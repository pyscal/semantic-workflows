"""
Pyiron wrapper for indentation workflow.
Imports core functions and decorates with @as_function_node.
"""

from pyiron_workflow import as_function_node
from workflows.indentation import (
    indentation_test as _indentation_test,
    read_final_structure as _read_final_structure,
    plot_force_depth as _plot_force_depth,
    plot_temperature as _plot_temperature,
    plot_centrosymmetry as _plot_centrosymmetry,
)

# Wrap all functions with @as_function_node decorator
indentation_test = as_function_node(_indentation_test)
read_final_structure = as_function_node(_read_final_structure)
plot_force_depth = as_function_node(_plot_force_depth)
plot_temperature = as_function_node(_plot_temperature)
plot_centrosymmetry = as_function_node(_plot_centrosymmetry)
