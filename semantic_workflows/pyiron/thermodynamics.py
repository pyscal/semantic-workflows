"""
Pyiron-decorated thermodynamics workflow functions.

These functions import core logic from the parent workflows module and
apply pyiron_workflow decorators.
"""

from pyiron_workflow import as_function_node
from .. import thermodynamics as _core


@as_function_node
def calculate_thermal_properties(
    structure,
    pair_style,
    pair_coeff,
    temperature,
    pressure=0,
    cores=1,
    n_equilibration_steps=15000,
    n_run_steps=25000,
    cdict=None,
    potential_type=None,
    potential_doi=None,
    autovalidate=False,
    cv_threshold=0.05,
    alpha_min=1e-6,
    alpha_max=1e-4,
    convergence_threshold=0.1,
):
    """Calculate thermal properties (pyiron-decorated)."""
    datadict = _core.calculate_thermal_properties(
        structure,
        pair_style,
        pair_coeff,
        temperature,
        pressure=pressure,
        cores=cores,
        n_equilibration_steps=n_equilibration_steps,
        n_run_steps=n_run_steps,
        cdict=cdict,
        potential_type=potential_type,
        potential_doi=potential_doi,
        autovalidate=autovalidate,
        cv_threshold=cv_threshold,
        alpha_min=alpha_min,
        alpha_max=alpha_max,
        convergence_threshold=convergence_threshold,
    )
    return datadict


@as_function_node
def calculate_free_energy(
    structure,
    pair_style,
    pair_coeff,
    temperature,
    pressure=0,
    n_switching_steps=10000,
    mode="fe",
    reference_phase="solid",
    cores=1,
    folder_prefix="calc",
    input_dict=None,
    cdict=None,
    potential_type=None,
    potential_doi=None,
):
    """Calculate free energy (pyiron-decorated)."""
    datadict = _core.calculate_free_energy(
        structure,
        pair_style,
        pair_coeff,
        temperature,
        pressure=pressure,
        n_switching_steps=n_switching_steps,
        mode=mode,
        reference_phase=reference_phase,
        cores=cores,
        folder_prefix=folder_prefix,
        input_dict=input_dict,
        cdict=cdict,
        potential_type=potential_type,
        potential_doi=potential_doi,
    )
    return datadict


# Re-export utility functions directly from core
calphy_default_input = _core.calphy_default_input
