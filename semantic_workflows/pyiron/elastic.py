"""
Pyiron-decorated elastic properties workflow functions.

These functions import core logic from the parent workflows module and
apply pyiron_workflow decorators.
"""

from pyiron_workflow import as_function_node
from .. import elastic as _core


@as_function_node
def calculate_elastic_constants(
    structure,
    pair_style,
    pair_coeff,
    cores=1,
    finite_deformation_size=1e-6,
    energy_tolerance=0.0,
    force_tolerance=1.0e-15,
    max_iterations=100,
    cdict=None,
    potential_type=None,
    potential_doi=None,
):
    """Calculate elastic constants (pyiron-decorated)."""
    datadict = _core.calculate_elastic_constants(
        structure,
        pair_style,
        pair_coeff,
        cores=cores,
        finite_deformation_size=finite_deformation_size,
        energy_tolerance=energy_tolerance,
        force_tolerance=force_tolerance,
        max_iterations=max_iterations,
        cdict=cdict,
        potential_type=potential_type,
        potential_doi=potential_doi,
    )
    return datadict


@as_function_node
def mechanical_response_test(
    structure,
    pair_style,
    pair_coeff,
    mode="hydrostatic",
    cores=1,
    temperature=10,
    annealing_temperature=600,
    n_equilibration_steps=5000,
    n_run_steps=10000,
    strain_rate=1e-5,
    dump_interval=10000,
    cdict=None,
    potential_type=None,
    potential_doi=None,
):
    """Perform mechanical response test (pyiron-decorated)."""
    datadict = _core.mechanical_response_test(
        structure,
        pair_style,
        pair_coeff,
        mode=mode,
        cores=cores,
        temperature=temperature,
        annealing_temperature=annealing_temperature,
        n_equilibration_steps=n_equilibration_steps,
        n_run_steps=n_run_steps,
        strain_rate=strain_rate,
        dump_interval=dump_interval,
        cdict=cdict,
        potential_type=potential_type,
        potential_doi=potential_doi,
    )
    return datadict
