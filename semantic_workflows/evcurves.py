"""
Core EV curves calculation functions.

These are framework-agnostic implementations that can be decorated
for use with pyiron_workflow or jobflow.
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from ase.io import read, write
from pylammpsmpi import LammpsLibrary


def calculate_ev_curves(
    structure,
    pair_style,
    pair_coeff,
    vol_range=0.3,
    num_of_points=5,
    cores=1,
    e_tol=0,
    f_tol=0.0001,
    n_energy_steps=1e5,
    n_force_steps=1e6,
    cdict=None,
    potential_type=None,
    potential_doi=None,
):
    """
    Calculate energy-volume curves for a structure.

    Parameters
    ----------
    structure : ase.Atoms
        Input atomic structure
    pair_style : str
        LAMMPS pair style
    pair_coeff : str
        LAMMPS pair coefficients
    vol_range : float, optional
        Volume range for scaling (default: 0.3)
    num_of_points : int, optional
        Number of points in the curve (default: 5)
    cores : int, optional
        Number of cores for LAMMPS (default: 1)
    e_tol : float, optional
        Energy tolerance for minimization (default: 0)
    f_tol : float, optional
        Force tolerance for minimization (default: 0.0001)
    n_energy_steps : float, optional
        Maximum energy minimization steps (default: 1e5)
    n_force_steps : float, optional
        Maximum force minimization steps (default: 1e6)
    cdict : dict, optional
        Knowledge graph dictionary for metadata (pyiron-specific, default: None)
    potential_type : str, optional
        Potential type for cdict (default: None)
    potential_doi : str, optional
        DOI for the potential (default: None)

    Returns
    -------
    dict
        Dictionary containing 'energy', 'volume', and 'bulk_modulus' arrays
    """
    volume_factors = np.linspace((1 - vol_range), (1.0 + vol_range), num_of_points)

    energies = []
    volumes = []
    initial_cell = structure.get_cell()
    initial_volume = structure.get_volume()

    for volume_factor in volume_factors:
        scaled_atoms = scale_atoms(structure, volume_factor)
        e, v = calculate_energy(
            scaled_atoms,
            pair_style,
            pair_coeff,
            cores=cores,
            e_tol=e_tol,
            f_tol=f_tol,
            n_energy_steps=n_energy_steps,
            n_force_steps=n_force_steps,
        )

        energies.append(e)
        volumes.append(v)
    V0, E0, B0, Bp = fit_bm(volumes, energies)
    v_fit = np.linspace(min(volumes), max(volumes), 100, endpoint=True)
    e_fit = birch_murnaghan_eval(v_fit, V0, E0, B0, Bp)
    bulk_modulus = B0 * 160.2176621
    datadict = {"energy": e_fit, "volume": v_fit, "bulk_modulus": bulk_modulus}

    # Knowledge graph integration (pyiron-specific)
    if cdict is not None:
        from .build import update_attributes
        from conceptual_dictionary import workflow_template

        final_volume = V0 * len(structure)
        scaling = final_volume / initial_volume
        final_structure = scale_atoms(structure, scaling)
        # A new structure is created and data is added!
        final_structure = update_attributes(final_structure, cdict, create_new=True)

        # Add workflow details
        workflow = workflow_template.copy()

        workflow["method"] = "MolecularStatics"
        workflow["algorithm"] = "EquationOfStateFit"
        workflow["input_sample"] = [structure.info["id"]]

        new_id = final_structure.info["id"]
        workflow["output_sample"] = [new_id]

        outputs = [
            {
                "label": "EquilibriumEnergy",
                "value": np.round(E0, decimals=4),
                "unit": "EV",
                "associate_to_sample": [new_id],
                "basename": "TotalEnergy",
            },
            {
                "label": "EquilibriumVolume",
                "value": np.round(V0, decimals=4),
                "unit": "ANGSTROM3",
                "associate_to_sample": [new_id],
                "basename": "Volume",
            },
            {
                "label": "BulkModulus",
                "basename": "BulkModulus",
                "value": np.round(bulk_modulus, decimals=2),
                "unit": "GigaPA",
                "associate_to_sample": [new_id],
            },
        ]

        workflow["degrees_of_freedom"] = ["AtomicPositionRelaxation"]
        workflow["calculated_property"] = outputs
        workflow["interatomic_potential"] = {
            "potential_type": potential_type,
            "uri": potential_doi,
        }
        workflow["software"] = [
            {
                "uri": "https://doi.org/10.1016/j.cpc.2021.108171",
                "label": "LAMMPS",
            }
        ]

        cdict["workflow"].append(workflow.copy())

    return datadict


def scale_atoms(structure, scale_factor):
    """
    Scale atomic structure uniformly.

    Parameters
    ----------
    structure : ase.Atoms
        Input structure
    scale_factor : float
        Scaling factor

    Returns
    -------
    ase.Atoms
        Scaled structure
    """
    scaled_atoms = structure.copy()
    scaled_atoms.set_cell(scaled_atoms.get_cell() * scale_factor, scale_atoms=True)
    return scaled_atoms


def fit_bm(vol, en):
    """
    Fit Birch-Murnaghan equation of state.

    Parameters
    ----------
    vol : array-like
        Volume values
    en : array-like
        Energy values

    Returns
    -------
    tuple
        (V0, E0, B0, Bp) - equilibrium volume, energy, bulk modulus, and derivative
    """
    a, b, c = np.polyfit(vol, en, 2)
    V0 = -b / (2 * a)
    E0 = a * V0**2 + b * V0 + c
    B0 = 2 * a * V0
    Bp = 4.0
    popt, pcov = curve_fit(birch_murnaghan_eval, vol, en, p0=[V0, E0, B0, Bp])
    return popt


def birch_murnaghan_eval(vol, V0, E0, B0, Bp):
    """
    Evaluate Birch-Murnaghan equation of state.

    Parameters
    ----------
    vol : array-like
        Volume values
    V0 : float
        Equilibrium volume
    E0 : float
        Equilibrium energy
    B0 : float
        Bulk modulus
    Bp : float
        Derivative of bulk modulus

    Returns
    -------
    array-like
        Energy values
    """
    eta = (vol / V0) ** (1.0 / 3.0)
    E = E0 + 9.0 * B0 * V0 / 16.0 * (eta**2 - 1.0) ** 2 * (
        6.0 + Bp * (eta**2 - 1.0) - 4.0 * eta**2
    )
    return E


def relax_structure(
    structure,
    pair_style,
    pair_coeff,
    cores=1,
    e_tol=0,
    f_tol=0.0001,
    n_energy_steps=1e5,
    n_force_steps=1e6,
    cdict=None,
    potential_type=None,
    potential_doi=None,
):
    """
    Relax structure with box relaxation.

    Parameters
    ----------
    structure : ase.Atoms
        Input structure
    pair_style : str
        LAMMPS pair style
    pair_coeff : str
        LAMMPS pair coefficients
    cores : int, optional
        Number of cores (default: 1)
    e_tol : float, optional
        Energy tolerance (default: 0)
    f_tol : float, optional
        Force tolerance (default: 0.0001)
    n_energy_steps : float, optional
        Max energy steps (default: 1e5)
    n_force_steps : float, optional
        Max force steps (default: 1e6)
    cdict : dict, optional
        Knowledge graph (default: None)
    potential_type : str, optional
        Potential type (default: None)
    potential_doi : str, optional
        Potential DOI (default: None)

    Returns
    -------
    tuple
        (final_structure, ecoh, vol) - relaxed structure, cohesive energy, volume
    """
    from mendeleev import element

    symbols, counts = list(
        np.unique(structure.get_chemical_symbols(), return_counts=True)
    )
    masses = [element(symbol).mass for symbol in symbols]
    write("tmp.data", structure, format="lammps-data")

    initial_cell = structure.get_cell()
    initial_volume = structure.get_volume()

    lmp = LammpsLibrary(cores=cores)
    lmp.command("units metal")
    lmp.command("dimension 3")
    lmp.command("boundary p p p")
    lmp.command("atom_style atomic")
    lmp.command("read_data tmp.data")
    for x in range(len(masses)):
        lmp.command(f"mass {x+1} {masses[x]}")

    lmp.command(f"pair_style {pair_style}")

    pair_coeff = np.atleast_1d(pair_coeff)
    for pair in pair_coeff:
        lmp.command(f"pair_coeff {pair}")

    lmp.command("fix ensemble all box/relax aniso 0.0")
    lmp.command(f"minimize {e_tol} {f_tol} {int(n_energy_steps)} {int(n_force_steps)}")

    lmp.command("run 0")
    ecoh = lmp.pe / lmp.natoms
    final_volume = lmp.vol
    vol = lmp.vol / lmp.natoms

    filename = "tmp.dump"
    lmp.command(
        "dump              2x all custom 1 %s id type mass x y z vx vy vz" % (filename)
    )
    lmp.command("run               0")
    lmp.command("undump            2x")

    lmp.close()

    final_structure = read("tmp.dump", format="lammps-dump-text")

    if cdict is not None:
        from .build import update_attributes
        from conceptual_dictionary import workflow_template

        final_structure.info["id"] = structure.info["id"]
        final_structure = update_attributes(final_structure, cdict, create_new=True)

        workflow = workflow_template.copy()

        workflow["method"] = "MolecularStatics"
        workflow["input_sample"] = [structure.info["id"]]

        new_id = final_structure.info["id"]
        workflow["output_sample"] = [new_id]

        outputs = [
            {
                "label": "EquilibriumEnergy",
                "value": np.round(ecoh, decimals=4),
                "unit": "EV",
                "associate_to_sample": [new_id],
                "basename": "TotalEnergy",
            },
            {
                "label": "EquilibriumVolume",
                "value": np.round(vol, decimals=4),
                "unit": "ANGSTROM3",
                "associate_to_sample": [new_id],
                "basename": "Volume",
            },
        ]

        workflow["degrees_of_freedom"] = [
            "AtomicPositionRelaxation",
            "CellVolumeRelaxation",
        ]
        workflow["calculated_property"] = outputs
        workflow["interatomic_potential"] = {
            "potential_type": potential_type,
            "uri": potential_doi,
        }
        workflow["software"] = [
            {
                "uri": "https://doi.org/10.1016/j.cpc.2021.108171",
                "label": "LAMMPS",
            }
        ]

        cdict["workflow"].append(workflow.copy())
    return final_structure, ecoh, vol


def calculate_energy(
    structure,
    pair_style,
    pair_coeff,
    cores=1,
    e_tol=0,
    f_tol=0.0001,
    n_energy_steps=1e5,
    n_force_steps=1e6,
):
    """
    Calculate energy of a structure.

    Parameters
    ----------
    structure : ase.Atoms
        Input structure
    pair_style : str
        LAMMPS pair style
    pair_coeff : str
        LAMMPS pair coefficients
    cores : int, optional
        Number of cores (default: 1)
    e_tol : float, optional
        Energy tolerance (default: 0)
    f_tol : float, optional
        Force tolerance (default: 0.0001)
    n_energy_steps : float, optional
        Max energy steps (default: 1e5)
    n_force_steps : float, optional
        Max force steps (default: 1e6)

    Returns
    -------
    tuple
        (ecoh, vol) - cohesive energy and volume per atom
    """
    from mendeleev import element

    symbols, counts = list(
        np.unique(structure.get_chemical_symbols(), return_counts=True)
    )
    masses = [element(symbol).mass for symbol in symbols]

    write("tmp.data", structure, format="lammps-data")
    lmp = LammpsLibrary(cores=cores)
    lmp.command("units metal")
    lmp.command("dimension 3")
    lmp.command("boundary p p p")
    lmp.command("atom_style atomic")
    lmp.command("read_data tmp.data")
    for x in range(len(masses)):
        lmp.command(f"mass {x+1} {masses[x]}")

    lmp.command(f"pair_style {pair_style}")
    lmp.command(f"pair_coeff {pair_coeff}")
    lmp.command("run 0")
    ecoh = lmp.pe / lmp.natoms
    vol = lmp.vol / lmp.natoms
    return ecoh, vol
