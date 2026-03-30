"""
Core thermodynamics workflow functions.

These functions contain the computational logic for thermal properties and
free energy calculations, without framework-specific decorators.
"""

import pandas as pd
from ase.io import read, write
from pylammpsmpi import LammpsLibrary
import scipy.constants as sc
import numpy as np


def calphy_default_input():
    """Return default input dictionary for calphy calculations."""
    return {
        "mode": None,
        "pressure": 0,
        "temperature": 0,
        "reference_phase": None,
        "npt": True,
        "n_equilibration_steps": 15000,
        "n_switching_steps": 25000,
        "n_print_steps": 1000,
        "n_iterations": 1,
        "spring_constants": None,
        "equilibration_control": "nose-hoover",
        "melting_cycle": True,
        "folder_prefix": "calc",
        "md": {
            "timestep": 0.001,
            "n_small_steps": 10000,
            "n_every_steps": 10,
            "n_repeat_steps": 10,
            "n_cycles": 100,
            "thermostat_damping": 0.5,
            "barostat_damping": 0.1,
        },
        "tolerance": {
            "lattice_constant": 0.0002,
            "spring_constant": 0.01,
            "solid_fraction": 0.7,
            "liquid_fraction": 0.05,
            "pressure": 0.5,
        },
        "nose_hoover": {
            "thermostat_damping": 0.1,
            "barostat_damping": 0.1,
        },
        "berendsen": {
            "thermostat_damping": 100.0,
            "barostat_damping": 100.0,
        },
        "queue": {
            "cores": 1,
        },
    }


def _write_structure(structure, filename, elements):
    """Write structure to LAMMPS data file format with specified element ordering."""
    from ase.io import write
    from mendeleev import element

    write(filename, structure, format="lammps-data", specorder=elements)
    masses = [element(symbol).mass for symbol in elements]
    return filename, masses


def _prepare_input(
    pair_style,
    pair_coeff,
    structure,
    temperature,
    pressure=0,
    n_switching_steps=10000,
    mode="fe",
    reference_phase="solid",
    cores=1,
    folder_prefix="calc",
    input_dict=None,
):
    """Prepare input for calphy calculation."""
    from calphy.input import Calculation
    import os
    import uuid
    import numpy as np

    if input_dict is None:
        input_dict = calphy_default_input()

    filename = f"tmp-{uuid.uuid4().hex[:8]}.data"
    elements = list(np.unique(structure.get_chemical_symbols()))
    # Ensure elements are present in the pair_coeff
    element_str = " ".join(elements)
    if element_str not in pair_coeff:
        raise ValueError(
            f"Elements {element_str} not present in pair_coeff {pair_coeff}"
        )

    filename, masses = _write_structure(structure, filename, elements)

    input_dict["pair_style"] = pair_style
    input_dict["pair_coeff"] = pair_coeff
    input_dict["element"] = elements
    input_dict["mass"] = masses
    input_dict["mode"] = mode
    input_dict["reference_phase"] = reference_phase
    input_dict["lattice"] = filename
    input_dict["temperature"] = temperature
    input_dict["pressure"] = pressure
    input_dict["n_switching_steps"] = n_switching_steps
    input_dict["n_equilibration_steps"] = n_switching_steps // 2

    if cores > 1:
        input_dict["queue"]["cores"] = cores

    calc = Calculation(**input_dict)
    calc.folder_prefix = folder_prefix
    return calc


def _run_cleanup(simfolder, lattice, delete_folder=False):
    """Clean up temporary files after calculation."""
    import shutil
    import os

    os.remove(lattice)
    if delete_folder:
        shutil.rmtree(simfolder)


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
    """
    Calculate thermal properties using NPT molecular dynamics.

    Computes specific heat, thermal expansion coefficient, and equilibrium volume
    at specified temperature and pressure conditions.

    Parameters
    ----------
    structure : ase.Atoms
        Input atomic structure
    pair_style : str
        LAMMPS pair style
    pair_coeff : str
        LAMMPS pair coefficients
    temperature : float
        Temperature in K
    pressure : float, optional
        Pressure in GPa (default: 0)
    cores : int, optional
        Number of cores for parallel execution (default: 1)
    n_equilibration_steps : int, optional
        Number of equilibration steps (default: 15000)
    n_run_steps : int, optional
        Number of production run steps (default: 25000)
    cdict : dict, optional
        Knowledge graph dictionary for metadata tracking (default: None)
    potential_type : str, optional
        Type of interatomic potential (default: None)
    potential_doi : str, optional
        DOI of interatomic potential (default: None)
    autovalidate : bool, optional
        If True, run automated validation checks on the trajectory (default: False)
    cv_threshold : float, optional
        Max allowed coefficient of variation for temperature stability and
        max relative half-half difference for pressure stability (default: 0.05)
    alpha_min : float, optional
        Minimum physically reasonable thermal expansion coefficient in K⁻¹ (default: 1e-6)
    alpha_max : float, optional
        Maximum physically reasonable thermal expansion coefficient in K⁻¹ (default: 1e-4)
    convergence_threshold : float, optional
        Max allowed relative difference between first and second half of energy
        trajectory for statistical convergence check (default: 0.1)

    Returns
    -------
    dict
        Dictionary containing specific_heat, thermal_expansion, and volume.
        If autovalidate=True, also contains is_valid (bool) and validation (dict)
        with keys temp_stable, press_stable, statistically_converged, physical_range.
    """
    ev_to_j = sc.physical_constants["electron volt-joule relationship"][0]
    Av = sc.physical_constants["Avogadro constant"][0]
    A_to_m = 1e-10
    kB = sc.physical_constants["Boltzmann constant"][0]

    write("tmp.data", structure, format="lammps-data")

    from mendeleev import element

    symbols, counts = list(
        np.unique(structure.get_chemical_symbols(), return_counts=True)
    )
    masses = [element(symbol).mass for symbol in symbols]

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

    T = temperature
    P = pressure

    lmp.velocity("all create", T, np.random.randint(0, 1000))
    lmp.fix("1 all npt temp", T, T, 0.1, "iso", 0.0, 0.0, 0.1)
    lmp.run(n_equilibration_steps)
    lmp.fix(
        'extra all print 100 " $(pe) $(ke) $(etotal) $(press) $(vol) $(temp)" file data.dat screen no'
    )

    lmp.run(n_run_steps)

    # Make it dump once
    filename = "tmp.dump"
    lmp.command(
        "dump              2x all custom 1 %s id type mass x y z vx vy vz" % (filename)
    )
    lmp.command("run               0")
    lmp.command("undump            2x")

    pe, ke, etotal, press, vol, temp = np.loadtxt("data.dat", unpack=True)
    efluct = etotal - np.mean(etotal)

    efluctsq = (efluct * ev_to_j) ** 2
    cp = np.mean(efluctsq) / (kB * T * T)
    w = 0
    for m, c in zip(masses, counts):
        w += (c / Av) * m * 1e3

    cp = (cp / w) / 1e-3

    crossfluct = (efluct * ev_to_j) * (vol - np.mean(vol))
    ap = np.mean(crossfluct) / (kB * T * T * np.mean(vol))  # 1/K

    # cdict-related code - only execute when cdict is not None
    if cdict is not None:
        from .build import update_attributes
        from conceptual_dictionary import workflow_template

        # Read dump system and update
        final_structure = read("tmp.dump", format="lammps-dump-text")
        final_structure.info["id"] = structure.info["id"]
        final_structure = update_attributes(final_structure, cdict, create_new=True)

        workflow = workflow_template.copy()

        workflow["method"] = "MolecularDynamics"
        workflow["input_sample"] = [structure.info["id"]]
        new_id = final_structure.info["id"]
        workflow["output_sample"] = [new_id]
        inputs = [
            {
                "basename": "Temperature",
                "label": "Temperature",
                "value": temperature,
                "unit": "K",
            },
            {
                "basename": "Pressure",
                "label": "Pressure",
                "value": pressure,
                "unit": "GigaPA",
            },
        ]
        workflow["input_parameter"] = inputs
        outputs = [
            {
                "label": "SpecificHeatCapacity",
                "value": cp,
                "unit": "J-PER-GM-K",
                "associate_to_sample": [new_id],
                "basename": "SpecificHeatCapacity",
            },
            {
                "label": "EquilibriumVolume",
                "value": np.round(np.mean(vol), decimals=4),
                "unit": "ANGSTROM3",
                "associate_to_sample": [new_id],
                "basename": "Volume",
            },
            {
                "label": "ThermalExpansionCoefficient",
                "basename": "ThermalExpansionCoefficient",
                "value": ap,
                "unit": "PER-K",
                "associate_to_sample": [new_id],
            },
        ]

        workflow["thermodynamic_ensemble"] = "IsothermalIsobaricEnsemble"
        workflow["degrees_of_freedom"] = [
            "AtomicPositionRelaxation",
            "CellVolumeRelaxation",
            "CellShapeRelaxation",
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

        # Append a *copy* to avoid overwriting in subsequent iterations
        cdict["workflow"].append(workflow.copy())

    results = {"specific_heat": cp, "thermal_expansion": ap, "volume": np.mean(vol)}

    if autovalidate:
        n = len(temp)
        mid = n // 2

        # 1. Temperature stability: trajectory mean close to target
        temp_stable = bool(
            abs(np.mean(temp) - temperature) / temperature < cv_threshold
        )

        # 2. Pressure stability: first and second half means agree
        press1, press2 = np.mean(press[:mid]), np.mean(press[mid:])
        press_denom = max(abs(np.mean(press)), 1.0)  # avoid div-by-zero at P=0
        press_stable = bool(abs(press1 - press2) / press_denom < cv_threshold)

        # 3. Statistical convergence: total energy mean converged between halves
        e_denom = abs(np.mean(etotal))
        if e_denom > 0:
            statistically_converged = bool(
                abs(np.mean(etotal[:mid]) - np.mean(etotal[mid:])) / e_denom
                < convergence_threshold
            )
        else:
            statistically_converged = False

        # 4. Physical range: positive Cp and reasonable thermal expansion
        physical_range = bool(cp > 0 and alpha_min <= ap <= alpha_max)

        is_valid = all(
            [temp_stable, press_stable, statistically_converged, physical_range]
        )
        results["is_valid"] = is_valid
        results["validation"] = {
            "temp_stable": temp_stable,
            "press_stable": press_stable,
            "statistically_converged": statistically_converged,
            "physical_range": physical_range,
        }

    return results


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
    """
    Calculate free energy using calphy thermodynamic integration.

    Supports different calculation modes (fe, ts, pscale) and reference phases
    (solid, liquid) for free energy calculations.

    Parameters
    ----------
    structure : ase.Atoms
        Input atomic structure
    pair_style : str
        LAMMPS pair style
    pair_coeff : str
        LAMMPS pair coefficients
    temperature : float
        Temperature in K
    pressure : float, optional
        Pressure in GPa (default: 0)
    n_switching_steps : int, optional
        Number of switching steps (default: 10000)
    mode : str, optional
        Calculation mode: 'fe', 'ts', or 'pscale' (default: 'fe')
    reference_phase : str, optional
        Reference phase: 'solid' or 'liquid' (default: 'solid')
    cores : int, optional
        Number of cores for parallel execution (default: 1)
    folder_prefix : str, optional
        Prefix for calculation folder (default: 'calc')
    input_dict : dict, optional
        Custom input dictionary for calphy (default: None)
    cdict : dict, optional
        Knowledge graph dictionary for metadata tracking (default: None)
    potential_type : str, optional
        Type of interatomic potential (default: None)
    potential_doi : str, optional
        DOI of interatomic potential (default: None)

    Returns
    -------
    dict
        Dictionary containing calculation results including free energy,
        temperature, and pressure
    """
    from calphy.solid import Solid
    from calphy.liquid import Liquid
    from calphy.routines import routine_fe, routine_ts, routine_pscale
    import os
    import numpy as np

    calc = _prepare_input(
        pair_style,
        pair_coeff,
        structure,
        temperature,
        pressure,
        n_switching_steps,
        mode,
        reference_phase,
        cores,
        folder_prefix,
        input_dict,
    )

    simfolder = calc.create_folders()
    if reference_phase == "solid":
        job = Solid(calculation=calc, simfolder=simfolder)
    elif reference_phase == "liquid":
        job = Liquid(calculation=calc, simfolder=simfolder)

    if mode == "fe":
        job = routine_fe(job)
    elif mode == "ts":
        job = routine_ts(job)
    elif mode == "pscale":
        job = routine_pscale(job)
    else:
        raise ValueError(f"Unknown mode {mode}, should be one of fe, ts, pscale")

    _run_cleanup(simfolder, calc.lattice)

    results = job.report["results"]
    results["temperature"] = temperature
    results["pressure"] = pressure

    if mode == "ts":
        outfile = os.path.join(simfolder, "temperature_sweep.dat")
        temp, fe = np.loadtxt(outfile, unpack=True, usecols=(0, 1))
        results["temperature"] = temp
        results["free_energy"] = fe
        pres = pressure
    elif mode == "pscale":
        outfile = os.path.join(simfolder, "pressure_sweep.dat")
        pres, fe = np.loadtxt(outfile, unpack=True, usecols=(0, 1))
        results["pressure"] = pres
        results["free_energy"] = fe
        temp = temperature
    else:
        # For mode == 'fe', use the values from results
        temp = temperature
        pres = pressure
        fe = results.get("free_energy", results.get("helmholtz_energy", None))

    # cdict-related code - only execute when cdict is not None
    if cdict is not None:
        from .build import update_attributes
        from conceptual_dictionary import workflow_template

        final_structure = read(
            os.path.join(simfolder, "conf.equilibration.data"),
            format="lammps-data",
            atom_style="atomic",
        )
        final_structure.info["id"] = structure.info["id"]
        final_structure = update_attributes(final_structure, cdict, create_new=True)

        workflow = workflow_template.copy()

        workflow["method"] = "MolecularDynamics"
        workflow["algorithm"] = "ThermodynamicIntegration"
        workflow["input_sample"] = [structure.info["id"]]
        new_id = final_structure.info["id"]
        workflow["output_sample"] = [new_id]
        inputs = [
            {
                "basename": "Temperature",
                "label": "Temperature",
                "value": temperature,
                "unit": "K",
            },
            {
                "basename": "Pressure",
                "label": "Pressure",
                "value": pressure,
                "unit": "GigaPA",
            },
        ]
        workflow["input_parameter"] = inputs
        outputs = []
        outputs.append(
            {
                "label": "FreeEnergy",
                "basename": "FreeEnergy",
                "value": fe,
                "unit": "EV",
                "associate_to_sample": [new_id],
            }
        )
        outputs.append(
            {
                "label": "VirialPressure",
                "basename": "VirialPressure",
                "value": pres,
                "unit": "GigaPA",
                "associate_to_sample": [new_id],
            }
        )
        outputs.append(
            {
                "label": "Temperature",
                "basename": "Temperature",
                "value": temp,
                "unit": "K",
                "associate_to_sample": [new_id],
            }
        )

        workflow["thermodynamic_ensemble"] = "IsothermalIsobaricEnsemble"
        workflow["degrees_of_freedom"] = [
            "AtomicPositionRelaxation",
            "CellVolumeRelaxation",
            "CellShapeRelaxation",
        ]
        workflow["calculated_property"] = outputs
        workflow["interatomic_potential"] = {
            "potential_type": potential_type,
            "uri": potential_doi,
        }
        software1 = {
            "uri": "https://doi.org/10.1016/j.cpc.2021.108171",
            "label": "LAMMPS",
        }
        software2 = {
            "uri": "https://doi.org/10.5281/zenodo.10527452",
            "label": "Calphy",
        }
        workflow["software"] = [software1, software2]

        # Append a *copy* to avoid overwriting in subsequent iterations
        cdict["workflow"].append(workflow.copy())

    return results
