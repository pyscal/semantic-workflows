"""
Core EV curves calculation functions.

These are framework-agnostic implementations that can be decorated
for use with pyiron_workflow or jobflow.

Supported engines: 'lammps', 'qe' (Quantum ESPRESSO).
"""

import os
import numpy as np
from scipy.optimize import curve_fit
from ase.io import read, write


# ---------------------------------------------------------------------------
# Shared math utilities
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# LAMMPS private helpers
# ---------------------------------------------------------------------------


def _calculate_energy_lammps(structure, pair_style, pair_coeff, cores=1):
    """Run a LAMMPS single-point and return (ecoh_per_atom, vol_per_atom)."""
    from mendeleev import element
    from pylammpsmpi import LammpsLibrary

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
    lmp.close()
    return ecoh, vol


def _relax_lammps(
    structure,
    pair_style,
    pair_coeff,
    cores=1,
    e_tol=0,
    f_tol=0.0001,
    n_energy_steps=1e5,
    n_force_steps=1e6,
):
    """Run a LAMMPS box relaxation and return (final_structure, ecoh_per_atom, vol_per_atom)."""
    from mendeleev import element
    from pylammpsmpi import LammpsLibrary

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
    pair_coeff = np.atleast_1d(pair_coeff)
    for pair in pair_coeff:
        lmp.command(f"pair_coeff {pair}")
    lmp.command("fix ensemble all box/relax aniso 0.0")
    lmp.command(f"minimize {e_tol} {f_tol} {int(n_energy_steps)} {int(n_force_steps)}")
    lmp.command("run 0")
    ecoh = lmp.pe / lmp.natoms
    vol = lmp.vol / lmp.natoms

    dump_file = "tmp.dump"
    lmp.command("dump 2x all custom 1 %s id type mass x y z vx vy vz" % dump_file)
    lmp.command("run 0")
    lmp.command("undump 2x")
    lmp.close()

    final_structure = read(dump_file, format="lammps-dump-text")
    return final_structure, ecoh, vol


# ---------------------------------------------------------------------------
# QE private helpers
# ---------------------------------------------------------------------------


def _parse_qe_template(template_file):
    """Parse a QE input template file to extract namelists and k-points.

    Returns
    -------
    tuple
        (input_data dict, kpts, koffset) where kpts is a 3-tuple of ints or
        None for Gamma-only, and koffset is a 3-tuple of ints.
    """
    from ase.io.espresso import read_fortran_namelist

    with open(template_file) as f:
        input_data, _ = read_fortran_namelist(f)

    kpts = (1, 1, 1)
    koffset = (0, 0, 0)
    with open(template_file) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "K_POINTS" in line.upper():
            mode = line.upper().replace("K_POINTS", "").strip().lower()
            if "gamma" in mode:
                kpts = None
                koffset = (0, 0, 0)
            else:  # automatic (default)
                for j in range(i + 1, len(lines)):
                    vals = lines[j].split()
                    if vals:
                        kpts = tuple(int(v) for v in vals[:3])
                        koffset = (
                            tuple(int(v) for v in vals[3:6])
                            if len(vals) >= 6
                            else (0, 0, 0)
                        )
                        break
            break
    return input_data, kpts, koffset


def _write_qe_input(filename, structure, input_data, pseudopotentials, kpts, koffset):
    """Write a QE pw.x input file via ASE."""
    write(
        filename,
        structure,
        format="espresso-in",
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kpts=kpts,
        koffset=koffset,
    )


def _run_qe(input_file, output_file, qe_command="pw.x"):
    """Run a QE calculation, piping stdout/stderr to output_file."""
    import subprocess

    with open(output_file, "w") as fout:
        subprocess.run(
            f"{qe_command} < {input_file}",
            shell=True,
            stdout=fout,
            stderr=subprocess.STDOUT,
            check=True,
        )


def _extract_xc_functional(outfile):
    """Extract the XC functional label from a QE output file."""
    with open(outfile) as f:
        for line in f:
            if "Exchange-correlation" in line and "=" in line:
                return line.split("=")[1].split("(")[0].strip()
    return None


def _extract_qe_input_params(input_data):
    """Build a cdict-compatible input parameter list from parsed QE namelists."""
    params = []
    ecutwfc = input_data.get("system", {}).get("ecutwfc")
    if ecutwfc is not None:
        params.append(
            {
                "label": "EnergyCutoff",
                "value": round(float(ecutwfc) * 13.6057039763, 4),  # Ry → eV
                "unit": "EV",
            }
        )
    return params


def _calculate_energy_qe(
    structure, input_template, pseudopotentials, working_dir, qe_command="pw.x"
):
    """Run a QE SCF single-point.

    Returns
    -------
    tuple
        (ecoh_per_atom, vol_per_atom, xc_functional, input_params)
    """
    os.makedirs(working_dir, exist_ok=True)
    input_data, kpts, koffset = _parse_qe_template(input_template)
    input_data.setdefault("control", {})["calculation"] = "scf"
    input_data.setdefault("control", {})["outdir"] = os.path.abspath(working_dir)

    infile = os.path.join(working_dir, "pw.in")
    outfile = os.path.join(working_dir, "pw.out")
    _write_qe_input(infile, structure, input_data, pseudopotentials, kpts, koffset)
    _run_qe(infile, outfile, qe_command)

    result = read(outfile, format="espresso-out")
    ecoh = result.get_total_energy() / len(structure)
    vol = result.get_volume() / len(structure)
    xc = _extract_xc_functional(outfile)
    input_params = _extract_qe_input_params(input_data)
    return ecoh, vol, xc, input_params


def _relax_qe(
    structure, input_template, pseudopotentials, working_dir, qe_command="pw.x"
):
    """Run a QE vc-relax calculation.

    Returns
    -------
    tuple
        (final_structure, ecoh_per_atom, vol_per_atom, xc_functional, input_params)
    """
    os.makedirs(working_dir, exist_ok=True)
    input_data, kpts, koffset = _parse_qe_template(input_template)
    input_data.setdefault("control", {})["calculation"] = "vc-relax"
    input_data.setdefault("control", {})["outdir"] = os.path.abspath(working_dir)

    infile = os.path.join(working_dir, "pw.in")
    outfile = os.path.join(working_dir, "pw.out")
    _write_qe_input(infile, structure, input_data, pseudopotentials, kpts, koffset)
    _run_qe(infile, outfile, qe_command)

    # index=-1 reads the last (converged) ionic step
    final_structure = read(outfile, format="espresso-out", index=-1)
    ecoh = final_structure.get_total_energy() / len(final_structure)
    vol = final_structure.get_volume() / len(final_structure)
    xc = _extract_xc_functional(outfile)
    input_params = _extract_qe_input_params(input_data)
    return final_structure, ecoh, vol, xc, input_params


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def calculate_ev_curves(
    structure,
    engine="lammps",
    vol_range=0.3,
    num_of_points=5,
    # LAMMPS params
    pair_style=None,
    pair_coeff=None,
    cores=1,
    e_tol=0,
    f_tol=0.0001,
    n_energy_steps=1e5,
    n_force_steps=1e6,
    # QE params
    input_template=None,
    pseudopotentials=None,
    working_dir="qe_ev_run",
    qe_command="pw.x",
    # cdict params
    cdict=None,
    potential_type=None,
    potential_doi=None,
):
    """
    Calculate energy-volume curves for a structure.

    The structure should be pre-relaxed before calling this function.
    Single-point energies are computed at each scaled volume.

    Parameters
    ----------
    structure : ase.Atoms
        Input atomic structure (pre-relaxed)
    engine : str, optional
        Simulation engine: 'lammps' or 'qe' (default: 'lammps')
    vol_range : float, optional
        Fractional volume range for scaling (default: 0.3)
    num_of_points : int, optional
        Number of volume points in the curve (default: 5)
    pair_style : str, optional
        LAMMPS pair style (engine='lammps' only)
    pair_coeff : str, optional
        LAMMPS pair coefficients (engine='lammps' only)
    cores : int, optional
        Number of cores for LAMMPS (default: 1)
    e_tol : float, optional
        Energy tolerance for LAMMPS (default: 0)
    f_tol : float, optional
        Force tolerance for LAMMPS (default: 0.0001)
    n_energy_steps : float, optional
        Maximum LAMMPS energy minimization steps (default: 1e5)
    n_force_steps : float, optional
        Maximum LAMMPS force minimization steps (default: 1e6)
    input_template : str, optional
        Path to QE input template file with namelists and K_POINTS (engine='qe' only)
    pseudopotentials : dict, optional
        Mapping of element symbol to UPF file path, e.g. ``{"Si": "Si.UPF"}``
        (engine='qe' only)
    working_dir : str, optional
        Working directory for QE calculations (default: 'qe_ev_run')
    qe_command : str, optional
        Command to invoke Quantum ESPRESSO pw.x (default: 'pw.x')
    cdict : dict, optional
        Conceptual dictionary for semantic metadata (default: None)
    potential_type : str, optional
        Interatomic potential type for cdict, LAMMPS only (default: None)
    potential_doi : str, optional
        DOI for the interatomic potential, LAMMPS only (default: None)

    Returns
    -------
    dict
        Dictionary containing 'energy', 'volume', and 'bulk_modulus' arrays
    """
    volume_factors = np.linspace((1 - vol_range), (1.0 + vol_range), num_of_points)
    initial_volume = structure.get_volume()

    energies = []
    volumes = []
    xc_functional = None
    input_params = []

    for i, volume_factor in enumerate(volume_factors):
        scaled_atoms = scale_atoms(structure, volume_factor)
        if engine == "lammps":
            e, v = _calculate_energy_lammps(
                scaled_atoms, pair_style, pair_coeff, cores=cores
            )
        elif engine == "qe":
            vol_dir = os.path.join(working_dir, f"vol_{i:03d}")
            e, v, xc_functional, input_params = _calculate_energy_qe(
                scaled_atoms, input_template, pseudopotentials, vol_dir, qe_command
            )
        else:
            raise ValueError(f"Unknown engine '{engine}'. Choose 'lammps' or 'qe'.")
        energies.append(e)
        volumes.append(v)

    V0, E0, B0, Bp = fit_bm(volumes, energies)
    v_fit = np.linspace(min(volumes), max(volumes), 100, endpoint=True)
    e_fit = birch_murnaghan_eval(v_fit, V0, E0, B0, Bp)
    bulk_modulus = B0 * 160.2176621
    datadict = {"energy": e_fit, "volume": v_fit, "bulk_modulus": bulk_modulus}

    if cdict is not None:
        from .build import update_attributes
        from conceptual_dictionary import workflow_template

        final_structure = scale_atoms(structure, (V0 * len(structure) / initial_volume))
        final_structure = update_attributes(final_structure, cdict, create_new=True)
        new_id = final_structure.info["id"]

        workflow = workflow_template.copy()
        workflow["algorithm"] = "EquationOfStateFit"
        workflow["input_sample"] = [structure.info["id"]]
        workflow["output_sample"] = [new_id]
        workflow["calculated_property"] = [
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
                "value": np.round(bulk_modulus, decimals=2),
                "unit": "GigaPA",
                "associate_to_sample": [new_id],
                "basename": "BulkModulus",
            },
        ]

        if engine == "lammps":
            workflow["method"] = "MolecularStatics"
            workflow["degrees_of_freedom"] = ["AtomicPositionRelaxation"]
            workflow["interatomic_potential"] = {
                "potential_type": potential_type,
                "uri": potential_doi,
            }
            workflow["software"] = [
                {"uri": "https://doi.org/10.1016/j.cpc.2021.108171", "label": "LAMMPS"}
            ]
        elif engine == "qe":
            workflow["method"] = "DensityFunctionalTheory"
            workflow["degrees_of_freedom"] = []
            workflow["xc_functional"] = xc_functional
            workflow["input_parameter"] = input_params
            workflow["software"] = [
                {"uri": "https://www.quantum-espresso.org/", "label": "QuantumEspresso"}
            ]

        cdict["workflow"].append(workflow.copy())

    return datadict


def relax_structure(
    structure,
    engine="lammps",
    # LAMMPS params
    pair_style=None,
    pair_coeff=None,
    cores=1,
    e_tol=0,
    f_tol=0.0001,
    n_energy_steps=1e5,
    n_force_steps=1e6,
    # QE params
    input_template=None,
    pseudopotentials=None,
    working_dir="qe_relax_run",
    qe_command="pw.x",
    # cdict params
    cdict=None,
    potential_type=None,
    potential_doi=None,
):
    """
    Relax atomic structure (full cell + ions).

    Parameters
    ----------
    structure : ase.Atoms
        Input structure
    engine : str, optional
        Simulation engine: 'lammps' or 'qe' (default: 'lammps')
    pair_style : str, optional
        LAMMPS pair style (engine='lammps' only)
    pair_coeff : str, optional
        LAMMPS pair coefficients (engine='lammps' only)
    cores : int, optional
        Number of cores for LAMMPS (default: 1)
    e_tol : float, optional
        Energy tolerance for LAMMPS (default: 0)
    f_tol : float, optional
        Force tolerance for LAMMPS (default: 0.0001)
    n_energy_steps : float, optional
        Maximum LAMMPS energy minimization steps (default: 1e5)
    n_force_steps : float, optional
        Maximum LAMMPS force minimization steps (default: 1e6)
    input_template : str, optional
        Path to QE input template file (engine='qe' only)
    pseudopotentials : dict, optional
        Mapping of element symbol to UPF file path, e.g. ``{"Si": "Si.UPF"}``
        (engine='qe' only)
    working_dir : str, optional
        Working directory for QE calculation (default: 'qe_relax_run')
    qe_command : str, optional
        Command to invoke Quantum ESPRESSO pw.x (default: 'pw.x')
    cdict : dict, optional
        Conceptual dictionary for semantic metadata (default: None)
    potential_type : str, optional
        Interatomic potential type for cdict, LAMMPS only (default: None)
    potential_doi : str, optional
        DOI for the interatomic potential, LAMMPS only (default: None)

    Returns
    -------
    tuple
        (final_structure, ecoh_per_atom, vol_per_atom)
    """
    if engine == "lammps":
        final_structure, ecoh, vol = _relax_lammps(
            structure,
            pair_style,
            pair_coeff,
            cores=cores,
            e_tol=e_tol,
            f_tol=f_tol,
            n_energy_steps=n_energy_steps,
            n_force_steps=n_force_steps,
        )
        xc_functional = None
        input_params = []
    elif engine == "qe":
        final_structure, ecoh, vol, xc_functional, input_params = _relax_qe(
            structure, input_template, pseudopotentials, working_dir, qe_command
        )
    else:
        raise ValueError(f"Unknown engine '{engine}'. Choose 'lammps' or 'qe'.")

    if cdict is not None:
        from .build import update_attributes, generate_id as _gen_id
        from conceptual_dictionary import workflow_template

        final_structure.info["id"] = structure.info["id"]
        final_structure = update_attributes(final_structure, cdict, create_new=True)

        workflow = workflow_template.copy()
        new_id = final_structure.info["id"]
        energy_id = _gen_id()
        final_structure.info["energy_id"] = energy_id

        workflow["input_sample"] = [structure.info["id"]]
        workflow["output_sample"] = [new_id]
        workflow["calculated_property"] = [
            {
                "id": energy_id,
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

        if engine == "lammps":
            workflow["method"] = "MolecularStatics"
            workflow["degrees_of_freedom"] = [
                "AtomicPositionRelaxation",
                "CellVolumeRelaxation",
            ]
            workflow["interatomic_potential"] = {
                "potential_type": potential_type,
                "uri": potential_doi,
            }
            workflow["software"] = [
                {"uri": "https://doi.org/10.1016/j.cpc.2021.108171", "label": "LAMMPS"}
            ]
        elif engine == "qe":
            workflow["method"] = "DensityFunctionalTheory"
            workflow["degrees_of_freedom"] = [
                "AtomicPositionRelaxation",
                "CellVolumeRelaxation",
                "CellShapeRelaxation",
            ]
            workflow["xc_functional"] = xc_functional
            workflow["input_parameter"] = input_params
            workflow["software"] = [
                {"uri": "https://www.quantum-espresso.org/", "label": "QuantumEspresso"}
            ]

        cdict["workflow"].append(workflow.copy())

    return final_structure, ecoh, vol
