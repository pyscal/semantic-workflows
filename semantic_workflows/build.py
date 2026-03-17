"""
Core structure building functions.

These are framework-agnostic implementations that can be decorated
for use with pyiron_workflow or jobflow.
"""

from typing import Optional
from ase.atoms import Atoms
from ase.build import bulk as ase_bulk
import numpy as np
from ase.data import atomic_numbers, chemical_symbols, reference_states
import spglib
from collections import Counter
import random
import string
import copy


def generate_id(length=7):
    """Generate a random alphanumeric ID of given length.

    Uses ``os.urandom`` rather than Python's ``random`` module so that
    third-party libraries that call ``random.seed()`` (e.g. PyIron/LAMMPS
    wrappers) cannot cause ID collisions across iterations.
    """
    import os
    chars = string.ascii_letters + string.digits  # A–Z, a–z, 0–9
    return "".join(chars[b % 62] for b in os.urandom(length))


def bulk(
    name: str,
    crystalstructure: Optional[str] = None,
    a: Optional[float | int] = None,
    c: Optional[float | int] = None,
    covera: Optional[float] | int = None,
    u: Optional[float | int] = None,
    orthorhombic: bool = False,
    cubic: bool = False,
    cdict=None,
):
    """
    Create a bulk crystal structure.

    Parameters
    ----------
    name : str
        Chemical symbol
    crystalstructure : str, optional
        Crystal structure type
    a : float, optional
        Lattice constant
    c : float, optional
        c lattice constant
    covera : float, optional
        c/a ratio
    u : float, optional
        Internal parameter
    orthorhombic : bool
        Create orthorhombic cell
    cubic : bool
        Create cubic cell
    cdict : dict, optional
        Knowledge graph for metadata (pyiron-specific)

    Returns
    -------
    ase.Atoms
        Bulk structure
    """
    sdict = _compute_structure_metadata(name, crystalstructure, a, a, c, covera)

    struct = ase_bulk(
        name,
        crystalstructure=crystalstructure,
        a=a,
        c=c,
        covera=covera,
        u=u,
        orthorhombic=orthorhombic,
        cubic=cubic,
    )
    sdict["spacegroup_symbol"] = get_spacegroup_symbol(struct)
    sdict["spacegroup_number"] = get_spacegroup_number(struct)

    if cdict is not None:
        from conceptual_dictionary import sample_template as template_dict

        data = _generate_atomic_sample_data(
            struct, sdict, repeat=None, template_dict=template_dict
        )
        id = generate_id()
        data["id"] = id
        cdict["computational_sample"].append(data)
        struct.info["id"] = id

    return struct


def repeat(
    structure: Atoms,
    repetitions,
    cdict=None,
) -> Atoms:
    """
    Repeat structure.

    Parameters
    ----------
    structure : ase.Atoms
        Input structure
    repetitions : tuple or int
        Repetition factors
    cdict : dict, optional
        Knowledge graph

    Returns
    -------
    ase.Atoms
        Repeated structure
    """
    structure = structure.repeat(repetitions)
    if cdict is not None:
        structure = update_attributes(structure, cdict, repeat=list(repetitions))
    return structure


def polycrystal(structure: Atoms, box_size, grain_size, cdict=None) -> Atoms:
    """
    Create polycrystalline structure.

    Parameters
    ----------
    structure : ase.Atoms
        Input structure
    box_size : tuple
        Box dimensions
    grain_size : float
        Grain size
    cdict : dict, optional
        Knowledge graph

    Returns
    -------
    ase.Atoms
        Polycrystalline structure
    """
    import os
    from ase.io import read
    from ase.io import write

    n_grains = int((box_size[0] * box_size[1] * box_size[2]) / (grain_size**3))

    # calculate repeats
    nx = np.ceil(structure.get_cell()[0] / box_size[0])
    ny = np.ceil(structure.get_cell()[1] / box_size[1])
    nz = np.ceil(structure.get_cell()[2] / box_size[2])

    with open("grain_sizes.txt", "w") as f:
        f.write(f"box {box_size[0]} {box_size[1]} {box_size[2]}\n")
        f.write(f"random {n_grains}\n")
    write("tmp.xsf", structure)

    os.system("atomsk --polycrystal tmp.xsf grain_sizes.txt final.cfg -ow")

    poly_struct = read("final.cfg")
    os.remove("grain_sizes.txt")
    os.remove("tmp.xsf")
    os.remove("final.cfg")

    if cdict is not None:
        poly_struct.info["id"] = structure.info["id"]
        # update system
        poly_struct = update_attributes(
            poly_struct,
            cdict,
            repeat=(nx, ny, nz),
            grain_size=grain_size,
            number_of_grains=n_grains,
            ignore_positions=True,
        )

    return poly_struct


# DATADICT properties
# ------------------------------------------
bravais_lattice_dict = {
    "b2": "https://www.wikidata.org/wiki/Q851536",
    "bcc": "https://www.wikidata.org/wiki/Q851536",
    "bct": "bct",
    "cubic": "https://www.wikidata.org/wiki/Q473227",
    "diamond": "https://www.wikidata.org/wiki/Q3006714",
    "diatom": "https://www.wikidata.org/wiki/Q1751859",
    "fcc": "https://www.wikidata.org/wiki/Q3006714",
    "hcp": "https://www.wikidata.org/wiki/Q663314",
    "monoclinic": "https://www.wikidata.org/wiki/Q624543",
    "orthorhombic": "https://www.wikidata.org/wiki/Q648961",
    "rhombohedral": "https://www.wikidata.org/wiki/Q13362463",
    "sc": "https://www.wikidata.org/wiki/Q2242450",
    "tetragonal": "https://www.wikidata.org/wiki/Q503601",
    "l12": "https://www.wikidata.org/wiki/Q3006714",
    "a15": "a15",
}


# SIMCELL properties
# --------------------------------------------
def get_chemical_composition(structure):
    """Get the chemical composition of the system."""
    symbols = structure.get_chemical_symbols()
    element_counts = Counter(symbols)
    total_atoms = len(symbols)
    composition = {
        element: count / total_atoms for element, count in element_counts.items()
    }
    return composition


def get_cell_volume(system):
    """Get the volume of the simulation cell."""
    return system.get_volume()


def get_number_of_atoms(system):
    """Get the number of atoms in the system."""
    return len(system)


def get_simulation_cell_length(system):
    """Get the length of the simulation cell."""
    return system.get_cell_lengths_and_angles()[:3]


def get_simulation_cell_vector(system):
    """Get the simulation cell vector of the given system."""
    return system.cell


def get_simulation_cell_angle(system):
    """Get the angles between the vectors of the simulation cell."""
    return system.get_cell_lengths_and_angles()[3:]


# LATTICE properties
# --------------------------------------------
def get_bravais_lattice(structure):
    """Get the Bravais lattice of a given system."""
    if structure in bravais_lattice_dict.keys():
        return bravais_lattice_dict[structure]
    return None


def get_spacegroup_symbol(system):
    """Get the symbol of the spacegroup for a given system."""
    try:
        results = _get_symmetry_dict(system)
        return results[0]
    except:
        return None


def get_spacegroup_number(system):
    """Get the spacegroup number of a given system."""
    try:
        results = _get_symmetry_dict(system)
        return results[1]
    except:
        return None


# SUPPORT functions
# --------------------------------------------
def _get_angle(vec1, vec2):
    """Get angle between two vectors in degrees."""
    return np.round(
        np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)))
        * 180
        / np.pi,
        decimals=2,
    )


def _get_symmetry_dict(system):
    """Get symmetry information using spglib."""
    box = system.get_cell()
    direct_coordinates = system.get_scaled_positions()
    symbols = system.get_chemical_symbols()
    atom_types = [list(dict.fromkeys(symbols).keys()).index(s) + 1 for s in symbols]

    results = spglib.get_symmetry_dataset((box, direct_coordinates, atom_types))
    return results.international, results.number


def _generate_atomic_sample_data(atoms, sdict=None, repeat=None, template_dict=None):
    """Generate atomic sample metadata."""
    data = template_dict.copy()
    data["material"]["element_ratio"] = get_chemical_composition(atoms)

    if sdict is not None:
        if "structure" in sdict.keys():
            data["material"]["crystal_structure"]["name"] = sdict["structure"]
        if "spacegroup_symbol" in sdict.keys():
            data["material"]["crystal_structure"]["spacegroup_symbol"] = sdict[
                "spacegroup_symbol"
            ]
        if "spacegroup_number" in sdict.keys():
            data["material"]["crystal_structure"]["spacegroup_number"] = sdict[
                "spacegroup_number"
            ]
        if "structure" in sdict.keys():
            data["material"]["crystal_structure"]["unit_cell"]["bravais_lattice"] = (
                get_bravais_lattice(sdict["structure"])
            )
        if "a" in sdict.keys():
            if "b" not in sdict.keys():
                sdict["b"] = sdict["a"]
            if "c" not in sdict.keys():
                sdict["c"] = sdict["a"]
            data["material"]["crystal_structure"]["unit_cell"]["lattice_parameter"] = [
                sdict["a"],
                sdict["b"],
                sdict["c"],
            ]

    data["material"]["crystal_structure"]["unit_cell"][
        "angle"
    ] = atoms.get_cell_lengths_and_angles()[3:]

    data["simulation_cell"]["volume"]["value"] = get_cell_volume(atoms)
    data["simulation_cell"]["number_of_atoms"] = get_number_of_atoms(atoms)
    data["simulation_cell"]["length"] = get_simulation_cell_length(atoms)
    data["simulation_cell"]["vector"] = get_simulation_cell_vector(atoms)
    data["simulation_cell"]["angle"] = get_simulation_cell_angle(atoms)

    if repeat is not None:
        if isinstance(repeat, int):
            data["simulation_cell"]["repetitions"] = (repeat, repeat, repeat)
        else:
            data["simulation_cell"]["repetitions"] = repeat

    data["atom_attribute"]["position"] = atoms.get_positions().tolist()
    data["atom_attribute"]["species"] = atoms.get_chemical_symbols()
    return data


def _compute_structure_metadata(name, crystalstructure, a, b, c, covera):
    """Compute structure metadata from parameters."""
    sdict = {"a": a, "b": b, "c": c, "covera": covera}
    atomic_number = atomic_numbers.get(name)
    ref = reference_states[atomic_number]

    xref = None
    if ref:
        xref = ref.get("symmetry")
        if xref and name in chemical_symbols:
            sdict["structure"] = xref

    if crystalstructure:
        sdict["structure"] = crystalstructure

    if a is None and ref and "a" in ref:
        sdict["a"] = ref["a"]

    if b is None and ref and (bovera := ref.get("b/a")) and a:
        sdict["b"] = bovera * a

    if crystalstructure in ["hcp", "wurtzite"]:
        if c:
            covera = c / a
        elif covera is None:
            covera = ref.get("c/a") if xref == crystalstructure else np.sqrt(8 / 3)

    if covera is None and ref and (ref_c_a := ref.get("c/a")):
        covera = ref_c_a
        if c is None and a:
            sdict["c"] = covera * a

    sdict["b"] = sdict["a"] if sdict["b"] is None else sdict["b"]
    sdict["c"] = sdict["a"] if sdict["c"] is None else sdict["c"]

    return sdict


def update_attributes(
    atoms,
    cdict,
    repeat=None,
    create_new=False,
    grain_size=None,
    number_of_grains=0,
    ignore_positions=False,
    interstitial=None,
    substitutional=None,
    vacancy=None,
):
    """
    Update the atom attributes based on the provided ASE Atoms object.
    This would also reset the id, since the structure has changed.
    """
    id = atoms.info["id"]

    # find data
    for d in cdict["computational_sample"]:
        if d["id"] == id:
            data = d

    if create_new:
        data = copy.deepcopy(data)
        data["id"] = generate_id()
        atoms = atoms.copy()
        atoms.info["id"] = data["id"]

    data["material"]["element_ratio"] = get_chemical_composition(atoms)
    data["simulation_cell"]["volume"]["value"] = get_cell_volume(atoms)
    data["simulation_cell"]["number_of_atoms"] = get_number_of_atoms(atoms)
    data["simulation_cell"]["length"] = get_simulation_cell_length(atoms)
    data["simulation_cell"]["vector"] = get_simulation_cell_vector(atoms)
    data["simulation_cell"]["angle"] = get_simulation_cell_angle(atoms)

    if grain_size is not None:
        data["simulation_cell"]["grain_size"] = grain_size
    if number_of_grains is not None:
        data["simulation_cell"]["number_of_grains"] = number_of_grains

    if repeat is not None:
        if isinstance(repeat, int):
            data["simulation_cell"]["repetitions"] = (repeat, repeat, repeat)
        else:
            data["simulation_cell"]["repetitions"] = repeat

    if not ignore_positions:
        data["atom_attribute"]["position"] = atoms.get_positions().tolist()
        data["atom_attribute"]["species"] = atoms.get_chemical_symbols()
    else:
        data["atom_attribute"]["position"] = []
        data["atom_attribute"]["species"] = []

    if interstitial is not None:
        data["interstitial"] = interstitial
    if substitutional is not None:
        data["substitutional"] = substitutional
    if vacancy is not None:
        data["vacancy"] = vacancy

    if create_new:
        cdict["computational_sample"].append(data)

    return atoms
