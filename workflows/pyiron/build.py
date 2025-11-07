from typing import Optional
from ase.atoms import Atoms
from ase.build import bulk as ase_bulk
from pyiron_workflow import as_function_node 
import numpy as np
from ase.data import atomic_numbers, chemical_symbols, reference_states
import numpy as np
import spglib
from collections import Counter
import random
import string
import copy

def generate_id(length=7):
    """Generate a random alphanumeric ID of given length."""
    chars = string.ascii_letters + string.digits  # A–Z, a–z, 0–9
    return ''.join(random.choices(chars, k=length))

template_dict = {
    "id": "sample1",
    "material": {
        "element_ratio": {},
        "crystal_structure": {
            "spacegroup_symbol": None,
            "spacegroup_number": None,
            "unit_cell": {
                "bravais_lattice": None,
                "lattice_parameter": None,
                "angle": []
            }
        }
    },
    "simulation_cell": {
        "volume": {'value': None},
        "number_of_atoms": None,
        "length": [],
        "vector": [],
        "angle": [],
        "repetitions": []
    },
    'atom_attribute': {
        'position': None,
        'species': None,
    }
}

@as_function_node
def bulk(
    name: str,
    crystalstructure: Optional[str] = None,
    a: Optional[float | int] = None,
    c: Optional[float | int] = None,
    covera: Optional[float] | int = None,
    u: Optional[float | int] = None,
    orthorhombic: bool = False,
    cubic: bool = False,
    kg = None,
):
    sdict = _compute_structure_metadata(name, crystalstructure, a, a, c, covera)
    
    struct = ase_bulk(
        name,
        crystalstructure=crystalstructure,
        a=a,
        c=c,
        covera = covera,
        u = u,
        orthorhombic = orthorhombic,
        cubic = cubic,
    )
    sdict["spacegroup_symbol"] = get_spacegroup_symbol(struct)
    sdict["spacegroup_number"] = get_spacegroup_number(struct)
    data = _generate_atomic_sample_data(struct, sdict, repeat)
    id = generate_id()
    data['id'] = id
    
    #add to kg if needed
    if kg is not None:
        kg['computational_sample'].append(data)
    struct.info['id'] = id
    return struct


#@as_function_node
#def change_composition(
#    structure: Atoms,
#    species: str,
#    fraction: float,
#):
#
#    # Copy structure
#    new_structure = structure.copy()
#    symbols = np.array(new_structure.get_chemical_symbols())
#
#    # Determine number of atoms to replace
#    n_atoms = len(symbols)
#    n_replace = int(round(fraction * n_atoms))
#    
#    # Randomly select unique indices
#    replace_indices = np.random.choice(n_atoms, n_replace, replace=False)
#
#    # Perform replacement
#    for i in replace_indices:
#        symbols[i] = species
#
#    new_structure.set_chemical_symbols(symbols)
#    return new_structure

@as_function_node
def repeat(
    structure: Atoms,
    repetitions: tuple[int, int, int],
    kg = None,
) -> Atoms:

    structure = structure.repeat(repetitions)
    if kg is not None:
        structure = update_attributes(structure, kg, repeat=repetitions)
    return structure

@as_function_node
def polycrystal(structure: Atoms, box_size: tuple[float, float, float], grain_size: float) -> Atoms:
    import numpy as np
    import os
    from ase.io import read
    from ase.io import write 

    n_grains = int((box_size[0]*box_size[1]*box_size[2])/(grain_size**3))
    with open('grain_sizes.txt', 'w') as f:
        f.write(f"box {box_size[0]} {box_size[1]} {box_size[2]}\n")
        f.write(f"random {n_grains}\n")
    write('tmp.xsf', structure)


    os.system('/u/smenon/ptmp/FeC/compression/atomsk --polycrystal tmp.xsf grain_sizes.txt final.cfg -ow')
    
    poly_struct = read('final.cfg')
    os.remove('grain_sizes.txt')
    os.remove('tmp.xsf')
    os.remove('final.cfg')
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
    """
    Get the chemical composition of the system.

    Parameters
    ----------
    system : object
        The system object.

    Returns
    -------
    composition : dict
        A dictionary containing the chemical elements as keys and their corresponding counts as values.

    """
    symbols = structure.get_chemical_symbols()
    element_counts = Counter(symbols)
    total_atoms = len(symbols)
    composition = {
        element: count / total_atoms for element, count in element_counts.items()
    }
    return composition


def get_cell_volume(system):
    """
    Get the volume of the simulation cell.

    Parameters
    ----------
    system : object
        The system object.

    Returns
    -------
    volume : float
        The volume of the simulation cell.

    """
    return system.get_volume()


def get_number_of_atoms(system):
    """
    Get the number of atoms in the system.

    Parameters
    ----------
    system : object
        The system object.

    Returns
    -------
    natoms : int
        The number of atoms in the system.

    """
    return len(system)


def get_simulation_cell_length(system):
    """
    Get the length of the simulation cell.

    Parameters
    ----------
    system : object
        The system object.

    Returns
    -------
    length : list
        A list containing the length of each dimension of the simulation cell.

    """
    return system.get_cell_lengths_and_angles()[:3]


def get_simulation_cell_vector(system):
    """
    Get the simulation cell vector of the given system.

    Parameters
    ----------
    system : object
        The system object containing the simulation cell information.

    Returns
    -------
    numpy.ndarray
        The simulation cell vector of the system.

    """
    return system.cell


def get_simulation_cell_angle(system):
    """
    Get the angles between the vectors of the simulation cell.

    Parameters
    ----------
    system : object
        The system object containing the simulation cell information.

    Returns
    -------
    angles : list
        A list containing the angles between the vectors of the simulation cell.

    """
    return system.get_cell_lengths_and_angles()[3:]


# LATTICE properties
# --------------------------------------------


def get_lattice_angle(system):
    """
    Calculate the lattice angles of a given system.

    Parameters
    ----------
    system : object
        The system object containing the structure information.

    Returns
    -------
    list
        A list of three lattice angles in degrees. If the structure information is not available, [None, None, None] is returned.

    """
    if system._structure_dict is None:
        return [None, None, None]
    if "box" in system._structure_dict.keys():
        return [
            _get_angle(
                system._structure_dict["box"][0], system._structure_dict["box"][1]
            ),
            _get_angle(
                system._structure_dict["box"][1], system._structure_dict["box"][2]
            ),
            _get_angle(
                system._structure_dict["box"][2], system._structure_dict["box"][0]
            ),
        ]
    else:
        return [None, None, None]


def get_lattice_parameter(system):
    """
    Calculate the lattice parameters of a system.

    Parameters
    ----------
    system : object
        The system object containing information about the atoms and structure.

    Returns
    -------
    list
        A list containing the lattice parameters of the system. If the lattice constant is not available,
        [None, None, None] is returned. If the system structure is available, the lattice parameters are
        calculated based on the box dimensions. Otherwise, the lattice constant is returned for all three
        dimensions.

    Examples
    --------
    >>> system = System()
    >>> system.atoms._lattice_constant = 3.5
    >>> system._structure_dict = {"box": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}
    >>> get_lattice_parameter(system)
    [3.5, 3.5, 3.5]

    >>> system.atoms._lattice_constant = None
    >>> get_lattice_parameter(system)
    [None, None, None]
    """

    if system.atoms._lattice_constant is None:
        return [None, None, None]
    else:
        if system._structure_dict is not None:
            if "box" in system._structure_dict.keys():
                return [
                    np.linalg.norm(system._structure_dict["box"][0])
                    * system.atoms._lattice_constant,
                    np.linalg.norm(system._structure_dict["box"][1])
                    * system.atoms._lattice_constant,
                    np.linalg.norm(system._structure_dict["box"][2])
                    * system.atoms._lattice_constant,
                ]
        return [
            system.atoms._lattice_constant,
            system.atoms._lattice_constant,
            system.atoms._lattice_constant,
        ]


def get_crystal_structure_name(system):
    """
    Get the name of the crystal structure for a given system.

    Parameters
    ----------
    system : object
        The system object containing the crystal structure information.

    Returns
    -------
    str or None
        The name of the crystal structure if available, otherwise None.

    """
    if system._structure_dict is None:
        return None
    return system.atoms._lattice


def get_repetitions(system):
    if system._structure_dict is None:
        return [None, None, None]
    if "repetitions" in system._structure_dict.keys():
        return system._structure_dict["repetitions"]
    return [None, None, None]


def get_bravais_lattice(structure):
    """
    Get the Bravais lattice of a given system.

    Parameters
    ----------
    system : object
        The system object for which the Bravais lattice is to be determined.

    Returns
    -------
    str or None
        The Bravais lattice of the system, or None if the system's structure dictionary is not available or the lattice is not found in the dictionary.

    """
    if structure in bravais_lattice_dict.keys():
        return bravais_lattice_dict[structure]
    return None


def get_basis_positions(system):
    """
    Get the basis positions from the given system.

    Parameters
    ----------
    system : object
        The system object containing the structure dictionary.

    Returns
    -------
    numpy.ndarray or None
        The basis positions if available, otherwise None.
    """
    if system._structure_dict is None:
        return None
    if "positions" in system._structure_dict.keys():
        return system._structure_dict["positions"]
    return None


# def get_basis_occupancy(system):
#    if system._structure_dict is None:
#        return None

#    if "species" in system._structure_dict.keys():
#        occ_numbers = system._structure_dict['species']
#        tdict = system.atoms._type_dict
#        vals = [val for key, val in tdict.items()]

#        if vals[0] is not None:
#            occ_numbers = [tdict[x] for x in occ_numbers]
#        return occ_numbers
#    return None


def get_lattice_vector(system):
    """
    Get the lattice vector of a system.

    Parameters
    ----------
    system : object
        The system object containing the structure information.

    Returns
    -------
    list
        A list representing the lattice vector of the system. If the structure
        dictionary is not available or the lattice vector is not defined, it
        returns [None, None, None].
    """
    if system._structure_dict is None:
        return [None, None, None]
    if "box" in system._structure_dict.keys():
        return system._structure_dict["box"]
    return [None, None, None]


def get_spacegroup_symbol(system):
    """
    Get the symbol of the spacegroup for a given system.

    Parameters:
        system (object): The system object for which to retrieve the spacegroup symbol.

    Returns:
        str: The symbol of the spacegroup if available, otherwise None.
    """
    try:
        results = _get_symmetry_dict(system)
        return results[0]
    except:
        return None


def get_spacegroup_number(system):
    """
    Get the spacegroup number of a given system.

    Parameters
    ----------
    system : object
        The system object for which the spacegroup number is to be determined.

    Returns
    -------
    int or None
        The spacegroup number of the system if it is available, otherwise None.
    """
    try:
        results = _get_symmetry_dict(system)
        return results[1]
    except:
        return None


# ATOM attributes
# --------------------------------------------
def get_position(system):
    """
    Get the positions of the atoms in the system.

    Parameters
    ----------
    system : object
        The system object containing the atom positions.

    Returns
    -------
    numpy.ndarray or None
        The positions of the atoms if available, otherwise None.

    """
    return system.atoms.positions


def get_species(system):
    """
    Get the species of atoms in the given system.

    Parameters
    ----------
    system : System
        The system object containing atoms.

    Returns
    -------
    list
        A list of species of atoms in the system.

    """
    return system.atoms.species


# SUPPORT functions
# --------------------------------------------
def _get_angle(vec1, vec2):
    """
    Get angle between two vectors in degrees

    Parameters
    ----------
    vec1: list
        first vector

    vec2: list
        second vector

    Returns
    -------
    angle: float
        angle in degrees

    Notes
    -----
    Angle is rounded to two decimal points

    """
    return np.round(
        np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)))
        * 180
        / np.pi,
        decimals=2,
    )


def _get_symmetry_dict(system):
    box = system.get_cell()
    direct_coordinates = system.get_scaled_positions()
    symbols = system.get_chemical_symbols()
    atom_types = [list(dict.fromkeys(symbols).keys()).index(s) + 1 for s in symbols]

    results = spglib.get_symmetry_dataset((box, direct_coordinates, atom_types))
    return results.international, results.number

def _generate_atomic_sample_data(atoms, sdict=None, repeat=None):
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

def update_attributes(atoms, kg, repeat=None, create_new=False):
    """
    Update the atom attributes based on the provided ASE Atoms object.
    This would also reset the id, since the structure has changed.
    """
    id = atoms.info['id']

    #find data
    for d in kg['computational_sample']:
        if d['id'] == id:
            data = d

    if create_new:
        data = copy.deepcopy(data)
        data['id'] = generate_id()
        atoms = atoms.copy()
        atoms.info['id'] = data['id']
    
    data["material"]["element_ratio"] = get_chemical_composition(atoms)
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

    if create_new:
        kg['computational_sample'].append(data)
        
    return atoms