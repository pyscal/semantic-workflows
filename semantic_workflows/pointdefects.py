from typing import Optional
from ase.atoms import Atoms
from ase.build import bulk as ase_bulk
import numpy as np
from ase.data import atomic_numbers, chemical_symbols, reference_states
import numpy as np
import spglib
from collections import Counter
import random
import string
import copy
from conceptual_dictionary import sample_template as template_dict
from .build import update_attributes
from pyscal3 import System


def create_interstitial(
    atoms,
    element,
    void_type,
    a=None,
    threshold=0.01,
    cdict=None,
):
    if void_type == "tetrahedral":
        sys = System(atoms, format="ase")
        element = np.atleast_1d(element)
        sys.find.neighbors(method="voronoi", cutoff=0.1)
        verts = sys.unique_vertices
        randindex = np.random.randint(0, len(verts), len(element))
        randpos = np.array(verts)[randindex]

    elif void_type == "octahedral":
        # we need lattice constant, we can find this if it exists
        if cdict is not None:
            data = None
            id = atoms.info["id"]
            for d in cdict["computational_sample"]:
                if d["id"] == id:
                    data = d
                    break
            if data is not None:
                a = data["material"]["crystal_structure"]["unit_cell"][
                    "lattice_parameter"
                ][0]
        if a is None:
            raise ValueError("please provide lattice constant")
        cutoff = a + threshold * 2
        sys.find.neighbors(method="cutoff", cutoff=cutoff)
        octa_pos = []
        for count, dist in enumerate(sys.atoms.neighbors.distance):
            diffs = np.abs(np.array(dist) - a)
            indices = np.where(diffs < 1e-2)[0]
            vector = np.array(sys.atoms["diff"][count])[indices]
            vector = sys.atoms.positions[count] + vector / 2
            for vect in vector:
                vect = sys.modify.remap_position_to_box(vect)
                octa_pos.append(vect)

        octa_pos = np.unique(octa_pos, axis=0)
        randindex = np.random.randint(0, len(octa_pos), len(element))
        randpos = octa_pos[randindex]

        if not len(randpos) == len(element):
            raise ValueError("not enough octahedral positions found!")
    else:
        raise ValueError("Unknown void type")

    no_of_impurities = len(randpos)
    conc_of_impurities = no_of_impurities / sys.natoms

    sys = System(source=sys.add_atoms({"positions": randpos, "species": element}))
    new_atoms = sys.write.ase()

    # ok now we need to update things
    if cdict is not None:
        new_atoms.info["id"] = atoms.info["id"]

        # add defect
        defect_record = {
            "label": void_type,
            "concentration": conc_of_impurities,
            "number": no_of_impurities,
        }
        new_atoms = update_attributes(
            new_atoms, cdict, create_new=False, interstitial=defect_record
        )

        operation = {
            "method": "AddAtom",
            "input_sample": atoms.info["id"],
            "output_sample": new_atoms.info["id"],
        }
        cdict["operation"].append(operation)
    return new_atoms


def create_substitutional(
    atoms,
    element,
    cdict=None,
):
    species = atoms.get_chemical_symbols()
    randint = np.random.randint(len(species))
    species[randint] = element
    atoms.set_chemical_symbols(species)

    no_of_impurities = 1
    conc_of_impurities = no_of_impurities / len(atoms)

    # ok now we need to update things
    if cdict is not None:
        # new_atoms.info['id'] = atoms.info['id']

        # add defect
        defect_record = {
            "concentration": conc_of_impurities,
            "number": no_of_impurities,
        }
        new_atoms = update_attributes(
            atoms, cdict, create_new=False, substitutional=defect_record
        )

        operation = {
            "method": "SubstituteAtom",
            "input_sample": atoms.info["id"],
            "output_sample": new_atoms.info["id"],
        }
        cdict["operation"].append(operation)
    return atoms


def create_vacancy(atoms, index=None, cdict=None):
    """
    Remove one atom to create a vacancy supercell.

    Parameters
    ----------
    atoms : ase.Atoms
        Input structure (N atoms). Must have ``atoms.info["id"]`` set when
        ``cdict`` is provided.
    index : int, optional
        Index of the atom to remove.  If None, a random atom is chosen.
    cdict : ConceptualDict, optional
        If provided, the sample record is updated and an operation entry is
        appended.

    Returns
    -------
    ase.Atoms
        Defect structure with N-1 atoms.
    """
    if index is None:
        index = np.random.randint(len(atoms))

    new_atoms = atoms.copy()
    del new_atoms[index]

    no_of_vacancies = 1
    conc_of_vacancies = no_of_vacancies / len(atoms)

    if cdict is not None:
        new_atoms.info["id"] = atoms.info["id"]
        vacancy_record = {
            "concentration": conc_of_vacancies,
            "number": no_of_vacancies,
        }
        new_atoms = update_attributes(
            new_atoms, cdict, create_new=False, vacancy=vacancy_record
        )
        operation = {
            "method": "DeleteAtom",
            "input_sample": atoms.info["id"],
            "output_sample": new_atoms.info["id"],
        }
        cdict["operation"].append(operation)
    return new_atoms


def calculate_vacancy_formation_energy(
    bulk_structure,
    defect_structure,
    e_bulk,
    e_defect,
    cdict=None,
):
    """
    Calculate the vacancy formation energy.

    E_vac = (N-1) * (e_defect - e_bulk)

    where ``e_bulk`` and ``e_defect`` are per-atom energies returned by
    ``relax_structure``, and N-1 = len(defect_structure).

    Parameters
    ----------
    bulk_structure : ase.Atoms
        Relaxed bulk structure (N atoms).
    defect_structure : ase.Atoms
        Relaxed defect structure (N-1 atoms).
    e_bulk : float
        Per-atom energy (eV/atom) of the relaxed bulk.
    e_defect : float
        Per-atom energy (eV/atom) of the relaxed defect cell.
    cdict : ConceptualDict, optional
        If provided, math_operation provenance entries are appended.

    Returns
    -------
    float
        Vacancy formation energy in eV.
    """
    n_defect = len(defect_structure)
    e_vac = n_defect * (e_defect - e_bulk)

    if cdict is not None:
        from conceptual_dictionary import math_operation_template

        uid = cdict.generate_id()
        delta_id = f"delta_e_vac_{uid}"
        sample_id = defect_structure.info.get("id")
        # Prefer explicit link to the defect structure's atom count; fall back to scalar
        n_defect_ref = (
            f"{sample_id}.simulation_cell.number_of_atoms" if sample_id else n_defect
        )

        # Step 1: delta_e = e_defect - e_bulk  (eV/atom)
        op1 = copy.deepcopy(math_operation_template)
        op1["type"] = "Subtraction"
        op1["minuend"] = float(e_defect)
        op1["subtrahend"] = float(e_bulk)
        op1["result"] = {
            "id": delta_id,
            "label": "CohesiveEnergyDifference",
            "basename": "TotalEnergy",
            "value": float(np.round(e_defect - e_bulk, decimals=6)),
            "unit": "EV",
            "associate_to_sample": [],
        }
        cdict["math_operation"].append(op1)

        # Step 2: E_vac = (N-1) * delta_e
        op2 = copy.deepcopy(math_operation_template)
        op2["type"] = "Multiplication"
        op2["factor"] = [n_defect_ref, delta_id]
        op2["result"] = {
            "id": None,
            "label": "VacancyFormationEnergy",
            "basename": "FormationEnergy",
            "value": float(np.round(e_vac, decimals=4)),
            "unit": "EV",
            "associate_to_sample": [sample_id] if sample_id else [],
        }
        cdict["math_operation"].append(op2)

    return e_vac


def calculate_substitutional_formation_energy(
    bulk_structure,
    defect_structure,
    e_bulk,
    e_defect,
    e_ref,
    cdict=None,
):
    """
    Calculate the substitutional defect formation energy.

    E_sub = N * e_defect - (N-1) * e_bulk - e_ref

    Requires three simulations: the defect supercell, the host bulk, and the
    pure bulk of the substituting species.  All energies are per-atom values
    returned by ``relax_structure``.

    Parameters
    ----------
    bulk_structure : ase.Atoms
        Relaxed bulk structure of the host element (N atoms).
    defect_structure : ase.Atoms
        Relaxed defect structure (N-1 host + 1 substituting atom, N atoms total).
    e_bulk : float
        Per-atom energy (eV/atom) of the relaxed host bulk.
    e_defect : float
        Per-atom energy (eV/atom) of the relaxed defect supercell.
    e_ref : float
        Per-atom energy (eV/atom) of the substituting species in its pure bulk.
    cdict : ConceptualDict, optional
        If provided, math_operation provenance entries are appended.

    Returns
    -------
    float
        Substitutional formation energy in eV.
    """
    n_atoms = len(defect_structure)  # same as len(bulk_structure)
    e_sub = n_atoms * e_defect - (n_atoms - 1) * e_bulk - e_ref

    if cdict is not None:
        from conceptual_dictionary import math_operation_template

        uid = cdict.generate_id()
        e_tot_defect_id = f"e_tot_def_sub_{uid}"
        e_ref_host_id = f"e_ref_host_sub_{uid}"
        e_diff_id = f"e_diff_sub_{uid}"
        sample_id = defect_structure.info.get("id")
        # N links directly to the defect supercell atom count; (N-1) has no
        # matching structure so it stays as a scalar.
        n_defect_ref = (
            f"{sample_id}.simulation_cell.number_of_atoms" if sample_id else n_atoms
        )

        # Step 1: e_tot_defect = N * e_defect
        op1 = copy.deepcopy(math_operation_template)
        op1["type"] = "Multiplication"
        op1["factor"] = [n_defect_ref, float(e_defect)]
        op1["result"] = {
            "id": e_tot_defect_id,
            "label": "TotalEnergyDefect",
            "basename": "TotalEnergy",
            "value": float(np.round(n_atoms * e_defect, decimals=4)),
            "unit": "EV",
            "associate_to_sample": [],
        }
        cdict["math_operation"].append(op1)

        # Step 2: e_ref_host = (N-1) * e_bulk
        # N-1 is not the atom count of any structure in this calculation, so
        # it remains a plain scalar.
        op2 = copy.deepcopy(math_operation_template)
        op2["type"] = "Multiplication"
        op2["factor"] = [n_atoms - 1, float(e_bulk)]
        op2["result"] = {
            "id": e_ref_host_id,
            "label": "TotalEnergyHostReference",
            "basename": "TotalEnergy",
            "value": float(np.round((n_atoms - 1) * e_bulk, decimals=4)),
            "unit": "EV",
            "associate_to_sample": [],
        }
        cdict["math_operation"].append(op2)

        # Step 3: e_diff = e_tot_defect - e_ref_host
        op3 = copy.deepcopy(math_operation_template)
        op3["type"] = "Subtraction"
        op3["minuend"] = e_tot_defect_id
        op3["subtrahend"] = e_ref_host_id
        op3["result"] = {
            "id": e_diff_id,
            "label": "EnergyDifference",
            "basename": "TotalEnergy",
            "value": float(
                np.round(n_atoms * e_defect - (n_atoms - 1) * e_bulk, decimals=4)
            ),
            "unit": "EV",
            "associate_to_sample": [],
        }
        cdict["math_operation"].append(op3)

        # Step 4: E_sub = e_diff - e_ref
        op4 = copy.deepcopy(math_operation_template)
        op4["type"] = "Subtraction"
        op4["minuend"] = e_diff_id
        op4["subtrahend"] = float(e_ref)
        op4["result"] = {
            "id": None,
            "label": "SubstitutionalFormationEnergy",
            "basename": "FormationEnergy",
            "value": float(np.round(e_sub, decimals=4)),
            "unit": "EV",
            "associate_to_sample": [sample_id] if sample_id else [],
        }
        cdict["math_operation"].append(op4)

    return e_sub


def calculate_interstitial_formation_energy(
    bulk_structure,
    defect_structure,
    e_bulk,
    e_defect,
    e_ref,
    cdict=None,
):
    """
    Calculate the interstitial defect formation energy.

    E_int = N_defect * e_defect - N_bulk * e_bulk - e_ref

    Requires three simulations: the defect supercell (N_bulk + 1 atoms), the
    host bulk, and the pure bulk of the interstitial species.  All energies are
    per-atom values returned by ``relax_structure``.

    Parameters
    ----------
    bulk_structure : ase.Atoms
        Relaxed bulk structure of the host element (N_bulk atoms).
    defect_structure : ase.Atoms
        Relaxed defect structure (N_bulk + 1 atoms).
    e_bulk : float
        Per-atom energy (eV/atom) of the relaxed host bulk.
    e_defect : float
        Per-atom energy (eV/atom) of the relaxed defect supercell.
    e_ref : float
        Per-atom energy (eV/atom) of the interstitial species in its pure bulk.
    cdict : ConceptualDict, optional
        If provided, math_operation provenance entries are appended.

    Returns
    -------
    float
        Interstitial formation energy in eV.
    """
    n_bulk = len(bulk_structure)
    n_defect = len(defect_structure)
    e_int = n_defect * e_defect - n_bulk * e_bulk - e_ref

    if cdict is not None:
        from conceptual_dictionary import math_operation_template

        uid = cdict.generate_id()
        e_tot_defect_id = f"e_tot_def_int_{uid}"
        e_tot_bulk_id = f"e_tot_bulk_int_{uid}"
        delta_id = f"delta_e_int_{uid}"
        sample_id = defect_structure.info.get("id")
        bulk_id = bulk_structure.info.get("id")
        # Prefer dotpath references to link atom counts to their structures explicitly
        n_defect_ref = (
            f"{sample_id}.simulation_cell.number_of_atoms" if sample_id else n_defect
        )
        n_bulk_ref = f"{bulk_id}.simulation_cell.number_of_atoms" if bulk_id else n_bulk

        # Step 1: e_tot_defect = N_defect * e_defect
        op1 = copy.deepcopy(math_operation_template)
        op1["type"] = "Multiplication"
        op1["factor"] = [n_defect_ref, float(e_defect)]
        op1["result"] = {
            "id": e_tot_defect_id,
            "label": "TotalEnergyDefect",
            "basename": "TotalEnergy",
            "value": float(np.round(n_defect * e_defect, decimals=4)),
            "unit": "EV",
            "associate_to_sample": [],
        }
        cdict["math_operation"].append(op1)

        # Step 2: e_tot_bulk = N_bulk * e_bulk
        op2 = copy.deepcopy(math_operation_template)
        op2["type"] = "Multiplication"
        op2["factor"] = [n_bulk_ref, float(e_bulk)]
        op2["result"] = {
            "id": e_tot_bulk_id,
            "label": "TotalEnergyBulkReference",
            "basename": "TotalEnergy",
            "value": float(np.round(n_bulk * e_bulk, decimals=4)),
            "unit": "EV",
            "associate_to_sample": [],
        }
        cdict["math_operation"].append(op2)

        # Step 3: delta_e = e_tot_defect - e_tot_bulk
        op3 = copy.deepcopy(math_operation_template)
        op3["type"] = "Subtraction"
        op3["minuend"] = e_tot_defect_id
        op3["subtrahend"] = e_tot_bulk_id
        op3["result"] = {
            "id": delta_id,
            "label": "TotalEnergyDifference",
            "basename": "TotalEnergy",
            "value": float(np.round(n_defect * e_defect - n_bulk * e_bulk, decimals=4)),
            "unit": "EV",
            "associate_to_sample": [],
        }
        cdict["math_operation"].append(op3)

        # Step 4: E_int = delta_e - e_ref
        op4 = copy.deepcopy(math_operation_template)
        op4["type"] = "Subtraction"
        op4["minuend"] = delta_id
        op4["subtrahend"] = float(e_ref)
        op4["result"] = {
            "id": None,
            "label": "InterstitialFormationEnergy",
            "basename": "FormationEnergy",
            "value": float(np.round(e_int, decimals=4)),
            "unit": "EV",
            "associate_to_sample": [sample_id] if sample_id else [],
        }
        cdict["math_operation"].append(op4)

    return e_int
