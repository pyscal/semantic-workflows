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
from .templates import sample_template as template_dict 
from .build import update_attributes
from pyscal3 import System

@as_function_node
def create_interstitial(
    atoms,
    element,
    void_type,
    a=None,
    threshold=0.01,
    kg=None,
):
    if void_type == 'tetrahedral':
        sys = System(atoms, format='ase')
        element = np.atleast_1d(element)
        sys.find.neighbors(method="voronoi", cutoff=0.1)
        verts = sys.unique_vertices
        randindex = np.random.randint(0, len(verts), len(element))
        randpos = np.array(verts)[randindex]
    
    elif void_type == 'octahedral':
        #we need lattice constant, we can find this if it exists
        if kg is not None:
            data = None
            id = atoms.info['id']
            for d in kg['computational_sample']:
                if d['id'] == id:
                    data = d
                    break
            if data is not None:
                a = data["material"]["crystal_structure"]["unit_cell"]["lattice_parameter"][0]
        if a is None:
            raise ValueError('please provide lattice constant')
        cutoff = a + threshold * 2
        self.find.neighbors(method="cutoff", cutoff=cutoff)
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
        raise ValueError('Unknown void type')

    no_of_impurities = len(randpos)
    conc_of_impurities = no_of_impurities/sys.natoms
    
    sys = System(source=sys.add_atoms({"positions": randpos, "species": element}))
    new_atoms = sys.write.ase()
    
    #ok now we need to update things
    if kg is not None:
        new_atoms.info['id'] = atoms.info['id']

        #add defect
        defect_record = {
            'label': void_type,
            'concentration': conc_of_impurities,
            'number': no_of_impurities,
        }
        new_atoms = update_attributes(new_atoms, kg, create_new=False,
                                     interstitial=defect_record)

        #ok good, now we need to add activity
        #for the moment we will just create defected structures, and thats it
    return new_atoms

@as_function_node
def create_substitutional(    
    atoms,
    element,
    kg = None,
):
    species = atoms.get_chemical_symbols()
    randint = np.random.randint(len(species))
    species[randint] = element
    atoms.set_chemical_symbols(species)

    no_of_impurities = 1
    conc_of_impurities = no_of_impurities/len(atoms)

    #ok now we need to update things
    if kg is not None:
        #new_atoms.info['id'] = atoms.info['id']

        #add defect
        defect_record = {
            'concentration': conc_of_impurities,
            'number': no_of_impurities,
        }
        new_atoms = update_attributes(atoms, kg, create_new=False,
                                     substitutional=defect_record)

        #ok good, now we need to add activity
        #for the moment we will just create defected structures, and thats it
    return atoms
    
    
    
        

    

