import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from ase.io import read, write
from pylammpsmpi import LammpsLibrary
from pyiron_workflow import as_function_node, as_macro_node

@as_function_node
def calculate_ev_curves(structure, pair_style, pair_coeff, vol_range=0.3, num_of_points=5, cores=1, 
                    e_tol=0, f_tol=0.0001, n_energy_steps=1e5, n_force_steps=1e6):
    
    #relax structure first
    #vol = relax_structure(structure, pair_style, pair_coeff, cores=cores, 
    #                e_tol=e_tol, f_tol=f_tol, 
    #                n_energy_steps=n_energy_steps, n_force_steps=n_force_steps)
    volume_factors = np.linspace((1 - vol_range), (1.0 + vol_range), num_of_points)
    
    energies = []
    volumes = []
    for volume_factor in volume_factors:
        scaled_atoms = scale_atoms(structure, volume_factor)
        e, v = calculate_energy(scaled_atoms, pair_style, pair_coeff, cores=cores, 
                    e_tol=e_tol, f_tol=f_tol, 
                    n_energy_steps=n_energy_steps, n_force_steps=n_force_steps)
        
        energies.append(e)
        volumes.append(v)
    V0, E0, B0, Bp = fit_bm(volumes, energies)
    v_fit = np.linspace(min(volumes), max(volumes), 100, endpoint=True)
    e_fit = birch_murnaghan_eval(v_fit, V0, E0, B0, Bp)
    bulk_modulus = B0*160.2176621
    return v_fit, e_fit, bulk_modulus


def scale_atoms(structure, scale_factor):
    scaled_atoms = structure.copy()
    scaled_atoms.set_cell(scaled_atoms.get_cell() * scale_factor, scale_atoms=True)
    return scaled_atoms
    

def fit_bm(vol, en):
    a, b, c = np.polyfit(vol, en, 2)
    V0 = -b/(2*a)
    E0 = a*V0**2 + b*V0 + c
    B0 = 2*a*V0
    Bp = 4.0
    popt, pcov = curve_fit(birch_murnaghan_eval, vol, en, p0=[V0, E0, B0, Bp])
    return popt


def birch_murnaghan_eval(vol, V0, E0, B0, Bp):    
    eta = (vol/V0)**(1.0/3.0)
    E = E0 + 9.0*B0*V0/16.0 * (eta**2-1.0)**2 * (6.0 + Bp*(eta**2-1.0) - 4.0*eta**2)
    return E
    

def relax_structure(structure, pair_style, pair_coeff, cores=1, 
                    e_tol=0, f_tol=0.0001, n_energy_steps=1e5, n_force_steps=1e6):
    write('tmp.data', structure, format='lammps-data')
    lmp = LammpsLibrary(cores=1)
    lmp.command("units metal") 
    lmp.command("dimension 3") 
    lmp.command("boundary p p p") 
    lmp.command("atom_style atomic")
    lmp.command("read_data tmp.data")
    lmp.command(f"pair_style {pair_style}")
    lmp.command(f"pair_coeff {pair_coeff}")
    lmp.command("fix ensemble all box/relax aniso 0.0")
    lmp.command("minimize {e_tol} {f_tol} {n_energy_steps} {n_force_steps}")    
    lmp.command("run 0")
    vol = lmp.vol/lmp.natoms
    return vol  


def calculate_energy(structure, pair_style, pair_coeff, cores=1, 
                    e_tol=0, f_tol=0.0001, n_energy_steps=1e5, n_force_steps=1e6):
    write('tmp.data', structure, format='lammps-data')
    lmp = LammpsLibrary(cores=1)
    lmp.command("units metal") 
    lmp.command("dimension 3") 
    lmp.command("boundary p p p") 
    lmp.command("atom_style atomic")
    lmp.command("read_data tmp.data")
    lmp.command(f"pair_style {pair_style}")
    lmp.command(f"pair_coeff {pair_coeff}")
    lmp.command("run 0")
    ecoh = lmp.pe/lmp.natoms
    vol = lmp.vol/lmp.natoms
    return ecoh, vol

