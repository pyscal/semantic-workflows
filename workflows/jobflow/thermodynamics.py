import pandas as pd
from ase.io import read, write
from pylammpsmpi import LammpsLibrary
import scipy.constants as sc
import numpy as np
from jobflow import job
from .utils import ase_pmg_bridge

@job
@ase_pmg_bridge
def calculate_thermal_properties(structure, pair_style, pair_coeff, 
                                temperature, pressure=0, cores=1,
                                n_equilibration_steps=15000,
                                n_run_steps=25000):
    
    ev_to_j = sc.physical_constants["electron volt-joule relationship"][0]
    Av =  sc.physical_constants["Avogadro constant"][0]
    A_to_m = 1E-10
    kB = sc.physical_constants["Boltzmann constant"][0]

    write('tmp.data', structure, format='lammps-data')

    from mendeleev import element

    symbols, counts = list(np.unique(structure.get_chemical_symbols(), return_counts=True))
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
    lmp.fix("extra all print 100 \" $(pe) $(ke) $(etotal) $(press) $(vol) $(temp)\" file data.dat screen no")

    lmp.run(n_run_steps)

    pe, ke, etotal, press, vol, temp = np.loadtxt("data.dat", unpack=True)
    efluct = etotal - np.mean(etotal)

    efluctsq = (efluct*ev_to_j)**2
    cp = np.mean(efluctsq)/(kB*T*T)
    w = 0
    for m, c in zip(masses, counts):
        w += (c/Av)*m*1E3

    cp = (cp/w)/1E-3
    
    crossfluct = (efluct*ev_to_j)*(vol - np.mean(vol))
    ap = np.mean(crossfluct)/(kB*T*T*np.mean(vol)) #1/K

    results = {'specific_heat': cp, 'thermal_expansion': ap, 'volume': np.mean(vol),}
    return results

def calphy_default_input():
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
        }
    }


def _write_structure(structure, filename, elements):
    from ase.io import write
    from mendeleev import element

    write(filename, structure, format="lammps-data", specorder=elements)
    masses = [element(symbol).mass for symbol in elements]
    return filename, masses


def _prepare_input(pair_style, pair_coeff, structure, 
                   temperature, pressure=0,
                   n_switching_steps=10000,
                   mode='fe', reference_phase='solid',
                   cores =1, 
                   input_dict=None):
    from calphy.input import Calculation
    import os
    import uuid
    import numpy as np

    if input_dict is None:
        input_dict = calphy_default_input()

    filename = f"tmp-{uuid.uuid4().hex[:8]}.data"
    elements = list(np.unique(structure.get_chemical_symbols()))
    #now ensure this is present in the pair_coeff
    element_str = " ".join(elements)
    if element_str not in pair_coeff:
        raise ValueError(f"Elements {element_str} not present in pair_coeff {pair_coeff}")
    
    filename, masses = _write_structure(structure, filename, elements)

    input_dict["pair_style"] = pair_style
    input_dict["pair_coeff"] = pair_coeff
    input_dict["element"] = elements
    input_dict["mass"] = masses
    input_dict['mode'] = mode
    input_dict['reference_phase'] = reference_phase
    input_dict['lattice'] = filename
    input_dict['temperature'] = temperature
    input_dict['pressure'] = pressure
    input_dict['n_switching_steps'] = n_switching_steps
    input_dict['n_equilibration_steps'] = n_switching_steps//2
    
    if cores > 1:
        input_dict["queue"]["cores"] = cores

    calc = Calculation(**input_dict)
    return calc


def _run_cleanup(simfolder, lattice, delete_folder=False):
    import shutil
    import os
    os.remove(lattice)
    if delete_folder:
        shutil.rmtree(simfolder)

@job
@ase_pmg_bridge
def calculate_free_energy(
    structure,
    pair_style, pair_coeff,
                   temperature, pressure=0,
                   n_switching_steps=10000,
                   mode='fe', reference_phase='solid',
                   cores =1, 
                   input_dict=None
):
    from calphy.solid import Solid
    from calphy.liquid import Liquid
    from calphy.routines import routine_fe, routine_ts, routine_pscale
    import os
    import numpy as np

    calc = _prepare_input(pair_style, pair_coeff, structure, temperature, pressure, n_switching_steps, mode, reference_phase, cores, input_dict)

    simfolder = calc.create_folders()
    if reference_phase == 'solid':
        job = Solid(calculation=calc, simfolder=simfolder)
    elif reference_phase == 'liquid':
        job = Liquid(calculation=calc, simfolder=simfolder)
    
    if mode == 'fe':
        job = routine_fe(job)
    elif mode == 'ts':
        job = routine_ts(job)
    elif mode == 'pscale':
        job = routine_pscale(job)
    else:
        raise ValueError(f"Unknown mode {mode}, should be one of fe, ts, pscale")
    
    _run_cleanup(simfolder, calc.lattice)

    results = job.report["results"]
    results['temperature'] = temperature
    results['pressure'] = pressure
    
    if mode == 'ts':
        outfile = os.path.join(simfolder, "temperature_sweep.dat")
        temp, fe = np.loadtxt(outfile, unpack=True, usecols=(0,1))
        results['temperature'] = temp
        results['free_energy'] = fe
    elif mode == 'pscale':
        outfile = os.path.join(simfolder, "pressure_sweep.dat")
        pres, fe = np.loadtxt(outfile, unpack=True, usecols=(0,1))
        results['pressure'] = pres
        results['free_energy'] = fe

    return results



