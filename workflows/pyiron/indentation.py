import pandas as pd
import numpy as np
from ase.io import write
from pylammpsmpi import LammpsLibrary
from pyiron_workflow import as_function_node

@as_function_node
def indentation_test(structure, pair_style, pair_coeff,
                     indenter_radius=80,
                     indenter_radius_factor=0.65,
                     indentation_velocity=1.0,
                     cores=1,
                     initial_temperature=10.0,
                     annealing_temperature=750,
                     n_equilibration_steps=5000,
                     dump_interval=1000,
                     print_interval=1000):
    """
    Perform nanoindentation simulation on a structure.
    
    Parameters
    ----------
    structure : ase.Atoms
        Input atomic structure
    pair_style : str
        LAMMPS pair_style command
    pair_coeff : str
        LAMMPS pair_coeff command
    indenter_radius : float
        Radius of spherical indenter in Angstrom (default: 80)
    indenter_radius_factor : float
        Factor for indentation depth (default: 0.65)
    indentation_velocity : float
        Velocity of indenter in Angstrom/ps (default: 1.0)
    cores : int
        Number of cores for parallel execution
    initial_temperature : float
        Starting temperature in K (default: 10.0)
    annealing_temperature : float
        Annealing temperature in K (default: 750)
    n_equilibration_steps : int
        Number of equilibration steps (default: 5000)
    dump_interval : int
        Interval for dump output (default: 1000)
    print_interval : int
        Interval for thermo output (default: 1000)
        
    Returns
    -------
    dict
        Dictionary containing indentation results
    """
    
    from mendeleev import element
    
    # Write structure to LAMMPS data file
    write('tmp.data', structure, format='lammps-data')
    
    # Get masses for elements
    symbols = np.unique(structure.get_chemical_symbols())
    masses = [element(symbol).mass for symbol in symbols]
    
    # Initialize LAMMPS
    lmp = LammpsLibrary(cores=cores)
    lmp.command("units metal")
    lmp.command("dimension 3")
    lmp.command("boundary p p p")
    lmp.command("atom_style atomic")
    lmp.command("atom_modify map array")
    
    # Read structure
    lmp.command("read_data tmp.data")
    
    # Set masses
    for i, mass in enumerate(masses):
        lmp.command(f"mass {i+1} {mass}")
    
    # Set potential
    lmp.command(f"pair_style {pair_style}")
    lmp.command(f"pair_coeff {pair_coeff}")
    
    # Setup neighbor style
    lmp.command("neighbor 1.0 nsq")
    lmp.command("neigh_modify once no every 1 delay 0 check yes")
    
    # Thermodynamics settings
    lmp.command("thermo 1000")
    lmp.command("thermo_style custom step temp press vol lx ly lz pe enthalpy")
    
    # Time settings
    ts = 0.001
    lmp.command(f"timestep {ts}")
    
    # Calculate vacuum and indentation parameters
    vacuum_lo = indenter_radius * 0.2
    vacuum_hi = 3 * indenter_radius
    inden_lo = indenter_radius * indenter_radius_factor
    inden_hi = 1.03 * indenter_radius
    
    # Equilibration phase
    lmp.command(f"velocity all create {initial_temperature} {np.random.randint(1000, 9999999)} mom yes rot yes dist gaussian")
    
    # Step 1: Heat to annealing temperature
    lmp.command(f"fix 2 all temp/rescale 1 {initial_temperature} {annealing_temperature} 1.0 1.0")
    lmp.command(f"fix 3 all nph x 0.0 0.0 {ts*1000} y 0.0 0.0 {ts*1000} z 0.0 0.0 {ts*1000} drag 1.0 nreset 1000")
    lmp.command(f"run {n_equilibration_steps}")
    lmp.command("unfix 2")
    lmp.command("unfix 3")
    
    # Step 2: Stabilize at annealing temperature
    lmp.command(f"fix 2 all temp/rescale 1 {annealing_temperature} {annealing_temperature} 1.0 1.0")
    lmp.command(f"fix 3 all nph x 0.0 0.0 {ts*1000} y 0.0 0.0 {ts*1000} z 0.0 0.0 {ts*1000} drag 1.0 nreset 1000")
    lmp.command(f"run {n_equilibration_steps * 5}")
    lmp.command("unfix 2")
    lmp.command("unfix 3")
    
    # Step 3: Cool to initial temperature
    lmp.command(f"fix 2 all temp/rescale 1 {annealing_temperature} {initial_temperature} 1.0 1.0")
    lmp.command(f"fix 3 all nph x 0.0 0.0 {ts*1000} y 0.0 0.0 {ts*1000} z 0.0 0.0 {ts*1000} drag 1.0 nreset 1000")
    lmp.command(f"run {n_equilibration_steps}")
    lmp.command("unfix 2")
    lmp.command("unfix 3")
    
    # Get box dimensions
    lmp.command("variable x1 equal xlo")
    lmp.command("variable X_lo equal ${x1}")
    lmp.command("variable x2 equal xhi")
    lmp.command("variable X_hi equal ${x2}")
    lmp.command("variable y1 equal ylo")
    lmp.command("variable Y_lo equal ${y1}")
    lmp.command("variable y2 equal yhi")
    lmp.command("variable Y_hi equal ${y2}")
    lmp.command("variable z1 equal zlo")
    lmp.command("variable Z_lo equal ${z1}")
    lmp.command("variable z2 equal zhi")
    lmp.command("variable Z_hi equal ${z2}")
    
    # Define layers
    lmp.command("variable Z1 equal ${Z_lo}")
    lmp.command("variable Z2 equal (${Z_lo}+10)")
    lmp.command("variable Z3 equal (${Z_lo}+20)")
    
    # Fixed atoms layer
    lmp.command("region b_fix block INF INF INF INF ${Z1} ${Z2} units box")
    lmp.command("group b_fix region b_fix")
    
    # Thermostatic atoms layer
    lmp.command("region b_the block INF INF INF INF ${Z2} ${Z3} units box")
    lmp.command("group b_the region b_the")
    
    # Substrate layer
    lmp.command("region b_sub block INF INF INF INF ${Z3} INF units box")
    lmp.command("group b_sub region b_sub")
    
    # Change box to add vacuum
    lmp.command(f"variable Z_hig equal (${{Z_hi}}+{vacuum_hi})")
    lmp.command(f"variable Z_low equal (${{Z_lo}}-{vacuum_lo})")
    lmp.command("change_box all z final ${Z_low} ${Z_hig}")
    
    # Setup run style
    lmp.command("run_style verlet")
    
    # Define computes
    lmp.command("compute Ene all pe/atom")
    lmp.command("compute Str all stress/atom NULL")
    lmp.command("compute f all property/atom fx fy fz")
    
    # Apply fixes to different layers
    lmp.command("fix fixed b_fix setforce 0.0 0.0 0.0")
    lmp.command(f"fix Therm b_the nvt temp {initial_temperature} {initial_temperature} {ts*100}")
    lmp.command("fix Subst b_sub nve")
    
    # Indenter position
    lmp.command("variable Ind_x equal 0.5*(${X_lo}+${X_hi})")
    lmp.command("variable Ind_y equal 0.5*(${Y_lo}+${Y_hi})")
    
    # Indentation parameters
    lmp.command(f"variable Iden_ini equal (${{Z_hi}}+{inden_hi})")
    lmp.command(f"variable Inde_fin equal (${{Z_hi}}+{inden_lo})")
    lmp.command("variable Inde_dis equal (${Iden_ini}-${Inde_fin})")
    lmp.command(f"variable vel equal {indentation_velocity}")
    lmp.command(f"variable Inde_tim equal (${{Inde_dis}}/(${{vel}}*{ts}))")
    lmp.command("variable Ind_z equal v_Iden_ini-step*dt*v_vel")
    
    # Apply indenter
    lmp.command(f"fix indent all indent 10.0 sphere v_Ind_x v_Ind_y v_Ind_z {indenter_radius} units box")
    
    # Thermo output
    lmp.command(f"thermo {print_interval}")
    lmp.command("thermo_style custom step temp ke pe press pxx pyy pzz pxy pxz pyz v_Ind_x v_Ind_y v_Ind_z f_indent f_indent[1] f_indent[2] f_indent[3]")
    
    # Dump output
    lmp.command(f"dump OUT1 all custom {dump_interval} Indent.* id type x y z vx vy vz c_Ene c_Str[1] c_Str[2] c_Str[3] c_Str[4] c_Str[5] c_Str[6] c_f[1] c_f[2] c_f[3]")
    
    # Print forces
    lmp.command("variable f0 equal f_indent")
    lmp.command("variable f1 equal f_indent[1]")
    lmp.command("variable f2 equal f_indent[2]")
    lmp.command("variable f3 equal f_indent[3]")
    lmp.command(f"fix def1 all print {print_interval} \"$(step), $(temp), ${{Ind_x}}, ${{Ind_y}}, ${{Ind_z}}, ${{f0}}, ${{f1}}, ${{f2}}, ${{f3}}\" append Result-indent.txt screen no")
    
    # Run indentation
    lmp.command("variable nos equal floor(${Inde_tim})")
    lmp.command("reset_timestep 0")  # Reset step counter to 0 before indentation
    lmp.command("run ${nos}")
    
    # Final dump
    lmp.command("undump OUT1")
    lmp.command("dump OUT1 all custom 1 Indent.* id type x y z vx vy vz c_Ene c_Str[1] c_Str[2] c_Str[3] c_Str[4] c_Str[5] c_Str[6] c_f[1] c_f[2] c_f[3]")
    lmp.command("run 1")
    lmp.command("write_restart restart.indent")
    
    # Cleanup
    lmp.command("undump OUT1")
    lmp.command("unfix def1")
    lmp.command("unfix indent")
    
    # Read results
    results_data = pd.read_csv("Result-indent.txt", names=["step", "temp", "Ind_x", "Ind_y", "Ind_z", "f0", "f1", "f2", "f3"])
    
    results = {
        'indentation_force': results_data['f0'].values,
        'indentation_depth': results_data['Ind_z'].values,
        'force_x': results_data['f1'].values,
        'force_y': results_data['f2'].values,
        'force_z': results_data['f3'].values,
        'step': results_data['step'].values,
        'temperature': results_data['temp'].values,
    }
    
    lmp.close()
    
    return results


@as_function_node
def read_final_structure(structure):
    """
    Read the final indented structure from the dump file.
    
    Parameters
    ----------
    structure : ase.Atoms
        Original structure (used to get atom types)
    results : dict
        Dictionary containing indentation results (used as dependency)
        
    Returns
    -------
    ase.Atoms
        Final structure after indentation
    """
    import glob
    from ase.io import read
    from ase import Atoms
    
    # Find the latest dump file (Indent.*)
    dump_files = sorted(glob.glob("Indent.*"))
    if not dump_files:
        raise FileNotFoundError("No dump files found matching 'Indent.*'")
    
    # Read the last dump file (which contains the final structure)
    final_dump = dump_files[-1]
    final_struct = read(final_dump, format='lammps-dump-text')

    return final_struct


@as_function_node
def plot_force_depth(results):
    """
    Plot force-depth curve and force components.
    
    Parameters
    ----------
    results : dict
        Dictionary containing indentation results
        
    Returns
    -------
    int
        Returns 1 as proxy value
    """
    import matplotlib.pyplot as plt
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot total indentation force vs depth
    ax1.plot(results['indentation_depth'], results['indentation_force'], 'b-', linewidth=2)
    ax1.set_xlabel('Indentation Depth (Å)', fontsize=12)
    ax1.set_ylabel('Indentation Force (eV/Å)', fontsize=12)
    ax1.set_title('Force-Depth Curve', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Plot force components
    ax2.plot(results['step'], results['force_x'], label='F_x', alpha=0.7)
    ax2.plot(results['step'], results['force_y'], label='F_y', alpha=0.7)
    ax2.plot(results['step'], results['force_z'], label='F_z', alpha=0.7)
    ax2.set_xlabel('Step', fontsize=12)
    ax2.set_ylabel('Force Components (eV/Å)', fontsize=12)
    ax2.set_title('Force Components vs Step', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return 1


@as_function_node
def plot_temperature(results):
    """
    Plot temperature evolution during indentation.
    
    Parameters
    ----------
    results : dict
        Dictionary containing indentation results
        
    Returns
    -------
    int
        Returns 1 as proxy value
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(results['step'], results['temperature'], 'r-', linewidth=2)
    ax.set_xlabel('Step', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Evolution During Indentation', fontsize=14)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    return 1


@as_function_node
def plot_centrosymmetry(final_struct, n_neighbors=8):
    """
    Calculate and plot centrosymmetry parameter for the final structure.
    
    Parameters
    ----------
    final_struct : ase.Atoms
        Final indented structure
    n_neighbors : int
        Number of neighbors for centrosymmetry calculation (default: 8)
        
    Returns
    -------
    int
        Returns 1 as proxy value
    """
    from pyscal3 import System
    from pyscal3.operations.visualize import plot_by_property
    
    # Create pyscal3 System from ASE structure
    sys = System(final_struct, format='ase')
    
    # Calculate centrosymmetry parameter
    sys.calculate.centrosymmetry(n_neighbors)
    
    # Plot centrosymmetry
    plot_by_property(sys, sys.atoms.centrosymmetry)
    
    return 1
