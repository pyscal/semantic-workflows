import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from ase.io import read, write
from pylammpsmpi import LammpsLibrary
from pyiron_workflow import as_function_node
from .templates import property_template, workflow_template

def _displace(lmp, dir, pair_style, pair_coeff):
    if dir == 1:
        lmp.command("variable len0 equal ${lx0}")
    if dir == 2:
        lmp.command("variable len0 equal ${ly0}")
    if dir == 3:
        lmp.command("variable len0 equal ${lz0}")
    if dir == 4:
        lmp.command("variable len0 equal ${lz0}")
    if dir == 5:
        lmp.command("variable len0 equal ${lz0}")
    if dir == 6:
        lmp.command("variable len0 equal ${ly0}")

    lmp.command("clear")
    lmp.command("box tilt large")
    lmp.command("read_restart restart.equil")
    lmp.command(f"pair_style {pair_style}")
    lmp.command(f"pair_coeff {pair_coeff}")

    lmp.command("variable delta equal -${up}*${len0}")
    lmp.command("variable deltaxy equal -${up}*xy")
    lmp.command("variable deltaxz equal -${up}*xz")
    lmp.command("variable deltayz equal -${up}*yz")
    if dir == 1:
        lmp.command("change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box")
    if dir == 2:
        lmp.command("change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box")
    if dir == 3:
        lmp.command("change_box all z delta 0 ${delta} remap units box")
    if dir == 4:
        lmp.command("change_box all yz delta ${delta} remap units box")
    if dir == 5:
        lmp.command("change_box all xz delta ${delta} remap units box")
    if dir == 6:
        lmp.command("change_box all xy delta ${delta} remap units box")

    lmp.command("minimize ${etol} ${ftol} ${maxiter} ${maxeval}")

    lmp.command("variable tmp equal pxx")
    lmp.command("variable pxx1 equal ${tmp}")
    lmp.command("variable tmp equal pyy")
    lmp.command("variable pyy1 equal ${tmp}")
    lmp.command("variable tmp equal pzz")
    lmp.command("variable pzz1 equal ${tmp}")
    lmp.command("variable tmp equal pxy")
    lmp.command("variable pxy1 equal ${tmp}")
    lmp.command("variable tmp equal pxz")
    lmp.command("variable pxz1 equal ${tmp}")
    lmp.command("variable tmp equal pyz")
    lmp.command("variable pyz1 equal ${tmp}")

    lmp.command("variable C1neg equal ${d1}")
    lmp.command("variable C2neg equal ${d2}")
    lmp.command("variable C3neg equal ${d3}")
    lmp.command("variable C4neg equal ${d4}")
    lmp.command("variable C5neg equal ${d5}")
    lmp.command("variable C6neg equal ${d6}")

    lmp.command("clear")
    lmp.command("box tilt large")
    lmp.command("read_restart restart.equil")
    lmp.command(f"pair_style {pair_style}")
    lmp.command(f"pair_coeff {pair_coeff}")

    lmp.command("variable delta equal ${up}*${len0}")
    lmp.command("variable deltaxy equal ${up}*xy")
    lmp.command("variable deltaxz equal ${up}*xz")
    lmp.command("variable deltayz equal ${up}*yz")
    if dir == 1:
        lmp.command("change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box")
    if dir == 2:
        lmp.command("change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box")
    if dir == 3:
        lmp.command("change_box all z delta 0 ${delta} remap units box")
    if dir == 4:
        lmp.command("change_box all yz delta ${delta} remap units box")
    if dir == 5:
        lmp.command("change_box all xz delta ${delta} remap units box")
    if dir == 6:
        lmp.command("change_box all xy delta ${delta} remap units box")

    lmp.command("minimize ${etol} ${ftol} ${maxiter} ${maxeval}")

    lmp.command("variable tmp equal pe")
    lmp.command("variable e1 equal ${tmp}")
    lmp.command("variable tmp equal press")
    lmp.command("variable p1 equal ${tmp}")
    lmp.command("variable tmp equal pxx")
    lmp.command("variable pxx1 equal ${tmp}")
    lmp.command("variable tmp equal pyy")
    lmp.command("variable pyy1 equal ${tmp}")
    lmp.command("variable tmp equal pzz")
    lmp.command("variable pzz1 equal ${tmp}")
    lmp.command("variable tmp equal pxy")
    lmp.command("variable pxy1 equal ${tmp}")
    lmp.command("variable tmp equal pxz")
    lmp.command("variable pxz1 equal ${tmp}")
    lmp.command("variable tmp equal pyz")
    lmp.command("variable pyz1 equal ${tmp}")

    # Compute elastic constant from pressure tensor

    lmp.command("variable C1pos equal ${d1}")
    lmp.command("variable C2pos equal ${d2}")
    lmp.command("variable C3pos equal ${d3}")
    lmp.command("variable C4pos equal ${d4}")
    lmp.command("variable C5pos equal ${d5}")
    lmp.command("variable C6pos equal ${d6}")

    # Combine positive and negative 

    lmp.command("variable C1${dir} equal 0.5*(${C1neg}+${C1pos})")
    lmp.command("variable C2${dir} equal 0.5*(${C2neg}+${C2pos})")
    lmp.command("variable C3${dir} equal 0.5*(${C3neg}+${C3pos})")
    lmp.command("variable C4${dir} equal 0.5*(${C4neg}+${C4pos})")
    lmp.command("variable C5${dir} equal 0.5*(${C5neg}+${C5pos})")
    lmp.command("variable C6${dir} equal 0.5*(${C6neg}+${C6pos})")

    # Delete dir to make sure it is not reused

    lmp.command("variable dir delete")
    return lmp

@as_function_node
def calculate_elastic_constants(structure, pair_style, pair_coeff, cores=1,
        finite_deformation_size=1e-6,
        energy_tolerance=0.0,
        force_tolerance=1.0e-15,
        max_iterations=100,
        kg = None,
        potential_type=None, potential_doi=None):

    write('tmp.data', structure, format='lammps-data')
    lmp = LammpsLibrary(cores=cores)

    lmp.command(f"variable up equal {finite_deformation_size}")
    lmp.command("variable atomjiggle equal 1.0e-5")

    lmp.command("units metal")
    lmp.command("variable cfac equal 1.0e-4")
    lmp.command("variable cunits string GPa")

    lmp.command(f"variable etol equal {energy_tolerance}")
    lmp.command(f"variable ftol equal {force_tolerance}")
    lmp.command(f"variable maxiter equal {max_iterations}")
    lmp.command("variable maxeval equal 1000")
    lmp.command("variable dmax equal 1.0e-2")

    #add read data command

    lmp.command("read_data tmp.data")
    lmp.command("change_box all triclinic")

    lmp.command(f"pair_style {pair_style}")
    lmp.command(f"pair_coeff {pair_coeff}")

    # Need to set mass to something, just to satisfy LAMMPS
    lmp.command("mass 1 55.85")

    # Setup neighbor style
    lmp.command("neighbor 1.0 nsq")
    lmp.command("neigh_modify once no every 1 delay 0 check yes")

    # Setup minimization style
    lmp.command("min_style cg")
    lmp.command("min_modify dmax ${dmax} line quadratic")

    # Setup output
    lmp.command("thermo 1")
    lmp.command("thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol")
    lmp.command("thermo_modify norm no")

    lmp.command("fix 3 all box/relax aniso 0.0")
    lmp.command("minimize ${etol} ${ftol} ${maxiter} ${maxeval}")

    lmp.command("variable tmp equal pxx")
    lmp.command("variable pxx0 equal ${tmp}")
    lmp.command("variable tmp equal pyy")
    lmp.command("variable pyy0 equal ${tmp}")
    lmp.command("variable tmp equal pzz")
    lmp.command("variable pzz0 equal ${tmp}")
    lmp.command("variable tmp equal pyz")
    lmp.command("variable pyz0 equal ${tmp}")
    lmp.command("variable tmp equal pxz")
    lmp.command("variable pxz0 equal ${tmp}")
    lmp.command("variable tmp equal pxy")
    lmp.command("variable pxy0 equal ${tmp}")

    lmp.command("variable tmp equal lx")
    lmp.command("variable lx0 equal ${tmp}")
    lmp.command("variable tmp equal ly")
    lmp.command("variable ly0 equal ${tmp}")
    lmp.command("variable tmp equal lz")
    lmp.command("variable lz0 equal ${tmp}")

    lmp.command("variable d1 equal -(v_pxx1-${pxx0})/(v_delta/v_len0)*${cfac}")
    lmp.command("variable d2 equal -(v_pyy1-${pyy0})/(v_delta/v_len0)*${cfac}")
    lmp.command("variable d3 equal -(v_pzz1-${pzz0})/(v_delta/v_len0)*${cfac}")
    lmp.command("variable d4 equal -(v_pyz1-${pyz0})/(v_delta/v_len0)*${cfac}")
    lmp.command("variable d5 equal -(v_pxz1-${pxz0})/(v_delta/v_len0)*${cfac}")
    lmp.command("variable d6 equal -(v_pxy1-${pxy0})/(v_delta/v_len0)*${cfac}")

    lmp.command("displace_atoms all random ${atomjiggle} ${atomjiggle} ${atomjiggle} 87287 units box")

    lmp.command("unfix 3")
    lmp.command("write_restart restart.equil")

    # uxx Perturbation

    lmp.command("variable dir equal 1")
    lmp = _displace(lmp, 1, pair_style, pair_coeff)

    # uyy Perturbation

    lmp.command("variable dir equal 2")
    lmp = _displace(lmp, 2, pair_style, pair_coeff)

    # uzz Perturbation

    lmp.command("variable dir equal 3")
    lmp = _displace(lmp, 3, pair_style, pair_coeff)

    # uyz Perturbation

    lmp.command("variable dir equal 4")
    lmp = _displace(lmp, 4, pair_style, pair_coeff)

    # uxz Perturbation

    lmp.command("variable dir equal 5")
    lmp = _displace(lmp, 5, pair_style, pair_coeff)

    # uxy Perturbation

    lmp.command("variable dir equal 6")
    lmp = _displace(lmp, 6, pair_style, pair_coeff)

    # Output final values

    lmp.command("variable C11all equal ${C11}")
    lmp.command("variable C22all equal ${C22}")
    lmp.command("variable C33all equal ${C33}")

    lmp.command("variable C12all equal 0.5*(${C12}+${C21})")
    lmp.command("variable C13all equal 0.5*(${C13}+${C31})")
    lmp.command("variable C23all equal 0.5*(${C23}+${C32})")

    lmp.command("variable C44all equal ${C44}")
    lmp.command("variable C55all equal ${C55}")
    lmp.command("variable C66all equal ${C66}")

    lmp.command("variable C14all equal 0.5*(${C14}+${C41})")
    lmp.command("variable C15all equal 0.5*(${C15}+${C51})")
    lmp.command("variable C16all equal 0.5*(${C16}+${C61})")

    lmp.command("variable C24all equal 0.5*(${C24}+${C42})")
    lmp.command("variable C25all equal 0.5*(${C25}+${C52})")
    lmp.command("variable C26all equal 0.5*(${C26}+${C62})")

    lmp.command("variable C34all equal 0.5*(${C34}+${C43})")
    lmp.command("variable C35all equal 0.5*(${C35}+${C53})")
    lmp.command("variable C36all equal 0.5*(${C36}+${C63})")

    lmp.command("variable C45all equal 0.5*(${C45}+${C54})")
    lmp.command("variable C46all equal 0.5*(${C46}+${C64})")
    lmp.command("variable C56all equal 0.5*(${C56}+${C65})")

    # Average moduli for cubic crystals
    c11 = lmp.extract_variable("C11all", None, 0)
    c12 = lmp.extract_variable("C12all", None, 0)
    c13 = lmp.extract_variable("C13all", None, 0)
    c14 = lmp.extract_variable("C14all", None, 0)
    c15 = lmp.extract_variable("C15all", None, 0)
    c16 = lmp.extract_variable("C16all", None, 0)

    c22 = lmp.extract_variable("C22all", None, 0)
    c23 = lmp.extract_variable("C23all", None, 0)
    c24 = lmp.extract_variable("C24all", None, 0)
    c25 = lmp.extract_variable("C25all", None, 0)
    c26 = lmp.extract_variable("C26all", None, 0)

    c33 = lmp.extract_variable("C33all", None, 0)
    c34 = lmp.extract_variable("C34all", None, 0)
    c35 = lmp.extract_variable("C35all", None, 0)
    c36 = lmp.extract_variable("C36all", None, 0)

    c44 = lmp.extract_variable("C44all", None, 0)
    c45 = lmp.extract_variable("C45all", None, 0)
    c46 = lmp.extract_variable("C46all", None, 0)

    c55 = lmp.extract_variable("C55all", None, 0)
    c56 = lmp.extract_variable("C56all", None, 0)
    c66 = lmp.extract_variable("C66all", None, 0)

    c_matrix = np.array([[c11, c12, c13, c14, c15, c16],
                         [c12, c22, c23, c24, c25, c26],
                         [c13, c23, c33, c34, c35, c36],
                         [c14, c24, c34, c44, c45, c46],
                         [c15, c25, c35, c45, c55, c56],
                         [c16, c26, c36, c46, c56, c66]])

    lmp.command("variable C11cubic equal (${C11all}+${C22all}+${C33all})/3.0")
    lmp.command("variable C12cubic equal (${C12all}+${C13all}+${C23all})/3.0")
    lmp.command("variable C44cubic equal (${C44all}+${C55all}+${C66all})/3.0")

    lmp.command("variable bulkmodulus equal (${C11cubic}+2*${C12cubic})/3.0")
    lmp.command("variable shearmodulus1 equal ${C44cubic}")
    lmp.command("variable shearmodulus2 equal (${C11cubic}-${C12cubic})/2.0")
    lmp.command("variable poissonratio equal 1.0/(1.0+${C11cubic}/${C12cubic})")

    B = lmp.extract_variable("bulkmodulus", None, 0)
    G1 = lmp.extract_variable("shearmodulus1", None, 0)
    G2 = lmp.extract_variable("shearmodulus2", None, 0)
    poisson = lmp.extract_variable("poissonratio", None, 0)

    c_matrix = np.array([[c11, c12, c13, c14, c15, c16],
                         [c12, c22, c23, c24, c25, c26],
                         [c13, c23, c33, c34, c35, c36],
                         [c14, c24, c34, c44, c45, c46],
                         [c15, c25, c35, c45, c55, c56],
                         [c16, c26, c36, c46, c56, c66]])
    
    if kg is not None:
        workflow = workflow_template.copy()
        workflow['method'] = 'MolecularStatics'
        workflow['input_sample'] = [structure.info['id']]
        sample_id = structure.info['id']
        
        outputs = [
            {"label": "C11", "value": c_matrix[0, 0], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C12", "value": c_matrix[0, 1], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C13", "value": c_matrix[0, 2], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C14", "value": c_matrix[0, 3], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C15", "value": c_matrix[0, 4], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C16", "value": c_matrix[0, 5], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
        
            {"label": "C22", "value": c_matrix[1, 1], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C23", "value": c_matrix[1, 2], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C24", "value": c_matrix[1, 3], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C25", "value": c_matrix[1, 4], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C26", "value": c_matrix[1, 5], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
        
            {"label": "C33", "value": c_matrix[2, 2], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C34", "value": c_matrix[2, 3], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C35", "value": c_matrix[2, 4], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C36", "value": c_matrix[2, 5], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
        
            {"label": "C44", "value": c_matrix[3, 3], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C45", "value": c_matrix[3, 4], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C46", "value": c_matrix[3, 5], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
        
            {"label": "C55", "value": c_matrix[4, 4], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "C56", "value": c_matrix[4, 5], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
        
            {"label": "C66", "value": c_matrix[5, 5], "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ElasticConstant"},
            {"label": "BulkModulus", "value": B, "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "BulkModulus"},
            {"label": "ShearModulus", "value": G1, "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ShearModulus"},
            {"label": "ShearModulus", "value": G2, "unit": "GigaPa", "associate_to_sample": [sample_id], "basename": "ShearModulus"},
            {"label": "PoissonRation", "value": poisson, "associate_to_sample": [sample_id], "basename": "PoissonRatio"},
        ]
        
        workflow['degrees_of_freedom'] = ["AtomicPositionRelaxation"]
        workflow['calculated_property'] = outputs
        workflow['interatomic_potential'] = {
            'potential_type': potential_type,
            'uri': potential_doi,
        }
        workflow['software'] = {
            'uri': "https://doi.org/10.1016/j.cpc.2021.108171",
            'label': 'LAMMPS',
        }
        
        # Append a *copy* to avoid overwriting in subsequent iterations
        kg['workflow'].append(workflow.copy())
    
    results = {'C_matrix': c_matrix,
               'Bulk_modulus': B,
               'Shear_modulus1': G1,
               'Shear_modulus2': G2,
               'Poisson_ratio': poisson}
    return results

@as_function_node
def compression_test(structure, pair_style, 
                     pair_coeff, cores=1,
                     temperature=10,
                     annealing_temperature=600,
                     n_equilibration_steps=5000,
                     n_run_steps=10000,
                     strain_rate=1e-5):

    import numpy as np
    from ase.io import write
    write('tmp.data', structure, format='lammps-data')

    lmp = LammpsLibrary(cores=cores)
    lmp.command("units metal") 
    lmp.command("dimension 3") 
    lmp.command("boundary p p p") 
    lmp.command("atom_style atomic")
    lmp.command("read_data tmp.data")
    lmp.command(f"pair_style {pair_style}")
    lmp.command(f"pair_coeff {pair_coeff}")

    #anneal the grains
    lmp.command(f"velocity        all create {temperature} {np.random.randint(10000)} mom yes rot yes dist gaussian")

    lmp.command(f"fix 2 all temp/rescale 1 {temperature} {annealing_temperature} 1.0 1.0")
    lmp.command(f"fix 3 all nph x 0.0 0.0 0.1 y 0.0 0.0 0.1 z 0.0 0.0 0.1 drag 1.0 nreset 1000")

    lmp.command(f"run {n_equilibration_steps}")
    lmp.command("unfix 2")
    lmp.command("unfix 3")

    lmp.command(f"fix 2 all temp/rescale 1 {annealing_temperature} {annealing_temperature} 1.0 1.0")
    lmp.command(f"fix 3 all nph x 0.0 0.0 0.1 y 0.0 0.0 0.1 z 0.0 0.0 0.1 drag 1.0 nreset 1000")
    lmp.command(f"run {n_equilibration_steps}")
    lmp.command("unfix 2")
    lmp.command("unfix 3")

    lmp.command(f"fix 2 all temp/rescale 1 {annealing_temperature} {temperature} 1.0 1.0")
    lmp.command(f"fix 3 all nph x 0.0 0.0 0.1 y 0.0 0.0 0.1 z 0.0 0.0 0.1 drag 1.0 nreset 1000")
    lmp.command(f"variable runing equal \"20/v_ts\"")
    lmp.command(f"run {n_equilibration_steps}")
    lmp.command("unfix 2")
    lmp.command("unfix 3")

    lmp.command("variable L0x equal lx")
    lmp.command("variable L0y equal ly")
    lmp.command("variable L0z equal lz")

    lmp.command("variable strainx equal \"(lx - v_L0x)/v_L0x\"")
    lmp.command("variable strainy equal \"(ly - v_L0y)/v_L0y\"")
    lmp.command("variable strainz equal \"(lz - v_L0z)/v_L0z\"")
    lmp.command("variable p1x equal \"-v_strainx*100\"")
    lmp.command("variable p1y equal \"-v_strainy*100\"")
    lmp.command("variable p1z equal \"-v_strainz*100\"")

    lmp.command("reset_timestep  0")
    lmp.command("compute         stress all stress/atom NULL")
    lmp.command("compute         peratom all pe/atom")

    lmp.command("dump myDump all custom 7500 dump.* id type x y z c_peratom c_stress[1] c_stress[2] c_stress[3] c_stress[4] c_stress[5] c_stress[6]")


    lmp.command("fix 2 all nve")
    lmp.command(f"fix 3 all deform 1 x erate {strain_rate} y erate {strain_rate} z erate {strain_rate}")

    lmp.command("fix def1 all print 100 \"$(temp), $(press), $(enthalpy), $(pe), ${p1x}, ${p1y}, ${p1z}\" append quantities.dat screen no")
    lmp.command("fix def2 all print 100 \"$(pxx), $(pyy), $(pzz), $(xy), $(xz), $(yz), ${p1x}, ${p1y}, ${p1z}\" append pressure.dat screen no")

    lmp.command(f"run {n_run_steps}")