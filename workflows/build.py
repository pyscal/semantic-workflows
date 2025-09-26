from typing import Optional
from ase.atoms import Atoms
from ase.build import bulk as ase_bulk


from pyiron_workflow import as_function_node, as_macro_node

@as_function_node("structure")
def bulk(
    name: str,
    crystalstructure: Optional[str] = None,
    a: Optional[float | int] = None,
    c: Optional[float | int] = None,
    c_over_a: Optional[float] | int = None,
    u: Optional[float | int] = None,
    orthorhombic: bool = False,
    cubic: bool = False,
):

    return ase_bulk(
        name,
        crystalstructure=crystalstructure,
        a=a,
        c=c,
        covera = c_over_a,
        u = u,
        orthorhombic = orthorhombic,
        cubic = cubic,
    )

@as_function_node("structure")
def repeat(
    structure: Atoms,
    repetitions: tuple[int, int, int],
) -> Atoms:
    return structure.repeat(repetitions)

@as_function_node("structure")
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


    os.system('atomsk --polycrystal tmp.xsf grain_sizes.txt final.cfg -ow')
    
    poly_struct = read('final.cfg')
    os.remove('grain_sizes.txt')
    os.remove('tmp.xsf')
    os.remove('final.cfg')
    return poly_struct