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
