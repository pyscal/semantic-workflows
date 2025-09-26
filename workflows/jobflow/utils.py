from functools import wraps
from collections.abc import Mapping, Iterable

from ase import Atoms as ASEAtoms
from pymatgen.core import Structure, Molecule
from pymatgen.io.ase import AseAtomsAdaptor


_adaptor = AseAtomsAdaptor()


def _to_ase(obj):
    """Convert a single object to ASE Atoms if it is a pymatgen Structure/Molecule."""
    if isinstance(obj, Structure) or isinstance(obj, Molecule):
        return _adaptor.get_atoms(obj)
    return obj


def _to_pmg(obj):
    """Convert a single object to pymatgen Structure or Molecule if it is ASE Atoms."""
    if isinstance(obj, ASEAtoms):
        return _adaptor.get_structure(obj)
    return obj


def _walk_convert(obj, conv_single):
    """Recursively convert within common containers without touching unknown types."""
    if isinstance(obj, Mapping):
        return obj.__class__((_walk_convert(k, conv_single), _walk_convert(v, conv_single)) for k, v in obj.items())
    # Treat strings as atomic, and avoid converting numpy arrays elementwise
    if isinstance(obj, (str, bytes)):
        return obj
    if isinstance(obj, tuple):
        return tuple(_walk_convert(x, conv_single) for x in obj)
    if isinstance(obj, list):
        return [_walk_convert(x, conv_single) for x in obj]
    if isinstance(obj, set):
        return {_walk_convert(x, conv_single) for x in obj}
    return conv_single(obj)


def ase_pmg_bridge(func):
    """
    Decorator for Jobflow tasks:
      - Convert any pymatgen Structure/Molecule in args/kwargs to ASE Atoms before calling `func`
      - Convert any ASE Atoms in the return value to pymatgen Structure/Molecule
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        args_ase = _walk_convert(args, _to_ase)
        kwargs_ase = _walk_convert(kwargs, _to_ase)
        result = func(*args_ase, **kwargs_ase)
        result_pmg = _walk_convert(result, _to_pmg)
        return result_pmg
    return wrapper