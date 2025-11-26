import yaml
from typing import Any, Dict
import numpy as np
import string
import random

class KnowledgeDict(dict):
    def __init__(self, *args, **kwargs):
        data = {'computational_sample': [], 'workflow': []}
        super().__init__(data, *args, **kwargs)


    
    @staticmethod
    def _clean_data(obj: Any) -> Any:
        if isinstance(obj, dict):
            return {k: KnowledgeDict._clean_data(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [KnowledgeDict._clean_data(v) for v in obj]
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.floating, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.integer, np.int32, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        elif obj is None or isinstance(obj, (str, float, int)):
            return obj
        else:
            return str(obj)

    def to_yaml(self, filepath: str, sort_keys: bool = False) -> None:
        clean_dict = KnowledgeDict._clean_data(dict(self))
        with open(filepath, 'w') as f:
            yaml.safe_dump(clean_dict, f, sort_keys=sort_keys, allow_unicode=True)

    @classmethod
    def from_yaml(cls, filepath: str) -> "KnowledgeDict":
        with open(filepath, 'r') as f:
            data = yaml.safe_load(f)
        kg = cls()
        kg.update(data)
        return kg