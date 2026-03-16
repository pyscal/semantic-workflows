sample_template = {
    "id": "sample1",
    "material": {
        "element_ratio": {},
        "crystal_structure": {
            "spacegroup_symbol": None,
            "spacegroup_number": None,
            "unit_cell": {
                "bravais_lattice": None,
                "lattice_parameter": None,
                "angle": [],
            },
        },
    },
    "simulation_cell": {
        "volume": {"value": None},
        "number_of_atoms": None,
        "length": [],
        "vector": [],
        "angle": [],
        "repetitions": [],
        "grain_size": None,
        "number_of_grains": 0,
    },
    "atom_attribute": {
        "position": None,
        "species": None,
        # file-reference mode (preferred for large structures)
        "file_path": None,
        "file_format": None,
        "file_species": None,
    },
    "calculated_property": [],
}

property_template = {
    "id": None,  # optional local ID for math_operation cross-referencing
    "label": None,  # primary name read by atomRDF WorkflowParser
    "basename": None,  # must match ASMO ontology term
    "value": None,
    "unit": None,
    "associate_to_sample": [],
}

workflow_template = {
    "algorithm": None,
    "method": None,
    "xc_functional": None,
    "input_parameter": [],
    "input_sample": [],
    "output_sample": [],
    "output_parameter": [],
    "calculated_property": [],
    "degrees_of_freedom": [],
    "interatomic_potential": {
        "potential_type": None,
        "uri": None,
    },
    "software": [],
    "workflow_manager": {
        "uri": None,
        "version": None,
        "label": None,
    },
    "thermodynamic_ensemble": None,
}
