### Knowledge-based workflows of mechanical behavior with atomistic simulations

Modular workflows for mechanical and thermodynamic properties of body-centered cubic iron as an example, covering tasks such as equation of state, elastic tensors, mechanical loading, defect energetics, and nanoindentation. Each workflow is implemented using pyiron workflow nodes that automatically generate semantic metadata aligned with the CMSO and ASMO ontologies, ensuring a consistent and machine-readable description of the computational sample, modelling method, and provenance.

Users can run atomistic simulations using the LAMMPS molecular dynamics software, with examples provided for embedded-atom method (EAM) potentials. The demonstrator also allows users to inspect the conceptual dictionary associated with each workflow, which captures workflow inputs, outputs, and semantic annotations. This interactive setup provides a lightweight environment for exploring how knowledge-based workflows are constructed, executed, and semantically enriched, without requiring local installation or configuration of simulation tools.

##### Highlights
- pyiron: used as a workflow manager
- atomRDF: Automated semantic annotation of atomistic simulation outputs.
- ASMO and CMSO: Standardized ontology framework for atomistic simulation and computational 
