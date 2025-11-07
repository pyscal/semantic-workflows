from workflows.pyiron.build import bulk, repeat
from workflows.pyiron.thermodynamics import calculate_free_energy
from pyiron_workflow import Workflow

pair_style = "eam/alloy"
pair_coeff = "* * workflows/potentials/Fe_Ack.eam Fe"

wf = Workflow('fe_bcc_0C')
wf.bulk = bulk('Fe', cubic=True)
wf.structure = repeat(wf.bulk, (10,10,10))
wf.free_energy = calculate_free_energy(wf.structure, 
                                  pair_style, pair_coeff,
                                  reference_phase='solid',
                                  temperature=100,
                                  mode='pscale',
                                  pressure = [120000, 150000],
                                  cores=24)
wf.run()

