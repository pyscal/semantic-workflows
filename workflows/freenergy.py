from pyiron_atomistics import Project
import pandas as pd


def calculate_free_energy(
    structure, potential, temperature, pressure=0, cores=4, n_switching_steps=10000
):
    pr = Project("FreeEnergy")
    pr.remove_jobs(recursive=True)

    job = pr.create.job.Calphy("c1", delete_existing_job=True)
    job.structure = structure
    job.potential = potential
    job.server.cores = 4

    job.calc_free_energy(
        temperature=temperature,
        pressure=pressure,
        n_switching_steps=n_switching_steps,
        reference_phase="solid",
    )
    job.run()
    return job.output.temperature, job.output.energy_free
