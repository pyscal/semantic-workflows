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


def calculate_quasiharmonic_free_energy(
    structure, potential, temperature, temperature_steps=2, cores=4
):
    pr = Project("FreeEnergyQH1")
    pr.remove_jobs(recursive=True)

    refjob = pr.create.job.Lammps("reflmp")
    refjob.structure = structure
    refjob.potential = potential

    phono = pr.create.job.PhonopyJob("phono")
    phono.ref_job = refjob

    quasi = pr.create.job.QuasiHarmonicJob("quasi")
    quasi.ref_job = phono
    quasi.server.cores = cores

    quasi.input["temperature_end"] = temperature
    quasi.input["temperature_steps"] = temperature_steps
    quasi.input["axes"] = ["x", "y", "z"]
    quasi.input["vol_range"] = 0.01
    quasi.input["strains"] = None

    quasi.run()

    ejob = pr.create.job.Lammps("lmp1", delete_existing_job=True)
    ejob.structure = structure
    ejob.potential = potential

    ejob.calc_minimize(pressure=0)
    ejob.run()

    int_energy = ejob.output.energy_pot[-1] / len(structure)
    temps = quasi["output"]["temperatures"][0]
    fes = quasi["output"]["free_energy"][0] + int_energy

    return temps, fes


def calculate_internal_energy(structure, potential):
    pr = Project("InternalEnergy")
    pr.remove_jobs(recursive=True)

    job = pr.create.job.Lammps("lmp1", delete_existing_job=True)
    job.structure = structure
    job.potential = potential

    job.calc_minimize(pressure=0)
    job.run()

    return job.output.energy_pot[-1] / len(structure)
