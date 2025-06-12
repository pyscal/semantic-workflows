from pyiron_atomistics import Project
import pandas as pd


def calculate_ev_curves(structure, potential, num_points=5):
    pr = Project("E-V")
    pr.remove_jobs(recursive=True)
    ref_job = pr.create.job.Lammps(
        "j1", delete_existing_job=True, delete_aborted_job=True
    )
    ref_job.structure = structure
    ref_job.potential = potential
    ref_job.calc_minimize()
    murn_job = ref_job.create_job(pr.job_type.Murnaghan, "murn_job5")
    murn_job.input["num_points"] = num_points
    murn_job.run(delete_existing_job=True)
    murn_job.plot()
