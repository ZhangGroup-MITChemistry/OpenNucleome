import MDAnalysis as mda
import numpy as np
import sys

def coor_transformation(traj_data, output_traj):
    """
    The function converts the trajectory with position in OpenMM default unit to reduced unit

    Parameters
    ----------
    traj_data (string, required):
        The trajectory file from a simulation

    output_traj (string, required):
        The trajectory with the reduced units
    """

    traj = mda.coordinates.LAMMPS.DCDReader(traj_data)

    def scale_down_by_10(ts):
        ts.positions /= 10
        return ts

    workflow = [scale_down_by_10]
    traj.trajectory.add_transformations(*workflow)

    with mda.Writer(output_traj, traj.n_atoms) as W:
        for ts in traj.trajectory:
            W.write(traj)
