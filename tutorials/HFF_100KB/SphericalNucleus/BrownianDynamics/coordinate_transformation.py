# This file converts the trajectory with position in OpenMM default unit to reduced unit
import MDAnalysis as mda
import numpy as np
import sys

traj_name = sys.argv[1]

traj = mda.coordinates.LAMMPS.DCDReader("%s"%traj_name)

def scale_down_by_10(ts):
    ts.positions /= 10
    return ts

workflow = [scale_down_by_10]
traj.trajectory.add_transformations(*workflow)

with mda.Writer("reduced_traj.dcd", traj.n_atoms) as W:
    for ts in traj.trajectory:
        W.write(traj)
