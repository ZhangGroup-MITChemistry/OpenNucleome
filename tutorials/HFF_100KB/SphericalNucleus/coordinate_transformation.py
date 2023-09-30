# This file converts the trajectory with position in OpenMM default unit to reduced unit
import MDAnalysis as mda
import numpy as np

traj = mda.coordinates.LAMMPS.DCDReader("HFF_3e6_every2000.dcd")

def scale_down_by_10(ts):
    ts.positions /= 10
    return ts

workflow = [scale_down_by_10]
traj.trajectory.add_transformations(*workflow)

with mda.Writer("DUMP_FILE.dcd", traj.n_atoms) as W:
    for ts in traj.trajectory:
        W.write(traj)
