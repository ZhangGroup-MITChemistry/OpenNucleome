from __future__ import print_function
import simtk.openmm as mm
import simtk.openmm.app as mmapp
import simtk.unit as u
import numpy as np
import pandas as pd
import os
import parmed as pmd
import json
import sys
from sys import platform
import mdtraj as md
import mdtraj.reporters
import random
from openNucleome import OpenNucleome

prob_P_dP = 0.2 # Transition probability from P to dP
prob_dP_P = 0.2 # Transition probability from dP to P
transition_freq = 4000
sampling_freq = 2000
dP_type = 6
P_type = 7
total_steps = 3000000

######################
## initialize the system
model = OpenNucleome(1.0, 0.1, 0.005, 1.0) # 1.0: temperature (LJ reduced unit); 0.1: damping coefficient (LJ reduced unit); 0.005: timestep (LJ reduced unit); 1.0: mass_scale
PDB_file = "human.pdb"
model.create_system(PDB_file, membrane_dynamics = False, membrane_bond = None) # Generate new elements and construct topology as well; membrane_dynamics: True for including lamina dynamics, False for excluding lamina dynamics; membrane_bond: A file contains the lamina bond when membrane_dynamics is on.

######################
## add force field

# add the default force field
model.load_default_settings()

# add the customized force field
#force_field = pd.read_csv('input.csv', sep=' ', header=0, index_col=0)
#model.load_customized_settings(force_field)

index_spec_spec_potential = 6
start_spec_index = model.N_chr_nuc+1
end_spec_index = model.N_chr_nuc_spec+1
N_spec = start_spec_index-end_spec_index

######################
## perform simulation, in this example, total step = 3,000,000, output to dcd every 2000 steps, and output the energy (similar to thermo in lammps) every 2000 steps

#model.save_system("model_before_simulation_0.xml")

simulation = model.create_simulation(platform_type = "CUDA")
simulation.context.setPositions(model.chr_positions) 

simulation.minimizeEnergy()

simulation.reporters.append(mdtraj.reporters.DCDReporter('HFF_3e6_every2000.dcd', sampling_freq))

def setVelocity(context):
    sigma = u.sqrt(1.0*u.kilojoule_per_mole / model.chr_system.getParticleMass(1)) 
    velocs = u.Quantity(1.0 * np.random.normal(size=(model.chr_system.getNumParticles(), 3)), u.meter) * (sigma / u.meter)
    context.setVelocities(velocs) 
setVelocity(simulation.context)

simulation.reporters.append(mmapp.statedatareporter.StateDataReporter(sys.stdout, sampling_freq, step=True, 
    potentialEnergy=True, kineticEnergy=True, temperature=True, progress=True, remainingTime=True, separator='\t', totalSteps = total_steps))

for i in range(total_steps//transition_freq):
    simulation.step(transition_freq)
    # Change the type of speckles every 4000 steps, non-equilibrium scheme.

    for j in np.random.randint(start_spec_index, end_spec_index, N_spec): # Do the chemical modification, and change the spec-spec potential on the fly

        if model.compart_type[j] == dP_type-1:
            model.compart_type[j] = P_type-1 if random.random() < prob_dP_P else dP_type-1
        else:
            model.compart_type[j] = dP_type-1 if random.random() < prob_P_dP else P_type-1

    for m in range(model.chr_system.getNumParticles()):
        model.chr_system.getForce(index_spec_spec_potential).setParticleParameters(m, [model.compart_type[m]])
    model.chr_system.getForce(index_spec_spec_potential).updateParametersInContext(simulation.context)

# Keep the final result of spec types in case constructing the configuration for the continuous simulation.
np.savetxt('compt_final_frame.txt', (np.array(model.compart_type)+1).reshape((-1,1)), fmt='%d')
