from __future__ import print_function
import simtk.openmm as mm
import simtk.openmm.app as mmapp
import simtk.unit as u
import numpy as np
import os
import parmed as pmd
import json
import sys
from sys import platform
import mdtraj as md
import mdtraj.reporters
import random
from openNucleome import OpenChrModel

prob_P_dP = 0.2 # Transition probability from P to dP
prob_dP_P = 0.2 # Transition probability from dP to P
transition_freq = 4000
sampling_freq = 2000
dP_type = 6
P_type = 7

######################
## initialize the system
model = OpenChrModel(1.0, 0.1, 0.005, 1.0) #1.0: temperature (LJ reduced unit); 0.1: damping coefficient (LJ reduced unit); 0.005: timestep (LJ reduced unit); 1.0: mass_scale#
PDB_file = "human.pdb"
ideal_file = "ideal_param.txt"
types_file = "types_param.txt"
chr_nuc_param = "chr_nuc_param.txt"
chr_spec_param = "chr_spec_param.txt"
chr_lam_param = "chr_lam_param.txt"
inter_file = 'inter_param.txt'
model.create_system(PDB_file) #generate new elements and construct topology as well

######################
## add force field
dict_chrom = {'bond':True, 'angle':True, 'softcore':True, 'ideal':True, 'compt':True, 'inter':True}
dict_spec = {'spec-spec':True, 'spec-chrom':True}
dict_nuc = {'nuc-nuc':True, 'nuc-spec':True, 'nuc-chrom':True}
dict_lam = {'lam-chrom':True, 'hard-wall':True}
model.add_chromosome_potential(dict_chrom, ideal_file, types_file, inter_file)
model.add_speckle_potential(dict_spec, chr_spec_param)
model.add_nucleolus_potential(dict_nuc, chr_nuc_param)
model.add_lamina_potential(dict_lam, chr_lam_param)

index_spec_spec_potential = 6
start_spec_index = model.N_chr_nuc+1
end_spec_index = model.N_chr_nuc_spec+1
N_index = start_spec_index-end_spec_index

######################
## perform simulation, in this example, total step = 3,000,000, output to dcd every 2000 steps, and output the energy (similar to thermo in lammps) every 2000 steps

#model.save_system("model_before_simulation_0.xml")

simulation = model.create_simulation(platform_type = "CUDA")
simulation.context.setPositions(model.chr_positions) 

#simulation.minimizeEnergy()

simulation.reporters.append(mdtraj.reporters.DCDReporter('HFF_3e6_every2000.dcd', sampling_freq))

def setVelocity(context):
    sigma = u.sqrt(1.0*u.kilojoule_per_mole / model.chr_system.getParticleMass(1)) 
    velocs = u.Quantity(1.0 * np.random.normal(size=(model.chr_system.getNumParticles(), 3)), u.meter) * (sigma / u.meter)
    context.setVelocities(velocs) 
setVelocity(simulation.context)

simulation.reporters.append(mmapp.statedatareporter.StateDataReporter(sys.stdout, sampling_freq, step=True, 
    potentialEnergy=True, kineticEnergy=True, temperature=True, progress=True, remainingTime=True, separator='\t', totalSteps = 3000000))

for i in range(3000000//transition_freq):
    simulation.step(transition_freq)
    # Change the type of speckles every 4000 steps, non-equilibrium scheme.

    for j in np.random.randint(start_spec_index, end_spec_index, N_index):

        if model.compart_type[j] == dP_type-1:
            model.compart_type[j] = P_type-1 if random.random() < prob_dP_P else dP_type-1
        else:
            model.compart_type[j] = dP_type-1 if random.random() < prob_P_dP else P_type-1

    for m in range(model.chr_system.getNumParticles()):
        model.chr_system.getForce(index_spec_spec_potential).setParticleParameters(m, [model.compart_type[m]])
    model.chr_system.getForce(index_spec_spec_potential).updateParametersInContext(simulation.context)

np.savetxt('type_final.txt', (np.array(model.compart_type)+1).reshape((-1,1)), fmt='%d')
