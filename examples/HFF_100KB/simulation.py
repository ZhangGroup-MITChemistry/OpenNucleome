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

######################
## initialize the system
model = OpenChrModel(1.0, 0.1, 0.005, 1.0)
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
model.add_class2_bond()
model.add_angle_force()
model.add_softcore()
model.add_ideal_potential(ideal_file)
model.add_type_type_potential(types_file)
model.add_LJ_plain()
model.add_LJ_nuc()
model.add_LJ_spec()
model.add_chr_nuc(chr_nuc_param)
model.add_chr_spec(chr_spec_param)
model.add_chr_lam(chr_lam_param)
model.add_hardwall()
model.add_inter_potential(inter_file)

######################
## perform simulation, in this example, total step = 3,000,000, output to dcd every 2000 steps, and output the energy (similar to thermo in lammps) every 2000 steps

#model.save_system("model_before_simulation_0.xml")

simulation = model.create_simulation(platform_type = "CUDA")
simulation.context.setPositions(model.chr_positions) 

#state = simulation.context.getState(getPositions=True)
#np.savetxt('bead_position.txt',np.array(state.getPositions().value_in_unit(u.nanometer)), fmt='%.6f')
#simulation.minimizeEnergy()

simulation.reporters.append(mdtraj.reporters.DCDReporter('HFF_3e6_every2000.dcd', 2000))
#simulation.context.setPositions(model.chr_positions) #assign a new configuration as the initial configuration, different from the pdb file.

def setVelocity(context):
    sigma = u.sqrt(1.0*u.kilojoule_per_mole / model.chr_system.getParticleMass(1)) 
    velocs = u.Quantity(1.0 * np.random.normal(size=(model.chr_system.getNumParticles(), 3)), u.meter) * (sigma / u.meter)
    context.setVelocities(velocs) 
setVelocity(simulation.context)

simulation.reporters.append(mmapp.statedatareporter.StateDataReporter(sys.stdout, 2000, step=True, 
    potentialEnergy=True, kineticEnergy=True, temperature=True, progress=True, remainingTime=True, separator='\t', totalSteps = 3000000))

for i in range(750):
    simulation.step(4000)
    # Change the type of speckles every 4000 steps, non-equilibrium scheme.

    print("Bonding potential:", simulation.context.getState(getEnergy=True, groups={10}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Softcore potential:", simulation.context.getState(getEnergy=True, groups={11}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Ideal(intra-chromosomal) potential:", simulation.context.getState(getEnergy=True, groups={12}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Compt-Compt potential:", simulation.context.getState(getEnergy=True, groups={13}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Hardwall:", simulation.context.getState(getEnergy=True, groups={14}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Nuc-Spec:", simulation.context.getState(getEnergy=True, groups={15}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Nuc-Nuc:", simulation.context.getState(getEnergy=True, groups={16}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Spec-Spec:", simulation.context.getState(getEnergy=True, groups={17}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Nuc-Chr:", simulation.context.getState(getEnergy=True, groups={18}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Spec-Chr:", simulation.context.getState(getEnergy=True, groups={19}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Chr-Lam:", simulation.context.getState(getEnergy=True, groups={20}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Inter-chromosomal potential:", simulation.context.getState(getEnergy=True, groups={21}).getPotentialEnergy() / u.kilojoule_per_mole)

    for j in np.random.randint(60942,62542,1600):

        if model.compart_type[j] == 5:
            model.compart_type[j] = 6 if random.random() < 0.2 else 5
        else:
            model.compart_type[j] = 5 if random.random() < 0.2 else 6

    for m in range(model.chr_system.getNumParticles()):
        model.chr_system.getForce(7).setParticleParameters(m, [model.compart_type[m]])
    model.chr_system.getForce(7).updateParametersInContext(simulation.context)

np.savetxt('type_final.txt', (np.array(model.compart_type)+1).reshape((-1,1)), fmt='%d')
