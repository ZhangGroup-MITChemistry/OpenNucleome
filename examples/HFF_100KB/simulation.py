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
from openChrModel import openChrModel #put openChrModel.py and this simulation.py in the same folder, otherwise, need to import sys and use sys.path to go the openChrModel.py folder


######################
## initialize the system
model = openChrModel(1.0, 0.1, 0.005, 1.0)
PDBFile = "human.pdb"
ICFile = "ideal_chromosome.txt"
TypesTable = "eij_compartment_uniform.txt"
rescalarFile = "nuc_rescaling.txt"
TSA = "TSA_8900.txt" #contains 300 nucleoli beads, 600 speckle beads, and 8000 lamina beads, 8900=300+600+8000
DAM_ID = "DamID_8900.txt"
InterFile = 'inter_chromosome.txt'
model.createSystem(PDBFile) #generate new elements and construct topology as well
radius_of_nucleus = 13.0

######################
## add force field
# model.addFeneBond()
model.addClass2Bond()
model.addAngleForce()
model.addSoftCore()
model.addIdealPotential(ICFile)
model.addType2TypePotential(TypesTable)
model.addLJplain()
model.addLJNuc()
model.addSpeLJ()
model.addNAD(rescalarFile)
model.addTSA(TSA)
model.add_p_DamID(DAM_ID)
model.addparticle_hw()
model.addinter_chrom(InterFile)

######################
## perform simulation, in this example, total step = 2,000,000, output to dcd every 100 steps, and output the energy (similar to thermo in lammps) every 2000 steps
simulation = model.createSimulation(platform_type = "CUDA")
simulation.context.setPositions(model.chrPositions) 

state = simulation.context.getState(getPositions=True)
np.savetxt('bead_position.txt',np.array(state.getPositions().value_in_unit(u.nanometer)), fmt='%.6f')
#simulation.minimizeEnergy()

simulation.reporters.append(mdtraj.reporters.DCDReporter('HFF_3e6_every2000.dcd', 2000))
#simulation.context.setPositions(model.chrPositions) #assign a new configuration as the initial configuration, different from the pdb file.

def setVelocity(context):
    sigma = u.sqrt(1.0*u.kilojoule_per_mole / model.chrSystem.getParticleMass(1)) 
    velocs = u.Quantity(1.0 * np.random.normal(size=(model.chrSystem.getNumParticles(), 3)), u.meter) * (sigma / u.meter)
    context.setVelocities(velocs) 
setVelocity(simulation.context)

simulation.reporters.append(mmapp.statedatareporter.StateDataReporter(sys.stdout, 2000, step=True, 
    potentialEnergy=True, kineticEnergy=True, temperature=True, progress=True, remainingTime=True, separator='\t', totalSteps = 3000000))

for i in range(1500):
    simulation.step(2000)
    state = simulation.context.getState(getPositions=True)
    if (np.amax(np.sqrt(np.sum(np.array(state.getPositions().value_in_unit(u.nanometer))**2, axis=1))) > 13.1):
        break # To ensure that no atom goes outside
    print("Bonding potential:", simulation.context.getState(getEnergy=True, groups={10}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Softcore potential:", simulation.context.getState(getEnergy=True, groups={11}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Ideal(intra) potential:", simulation.context.getState(getEnergy=True, groups={12}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Compt-Compt potential:", simulation.context.getState(getEnergy=True, groups={13}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Hardwall:", simulation.context.getState(getEnergy=True, groups={14}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Nuc-Spec:", simulation.context.getState(getEnergy=True, groups={15}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Nuc-Nuc:", simulation.context.getState(getEnergy=True, groups={16}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Spec-Spec:", simulation.context.getState(getEnergy=True, groups={17}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Nuc-Chr:", simulation.context.getState(getEnergy=True, groups={18}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Spec-Chr:", simulation.context.getState(getEnergy=True, groups={19}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Chr-Lamina:", simulation.context.getState(getEnergy=True, groups={20}).getPotentialEnergy() / u.kilojoule_per_mole)
    print("Inter potential:", simulation.context.getState(getEnergy=True, groups={21}).getPotentialEnergy() / u.kilojoule_per_mole)
