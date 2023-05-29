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
from openNucleome import Chromosome, Speckle, Nucleolus, Lamina

class OpenChrModel:
    '''
    The OpenChrModel class performes the whole nucleus dynamics based on the compartment annotations sequence of chromosomes.

    The simulations can be performed using customed values for the type-to-type, intra, inter-chromosomal, chromosomes-nucleoli, and chromosomes-speckles parameters.

    In our model, each component has specific type. Compartment A: 1, Compartment B: 2, Pericentromeres: 3, Centromeres: 4, Nucleoli: 5, dP Speckles: 6, P Speckles: 7, Lamina: 8
    '''
    def __init__(self, temperature = 1.0, gamma = 0.1, timestep = 0.005, mass_scale=1.0):
        r'''
        Initialize the whole nucleus model.

        Parameters
        ----------
        temperature (float, required) :
            Temperature in reduced units. (Default value = 1.0).

        gamma (float, required) :
            Friction/Damping constant in units of reciprocal time (:math:`1/\tau`) (Default value = 0.1).

        timestep (float, required):
            Simulation time step in units of :math:`\tau` (Default value = 0.005).

        mass_scale (float, required):
            Mass scale used in units of :math:`\mu` (Default value = 1.0).

        '''
        self.temperature = temperature / 0.008314 # temperature in reduced unit
        self.gamma = gamma # fiction coefficient (1/time)
        self.timestep = timestep  # timestep in reduced unit (time)
        self.kB = u.BOLTZMANN_CONSTANT_kB * u.AVOGADRO_CONSTANT_NA
        self.kT = self.kB * self.temperature
        self.chr_mass = 10.0 * u.amu * mass_scale

    def create_system(self, PDB_file):
        '''
        Create the system with initial configurations.

        Parameters
        ----------

        PDB_file (string, required):
            The path to the chromosome PDB file, including nucleoli, speckles, and lamina.
        '''
        chr_PDB_path = PDB_file
        self.chr_PDB = mmapp.PDBFile(chr_PDB_path)

        self.chr_system = mm.System()
        self.chr_positions = []
        self.N_type = 8

        #record the position information of each bead in self.chr_positions
        for i in range(len(self.chr_PDB.positions)):
        #unit conversion, the default imported PDB position information is angstrom, which will be 10 times smaller for the normal nm calculation.
            self.chr_PDB.positions[i] *= 10
            self.chr_positions.append(self.chr_PDB.positions[i])

        self.N_total = len(self.chr_positions) #total number of chr beads
        self.compart_type = [] #compartment type information
        self.mol_type = [] #mol type information
        self.bead_groups = [[] for _ in range(8)] # group the indexes of atoms by their types; if you do not want to add some nuclear landmarks, then use the different value. compartment A: 1, compartment B: 2, regions around the centromeres: 3, centromeres: 4, nucleoli: 5, dP speckles: 6, P speckles: 7, lamina: 8
        self.chr_groups = [[] for _ in range(46)] # represent the 46 chromosomes of human
        prev_res_ID = -1
        self.chr_residues = [] #(start, end) index of the atoms for each residue
        start = 0

        for a in self.chr_PDB.topology.atoms():
            #record compartment type info
            self.compart_type.append(int(a.name) - 1)
            self.mol_type.append(int(a.residue.index))

            #add atom into the system with mass = chr_mass
            m = 1. if self.compart_type[-1] != 7 else 0. # set mass of Lamina beads to zero to freeze them
            self.chr_system.addParticle(self.chr_mass * m)

            #find the start and end index in each residue, self.chrResidue will be (start, end, residue) for each residue after self.construct_topology
            if a.residue.index != prev_res_ID and prev_res_ID != -1:
                self.chr_residues.append([start, a.index - 1])
                start = a.index
            prev_res_ID = a.residue.index

            #group the indexes of atoms by their types
            self.bead_groups[int(a.name) - 1].append(a.index)
            if a.residue.index < 46:
                self.chr_groups[int(a.residue.index)].append(a.index)

        self.N_chr = 0
        self.N_chr_nuc = 0
        self.N_chr_nuc_spec = 0
        for i in range(7):
            if i<=3:
                self.N_chr += len(self.bead_groups[i])
                self.N_chr_nuc += len(self.bead_groups[i])
                self.N_chr_nuc_spec += len(self.bead_groups[i])
            elif i<=4:
                self.N_chr_nuc += len(self.bead_groups[i])
                self.N_chr_nuc_spec += len(self.bead_groups[i])
            else:
                self.N_chr_nuc_spec += len(self.bead_groups[i])

        self.chr_residues.append([start, a.index])

        self.generate_element()

        self.construct_topology()

        chromosome_model = Chromosome(self.chr_system.getNumParticles(), self.N_type, self.all_bonds, self.compart_type, self.chr_groups, self.bead_groups, self.mol_type)
        speckle_model = Speckle(self.chr_system.getNumParticles(), self.N_type, self.all_bonds, self.compart_type, self.chr_groups, self.bead_groups, self.mol_type)
        nucleolus_model = Nucleolus(self.chr_system.getNumParticles(), self.N_type, self.all_bonds, self.compart_type, self.chr_groups, self.bead_groups, self.mol_type)
        lamina_model = Lamina(self.chr_system.getNumParticles(), self.N_type, self.all_bonds, self.compart_type, self.chr_groups, self.bead_groups, self.mol_type)

    def generate_element(self):
        '''
        Generate elements for each coarse-grained bead.

        This function is called automatically when creating a system.
        '''
        name_to_element = ['ASP', 'GLU', 'HIS', 'LYS', 'ARG', 'ASN', 'GLN', 'PRO']
        self.Elements = []
        for i in range(8):
            m = 1. if i!= 7 else 0.
            self.Elements.append(mmapp.Element(1000+i, name_to_element[i], name_to_element[i], self.chr_mass * m))

    def construct_topology(self):
        '''
        Construct the topology for the system.

        This function is called automatically when creating a system.
        '''
        # Construct topology, add chain, residues, atoms, bonds
        self.chr_topo = mmapp.topology.Topology()
        chr_chain = self.chr_topo.addChain('0')
        for i in range(len(self.chr_residues)):
            curr_res_idx = self.chr_topo.addResidue(str(i), chr_chain)
            self.chr_residues[i].append(curr_res_idx)
            for j in range(self.chr_residues[i][0], self.chr_residues[i][1] + 1):
                self.chr_topo.addAtom(str(j), self.Elements[self.compart_type[j]], curr_res_idx)
                if j != self.chr_residues[i][0]:
                    self.chr_topo.addBond(j - 1, j) #added atom index, instead of atom item here

        #store all bonds in all_bonds
        self.all_bonds = []
        for b in self.chr_topo.bonds():
            self.all_bonds.append((b.atom1, b.atom2))

    def save_system(self, system_xml):
        '''
        Save the simulation state to readable xml format.
        '''
        with open(system_xml, 'w') as output_writer:
            output_writer.write(mm.XmlSerializer.serialize(self.chr_system))

    def create_simulation(self, platform_type = 'CPU'):
        '''
        Specify a platform when creating the simulation.

        Parameters
        ----------
        platform_type (string, required) :
            The CPU platform is set as the default, while CUDA, OpenCL, and Reference are available as optional platforms.
        '''
        integrator = mm.LangevinIntegrator(self.temperature, self.gamma, self.timestep)
#       integrator = mm.VerletIntegrator(self.timestep) # NVE

        if platform_type == 'CUDA':
            platform = mm.Platform.getPlatformByName('CUDA')
            properties = {'CudaPrecision': 'double'}
            simulation = mmapp.Simulation(self.chr_topo, self.chr_system, integrator, platform, properties)
        elif platform_type == 'OpenCL':
            platform = mm.Platform.getPlatformByName('OpenCL')
            properties = {'OpenCLPrecision': 'double'}
            simulation = mmapp.Simulation(self.chr_topo, self.chr_system, integrator, platform, properties)
        elif platform_type == 'Reference':
            platform = mm.Platform.getPlatformByName('Reference')
            simulation = mmapp.Simulation(self.chr_topo, self.chr_system, integrator, platform)
        elif platform_type == 'CPU':
            platform = mm.Platform.getPlatformByName('CPU')
            simulation = mmapp.Simulation(self.chr_topo, self.chr_system, integrator, platform)
        else:
            print("platform_type can be either CUDA, OpenCL, or CPU")
        return simulation
