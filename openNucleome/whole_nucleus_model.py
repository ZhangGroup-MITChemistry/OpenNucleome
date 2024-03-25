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
from .chromosome import Chromosome
from .nucleolus import Nucleolus
from .speckle import Speckle
from .lamina import Lamina

current_path = os.path.abspath(__file__)
father_path = os.path.abspath(os.path.dirname(current_path) + os.path.sep + '.')
config_path = os.path.join(father_path, 'parameters', 'HFF_100KB')

class OpenNucleome:
    '''
    The OpenNucleome class performes the whole nucleus dynamics based on the compartment annotations sequence of chromosomes.

    The simulations can be performed using customed values for the type-to-type, intra, inter-chromosomal, chromosomes-nucleoli, and chromosomes-speckles parameters.

    In our model, each component has specific type. Compartment A: 1, Compartment B: 2, Pericentromeres: 3, Centromeres: 4, Nucleoli: 5, dP Speckles: 6, P Speckles: 7, Lamina: 8
    '''
    def __init__(self, temperature = 1.0, gamma = 0.1, timestep = 0.005, mass_scale = 1.0):
        r'''
        Initialize the whole nucleus model.

        Parameters
        ----------
        temperature (float, required):
            Temperature in reduced units. (Default value = 1.0).

        gamma (float, required):
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

    def create_system(self, PDB_file, flag_membrane = False, lam_bond = None):
        '''
        Create the system with initial configurations.

        Parameters
        ----------

        PDB_file (string, required):
            The path to the chromosome PDB file, including nucleoli, speckles, and lamina.

        flag_membrane (bool, required):
            A flag to control the dynamics of the lamina membrane, True: include the dynamics; False: exclude the dynamics (Default: False)

        lam_bond (string, required):
            The file containing the bond between specific lamina beads. (Default: None)

        '''
        chr_PDB_path = PDB_file
        self.flag_membrane = flag_membrane
        if flag_membrane: self.membrane_bond = np.loadtxt(lam_bond, dtype='int')

        self.chr_PDB = mmapp.PDBFile(chr_PDB_path)

        self.chr_system = mm.System()
        self.chr_positions = []
        self.N_type = 8
        self.chrbead_type = 4

        #record the position information of each bead in self.chr_positions
        for i in range(len(self.chr_PDB.positions)):
        #unit conversion, the default imported PDB position information is angstrom, which will be 10 times smaller for the normal nm calculation.
            self.chr_PDB.positions[i] *= 10
            self.chr_positions.append(self.chr_PDB.positions[i])

        self.N_total = len(self.chr_positions) #total number of chr beads
        self.compart_type = [] #compartment type information
        self.mol_type = [] #mol type information
        self.bead_groups = [[] for _ in range(self.N_type)] # group the indexes of atoms by their types; if you do not want to add some nuclear landmarks, then use the different value. compartment A: 1, compartment B: 2, regions around the centromeres: 3, centromeres: 4, nucleoli: 5, dP speckles: 6, P speckles: 7, lamina: 8
        self.chr_groups = [[] for _ in range(46)] # represent the 46 chromosomes of human
        prev_res_ID = -1
        self.chr_residues = [] #(start, end) index of the atoms for each residue
        start = 0

        for a in self.chr_PDB.topology.atoms():
            #record compartment type info
            self.compart_type.append(int(a.name) - 1)
            self.mol_type.append(int(a.residue.index))

            #add atom into the system with mass = chr_mass
            if not flag_membrane:
                m = 1. if self.compart_type[-1] != self.N_type-1 else 0. # set mass of Lamina beads to zero to freeze them
            else:
                m = 1.
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
        for i in range(self.N_type-1):
            if i<=self.chrbead_type-1:
                self.N_chr += len(self.bead_groups[i])
                self.N_chr_nuc += len(self.bead_groups[i])
                self.N_chr_nuc_spec += len(self.bead_groups[i])
            elif i<=self.chrbead_type:
                self.N_chr_nuc += len(self.bead_groups[i])
                self.N_chr_nuc_spec += len(self.bead_groups[i])
            else:
                self.N_chr_nuc_spec += len(self.bead_groups[i])

        self.chr_residues.append([start, a.index])

        self.generate_element(flag_membrane)

        self.construct_topology(flag_membrane)

        self.chromosome_model = Chromosome(self.chr_system.getNumParticles(), self.N_chr, self.N_type, self.all_bonds, self.compart_type, self.chr_groups, self.bead_groups, self.mol_type)
        self.speckle_model = Speckle(self.chr_system.getNumParticles(), self.N_type, self.all_bonds, self.compart_type, self.chr_groups, self.bead_groups, self.mol_type)
        self.nucleolus_model = Nucleolus(self.chr_system.getNumParticles(), self.N_type, self.N_chr, self.N_chr_nuc, self.all_bonds, self.compart_type, self.chr_groups, self.bead_groups, self.mol_type)
        self.lamina_model = Lamina(self.chr_system.getNumParticles(), self.N_type, self.N_chr_nuc_spec, self.all_bonds, self.compart_type, self.chr_groups, self.bead_groups, self.mol_type)

    def generate_element(self, flag_membrane):
        '''
        Generate elements for each coarse-grained bead.

        This function is called automatically when creating a system.

        Parameters
        ----------

        flag_membrane (bool, required):
            A flag to control the dynamics of the lamina membrane, True: include the dynamics; False: exclude the dynamics (Default: False)
        '''
        name_to_element = ['ASP', 'GLU', 'HIS', 'LYS', 'ARG', 'ASN', 'GLN', 'PRO']
        self.Elements = []
        for i in range(self.N_type):
            if not flag_membrane:
                m = 1. if i!= self.N_type-1 else 0.
            else:
                m = 1.
            self.Elements.append(mmapp.Element(1000+i, name_to_element[i], name_to_element[i], self.chr_mass * m))

    def construct_topology(self, flag_membrane):
        '''
        Construct the topology for the system.

        This function is called automatically when creating a system.

        Parameters
        ----------

        flag_membrane (bool, required):
            A flag to control the dynamics of the lamina membrane, True: include the dynamics; False: exclude the dynamics (Default: False)
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

        if flag_membrane:

            self.membrane_bond += self.N_chr_nuc_spec

            for i, j in self.membrane_bond:
                self.chr_topo.addBond(i, j)

        #store all bonds in all_bonds
        self.all_bonds = []
        for b in self.chr_topo.bonds():
            self.all_bonds.append((b.atom1, b.atom2))

    def add_chromosome_potential(self, flag, ideal_file, types_file, inter_file):
        '''
        Add the potential related to the genome in the system.

        Parameters
        ----------

        flag (dict, required):
            A dict contains the switch of each chromosome-related potential.

        ideal_file (string, required):
            The path to the ideal potential scaling factor txt file.

        types_file (string, required):
            The path to the type-to-type potential scaling factor txt file.

        inter_file (string, required):
            The path to the inter-chromosomal potential scaling factor txt file.
        '''
        if flag['bond']:
            self.chr_system.addForce(self.chromosome_model.add_class2_bond(10))
        if flag['angle']:
            self.chr_system.addForce(self.chromosome_model.add_angle_force())
        if flag['softcore']:
            self.chr_system.addForce(self.chromosome_model.add_softcore(11))
        if flag['ideal']:
            self.chr_system.addForce(self.chromosome_model.add_ideal_potential(ideal_file, 12))
        if flag['compt']:
            self.chr_system.addForce(self.chromosome_model.add_type_type_potential(types_file, 13))
        if flag['inter']:
            self.chr_system.addForce(self.chromosome_model.add_inter_potential(inter_file, 14))

    def add_speckle_potential(self, flag, chr_spec_param):
        '''
        Add the potential related to the nuclear speckles in the system.

        Parameters
        ----------

        flag (dict, required):
            A dict contains the switch of each speckle-related potential

        chr_spec_param (string, required):
            The path to the potential scaling factor txt file.
        '''
        if flag['spec-spec']:
            self.chr_system.addForce(self.speckle_model.add_spec_spec(15))
        if flag['spec-chrom']:
            self.chr_system.addForce(self.speckle_model.add_chr_spec(chr_spec_param, 16))

    def add_nucleolus_potential(self, flag, rescalar_file):
        '''
        Add the potential related to the nucleoli in the system.

        Parameters
        ----------

        flag (dict, required):
            A dict contains the switch of each nucleolus-related potential

        rescalar_file (string, required):
            The path to the txt file of rescaling factors (SPIN state probabilities).
        '''
        if flag['nuc-nuc']:
            self.chr_system.addForce(self.nucleolus_model.add_nuc_nuc(17))
        if flag['nuc-spec']:
            self.chr_system.addForce(self.nucleolus_model.add_nuc_spec(18))
        if flag['nuc-chrom']:
            self.chr_system.addForce(self.nucleolus_model.add_chr_nuc(rescalar_file, 19))

    def add_lamina_potential(self, flag, chr_lam_param):
        '''
        Add the potential related to the nuclear lamina in the system.

        Parameters
        ----------

        flag (dict, required):
            A dict contains the switch of each lamina-related potential

        chr_lam_param (string, required):
            The path to the potential scaling factor txt file.
        '''
        if flag['lam-chrom']:
            self.chr_system.addForce(self.lamina_model.add_chr_lam(chr_lam_param, 20))
        if flag['hard-wall']:
            self.chr_system.addForce(self.lamina_model.add_hardwall(21))
        if flag['lam-lam']:
            self.chr_system.addForce(self.lamina_model.add_lam_lam(22))
        if flag['squeeze_nucleus']:
            self.chr_system.addForce(self.lamina_model.add_squeeze_nucleus(23, k = flag['strength']))

    def load_default_settings(self):
        '''
        Load default force field settings
        '''
        dict_chrom_ff = {'bond':True, 'angle':True, 'softcore':True, 'ideal':True, 'compt':True, 'inter':True}
        dict_spec_ff = {'spec-spec':True, 'spec-chrom':True}
        dict_nuc_ff = {'nuc-nuc':True, 'nuc-spec':True, 'nuc-chrom':True}
        dict_lam_ff = {'lam-chrom':True, 'hard-wall':True, 'lam-lam':self.flag_membrane, 'squeeze_nucleus':self.flag_membrane}

        ideal_param_file = os.path.join(config_path, "ideal_param_file.txt")
        compt_param_file = os.path.join(config_path, "compt_param_file.txt")
        interchr_param_file = os.path.join(config_path, "interchr_param_file.txt")
        chr_nuc_param = os.path.join(config_path, "chr_nuc_param.txt")
        chr_spec_param = os.path.join(config_path, "chr_spec_param.txt")
        chr_lam_param = os.path.join(config_path, "chr_lam_param.txt")

        self.add_chromosome_potential(dict_chrom_ff, ideal_param_file, compt_param_file, interchr_param_file)
        self.add_speckle_potential(dict_spec_ff, chr_spec_param)
        self.add_nucleolus_potential(dict_nuc_ff, chr_nuc_param)
        self.add_lamina_potential(dict_lam_ff, chr_lam_param)

    def load_customized_settings(self, force_field, param_folder, k = 1.0):
        '''
        Load customized force field settings

        Parameters
        ----------
        force_field (Pandas Dataframe, required):
            Store the flag of all the potentials and the parameter file names

        param_folder (string, required):
            Save the customized parameters used in the simulations

        k (float, required):
            The strength of force squeezing the nucleus (Default: 1.0)
        '''
        dict_chrom_ff = {'bond':force_field.loc['chromosome']['bond'], 'angle':force_field.loc['chromosome']['angle'], 'softcore':force_field.loc['chromosome']['softcore'], 'ideal':force_field.loc['chromosome']['ideal'], 'compt':force_field.loc['chromosome']['compt'], 'inter':force_field.loc['chromosome']['inter']}
        dict_spec_ff = {'spec-spec':force_field.loc['speckle']['spec-spec'], 'spec-chrom':force_field.loc['speckle']['spec-chrom']}
        dict_nuc_ff = {'nuc-nuc':force_field.loc['nucleolus']['nuc-nuc'], 'nuc-spec':force_field.loc['nucleolus']['nuc-spec'], 'nuc-chrom':force_field.loc['nucleolus']['nuc-chrom']}
        dict_lam_ff = {'lam-chrom':force_field.loc['lamina']['lam-chrom'], 'hard-wall':force_field.loc['lamina']['hard-wall'], 'lam-lam':self.flag_membrane, 'squeeze_nucleus':self.flag_membrane, 'strength': k}
        ideal_param_file = force_field.loc['chromosome']['ideal_param_file']
        if ideal_param_file:
            ideal_param_file = os.path.join('.', param_folder, ideal_param_file)
        compt_param_file = force_field.loc['chromosome']['compt_param_file']
        if compt_param_file:
            compt_param_file = os.path.join('.', param_folder, compt_param_file)
        interchr_param_file = force_field.loc['chromosome']['interchr_param_file']
        if interchr_param_file:
            interchr_param_file = os.path.join('.', param_folder, interchr_param_file)
        chr_nuc_param = force_field.loc['nucleolus']['chr_nuc_param']
        if chr_nuc_param:
            chr_nuc_param = os.path.join('.', param_folder, chr_nuc_param)
        chr_spec_param = force_field.loc['speckle']['chr_spec_param']
        if chr_spec_param:
            chr_spec_param = os.path.join('.', param_folder, chr_spec_param)
        chr_lam_param = force_field.loc['lamina']['chr_lam_param']
        if chr_lam_param:
            chr_lam_param = os.path.join('.', param_folder, chr_lam_param)
        self.add_chromosome_potential(dict_chrom_ff, ideal_param_file, compt_param_file, interchr_param_file)
        self.add_speckle_potential(dict_spec_ff, chr_spec_param)
        self.add_nucleolus_potential(dict_nuc_ff, chr_nuc_param)
        self.add_lamina_potential(dict_lam_ff, chr_lam_param)

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
        platform_type (string, required):
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
