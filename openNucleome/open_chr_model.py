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
        self.chr_PDB = mmapp.PDB_file(chr_PDB_path)

        self.chr_system = mm.System()
        self.chr_positions = []

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

    def generate_element(self):
        '''
        Generate elements for each coarse-grained bead.

        This function is called automatically when creating a system.
        '''
        name_to_element = ['ASP', 'GLU', 'HIS', 'LYS', 'ARG', 'ASN', 'GLN', 'PRO']
        self.Elements = []
        for i in range(8):
            m = 1. if i!= 7 else 0.
#           self.Elements.append(mmapp.element.Element.getBySymbol(name_to_element[i]))
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

    # The diameter of 100KB bead is 0.5; 1.0 for 1MB bead. For example, I set sigma_fene as 1.0, so currently, these parameters are for 1MB bead.
    def add_fene_bond(self, kfb = 30.0, r0 = 1.5, epsilon_fene = 1.0, sigma_fene = 1.0, cutoff_fene = 1.12):
        '''
        Internal function that inits FENE bond force. The default parameters value is for the 1MB resolution.
        '''
        bond = mm.CustomBondForce("- 0.5 * kfb * r0 * r0 * log(1-(r/r0)*(r/r0)) + (4 * epsilon_fene * ((sigma_fene/r)^12 - (sigma_fene/r)^6) + epsilon_fene) * step(cutoff_fene - r)")
        bond.addGlobalParameter("kfb", kfb) 
        bond.addGlobalParameter("r0", r0) 
        bond.addGlobalParameter("epsilon_fene", epsilon_fene) 
        bond.addGlobalParameter("sigma_fene", sigma_fene) 
        bond.addGlobalParameter("cutoff_fene", cutoff_fene) 
        bond.setForceGroup(10)

        bond.setUsesPeriodicBoundaryConditions(False)#avoid periodic boundary
        for b in self.chr_topo.bonds():
            bond.addBond(b.atom1, b.atom2)
        self.chr_system.addForce(bond)

    def add_class2_bond(self, r0c = 0.5, kc2 = 20.0, kc3 = 20.0, kc4 = 20.0):
        '''
        Internal function that inits Class2 bond force. The default parameters value is for the 100KB resolution. (Recommended)
        '''
        class2_bond = mm.CustomBondForce("kc2 * (r-r0c) ^ 2 + kc3 * (r-r0c) ^ 3 + kc4 * (r-r0c) ^ 4")
        class2_bond.addGlobalParameter("r0c", r0c) 
        class2_bond.addGlobalParameter("kc2", kc2) 
        class2_bond.addGlobalParameter("kc3", kc3) 
        class2_bond.addGlobalParameter("kc4", kc4) 
        class2_bond.setForceGroup(10)

        class2_bond.setUsesPeriodicBoundaryConditions(False)#avoid periodic boundary
        for b in self.chr_topo.bonds():
            class2_bond.addBond(b.atom1, b.atom2)
        self.chr_system.addForce(class2_bond)

    def add_angle_force(self, ka = 2.0):
        '''
        Add an angular potential between bonds connecting beads.

        Parameters
        ----------
        ka (float, required) :
            angle coefficients, default = 2.0
        '''
        angle = mm.CustomAngleForce("ka*(1 + cos(theta))")
        angle.addGlobalParameter("ka", ka) #k_a in reduced unit

        #Apply to all angle triplets
        bond1 = None
        for bond2 in self.chr_topo.bonds():
            if bond1 is None:
                bond1 = bond2
                continue
            if bond1.atom2 == bond2.atom1:
                angle.addAngle(bond1.atom1, bond1.atom2, bond2.atom2)
            bond1 = bond2

        angle.setUsesPeriodicBoundaryConditions(False) #avoid periodic boundary
        self.chr_system.addForce(angle)

    def add_softcore(self, epsilon_sc = 1.0, sigma_sc = 0.5, cutoff_sc = 0.56, Ecut_sc = 4.0):
        r'''
        We plug the tanh potential into the WCA potential to represent the excluded volume effect;

        The softcore potential will only be applied to the chromosomes (type 1, 2, 3, 4).

        When :math:`r_{sc} < r < r_{cutoff}`, use :math:`E_{LJ} = 4\epsilon_{sc}\left(\left(\frac{\sigma_{sc}}{r}\right)^{12} -\left(\frac{\sigma_{sc}}{r}\right)^{6}   \right) + \epsilon_{rc}`; :math:`r < r_{sc}`, use :math:`E = \frac{1}{2}E_{cut}\left(1+\text{tanh}(\frac{2E_{LJ}}{E_{cut}}-1) \right)`, where :math:`r_{sc} \approx 0.485`.

        Parameters
        ----------
        epsilon_sc (float, required) :
            energy units (Default value = 1.0).

        sigma_sc (float, required) :
            distance units, which represents the diameter of a 100KB-resolution bead (Default value = 0.5).

        cutoff_sc (float, required) :
            distance units, the cutoff (Default value = 0.5*1.12).

        Ecut_sc (float, required) :
            energy units, the energy cost for the chain passing (Default value = 4.0).
        '''
        Ecut_sc = Ecut_sc * epsilon_sc
        r0_sc = sigma_sc*(((0.5*Ecut_sc)/(4.0*epsilon_sc) - 0.25 +((0.5)**(2.0)))**(1.0/2.0) +0.5)**(-1.0/6.0)
 
        repul_energy = ("LJ*step(r-r0_sc)*step(cutoff_sc-r)*step(3.5-max(c1,c2))"
                              " + step(3.5-max(c1,c2))*step(r0_sc-r)*0.5*Ecut_sc*(1.0+tanh((2.0*LJ/Ecut_sc)-1.0));"
                              "LJ = 4.0 * epsilon_sc * ((sigma_sc/r)^12 - (sigma_sc/r)^6) + epsilon_sc")

        soft_core = mm.CustomNonbondedForce(repul_energy)
        soft_core.addGlobalParameter('epsilon_sc', epsilon_sc)
        soft_core.addGlobalParameter('sigma_sc', sigma_sc)
        soft_core.addGlobalParameter('Ecut_sc', Ecut_sc)
        soft_core.addGlobalParameter('r0_sc', r0_sc)
        soft_core.addGlobalParameter('cutoff_sc', cutoff_sc)
        soft_core.setCutoffDistance(cutoff_sc)

        soft_core.setForceGroup(11)
        soft_core.addPerParticleParameter("c")

        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            soft_core.addExclusion(idx1, idx2)

        #avoid periodic boundary
        soft_core.setNonbondedMethod(soft_core.CutoffNonPeriodic)

        #apply the potential to all beads
        for i in range(self.chr_system.getNumParticles()):
            soft_core.addParticle([self.compart_type[i]])

        self.chr_system.addForce(soft_core)

    def add_ideal_potential(self, ideal_file, rc_tanh_ideal = 0.54, dend_ideal = 1000, cutoff_ideal = 2.0):
        r'''
        Add the Intra-chromosomal ideal potential using custom values for interactions between beads separated by a genomic distance :math:`d`. 

        The parameter rc_tanh_ideal is part of the probability of crosslink function 

        .. math::

            f(r_{i,j})  & = \frac{1}{2}\left(1+\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right) \\
                        & + \frac{1}{2}\left(1-\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right)\left(\frac{r_c}{r_{i,j}}\right)^4,

        where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        ideal_file (string, required) :
            The path to the ideal potential scaling factor txt file.

        dend_ideal (int, required) : 
            Cutoff of genomic separation (Default value = 1000 for the 100KB model).

        rc_tanh_ideal (float, required) : 
            :math:`r_c` in the equation, a good estimation of the bond length during the simulation (Default value = 0.54 for the 100KB model).

        cutoff_ideal (float, required) : 
            Cutoff of Ideal potential (Default value = 2.0 for the 100KB model).
        '''
        energy_ideal = ("ideal_list(d)*step(dend_ideal-d)*step(cutoff_ideal-r)*(tanh1_ideal+tanh2_ideal*(rc_tanh_ideal/r)^4-1/2*(1.+tanh((rc_tanh_ideal-cutoff_ideal)^(-5)+5*(rc_tanh_ideal-cutoff_ideal)))-1/2*(1.-tanh((rc_tanh_ideal-cutoff_ideal)^(-5)+5*(rc_tanh_ideal-cutoff_ideal)))*(rc_tanh_ideal/cutoff_ideal)^4);"
                        "tanh1_ideal = 1/2*(1.+tanh((rc_tanh_ideal-r)^(-5)+5*(rc_tanh_ideal-r)));"
                        "tanh2_ideal = 1/2*(1.-tanh((rc_tanh_ideal-r)^(-5)+5*(rc_tanh_ideal-r)));"
                        "d = abs(idx2-idx1)")

        ideal = mm.CustomNonbondedForce(energy_ideal)
        ideal_alpha = np.loadtxt(ideal_file)[:,1] 

        tab_ideal_list = mm.Discrete1DFunction(ideal_alpha)
        ideal.addTabulatedFunction('ideal_list', tab_ideal_list) 
        ideal.addGlobalParameter('dend_ideal', dend_ideal)
        ideal.addGlobalParameter('rc_tanh_ideal', rc_tanh_ideal) 
        ideal.addGlobalParameter('cutoff_ideal', cutoff_ideal) 
        ideal.setForceGroup(12)

        ideal.setCutoffDistance(cutoff_ideal)

        ideal.addPerParticleParameter("idx")

        #apply within each chromosome
        for i in range(len(self.chr_groups)):
            ideal.addInteractionGroup(self.chr_groups[i], self.chr_groups[i])

        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            ideal.addExclusion(idx1, idx2)

        #avoid periodic boundary
        ideal.setNonbondedMethod(ideal.CutoffNonPeriodic)

        #apply the potential to all beads
        for i in range(self.chr_system.getNumParticles()):
            ideal.addParticle([i])

        self.chr_system.addForce(ideal)

    def add_type_type_potential(self, types_file, rc_tanh_type = 0.54, cutoff_type = 2.0):
        r'''
        Add the type-to-type potential using custom values for interactions between the chromatin types. 

        The parameter rc_tanh_type is part of the probability of crosslink function

        .. math::

            f(r_{i,j})  & = \frac{1}{2}\left(1+\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right) \\
                        & + \frac{1}{2}\left(1-\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right)\left(\frac{r_c}{r_{i,j}}\right)^4,

        where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        types_file (string, required) :
            The path to the type-to-type potential scaling factor txt file.

        rc_tanh_type (float, required) : 
            :math:`r_c` in the equation, a good estimation of the bond length during the simulation (Default value = 0.54 for the 100KB model).
        
        cutoff_type (float, required) :
            Cutoff of Type-Type potential (Default value = 2.0 for the 100KB model).
        '''

        energy_type = ("map_type(t1,t2)*step(cutoff_type-r)*(tanh1_LP+tanh2_LP*(rc_tanh_type/r)^4-1/2*(1.+tanh((rc_tanh_type-cutoff_type)^(-5)+5*(rc_tanh_type-cutoff_type)))-1/2*(1.-tanh((rc_tanh_type-cutoff_type)^(-5)+5*(rc_tanh_type-cutoff_type)))*(rc_tanh_type/cutoff_type)^4);"
                        "tanh1_LP = 1/2*(1.+tanh((rc_tanh_type-r)^(-5)+5*(rc_tanh_type-r)));"
                        "tanh2_LP = 1/2*(1.-tanh((rc_tanh_type-r)^(-5)+5*(rc_tanh_type-r)));")

        cross_type = mm.CustomNonbondedForce(energy_type)

        cross_type.addGlobalParameter('rc_tanh_type', rc_tanh_type)
        cross_type.addGlobalParameter('cutoff_type', cutoff_type)
        cross_type.setCutoffDistance(cutoff_type)
        
        cross_type.setForceGroup(13)
        #generate the above shown 8*8 matrix, 8 is the number of types
        tab = np.zeros((8, 8))
        with open(types_file, "r") as fin:
            for line in fin.readlines():
                line = line.split()
                i = int(line[1]) - 1
                j = int(line[2]) - 1
                tab[i, j] = tab[j, i] = float(line[1])

        tab = list(np.ravel(tab))
        diff_types_size = 8

        f_types = mm.Discrete2DFunction(diff_types_size, diff_types_size, tab)
        cross_type.addTabulatedFunction('map_type', f_types) 

        cross_type.addPerParticleParameter("t")

        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            cross_type.addExclusion(idx1, idx2)

        cross_type.setNonbondedMethod(cross_type.CutoffNonPeriodic) #avoid periodic boundary

        for i in range(self.chr_system.getNumParticles()):
            cross_type.addParticle([self.compart_type[i]])

        self.chr_system.addForce(cross_type)

    def add_hardwall(self, epsilon_hw = 1.0, sigma_hw = 0.75, cutoff_hw = 0.84):
        '''
        Add nonbonded hard-wall potential between (chromosomes, speckles, nucleoli) and lamina

        The size of lamina is 1.0 (diameter) and the size of other beads is 0.5, so the sigma is 1.0/2 + 0.5/2 = 0.75 for the 100KB model

        Parameters
        ----------
        epsilon_hw (float, required) : 
            energy units (Default value = 1.0).

        sigma_hw (float, required) : 
            distance units (Default value = 0.75).

        cutoff_hw (float, required) : 
            The cutoff, distance units (Default value = 0.75*1.12).
        '''

        particle_hw_energy = mm.CustomNonbondedForce("step(max(phw1,phw2)-N_hw)*step(N_hw-min(phw1,phw2))*4*epsilon_hw*((sigma_hw/r)^12-(sigma_hw/r)^6-(sigma_hw/cutoff_hw)^12+(sigma_hw/cutoff_hw)^6) * step(cutoff_hw-r)")

        M = self.chr_system.getNumParticles()

        particle_hw_energy.addGlobalParameter("sigma_hw", sigma_hw)
        particle_hw_energy.addGlobalParameter("epsilon_hw", epsilon_hw)
        particle_hw_energy.addGlobalParameter("cutoff_hw", cutoff_hw)
        particle_hw_energy.addGlobalParameter("N_hw", self.N_chr_nuc_spec-0.5)
        particle_hw_energy.setCutoffDistance(cutoff_hw)

        particle_hw_energy.addPerParticleParameter("phw")

        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            particle_hw_energy.addExclusion(idx1, idx2)

        particle_hw_energy.setNonbondedMethod(particle_hw_energy.CutoffNonPeriodic) #avoid periodic boundary

        particle_hw_energy.setForceGroup(14)

        for i in range(self.chr_system.getNumParticles()):
            particle_hw_energy.addParticle([i])

        self.chr_system.addForce(particle_hw_energy)

    def add_LJ_plain(self, epsilon_LJ_plain = 1.0, sigma_LJ_plain = 0.5, cutoff_LJ_plain = 0.56):
        '''
        Add nonbonded plain LJpotential between Nucleoli and Speckles

        Parameters
        ----------
        epsilon_LJ_plain (float, required) :
            energy units (Default value = 1.0).

        sigma_LJ_plain (float, required) :
            distance units (Default value = 0.5).

        cutoff_LJ_plain (float, required) :
            The cutoff, distance units (Default value = 0.5*1.12).
        '''

        LJ_plain = mm.CustomNonbondedForce("4.*epsilon_LJ_plain * ((sigma_LJ_plain/r)^12-(sigma_LJ_plain/r)^6-(sigma_LJ_plain/cutoff_LJ_plain)^12+(sigma_LJ_plain/cutoff_LJ_plain)^6) * step(cutoff_LJ_plain - r)")
        LJ_plain.setCutoffDistance(cutoff_LJ_plain)
        LJ_plain.addGlobalParameter("epsilon_LJ_plain", 1.0)
        LJ_plain.addGlobalParameter("sigma_LJ_plain", 0.5)
        LJ_plain.addGlobalParameter('cutoff_LJ_plain', cutoff_LJ_plain)

        speckle_group = self.bead_groups[5]+self.bead_groups[6]

        LJ_plain.addInteractionGroup(self.bead_groups[4], speckle_group)
        LJ_plain.setNonbondedMethod(LJ_plain.CutoffNonPeriodic) #avoid periodic boundary

        LJ_plain.setForceGroup(15)
        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            LJ_plain.addExclusion(idx1, idx2)

        for i in range(self.chr_system.getNumParticles()):
            LJ_plain.addParticle([])

        self.chr_system.addForce(LJ_plain)

    def add_LJ_nuc(self, epsilon_nuc = 3.0, sigma_nuc = 0.5, cutoff_nuc = 1.5):
        '''
        Add nonbonded rescaled LJpotential between Nucleoli-Nucleoli.

        Parameters
        ----------
        epsilon_nuc (float, required) :
            energy units (Default value = 3.0).

        sigma_nuc (float, required) :
            distance units (Default value = 0.5).

        cutoff_nuc (float, required) : 
            The cutoff, distance units (Default value = 1.5).
        '''
        LJ_nuc = mm.CustomNonbondedForce("4.*epsilon_nuc*((sigNuc/r)^12-(sigNuc/r)^6 - (sigNuc/cutoff_nuc)^12 + (sigNuc/cutoff_nuc)^6) * step(cutoff_nuc - r)")
        LJ_nuc.addGlobalParameter("epsilon_nuc", epsilon_nuc)
        LJ_nuc.addGlobalParameter("sigNuc", sigma_nuc)
        LJ_nuc.addGlobalParameter("cutoff_nuc", cutoff_nuc)
        LJ_nuc.setCutoffDistance(cutoff_nuc)

        LJ_nuc.addInteractionGroup(self.bead_groups[4], self.bead_groups[4])
        LJ_nuc.setNonbondedMethod(LJ_nuc.CutoffNonPeriodic) #avoid periodic boundary

        LJ_nuc.setForceGroup(16)

        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            LJ_nuc.addExclusion(idx1, idx2)

        for i in range(self.chr_system.getNumParticles()):
            LJ_nuc.addParticle([])

        self.chr_system.addForce(LJ_nuc)

    def add_LJ_spec(self, epsilon_dP_dP = 3.0, sigma_spec = 0.5, cutoff_dP_dP = 1.5, epsilon_P = 1., cutoff_P = 0.56):
        """
        Add nonbonded rescaled LJpotential between Speckles-Speckles

        Notice that the system has two different speckles dP (type 6) and P (type 7).

        Parameters
        ----------
        epsilon_dP_dP (float, required) : 
            Epsilon of LJ between (dP,dP), energy units (Default value = 3.0).

        sigma_spec (float, required) : 
            Sigma of LJ between speckles and speckles, distance units (Default value = 0.5).

        cutoff_dP_dP (float, required) : 
            Cutoff of LJ between (dP,dP), distance units (Default value = 1.5).

        epsilon_P (float, required) :
            Epsilon of LJ between (dP,P) and (P,P), energy units (Default value = 1.0).

        cutoff_P (float, required) :
            Cutoff of LJ between (dP,P) and (P,P), distance units (Default value = 0.5*1.12).
        """

        LJ_spec = mm.CustomNonbondedForce("spec_dP_dP + spec_P;"
                "spec_dP_dP = delta(spec1-5)*delta(spec2-5)*4.*epsilon_dP_dP*((sigma_spec/r)^12-(sigma_spec/r)^6-(sigma_spec/cutoff_dP_dP)^12+(sigma_spec/cutoff_dP_dP)^6)*step(cutoff_dP_dP-r);"
                "spec_P = delta(6-max(spec1,spec2))*step(min(spec1,spec2)-4.5)*4.*epsilon_P*((sigma_spec/r)^12-(sigma_spec/r)^6-(sigma_spec/cutoff_P)^12+(sigma_spec/cutoff_P)^6)*step(cutoff_P-r);")
        LJ_spec.setCutoffDistance(cutoff_dP_dP)
        LJ_spec.addGlobalParameter("epsilon_dP_dP", epsilon_dP_dP)
        LJ_spec.addGlobalParameter("sigma_spec", sigma_spec)
        LJ_spec.addGlobalParameter('cutoff_dP_dP', cutoff_dP_dP)
        LJ_spec.addGlobalParameter('epsilon_P', epsilon_P)
        LJ_spec.addGlobalParameter('cutoff_P', cutoff_P)

        LJ_spec.setNonbondedMethod(LJ_spec.CutoffNonPeriodic) #avoid periodic boundary
        LJ_spec.addPerParticleParameter("spec")

        LJ_spec.setForceGroup(17)

        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            LJ_spec.addExclusion(idx1, idx2)

        for i in range(self.chr_system.getNumParticles()):
            LJ_spec.addParticle([self.compart_type[i]])
        self.chr_system.addForce(LJ_spec)
			
    def add_chr_nuc(self, rescalar_file, epsilon_chr_nuc = 1.0, sigma_chr_nuc = 0.5, cutoff_chr_nuc = 1.5):
        '''
        Add nonbonded rescaled LJpotential between chromosomes and nucleoli

        Parameters
        ----------
        rescalar_file (string, required) : 
            The path to the txt file of rescaling factors (SPIN state probabilities).

        epsilon_chr_nuc (float, required) : 
            Epsilon of LJ between nucleolus beads (Default value = 1.0).

        sigma_chr_nuc (float, required) :
            Sigma of LJ between nucleolus beads (Default value = 0.5).

        cutoff_chr_nuc (float, required):
            Cutoff of LJ between nucleolus beads (Default value = 1.5).
        '''

        chr_nuc_energy = ("chr_nuc(min(f1,f2))*step(max(f1,f2)-N_chr)*step(N_chr_nuc-max(f1,f2))*epsilon_chr_nuc*4.*((sigma_chr_nuc/r)^12-(sigma_chr_nuc/r)^6-(sigma_chr_nuc/cutoff_chr_nuc)^12+(sigma_chr_nuc/cutoff_chr_nuc)^6)*step(cutoff_chr_nuc-r);")

        rescalar_mat = np.loadtxt(rescalar_file)
        M = self.chr_system.getNumParticles()
        
        NAD_energy = mm.CustomNonbondedForce(chr_nuc_energy)

        rescalar = mm.Discrete1DFunction(rescalar_mat)
        NAD_energy.addTabulatedFunction('chr_nuc', rescalar)
        NAD_energy.addGlobalParameter("sigma_chr_nuc", sigma_chr_nuc)
        NAD_energy.addGlobalParameter("epsilon_chr_nuc", epsilon_chr_nuc)
        NAD_energy.addGlobalParameter("cutoff_chr_nuc", cutoff_chr_nuc)
        NAD_energy.addGlobalParameter("N_chr", self.N_chr-0.5)
        NAD_energy.addGlobalParameter("N_chr_nuc", self.N_chr_nuc-0.5)
        NAD_energy.setCutoffDistance(cutoff_chr_nuc)

        NAD_energy.addPerParticleParameter("f")

        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            NAD_energy.addExclusion(idx1, idx2)

        NAD_energy.setNonbondedMethod(NAD_energy.CutoffNonPeriodic) #avoid periodic boundary

        NAD_energy.setForceGroup(18)

        for i in range(self.chr_system.getNumParticles()):
            NAD_energy.addParticle([i])

        self.chr_system.addForce(NAD_energy)

    def add_chr_spec(self, chr_spec_param, sigma_tanh_chr_spec = 4., rc_tanh_chr_spec = 0.75, cutoff_chr_spec = 2.0):
        r'''
        Add nonbonded potential using custom values for interactions between chromosomes and speckles

        The parameters are part of the probability of function :math:`g(r_{i,j}) = \frac{1}{2}\left(1 + \text{tanh}\left[\sigma(r_c - r_{i,j})\right] \right)`, where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        chr_spec_param (string, required) : 
            The path to the potential scaling factor txt file.

        sigma_tanh_chr_spec (float, required) :
            Distance units (Default value = 4.0).

        rc_tanh_chr_spec (float, required) : 
            Distance units, where :math:`g(r_{i,j}) = 0.5` (Default value = 0.75).

        cutoff_chr_spec (float, required) : 
            The cutoff (Default value = 2.0)
        '''
        chr_spec_energy = ("chr_dP + chr_P;"
                "chr_dP = tsaseq(min(p1,p2))*delta(max(v1,v2)-5)*(0.5*(1.+tanh(sigma_tanh_chr_spec*(rc_tanh_chr_spec-r)))-0.5*(1.+tanh(sigma_tanh_chr_spec*(rc_tanh_chr_spec-cutoff_chr_spec))))*step(cutoff_chr_spec-r);"
                "chr_P = tsaseq(min(p1,p2))*delta(max(v1,v2)-6)*(0.5*(1.+tanh(sigma_tanh_chr_spec*(rc_tanh_chr_spec-r)))-0.5*(1.+tanh(sigma_tanh_chr_spec*(rc_tanh_chr_spec-cutoff_chr_spec))))*step(cutoff_chr_spec-r);")
        chr_spec_mat = np.loadtxt(chr_spec_param)
        M = self.chr_system.getNumParticles()

        chr_spec_energy = mm.CustomNonbondedForce(chr_spec_energy)

        tsaseq = mm.Discrete1DFunction(chr_spec_mat)
        chr_spec_energy.addTabulatedFunction('tsaseq', tsaseq)
        chr_spec_energy.addGlobalParameter("sigma_tanh_chr_spec", sigma_tanh_chr_spec)
        chr_spec_energy.addGlobalParameter("rc_tanh_chr_spec", rc_tanh_chr_spec)
        chr_spec_energy.addGlobalParameter("cutoff_chr_spec", cutoff_chr_spec)
        chr_spec_energy.setCutoffDistance(cutoff_chr_spec)

        chr_spec_energy.addPerParticleParameter("p")
        chr_spec_energy.addPerParticleParameter("v")

        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            chr_spec_energy.addExclusion(idx1, idx2)

        chr_spec_energy.setNonbondedMethod(chr_spec_energy.CutoffNonPeriodic) #avoid periodic boundary

        chr_spec_energy.setForceGroup(19)

        for i in range(self.chr_system.getNumParticles()):
            chr_spec_energy.addParticle([i,self.compart_type[i]])

        self.chr_system.addForce(chr_spec_energy)

    def add_chr_lam(self, chr_lam_param, sigma_tanh_chr_lam = 4., rc_tanh_chr_lam = 0.75, cutoff_chr_lam = 2.0):
        r'''
        Add nonbonded potential using custom values for interactions between chromosomes and lamina

        The parameters are part of the probability of function :math:`g(r_{i,j}) = \frac{1}{2}\left(1 + \text{tanh}\left[\sigma(r_c - r_{i,j})\right] \right)`, where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        chr_lam_param (string, required) :
            The path to the potential scaling factor txt file.

        sigma_tanh_chr_lam (float, required) :
            Distance units (Default value = 4.0).

        rc_tanh_chr_lam (float, required) : 
            Distance units, where :math:`g(r_{i,j}) = 0.5` (Default value = 0.75).

        cutoff_chr_lam (float, required) :
            The cutoff (Default value = 2.0)
        '''
        chr_lam_energy = ("damid(min(idx_damid1,idx_damid2))*step(max(idx_damid1,idx_damid2)-N_chr_nuc_spec)*(f-0.5*(1.+tanh(sigma_tanh_chr_lam*(rc_tanh_chr_lam-cutoff_chr_lam))))*step(cutoff_chr_lam-r);"
                         "f = 0.5*(1.+tanh(sigma_tanh_chr_lam*(rc_tanh_chr_lam-r)));")

        chr_lam_mat = np.loadtxt(chr_lam_param)
        M = self.chr_system.getNumParticles()

        chr_lam_energy = mm.CustomNonbondedForce(chr_lam_energy)

        damid = mm.Discrete1DFunction(chr_lam_mat)
        chr_lam_energy.addTabulatedFunction('damid', damid)
        chr_lam_energy.addGlobalParameter("sigma_tanh_chr_lam", sigma_tanh_chr_lam)
        chr_lam_energy.addGlobalParameter("rc_tanh_chr_lam", rc_tanh_chr_lam)
        chr_lam_energy.addGlobalParameter("cutoff_chr_lam", cutoff_chr_lam)
        chr_lam_energy.addGlobalParameter("N_chr_nuc_spec", self.N_chr_nuc_spec-0.5)
        chr_lam_energy.setCutoffDistance(cutoff_chr_lam)

        chr_lam_energy.addPerParticleParameter("idx_damid")

        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            chr_lam_energy.addExclusion(idx1, idx2)

        chr_lam_energy.setNonbondedMethod(chr_lam_energy.CutoffNonPeriodic) #avoid periodic boundary

        chr_lam_energy.setForceGroup(20)

        for i in range(self.chr_system.getNumParticles()):
            chr_lam_energy.addParticle([i])

        self.chr_system.addForce(chr_lam_energy)

    def add_inter_potential(self, inter_file, rc_tanh_inter=0.54, cutoff_inter = 2.0):
        r'''
        Add the Inter-chromosome Chromosome potential using custom values for interactions between beads separated by a genomic distance :math:`d`.

        The parameter rc_tanh_inter is part of the probability of crosslink function

        .. math::

            f(r_{i,j})  & = \frac{1}{2}\left(1+\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right) \\
                        & + \frac{1}{2}\left(1-\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right)\left(\frac{r_c}{r_{i,j}}\right)^4,

        where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        inter_file (string, required) :
            The path to the inter-chromosomal potential scaling factor txt file.

        rc_tanh_inter (float, required) :
            :math:`r_c` in the equation, a good estimation of the bond length during the simulation (Default value = 0.54 for the 100KB model).

        cutoff_inter (float, required) :
            Cutoff of inter-chromosomal potential (Default value = 2.0 for the 100KB model).
        '''
        inter_energy = ("mapid(m1,m2)*(tanh1_inter+tanh2_inter*(rc_tanh_inter/r)^4-1/2*(1.+tanh((rc_tanh_inter-cutoff_inter)^(-5)+5*(rc_tanh_inter-cutoff_inter)))-1/2*(1.-tanh((rc_tanh_inter-cutoff_inter)^(-5)+5*(rc_tanh_inter-cutoff_inter)))*(rc_tanh_inter/cutoff_inter)^4)*step(cutoff_inter-r);"
                        "tanh1_inter = 1/2*(1.+tanh((rc_tanh_inter-r)^(-5)+5*(rc_tanh_inter-r)));"
                        "tanh2_inter = 1/2*(1.-tanh((rc_tanh_inter-r)^(-5)+5*(rc_tanh_inter-r)));")

        cross_inter = mm.CustomNonbondedForce(inter_energy)
        
        cross_inter.addGlobalParameter('rc_tanh_inter', rc_tanh_inter)
        cross_inter.addGlobalParameter('cutoff_inter', cutoff_inter)
        cross_inter.setCutoffDistance(cutoff_inter)
        cross_inter.setForceGroup(21)

        tab_inter = np.zeros((23, 23))
        with open(inter_file, "r") as fin_inter:
            for line in fin_inter.readlines():
                line = line.split()
                i = int(line[0]) - 1
                j = int(line[1]) - 1
                tab_inter[i, j] = tab_inter[j, i] = float(line[2])

        tab_inter = list(np.ravel(tab_inter))
        diff_chrom_size = 23

        f_types_inter = mm.Discrete2DFunction(diff_chrom_size, diff_chrom_size, tab_inter)
        cross_inter.addTabulatedFunction('mapid', f_types_inter)

        cross_inter.addPerParticleParameter("m")

        for idx1, idx2 in self.all_bonds: #exclude all nearest neighbor interactions
            cross_inter.addExclusion(idx1, idx2)

        cross_inter.setNonbondedMethod(cross_inter.CutoffNonPeriodic) #avoid periodic boundary

        for i in range(self.chr_system.getNumParticles()):
            cross_inter.addParticle([min(int(np.floor(self.mol_type[i]/2.)),22)])

        self.chr_system.addForce(cross_inter)

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
          
