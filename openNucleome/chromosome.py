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

class Chromosome(object):
    '''
    The Chromosome class stores all the useful potentials amoung chromatin beads used in the whole nucleus model.
    '''
    def __init__(self, N_bead, N_chr, N_type, bond_list, compart_type, chr_groups, bead_groups, mol_type):
        '''
        Initialize the useful parameters

        Parameters
        ----------
        N_bead (int, required):
            the number of beads in the system.

        N_chr (int, required):
            the number of chromosome beads in the system

        N_type (int, required):
            the number of types in the system.

        bond_list (list, required):
            contains the corresponding atom indexes of all bonds.

        compart_type (list, required):
            contains the corresponding compartment types of all beads.

        chr_groups (list, required):
            contains the indexes of all chromatin beads.

        bead_groups (list, required):
            contains the indexes of beads of each type

        mol_type (list, required):
            contains the indexes of all molecules (each chromosome is one molecule, but each neclear body bead is also one molecule).
        '''
        self.N_bead = N_bead
        self.N_chr = N_chr
        self.N_type = N_type
        self.N_hapchrom = 23 # the number of haploid chromsomes, always 23 for human cells
        self.bond_list = bond_list
        self.compart_type = compart_type
        self.chr_groups = chr_groups
        self.bead_groups = bead_groups
        self.mol_type = mol_type

    def add_fene_bond(self, force_group, kfb = 30.0, r0 = 1.5, epsilon_fene = 1.0, sigma_fene = 1.0, cutoff_fene = 1.12):
        '''
        Internal function that inits FENE bond force. The default parameters value is for the 1MB resolution.

        Parameters
        ----------
        force_group (int, required):
            labels the index of the current force.

        r0c, kfb, r0, epsilon_fene, sigma_fene, cutoff_fene (float, required) :
            fene bond coefficients
        '''
        bond = mm.CustomBondForce("- 0.5 * kfb * r0 * r0 * log(1-(r/r0)*(r/r0)) + (4 * epsilon_fene * ((sigma_fene/r)^12 - (sigma_fene/r)^6) + epsilon_fene) * step(cutoff_fene - r)")
        bond.addGlobalParameter("kfb", kfb)
        bond.addGlobalParameter("r0", r0)
        bond.addGlobalParameter("epsilon_fene", epsilon_fene)
        bond.addGlobalParameter("sigma_fene", sigma_fene)
        bond.addGlobalParameter("cutoff_fene", cutoff_fene)
        bond.setForceGroup(force_group)

        bond.setUsesPeriodicBoundaryConditions(False)#avoid periodic boundary
        for idx1, idx2 in self.bond_list:
            bond.addBond(int(idx1), int(idx2))

        return bond

    def add_class2_bond(self, force_group, r0c = 0.5, kc2 = 20.0, kc3 = 20.0, kc4 = 20.0):
        '''
        Internal function that inits Class2 bond force. The default parameters value is for the 100KB resolution. (Recommended)

        Parameters
        ----------
        force_group (int, required):
            labels the index of the current force.

        r0c, kc2, kc3, kc4 (float, required) :
            class2 bond coefficients
        '''
        class2_bond = mm.CustomBondForce("kc2 * (r-r0c) ^ 2 + kc3 * (r-r0c) ^ 3 + kc4 * (r-r0c) ^ 4")
        class2_bond.addGlobalParameter("r0c", r0c)
        class2_bond.addGlobalParameter("kc2", kc2)
        class2_bond.addGlobalParameter("kc3", kc3)
        class2_bond.addGlobalParameter("kc4", kc4)
        class2_bond.setForceGroup(force_group)

        class2_bond.setUsesPeriodicBoundaryConditions(False)#avoid periodic boundary
        for idx1, idx2 in self.bond_list:
            class2_bond.addBond(int(idx1), int(idx2))

        return class2_bond

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
        for bond2 in self.bond_list:
            if bond1 is None:
                bond1 = bond2
                continue
            if bond1[1] == bond2[0] and bond2[1] <= self.N_chr:
                angle.addAngle(bond1[0], bond1[1], bond2[1])
            bond1 = bond2

        angle.setUsesPeriodicBoundaryConditions(False) #avoid periodic boundary
        
        return angle

    def add_softcore(self, force_group, epsilon_sc = 1.0, sigma_sc = 0.5, cutoff_sc = 0.56, Ecut_sc = 4.0):
        r'''
        We plug the tanh potential into the WCA potential to represent the excluded volume effect;

        The softcore potential will only be applied to the chromosomes (type 1, 2, 3, 4).

        When :math:`r_{sc} < r < r_{cutoff}`, use :math:`E_{LJ} = 4\epsilon_{sc}\left(\left(\frac{\sigma_{sc}}{r}\right)^{12} -\left(\frac{\sigma_{sc}}{r}\right)^{6}   \right) + \epsilon_{rc}`; :math:`r < r_{sc}`, use :math:`E = \frac{1}{2}E_{cut}\left(1+\text{tanh}(\frac{2E_{LJ}}{E_{cut}}-1) \right)`, where :math:`r_{sc} \approx 0.485`.

        Parameters
        ----------
        force_group (int, required):
            labels the index of the current force.

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

        soft_core.setForceGroup(force_group)
        soft_core.addPerParticleParameter("c")

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            soft_core.addExclusion(int(idx1), int(idx2))

        #avoid periodic boundary
        soft_core.setNonbondedMethod(soft_core.CutoffNonPeriodic)

        #apply the potential to all beads
        for i in range(self.N_bead):
            soft_core.addParticle([self.compart_type[i]])

        return soft_core

    def add_ideal_potential(self, ideal_file, force_group, rc_tanh_ideal = 0.54, dend_ideal = 1000, cutoff_ideal = 2.0):
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

        force_group (int, required):
            labels the index of the current force.

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
                        "d = abs(id2-id1)")

        ideal = mm.CustomNonbondedForce(energy_ideal)
        ideal_alpha = np.loadtxt(ideal_file)[:,1]

        tab_ideal_list = mm.Discrete1DFunction(ideal_alpha)
        ideal.addTabulatedFunction('ideal_list', tab_ideal_list)
        ideal.addGlobalParameter('dend_ideal', dend_ideal)
        ideal.addGlobalParameter('rc_tanh_ideal', rc_tanh_ideal)
        ideal.addGlobalParameter('cutoff_ideal', cutoff_ideal)
        ideal.setForceGroup(force_group)

        ideal.setCutoffDistance(cutoff_ideal)

        ideal.addPerParticleParameter("id")

        #apply within each chromosome
        for i in range(len(self.chr_groups)):
            ideal.addInteractionGroup(self.chr_groups[i], self.chr_groups[i])

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            ideal.addExclusion(int(idx1), int(idx2))

        #avoid periodic boundary
        ideal.setNonbondedMethod(ideal.CutoffNonPeriodic)

        #apply the potential to all beads
        for i in range(self.N_bead):
            ideal.addParticle([i])

        return ideal

    def add_type_type_potential(self, types_file, force_group, rc_tanh_type = 0.54, cutoff_type = 2.0):
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

        force_group (int, required):
            labels the index of the current force.        

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

        cross_type.setForceGroup(force_group)
        #generate the above shown 8*8 matrix, 8 is the number of types
        tab = np.zeros((self.N_type, self.N_type))
        with open(types_file, "r") as fin:
            for line in fin.readlines():
                line = line.split()
                i = int(line[0]) - 1
                j = int(line[1]) - 1
                tab[i, j] = tab[j, i] = float(line[2])

        tab = list(np.ravel(tab))

        diff_types_size = self.N_type

        f_types = mm.Discrete2DFunction(diff_types_size, diff_types_size, tab)
        cross_type.addTabulatedFunction('map_type', f_types)

        cross_type.addPerParticleParameter("t")

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            cross_type.addExclusion(int(idx1), int(idx2))

        cross_type.setNonbondedMethod(cross_type.CutoffNonPeriodic) #avoid periodic boundary

        for i in range(self.N_bead):
            cross_type.addParticle([self.compart_type[i]])

        return cross_type

    def add_inter_potential(self, inter_file, force_group, rc_tanh_inter=0.54, cutoff_inter = 2.0):
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

        force_group (int, required):
            labels the index of the current force.

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
        cross_inter.setForceGroup(force_group)

        tab_inter = np.zeros((self.N_hapchrom, self.N_hapchrom))
        with open(inter_file, "r") as fin_inter:
            for line in fin_inter.readlines():
                line = line.split()
                i = int(line[0]) - 1
                j = int(line[1]) - 1
                tab_inter[i, j] = tab_inter[j, i] = float(line[2])

        tab_inter = list(np.ravel(tab_inter))
        diff_chrom_size = self.N_hapchrom

        f_types_inter = mm.Discrete2DFunction(diff_chrom_size, diff_chrom_size, tab_inter)
        cross_inter.addTabulatedFunction('mapid', f_types_inter)

        cross_inter.addPerParticleParameter("m")

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            cross_inter.addExclusion(int(idx1), int(idx2))

        cross_inter.setNonbondedMethod(cross_inter.CutoffNonPeriodic) #avoid periodic boundary

        for i in range(self.N_bead):
            cross_inter.addParticle([min(int(np.floor(self.mol_type[i]/2.)),22)])

        return cross_inter
