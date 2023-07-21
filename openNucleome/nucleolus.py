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

class Nucleolus(object):
    '''
    The Nucleolus class stores all the useful potentials related to nucleoli used in the whole nucleus model.
    '''
    def __init__(self, N_bead, N_type, N_chr, N_chr_nuc, bond_list, compart_type, chr_groups, bead_groups, mol_type):
        '''
        Initialize the useful parameters

        Parameters
        ----------
        N_bead (int, required):
            the number of beads in the system.

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
        self.N_type = N_type
        self.N_chr = N_chr
        self.N_chr_nuc = N_chr_nuc
        self.N_hapchrom = 23 # the number of haploid chromsomes, always 23 for human cells
        self.bond_list = bond_list
        self.compart_type = compart_type
        self.chr_groups = chr_groups
        self.bead_groups = bead_groups
        self.mol_type = mol_type

    def add_LJ_plain(self, force_group, epsilon_LJ_plain = 1.0, sigma_LJ_plain = 0.5, cutoff_LJ_plain = 0.56):
        '''
        Add nonbonded plain LJpotential between Nucleoli and Speckles

        Parameters
        ----------
        force_group (int, required):
            labels the index of the current force.

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

        LJ_plain.setForceGroup(force_group)
        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            LJ_plain.addExclusion(int(idx1), int(idx2))

        for i in range(self.N_bead):
            LJ_plain.addParticle([])

        return LJ_plain

    def add_LJ_nuc(self, force_group, epsilon_nuc = 3.0, sigma_nuc = 0.5, cutoff_nuc = 1.5):
        '''
        Add nonbonded LJpotential between Nucleoli-Nucleoli.

        Parameters
        ----------
        force_group (int, required):
            labels the index of the current force.

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

        LJ_nuc.setForceGroup(force_group)

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            LJ_nuc.addExclusion(int(idx1), int(idx2))

        for i in range(self.N_bead):
            LJ_nuc.addParticle([])

        return LJ_nuc

    def add_chr_nuc(self, rescalar_file, force_group, sigma_tanh_chr_nuc = 4.0, rc_tanh_chr_nuc = 0.75, cutoff_chr_nuc = 2.0):
        r'''
        Add nonbonded potential using custom values for interactions between chromosomes and nucleolus.

        The parameters are part of the probability of function :math:`g(r_{i,j}) = \frac{1}{2}\left(1 + \text{tanh}\left[\sigma(r_c - r_{i,j})\right] \right)`, where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        rescalar_file (string, required) :
            The path to the txt file of rescaling factors (SPIN state probabilities).

        force_group (int, required):
            labels the index of the current force.

        sigma_tanh_chr_nuc (float, required) :
            Distance units (Default value = 4.0).

        rc_tanh_chr_nuc (float, required) :
            Distance units, where :math:`g(r_{i,j}) = 0.5` (Default value = 0.75).

        cutoff_chr_nuc (float, required):
            The cutoff (Default value = 2.0).
        '''

        chr_nuc_energy =("chr_nuc(min(f1,f2))*step(max(f1,f2)-N_chr)*step(N_chr_nuc-max(f1,f2))*(f-0.5*(1.+tanh(sigma_tanh_chr_nuc*(rc_tanh_chr_nuc-cutoff_chr_nuc))))*step(cutoff_chr_nuc-r);"
                         "f = 0.5*(1.+tanh(sigma_tanh_chr_nuc*(rc_tanh_chr_nuc-r)));")
        rescalar_mat = np.loadtxt(rescalar_file)

        NAD_energy = mm.CustomNonbondedForce(chr_nuc_energy)

        rescalar = mm.Discrete1DFunction(rescalar_mat)
        NAD_energy.addTabulatedFunction('chr_nuc', rescalar)
        NAD_energy.addGlobalParameter("sigma_tanh_chr_nuc", sigma_tanh_chr_nuc)
        NAD_energy.addGlobalParameter("rc_tanh_chr_nuc", rc_tanh_chr_nuc)
        NAD_energy.addGlobalParameter("cutoff_chr_nuc", cutoff_chr_nuc)
        NAD_energy.addGlobalParameter("N_chr", self.N_chr-0.5)
        NAD_energy.addGlobalParameter("N_chr_nuc", self.N_chr_nuc-0.5)
        NAD_energy.setCutoffDistance(cutoff_chr_nuc)

        NAD_energy.addPerParticleParameter("f")

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            NAD_energy.addExclusion(int(idx1), int(idx2))

        NAD_energy.setNonbondedMethod(NAD_energy.CutoffNonPeriodic) #avoid periodic boundary

        NAD_energy.setForceGroup(force_group)

        for i in range(self.N_bead):
            NAD_energy.addParticle([i])

        return NAD_energy
