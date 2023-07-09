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

class Speckle(object):
    '''
    The Speckle class stores all the useful potentials related to speckles used in the whole nucleus model.
    '''
    def __init__(self, N_bead, N_type, bond_list, compart_type, chr_groups, bead_groups, mol_type):
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
        self.N_hapchrom = 23 # the number of haploid chromsomes, always 23 for human cells
        self.bond_list = bond_list
        self.compart_type = compart_type
        self.chr_groups = chr_groups
        self.bead_groups = bead_groups
        self.mol_type = mol_type

    def add_LJ_spec(self, force_group, epsilon_dP_dP = 3.0, sigma_spec = 0.5, cutoff_dP_dP = 1.5, epsilon_P = 1., cutoff_P = 0.56):
        """
        Add nonbonded rescaled LJpotential between Speckles-Speckles

        Notice that the system has two different speckles dP (type 6) and P (type 7).

        Parameters
        ----------
        force_group (int, required):
            labels the index of the current force.

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

        LJ_spec.setForceGroup(force_group)

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            LJ_spec.addExclusion(int(idx1), int(idx2))

        for i in range(self.N_bead):
            LJ_spec.addParticle([self.compart_type[i]])
        
        return LJ_spec

    def add_chr_spec(self, chr_spec_param, force_group, sigma_tanh_chr_spec = 4., rc_tanh_chr_spec = 0.75, cutoff_chr_spec = 2.0):
        r'''
        Add nonbonded potential using custom values for interactions between chromosomes and speckles

        The parameters are part of the probability of function :math:`g(r_{i,j}) = \frac{1}{2}\left(1 + \text{tanh}\left[\sigma(r_c - r_{i,j})\right] \right)`, where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        chr_spec_param (string, required) :
            The path to the potential scaling factor txt file.

        force_group (int, required):
            labels the index of the current force.

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

        chr_spec_energy = mm.CustomNonbondedForce(chr_spec_energy)

        tsaseq = mm.Discrete1DFunction(chr_spec_mat)
        chr_spec_energy.addTabulatedFunction('tsaseq', tsaseq)
        chr_spec_energy.addGlobalParameter("sigma_tanh_chr_spec", sigma_tanh_chr_spec)
        chr_spec_energy.addGlobalParameter("rc_tanh_chr_spec", rc_tanh_chr_spec)
        chr_spec_energy.addGlobalParameter("cutoff_chr_spec", cutoff_chr_spec)
        chr_spec_energy.setCutoffDistance(cutoff_chr_spec)

        chr_spec_energy.addPerParticleParameter("p")
        chr_spec_energy.addPerParticleParameter("v")

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            chr_spec_energy.addExclusion(int(idx1), int(idx2))

        chr_spec_energy.setNonbondedMethod(chr_spec_energy.CutoffNonPeriodic) #avoid periodic boundary

        chr_spec_energy.setForceGroup(force_group)

        for i in range(self.N_bead):
            chr_spec_energy.addParticle([i,self.compart_type[i]])

        return chr_spec_energy
