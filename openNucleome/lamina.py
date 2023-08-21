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

class Lamina(object):
    '''
    The Lamina class stores all the useful potentials related to lamina used in the whole nucleus model.
    '''
    def __init__(self, N_bead, N_type, N_chr_nuc_spec, bond_list, compart_type, chr_groups, bead_groups, mol_type):
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
        self.N_chr_nuc_spec = N_chr_nuc_spec
        self.N_hapchrom = 23 # the number of haploid chromsomes, always 23 for human cells
        self.bond_list = bond_list
        self.compart_type = compart_type
        self.chr_groups = chr_groups
        self.bead_groups = bead_groups
        self.mol_type = mol_type

    def add_hardwall(self, force_group, epsilon_hw = 1.0, sigma_hw = 0.75, cutoff_hw = 0.84):
        '''
        Add nonbonded hard-wall potential between (chromosomes, speckles, nucleoli) and lamina

        The size of lamina is 1.0 (diameter) and the size of other beads is 0.5, so the sigma is 1.0/2 + 0.5/2 = 0.75 for the 100KB model

        Parameters
        ----------
        force_group (int, required):
            labels the index of the current force.

        epsilon_hw (float, required) :
            energy units (Default value = 1.0).

        sigma_hw (float, required) :
            distance units (Default value = 0.75).

        cutoff_hw (float, required) :
            The cutoff, distance units (Default value = 0.75*1.12).
        '''

        particle_hw_energy = mm.CustomNonbondedForce("step(max(phw1,phw2)-N_hw)*step(N_hw-min(phw1,phw2))*4*epsilon_hw*((sigma_hw/r)^12-(sigma_hw/r)^6-(sigma_hw/cutoff_hw)^12+(sigma_hw/cutoff_hw)^6) * step(cutoff_hw-r)")

        particle_hw_energy.addGlobalParameter("sigma_hw", sigma_hw)
        particle_hw_energy.addGlobalParameter("epsilon_hw", epsilon_hw)
        particle_hw_energy.addGlobalParameter("cutoff_hw", cutoff_hw)
        particle_hw_energy.addGlobalParameter("N_hw", self.N_chr_nuc_spec-0.5)
        particle_hw_energy.setCutoffDistance(cutoff_hw)

        particle_hw_energy.addPerParticleParameter("phw")

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            particle_hw_energy.addExclusion(int(idx1), int(idx2))

        particle_hw_energy.setNonbondedMethod(particle_hw_energy.CutoffNonPeriodic) #avoid periodic boundary

        particle_hw_energy.setForceGroup(force_group)

        for i in range(self.N_bead):
            particle_hw_energy.addParticle([i])

        return particle_hw_energy

    def add_chr_lam(self, chr_lam_param, force_group, sigma_tanh_chr_lam = 4., rc_tanh_chr_lam = 0.75, cutoff_chr_lam = 2.0):
        r'''
        Add nonbonded potential using custom values for interactions between chromosomes and lamina

        The parameters are part of the probability of function :math:`g(r_{i,j}) = \frac{1}{2}\left(1 + \text{tanh}\left[\sigma(r_c - r_{i,j})\right] \right)`, where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------

        chr_lam_param (string, required) :
            The path to the potential scaling factor txt file.

        force_group (int, required):
            labels the index of the current force.

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

        chr_lam_energy = mm.CustomNonbondedForce(chr_lam_energy)

        damid = mm.Discrete1DFunction(chr_lam_mat)
        chr_lam_energy.addTabulatedFunction('damid', damid)
        chr_lam_energy.addGlobalParameter("sigma_tanh_chr_lam", sigma_tanh_chr_lam)
        chr_lam_energy.addGlobalParameter("rc_tanh_chr_lam", rc_tanh_chr_lam)
        chr_lam_energy.addGlobalParameter("cutoff_chr_lam", cutoff_chr_lam)
        chr_lam_energy.addGlobalParameter("N_chr_nuc_spec", self.N_chr_nuc_spec-0.5)
        chr_lam_energy.setCutoffDistance(cutoff_chr_lam)

        chr_lam_energy.addPerParticleParameter("idx_damid")

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            chr_lam_energy.addExclusion(int(idx1), int(idx2))

        chr_lam_energy.setNonbondedMethod(chr_lam_energy.CutoffNonPeriodic) #avoid periodic boundary

        chr_lam_energy.setForceGroup(force_group)

        for i in range(self.N_bead):
            chr_lam_energy.addParticle([i])

        return chr_lam_energy

    def add_lam_lam(self, force_group, epsilon_lam = 1., sigma_lam = 0.5, cutoff_lam = 0.5*(2.**(1./6.))):
        '''
        Add nonbonded plain LJpotential between Lamina and Lamina

        Parameters
        ----------

        force_group (int, required) :
            labels the index of the current force.

        epsilon_lam (float, required) :
            energy units (Default value = 1.0).

        sigma_lam (float, required) :
            Distance units (Default value = 0.5).

        cutoff_lam (float, required) :
            The cutoff (Default value = 0.5*1.12)
        '''
        #load the bead-specific rescaling factors

        LJ_lam = mm.CustomNonbondedForce("step(min(h1,h2)-N_lam)*4*epsilon_lam*((sigma_lam/r)^12-(sigma_lam/r)^6-(sigma_lam/cutoff_lam)^12+(sigma_lam/cutoff_lam)^6) * step(cutoff_lam-r)")

        LJ_lam.addGlobalParameter("epsilon_lam", epsilon_lam)
        LJ_lam.addGlobalParameter("sigma_lam", sigma_lam)
        LJ_lam.addGlobalParameter("cutoff_lam", cutoff_lam)
        LJ_lam.addGlobalParameter("N_lam", self.N_chr_nuc_spec-0.5)
        LJ_lam.setCutoffDistance(cutoff_lam)

        LJ_lam.addPerParticleParameter("h")

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            LJ_lam.addExclusion(int(idx1), int(idx2))

        LJ_lam.setNonbondedMethod(LJ_lam.CutoffNonPeriodic) #avoid periodic boundary

        LJ_lam.setForceGroup(force_group)

        for i in range(self.N_bead):
            LJ_lam.addParticle([i])

        return LJ_lam

    def add_squeeze_nucleus(self, force_group, k = 1., R = 13.):
        '''
        Add external force to squeeze the cell nucleus (Default: spring force)

        Parameters
        ----------

        force_group (int, required) :
            labels the index of the current force.

        k (float, required) :
            energy units (Default value = 1.0).

        R (float, required) :
            perfect sphere nucleus radius (Default value = 13.0)
        '''
        #load the bead-specific rescaling factors

        squeeze = mm.CustomExternalForce("step(x-N_lam)*k*(2*step(z)-1)*z*abs(z)/R")

        squeeze.addGlobalParameter("k", k)
        squeeze.addGlobalParameter("R", R)
        squeeze.addGlobalParameter("N_lam", self.N_chr_nuc_spec-0.5)

        squeeze.addPerParticleParameter("x")

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            squeeze.addExclusion(int(idx1), int(idx2))

        squeeze.setForceGroup(force_group)

        for i in range(self.N_bead):
            squeeze.addParticle(i, [i])

        return squeeze
