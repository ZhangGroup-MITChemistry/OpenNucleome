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

class NewFeature(object):
    '''
    The NewFeature class stores all the useful potentials related to new beads used in the whole nucleus model.
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

    def add_new(self, force_group, epsilon_new = 3.0, sigma_new = 0.5, cutoff_new = 1.5):
        '''
        Add nonbonded LJpotential between new beads and new beads, and new beads and centromeres.

        Parameters
        ----------
        force_group (int, required):
            labels the index of the current force.

        epsilon_new (float, required):
            energy units (Default value = 3.0).

        sigma_new (float, required):
            distance units (Default value = 0.5).

        cutoff_new (float, required):
            The cutoff, distance units (Default value = 1.5).
        '''
        LJ_new = mm.CustomNonbondedForce("4.*epsilon_new*((sig_new/r)^12-(sig_new/r)^6 - (sig_new/cutoff_new)^12 + (sig_new/cutoff_new)^6) * step(cutoff_new - r)")
        LJ_new.addGlobalParameter("epsilon_new", epsilon_new)
        LJ_new.addGlobalParameter("sig_new", sigma_new)
        LJ_new.addGlobalParameter("cutoff_new", cutoff_new)
        LJ_new.setCutoffDistance(cutoff_new)

        LJ_new.addInteractionGroup(self.bead_groups[7]+self.bead_groups[3], self.bead_groups[7])
        LJ_new.setNonbondedMethod(LJ_new.CutoffNonPeriodic) #avoid periodic boundary

        LJ_new.setForceGroup(force_group)

        for idx1, idx2 in self.bond_list: #exclude all nearest neighbor interactions
            LJ_new.addExclusion(int(idx1), int(idx2))

        for i in range(self.N_bead):
            LJ_new.addParticle([])

        return LJ_new
