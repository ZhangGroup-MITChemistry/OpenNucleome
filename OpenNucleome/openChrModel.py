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

class openChrModel:
    '''
    The openChrModel class performes the whole genome dynamics based on the compartment annotations sequence of chromosomes.
    
    The simulations can be performed using customed values for the type-to-type, intra-chromosomal, and inter-chromosomal parameters.
    '''
    def __init__(self, temperature = 1.0, gamma = 0.1, timestep = 0.01, mass_scale=1.0):
        r'''
        Initialize the whole genome model.

        Parameters
        ----------
        temperature (float, required) :
            Temperature in reduced units. (Default value = 1.0).

        gamma (float, required) :
            Friction/Damping constant in units of reciprocal time (:math:`1/\tau`) (Default value = 0.1).

        timestep (float, required):
            Simulation time step in units of :math:`\tau` (Default value = 0.01).

        mass_scale (float, required):
            Mass scale used in units of :math:`\mu` (Default value = 1.0).

        '''
        self.temperature = temperature / 0.008314 # temperature in reduced unit
        self.gamma = gamma # fiction coefficient (1/time)
        self.timestep = timestep  # timestep in reduced unit (time)
        self.kB = u.BOLTZMANN_CONSTANT_kB * u.AVOGADRO_CONSTANT_NA
        self.kT = self.kB * self.temperature
        self.chrMass = 10.0 * u.amu * mass_scale

    def createSystem(self, PDBFile):
        '''
        Create the system with initial configurations.

        Parameters
        ----------

        PDBFile (string, required):
            The path to the chromosome PDB file, including lamina, speckles, and nucleoli.
        '''
        chrPDBpath = PDBFile
        self.chrPDB = mmapp.PDBFile(chrPDBpath)

        self.chrSystem = mm.System()
        self.chrPositions = []

        #record the position information of each bead in self.chrPositions
        for i in range(len(self.chrPDB.positions)):
        #unit conversion, the default imported PDB position information is angstrom, which will be 10 times smaller for the normal nm calculation.
            self.chrPDB.positions[i] *= 10 
            self.chrPositions.append(self.chrPDB.positions[i])

        self.Ntotal = len(self.chrPositions) #total number of chr beads
        self.compartType = [] #compartment type information
        self.molType = [] #mol type information
        self.beadGroups = [[] for _ in range(8)] # group the indexes of atoms by their types; if you do not want to add some nuclear landmarks, then use the different value. compartment A: 1, compartment B: 2, regions around the centromeres: 3, centromeres: 4, nucleoli: 5, speckles: 6, lamina: 7
        self.chrGroups = [[] for _ in range(46)] # represent the 46 chromosomes of human
        prevResID = -1
        self.chrResidues = [] #(start, end) index of the atoms for each residue
        start = 0

        for a in self.chrPDB.topology.atoms():
            #record compartment type info
            self.compartType.append(int(a.name) - 1)
            self.molType.append(int(a.residue.index))

            #add atom into the system with mass = chrMass
            m = 1. if self.compartType[-1] != 7 else 0. # set mass of Lamina beads to zero to freeze them
            self.chrSystem.addParticle(self.chrMass * m)

            #find the start and end index in each residue, self.chrResidue will be (start, end, residue) for each residue after self.constructTopology
            if a.residue.index != prevResID and prevResID != -1:
                self.chrResidues.append([start, a.index - 1])
                start = a.index
            prevResID = a.residue.index

            #group the indexes of atoms by their types
            self.beadGroups[int(a.name) - 1].append(a.index)
            if a.residue.index < 46:
                self.chrGroups[int(a.residue.index)].append(a.index)
        
        self.N_chr = 0
        self.N_chr_nuc = 0
        self.N_chr_nuc_spec = 0
        for i in range(7):
            if i<=3:
                self.N_chr += len(self.beadGroups[i])
                self.N_chr_nuc += len(self.beadGroups[i])
                self.N_chr_nuc_spec += len(self.beadGroups[i])
            elif i<=4:
                self.N_chr_nuc += len(self.beadGroups[i])
                self.N_chr_nuc_spec += len(self.beadGroups[i])
            else:
                self.N_chr_nuc_spec += len(self.beadGroups[i])

        self.chrResidues.append([start, a.index])

        self.generateElement()

        self.constructTopology()

    def generateElement(self):
        '''
        Generate elements for each coarse-grained bead.

        This function is called automatically when creating a system.
        '''
        name2Element = ['ASP', 'GLU', 'HIS', 'LYS', 'ARG', 'ASN', 'GLN', 'PRO']
        self.Elements = []
        for i in range(8):
            m = 1. if i!= 7 else 0.
#           self.Elements.append(mmapp.element.Element.getBySymbol(name2Element[i]))
            self.Elements.append(mmapp.Element(1000+i, name2Element[i], name2Element[i], self.chrMass * m))

    def constructTopology(self):
        '''
        Construct the topology for the system.

        This function is called automatically when creating a system.
        '''
        # Construct topology, add chain, residues, atoms, bonds
        self.chrTopo = mmapp.topology.Topology()
        chrChain = self.chrTopo.addChain('0')
        for i in range(len(self.chrResidues)):
            currResIdx = self.chrTopo.addResidue(str(i), chrChain)
            self.chrResidues[i].append(currResIdx)
            for j in range(self.chrResidues[i][0], self.chrResidues[i][1] + 1):
                self.chrTopo.addAtom(str(j), self.Elements[self.compartType[j]], currResIdx)
                if j != self.chrResidues[i][0]:
                    self.chrTopo.addBond(j - 1, j) #added atom index, instead of atom item here

        #store all bonds in Allbonds
        self.Allbonds = []
        for b in self.chrTopo.bonds():
            self.Allbonds.append((b.atom1, b.atom2))

    # The diameter of 100KB bead is 0.5; 1.0 for 1MB bead. For example, I set sigmaFene as 1.0, so currently, these parameters are for 1MB bead.
    def addFeneBond(self, kfb = 30.0, r0 = 1.5, epsFene = 1.0, sigmaFene = 1.0, cutFene = 1.12):
        '''
        Internal function that inits FENE bond force. The default parameters value is for the 1MB resolution.
        '''
        bond = mm.CustomBondForce("- 0.5 * kfb * r0 * r0 * log(1-(r/r0)*(r/r0)) + (4 * epsFene * ((sigmaFene/r)^12 - (sigmaFene/r)^6) + epsFene) * step(cutFene - r)")
        bond.addGlobalParameter("kfb", kfb) 
        bond.addGlobalParameter("r0", r0) 
        bond.addGlobalParameter("epsFene", epsFene) 
        bond.addGlobalParameter("sigmaFene", sigmaFene) 
        bond.addGlobalParameter("cutFene", cutFene) 
        bond.setForceGroup(10)

        bond.setUsesPeriodicBoundaryConditions(False)#avoid periodic boundary
        for b in self.chrTopo.bonds():
            bond.addBond(b.atom1, b.atom2)
        self.chrSystem.addForce(bond)

    def addClass2Bond(self, r0c = 0.5, kc2 = 20.0, kc3 = 20.0, kc4 = 20.0):
        '''
        Internal function that inits Class2 bond force. The default parameters value is for the 100KB resolution.
        '''
        Hbond = mm.CustomBondForce("kc2 * (r-r0c) ^ 2 + kc3 * (r-r0c) ^ 3 + kc4 * (r-r0c) ^ 4")
        Hbond.addGlobalParameter("r0c", r0c) 
        Hbond.addGlobalParameter("kc2", kc2) 
        Hbond.addGlobalParameter("kc3", kc3) 
        Hbond.addGlobalParameter("kc4", kc4) 
        Hbond.setForceGroup(10)


        Hbond.setUsesPeriodicBoundaryConditions(False)#avoid periodic boundary
        for b in self.chrTopo.bonds():
            Hbond.addBond(b.atom1, b.atom2)
        self.chrSystem.addForce(Hbond)

    def addAngleForce(self, k_a = 2.0):
        '''
        Add an angular potential between bonds connecting beads.

        Parameters
        ----------
        k_a (float, required) :
            angle coefficients, default = 2.0
        '''
        angle = mm.CustomAngleForce("k_a*(1 + cos(theta))")
        angle.addGlobalParameter("k_a", k_a) #k_a in reduced unit

        #Apply to all angle triplets
        bond1 = None
        for bond2 in self.chrTopo.bonds():
            if bond1 is None:
                bond1 = bond2
                continue
            if bond1.atom2 == bond2.atom1:
                angle.addAngle(bond1.atom1, bond1.atom2, bond2.atom2)
            bond1 = bond2

        angle.setUsesPeriodicBoundaryConditions(False) #avoid periodic boundary
        self.chrSystem.addForce(angle)

    def addSoftCore(self, epsSC = 1.0, sigmaSC = 0.5, cutOffSC = 0.56, EcutSC = 4.0):
        r'''
        We plug the tanh potential into the WCA potential to represent the excluded volume effect;

        The softcore potential will only be applied to the chromosomes (type 1, 2, 3, 4).

        Here, we use step function to select the atoms.

        When :math:`r < r_{sc}`, use tanh potential; :math:`r_{sc} < r < r_{cutoff}`, use WCA potential, where :math:`r_{sc} \approx 0.485`.

        Parameters
        ----------
        epsSC (float, required) :
            energy units

        sigmaSC (float, required) :
            distance units, which represents the diameter of a 100KB-resolution bead.

        cutOffSC (float, required) :
            distance units, the cutoff.

        EcutSC (float, required) :
            energy units, the energy cost for the chain passing (Default value = 4.0).
        '''
        EcutSC = EcutSC * epsSC
        r_0SC = sigmaSC*(((0.5*EcutSC)/(4.0*epsSC) - 0.25 +((0.5)**(2.0)))**(1.0/2.0) +0.5)**(-1.0/6.0)
 
        repul_energy = ("LJ*step(r-r_0SC)*step(cutOffSC-r)*step(3.5-max(c1,c2))"
                              " + step(3.5-max(c1,c2))*step(r_0SC-r)*0.5*EcutSC*(1.0+tanh((2.0*LJ/EcutSC)-1.0));"
                              "LJ = 4.0 * epsSC * ((sigmaSC/r)^12 - (sigmaSC/r)^6) + epsSC")

        softCore = mm.CustomNonbondedForce(repul_energy)
        softCore.addGlobalParameter('epsSC', epsSC)
        softCore.addGlobalParameter('sigmaSC', sigmaSC)
        softCore.addGlobalParameter('EcutSC', EcutSC)
        softCore.addGlobalParameter('r_0SC', r_0SC)
        softCore.addGlobalParameter('cutOffSC', cutOffSC)
        softCore.setCutoffDistance(cutOffSC)

        softCore.setForceGroup(11)
        softCore.addPerParticleParameter("c")

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            softCore.addExclusion(idx1, idx2)

        #avoid periodic boundary
        softCore.setNonbondedMethod(softCore.CutoffNonPeriodic)

        #apply the potential to all beads
        for i in range(self.chrSystem.getNumParticles()):
            softCore.addParticle([self.compartType[i]])

        self.chrSystem.addForce(softCore)

    def addIdealPotential(self, ICFile, rctanhIC = 0.54, dendIC = 1000, cutOffIC = 2.0):
        r'''
        Add the Intra-chromosomal ideal potential using custom values for interactions between beads separated by a genomic distance :math:`d`. 

        The parameter rctanhIC is part of the probability of crosslink function :math:`f(r_{i,j})=\frac{1}{2}\left(1+\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right)+\frac{1}{2}\left(1-\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right)\left(\frac{r_c}{r_{i,j}}\right)^4`, 

        where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        ICFile (string, required) :
            The path to the ideal potential scaling factor txt file.

        dendIC (int, required) : 
            Cutoff of genomic separation (Default value = 1000 for the 100KB model).

        rctanhIC (float, required) : 
            :math:`r_c` in the equation, a good estimation of the bond length during the simulation (Default value = 0.54 for the 100KB model).

        cutOffIC (float, required) : 
            Cutoff of Ideal potential (Default value = 2.0 for the 100KB model).
        '''
        energyIC = ("IClist(d)*step(dendIC-d)*step(cutOffIC-r)*(tanh1_IC+tanh2_IC*(rctanhIC/r)^4-1/2*(1.+tanh((rctanhIC-cutOffIC)^(-5)+5*(rctanhIC-cutOffIC)))-1/2*(1.-tanh((rctanhIC-cutOffIC)^(-5)+5*(rctanhIC-cutOffIC)))*(rctanhIC/cutOffIC)^4);"
                        "tanh1_IC = 1/2*(1.+tanh((rctanhIC-r)^(-5)+5*(rctanhIC-r)));"
                        "tanh2_IC = 1/2*(1.-tanh((rctanhIC-r)^(-5)+5*(rctanhIC-r)));"
                        "d = abs(idx2-idx1)")

        IC = mm.CustomNonbondedForce(energyIC)
        IClist_listfromfile = np.loadtxt(ICFile)[:,1] 

        tabIClist = mm.Discrete1DFunction(IClist_listfromfile)
        IC.addTabulatedFunction('IClist', tabIClist) 
        IC.addGlobalParameter('dendIC', dendIC)
        IC.addGlobalParameter('rctanhIC', rctanhIC) 
        IC.addGlobalParameter('cutOffIC', cutOffIC) 
        IC.setForceGroup(12)

        IC.setCutoffDistance(cutOffIC)

        IC.addPerParticleParameter("idx")

        #apply within each chromosome
        for i in range(len(self.chrGroups)):
            IC.addInteractionGroup(self.chrGroups[i], self.chrGroups[i])

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            IC.addExclusion(idx1, idx2)

        #avoid periodic boundary
        IC.setNonbondedMethod(IC.CutoffNonPeriodic)

        #apply the potential to all beads
        for i in range(self.chrSystem.getNumParticles()):
            IC.addParticle([i])

        self.chrSystem.addForce(IC)

    def addType2TypePotential(self, TypesTable, rctanhLP = 0.54, crossLPcutOff = 2.0):
        r'''
        Add the type-to-type potential using custom values for interactions between the chromatin types. 

        The parameter rctanhLP is part of the probability of crosslink function :math:`f(r_{i,j})=\frac{1}{2}\left(1+\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right)+\frac{1}{2}\left(1-\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right)\left(\frac{r_c}{r_{i,j}}\right)^4`,

        where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        TypesTable (string, required) :
            The path to the type-to-type potential scaling factor txt file.

        rctanhLP (float, required) : 
            :math:`r_c` in the equation, a good estimation of the bond length during the simulation (Default value = 0.54 for the 100KB model).
        
        crossLPcutOff (float, required) :
            Cutoff of Type-Type potential (Default value = 2.0 for the 100KB model).
        '''

        energy = ("mapType(t1,t2)*step(crossLPcutOff-r)*(tanh1_LP+tanh2_LP*(rctanhLP/r)^4-1/2*(1.+tanh((rctanhLP-crossLPcutOff)^(-5)+5*(rctanhLP-crossLPcutOff)))-1/2*(1.-tanh((rctanhLP-crossLPcutOff)^(-5)+5*(rctanhLP-crossLPcutOff)))*(rctanhLP/crossLPcutOff)^4);"
                        "tanh1_LP = 1/2*(1.+tanh((rctanhLP-r)^(-5)+5*(rctanhLP-r)));"
                        "tanh2_LP = 1/2*(1.-tanh((rctanhLP-r)^(-5)+5*(rctanhLP-r)));")

        crossLP = mm.CustomNonbondedForce(energy)

        crossLP.addGlobalParameter('rctanhLP', rctanhLP)
        crossLP.addGlobalParameter('crossLPcutOff', crossLPcutOff)
        crossLP.setCutoffDistance(crossLPcutOff)
        
        crossLP.setForceGroup(13)
        #generate the above shown 8*8 matrix, 8 is the number of types
        tab = np.zeros((8, 8))
        with open(TypesTable, "r") as fin:
            for line in fin.readlines():
                line = line.split()
                i = int(line[1]) - 1
                j = int(line[2]) - 1
                tab[i, j] = tab[j, i] = float(line[5])

        tab = list(np.ravel(tab))
        diff_types_size = 8

        fTypes = mm.Discrete2DFunction(diff_types_size, diff_types_size, tab)
        crossLP.addTabulatedFunction('mapType', fTypes) 

        crossLP.addPerParticleParameter("t")

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            crossLP.addExclusion(idx1, idx2)

        crossLP.setNonbondedMethod(crossLP.CutoffNonPeriodic) #avoid periodic boundary

        for i in range(self.chrSystem.getNumParticles()):
            crossLP.addParticle([self.compartType[i]])

        self.chrSystem.addForce(crossLP)   

    def addparticle_hw(self, eps_p_hw = 1.0, sigma_p_hw = 0.75, cutoff_p_hw = 0.84):
        '''
        Add nonbonded hard-wall potential between (chromosomes, speckles, nucleoli) and lamina

        The size of lamina is 1.0 (diameter) and the size of other beads is 0.5, so the sigma is 1.0/2 + 0.5/2 = 0.75 for the 100KB model

        Parameters
        ----------
        eps_p_hw (float, required) : 
            energy units (Default value = 1.0).

        sigma_p_hw (float, required) : 
            distance units (Default value = 0.75).

        cutoff_p_hw (float, required) : 
            The cutoff, distance units (Default value = 0.75*1.12).
        '''

        particle_hw_energy = mm.CustomNonbondedForce("step(max(phw1,phw2)-N_hw)*step(N_hw-min(phw1,phw2))*4*eps_p_hw*((sigma_p_hw/r)^12-(sigma_p_hw/r)^6-(sigma_p_hw/cutoff_p_hw)^12+(sigma_p_hw/cutoff_p_hw)^6) * step(cutoff_p_hw-r)")

        M = self.chrSystem.getNumParticles()

        particle_hw_energy.addGlobalParameter("sigma_p_hw", sigma_p_hw)
        particle_hw_energy.addGlobalParameter("eps_p_hw", eps_p_hw)
        particle_hw_energy.addGlobalParameter("cutoff_p_hw", cutoff_p_hw)
        particle_hw_energy.addGlobalParameter("N_hw", self.N_chr_nuc_spec-0.5)
        particle_hw_energy.setCutoffDistance(cutoff_p_hw)

        particle_hw_energy.addPerParticleParameter("phw")

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            particle_hw_energy.addExclusion(idx1, idx2)

        particle_hw_energy.setNonbondedMethod(particle_hw_energy.CutoffNonPeriodic) #avoid periodic boundary

        particle_hw_energy.setForceGroup(14)

        for i in range(self.chrSystem.getNumParticles()):
            particle_hw_energy.addParticle([i])

        self.chrSystem.addForce(particle_hw_energy)

    def addLJplain(self):
        '''
        Add nonbonded plain LJpotential between Nucleoli and Speckles
        '''
        CutOffLJplain = 0.5 * 2. ** (1 / 6) 

        LJplain = mm.CustomNonbondedForce("4.*epsLJplain * ((sigLJplain/r)^12-(sigLJplain/r)^6-(sigLJplain/CutOffLJplain)^12+(sigLJplain/CutOffLJplain)^6) * step(CutOffLJplain - r)")
        LJplain.setCutoffDistance(CutOffLJplain)
        LJplain.addGlobalParameter("epsLJplain", 1.0)
        LJplain.addGlobalParameter("sigLJplain", 0.5)
        LJplain.addGlobalParameter('CutOffLJplain', CutOffLJplain)

        speckle_group = self.beadGroups[5]+self.beadGroups[6]

        LJplain.addInteractionGroup(self.beadGroups[4], speckle_group)
        LJplain.setNonbondedMethod(LJplain.CutoffNonPeriodic) #avoid periodic boundary

        LJplain.setForceGroup(15)
        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            LJplain.addExclusion(idx1, idx2)

        for i in range(self.chrSystem.getNumParticles()):
            LJplain.addParticle([])

        self.chrSystem.addForce(LJplain)

    def addLJNuc(self, epsNuc = 3.0, sigmaNuc = 0.5, cutOffNuc = 1.5):
        '''
        Add nonbonded rescaled LJpotential between Nucleoli-Nucleoli.

        Parameters
        ----------
        epsNuc (float, required) :
            energy units (Default value = 3.0).

        sigmaNuc (float, required) :
            distance units (Default value = 0.5).

        cutOffNuc (float, required) : 
            The cutoff, distance units (Default value = 1.5).
        '''
        LJNuc = mm.CustomNonbondedForce("4.*epsNuc*((sigNuc/r)^12-(sigNuc/r)^6 - (sigNuc/cutOffNuc)^12 + (sigNuc/cutOffNuc)^6) * step(cutOffNuc - r)")
        LJNuc.addGlobalParameter("epsNuc", epsNuc)
        LJNuc.addGlobalParameter("sigNuc", sigmaNuc)
        LJNuc.addGlobalParameter("cutOffNuc", cutOffNuc)
        LJNuc.setCutoffDistance(cutOffNuc)

        LJNuc.addInteractionGroup(self.beadGroups[4], self.beadGroups[4])
        LJNuc.setNonbondedMethod(LJNuc.CutoffNonPeriodic) #avoid periodic boundary

        LJNuc.setForceGroup(16)

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            LJNuc.addExclusion(idx1, idx2)

        for i in range(self.chrSystem.getNumParticles()):
            LJNuc.addParticle([])

        self.chrSystem.addForce(LJNuc)

    def addSpeLJ(self, epsSpe = 3.0, sigmaSpe66 = 0.5, cutOffSpe = 1.5, eps67 = 1., eps77 = 1., sigmaSpe67_77 = 0.5, cutOff67_77 = 0.56):
        """
        Add nonbonded rescaled LJpotential between Speckles-Speckles

        Notice that the system has two different speckles dP (type 6) and P (type 7).

        Parameters
        ----------
        epsSpe (float, required) : 
            Epsilon of LJ between dP and dP, energy units (Default value = 3.0).

        sigmaSpe66 (float, required) : 
            Sigma of LJ between dP and dP, distance units (Default value = 0.5).

        cutOffSpe (float, required) : 
            Cutoff of LJ between dP and dP, distance units (Default value = 1.5).

        eps67, eps77 (float, required) :
            Epsilon of LJ between (dP,P) and (P,P), energy units (Default value = 1.0).

        sigmaSpe67_77 (float, required) : 
            Sigma of LJ between (dP,P) and (P,P), distance units (Default value = 0.5).

        curOff67_77 (float, required) :
            Cutoff of LJ between (dP,P) and (P,P), distance units (Default value = 0.5*1.12).
        """

        LJSpe = mm.CustomNonbondedForce("spec66 + spec77 + spec67;"
                "spec66 = delta(spe1-5)*delta(spe2-5)*4.*epsSpe*((sigmaSpe66/r)^12-(sigmaSpe66/r)^6-(sigmaSpe66/CutOffSpe)^12+(sigmaSpe66/CutOffSpe)^6)*step(CutOffSpe-r);"
                "spec67 = delta(max(spe1,spe2)-6)*delta(min(spe1,spe2)-5)*4.*eps67*((sigmaSpe67_77/r)^12-(sigmaSpe67_77/r)^6-(sigmaSpe67_77/CutOff67_77)^12+(sigmaSpe67_77/CutOff67_77)^6)*step(CutOff67_77-r);"
                "spec77 = delta(spe1-6)*delta(spe2-6)*4.*eps77*((sigmaSpe67_77/r)^12-(sigmaSpe67_77/r)^6-(sigmaSpe67_77/CutOff67_77)^12+(sigmaSpe67_77/CutOff67_77)^6)*step(CutOff67_77-r);")
        LJSpe.setCutoffDistance(cutOffSpe)
        LJSpe.addGlobalParameter("epsSpe", epsSpe)
        LJSpe.addGlobalParameter("sigmaSpe66", sigmaSpe66)
        LJSpe.addGlobalParameter('CutOffSpe', cutOffSpe)
        LJSpe.addGlobalParameter('eps67', eps67)
        LJSpe.addGlobalParameter('eps77', eps77)
        LJSpe.addGlobalParameter('sigmaSpe67_77', sigmaSpe67_77)
        LJSpe.addGlobalParameter('CutOff67_77', cutOff67_77)

        LJSpe.setNonbondedMethod(LJSpe.CutoffNonPeriodic) #avoid periodic boundary
        LJSpe.addPerParticleParameter("spe")

        LJSpe.setForceGroup(17)

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            LJSpe.addExclusion(idx1, idx2)

        for i in range(self.chrSystem.getNumParticles()):
            LJSpe.addParticle([self.compartType[i]])
        self.chrSystem.addForce(LJSpe)
			
    def addNAD(self, rescalarFile, epsilon_nad = 1.0, sigma_nad = 0.5, cutOff_nad = 1.5):
        '''
        Add nonbonded rescaled LJpotential between chromosomes and Nucleoli

        Parameters
        ----------
        rescalarFile (string, required) : 
            The path to the txt file of rescaling factors (SPIN state probabilities).

        epsilon_nad (float, required) : 
            Epsilon of LJ between nucleolus beads (Default value = 1.0).

        sigma_nad (float, required) :
            Sigma of LJ between nucleolus beads (Default value = 0.5).

        cutOff_nad (float, required):
            Cutoff of LJ between nucleolus beads (Default value = 1.5).
        '''

        nad_energy = ("nad(min(f1,f2))*step(max(f1,f2)-N_nad_chr)*step(N_nad_chr_nuc-max(f1,f2))*epsilon_nad*4.*((sigma_nad/r)^12-(sigma_nad/r)^6-(sigma_nad/cutOff_nad)^12+(sigma_nad/cutOff_nad)^6)*step(cutOff_nad-r);")

        rescalarMat = np.loadtxt(rescalarFile)
        M = self.chrSystem.getNumParticles()
        
        NAD_energy = mm.CustomNonbondedForce(nad_energy)

        rescalar = mm.Discrete1DFunction(rescalarMat)
        NAD_energy.addTabulatedFunction('nad', rescalar)
        NAD_energy.addGlobalParameter("sigma_nad", sigma_nad)
        NAD_energy.addGlobalParameter("epsilon_nad", epsilon_nad)
        NAD_energy.addGlobalParameter("cutOff_nad", cutOff_nad)
        NAD_energy.addGlobalParameter("N_nad_chr", self.N_chr-0.5)
        NAD_energy.addGlobalParameter("N_nad_chr_nuc", self.N_chr_nuc-0.5)
        NAD_energy.setCutoffDistance(cutOff_nad)

        NAD_energy.addPerParticleParameter("f")

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            NAD_energy.addExclusion(idx1, idx2)

        NAD_energy.setNonbondedMethod(NAD_energy.CutoffNonPeriodic) #avoid periodic boundary

        NAD_energy.setForceGroup(18)

        for i in range(self.chrSystem.getNumParticles()):
            NAD_energy.addParticle([i])

        self.chrSystem.addForce(NAD_energy)

    def addTSA(self, tsaTable, sigmatanhtsa = 4., rctanhtsa = 0.75, cutOfftsa = 2.0):
        r'''
        Add nonbonded potential using custom values for interactions between chromosomes and Speckles

        The parameters are part of the probability of function :math:`g(r_{i,j}) = \frac{1}{2}\left(1 + \text{tanh}\left[\sigma(r_c - r_{i,j})\right] \right)`, where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        tsaTable (string, required) : 
            The path to the potential scaling factor txt file.

        sigmatanhtsa (float, required) :
            Distance units (Default value = 4.0).

        rctanhtsa (float, required) : 
            Distance units, where :math:`g(r_{i,j}) = 0.5` (Default value = 0.75).

        cutOfftsa (float, required) : 
            The cutoff (Default value = 2.0)
        '''
        tsa_energy = ("chr_s6 + chr_s7;"
                "chr_s6 = tsa(min(p1,p2))*delta(max(v1,v2)-5)*(0.5*(1.+tanh(sigmatanhtsa*(rctanhtsa-r)))-0.5*(1.+tanh(sigmatanhtsa*(rctanhtsa-cutOfftsa))))*step(cutOfftsa-r);"
                "chr_s7 = tsa(min(p1,p2))*delta(max(v1,v2)-6)*(0.5*(1.+tanh(sigmatanhtsa*(rctanhtsa-r)))-0.5*(1.+tanh(sigmatanhtsa*(rctanhtsa-cutOfftsa))))*step(cutOfftsa-r);")
        tsaMat = np.loadtxt(tsaTable)
        M = self.chrSystem.getNumParticles()

        tsa_energy = mm.CustomNonbondedForce(tsa_energy)

        tsa = mm.Discrete1DFunction(tsaMat)
        tsa_energy.addTabulatedFunction('tsa', tsa)
        tsa_energy.addGlobalParameter("sigmatanhtsa", sigmatanhtsa)
        tsa_energy.addGlobalParameter("rctanhtsa", rctanhtsa)
        tsa_energy.addGlobalParameter("cutOfftsa", cutOfftsa)
        tsa_energy.setCutoffDistance(cutOfftsa)

        tsa_energy.addPerParticleParameter("p")
        tsa_energy.addPerParticleParameter("v")

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            tsa_energy.addExclusion(idx1, idx2)

        tsa_energy.setNonbondedMethod(tsa_energy.CutoffNonPeriodic) #avoid periodic boundary

        tsa_energy.setForceGroup(19)

        for i in range(self.chrSystem.getNumParticles()):
            tsa_energy.addParticle([i,self.compartType[i]])

        self.chrSystem.addForce(tsa_energy)

    def add_p_DamID(self, DamIDTable, sigmatanh_p_damid = 4., rctanh_p_damid = 0.75, cutOff_p_damid = 2.0):
        r'''
        Add nonbonded potential using custom values for interactions between chromosomes and Lamina

        The parameters are part of the probability of function :math:`g(r_{i,j}) = \frac{1}{2}\left(1 + \text{tanh}\left[\sigma(r_c - r_{i,j})\right] \right)`, where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        DamIDTable (string, required) :
            The path to the potential scaling factor txt file.

        sigmatanh_p_damid (float, required) :
            Distance units (Default value = 4.0).

        rctanh_p_damid (float, required) : 
            Distance units, where :math:`g(r_{i,j}) = 0.5` (Default value = 0.75).

        cutOff_p_damid (float, required) :
            The cutoff (Default value = 2.0)
        '''
        p_dam_energy = ("damid(min(pdam1,pdam2))*step(max(pdam1,pdam2)-N_damid)*(f-0.5*(1.+tanh(sigmatanh_p_damid*(rctanh_p_damid-cutOff_p_damid))))*step(cutOff_p_damid-r);"
                         "f = 0.5*(1.+tanh(sigmatanh_p_damid*(rctanh_p_damid-r)));")

        damMat = np.loadtxt(DamIDTable)
        M = self.chrSystem.getNumParticles()

        p_dam_energy = mm.CustomNonbondedForce(p_dam_energy)

        dam = mm.Discrete1DFunction(damMat)
        p_dam_energy.addTabulatedFunction('damid', dam)
        p_dam_energy.addGlobalParameter("sigmatanh_p_damid", sigmatanh_p_damid)
        p_dam_energy.addGlobalParameter("rctanh_p_damid", rctanh_p_damid)
        p_dam_energy.addGlobalParameter("cutOff_p_damid", cutOff_p_damid)
        p_dam_energy.addGlobalParameter("N_damid", self.N_chr_nuc_spec-0.5)
        p_dam_energy.setCutoffDistance(cutOff_p_damid)

        p_dam_energy.addPerParticleParameter("pdam")

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            p_dam_energy.addExclusion(idx1, idx2)

        p_dam_energy.setNonbondedMethod(p_dam_energy.CutoffNonPeriodic) #avoid periodic boundary

        p_dam_energy.setForceGroup(20)

        for i in range(self.chrSystem.getNumParticles()):
            p_dam_energy.addParticle([i])

        self.chrSystem.addForce(p_dam_energy)

    def addinter_chrom(self, InterFile, rc_tanh_inter=0.54, cutOff_inter = 2.0):
        r'''
        Add the Inter-chromosome Chromosome potential using custom values for interactions between beads separated by a genomic distance :math:`d`.

        The parameter rc_tanh_inter is part of the probability of crosslink function :math:`f(r_{i,j})=\frac{1}{2}\left(1+\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right)+\frac{1}{2}\left(1-\text{tanh}\left[(r_c - r_{i,j})^{-5}+5(r_c - r_{i,j})\right] \right)\left(\frac{r_c}{r_{i,j}}\right)^4`,

        where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        Parameters
        ----------
        InterFile (string, required) :
            The path to the inter-chromosomal potential scaling factor txt file.

        rc_tanh_inter (float, required) :
            :math:`r_c` in the equation, a good estimation of the bond length during the simulation (Default value = 0.54 for the 100KB model).

        cutOff_inter (float, required) :
            Cutoff of inter-chromosomal potential (Default value = 2.0 for the 100KB model).
        '''
        inter_energy = ("mapid(m1,m2)*(tanh1_inter+tanh2_inter*(rc_tanh_inter/r)^4-1/2*(1.+tanh((rc_tanh_inter-cutOff_inter)^(-5)+5*(rc_tanh_inter-cutOff_inter)))-1/2*(1.-tanh((rc_tanh_inter-cutOff_inter)^(-5)+5*(rc_tanh_inter-cutOff_inter)))*(rc_tanh_inter/cutOff_inter)^4)*step(cutOff_inter-r);"
                        "tanh1_inter = 1/2*(1.+tanh((rc_tanh_inter-r)^(-5)+5*(rc_tanh_inter-r)));"
                        "tanh2_inter = 1/2*(1.-tanh((rc_tanh_inter-r)^(-5)+5*(rc_tanh_inter-r)));")

        cross_inter = mm.CustomNonbondedForce(inter_energy)
        
        cross_inter.addGlobalParameter('rc_tanh_inter', rc_tanh_inter)
        cross_inter.addGlobalParameter('cutOff_inter', cutOff_inter)
        cross_inter.setCutoffDistance(cutOff_inter)
        cross_inter.setForceGroup(21)

        tab_inter = np.zeros((23, 23))
        with open(InterFile, "r") as fin_inter:
            for line in fin_inter.readlines():
                line = line.split()
                i = int(line[0]) - 1
                j = int(line[1]) - 1
                tab_inter[i, j] = tab_inter[j, i] = float(line[2])

        tab_inter = list(np.ravel(tab_inter))
        diff_chrom_size = 23

        fTypes_inter = mm.Discrete2DFunction(diff_chrom_size, diff_chrom_size, tab_inter)
        cross_inter.addTabulatedFunction('mapid', fTypes_inter)

        cross_inter.addPerParticleParameter("m")

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            cross_inter.addExclusion(idx1, idx2)

        cross_inter.setNonbondedMethod(cross_inter.CutoffNonPeriodic) #avoid periodic boundary

        for i in range(self.chrSystem.getNumParticles()):
            cross_inter.addParticle([min(int(np.floor(self.molType[i]/2.)),22)])

        self.chrSystem.addForce(cross_inter)

    def save_system(self, system_xml):
        '''
        Save the simulation state to readable xml format.
        '''
        with open(system_xml, 'w') as output_writer:
            output_writer.write(mm.XmlSerializer.serialize(self.chrSystem))

    def createSimulation(self, platform_type = 'CPU'):
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
            simulation = mmapp.Simulation(self.chrTopo, self.chrSystem, integrator, platform, properties)
        elif platform_type == 'OpenCL':
            platform = mm.Platform.getPlatformByName('OpenCL')
            properties = {'OpenCLPrecision': 'double'}
            simulation = mmapp.Simulation(self.chrTopo, self.chrSystem, integrator, platform, properties)
        elif platform_type == 'Reference':
            platform = mm.Platform.getPlatformByName('Reference')
            simulation = mmapp.Simulation(self.chrTopo, self.chrSystem, integrator, platform)
        elif platform_type == 'CPU':
            platform = mm.Platform.getPlatformByName('CPU')
            simulation = mmapp.Simulation(self.chrTopo, self.chrSystem, integrator, platform)
        else:
            print("platform_type can be either CUDA, OpenCL, or CPU")
        return simulation
          
