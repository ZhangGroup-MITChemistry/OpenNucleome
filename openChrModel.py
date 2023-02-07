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
    def __init__(self, temperature = 1.0, gamma = 0.1, timestep = 0.01, mass_scale=1.0):
        self.temperature = temperature / 0.008314 # temperature in reduced unit
        self.gamma = gamma # fiction coefficient (1/time)
        self.timestep = timestep  # timestep in reduced unit (time)
        self.kB = u.BOLTZMANN_CONSTANT_kB * u.AVOGADRO_CONSTANT_NA
        self.kT = self.kB * self.temperature
        self.chrMass = 10.0 * u.amu * mass_scale

    def createSystem(self, PDBFile):
        """
        PDBFile: the path to the chromosome PDB file
        """
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
        self.beadGroups = [[] for _ in range(7)] # group the indexes of atoms by their types; if you do not want to add some nuclear landmarks, then use the different value. compartment A: 1, compartment B: 2, regions around the centromeres: 3, centromeres: 4, nucleoli: 5, speckles: 6, lamina: 7
        self.chrGroups = [[] for _ in range(46)] # represent the 46 chromosomes of human
        prevResID = -1
        self.chrResidues = [] #(start, end) index of the atoms for each residue
        start = 0

        for a in self.chrPDB.topology.atoms():
            #record compartment type info
            self.compartType.append(int(a.name) - 1)
            self.molType.append(int(a.residue.index))

            #add atom into the system with mass = chrMass
            m = 1. if self.compartType[-1] != 6 else 0. # set mass of Lamina beads to zero to freeze them
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
        for i in range(6):
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
        name2Element = ['ASP', 'GLU', 'HIS', 'LYS', 'ARG', 'ASN', 'GLN']
        self.Elements = []
        for i in range(7):
            m = 1. if i!= 6 else 0.
#           self.Elements.append(mmapp.element.Element.getBySymbol(name2Element[i]))
            self.Elements.append(mmapp.Element(1000+i, name2Element[i], name2Element[i], self.chrMass * m))

    def constructTopology(self):
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
    def addFeneBond(self, kfb = 30.0, r0 = 1.5, epsFene = 1.0, sigmaFene = 1.0, cutFene = 2.**(1./6.)):
        '''
        epsFene: epsilon
        sigmaFene: sigma
        cutFene: cutoff distance
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
        kc2, kc3, kc4: bond coefficients
        r0c: distance
        Currently, the Class2 bond is for the 100KB resolution
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
        """
        k_a: angle coefficients, default = 2.0
        """
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

    def addSoftCore(self, epsSC = 1.0, sigmaSC = 0.5, cutOffSC = 0.5*(2.**(1./6.)), EcutSC = 4.0):
        """
        We plug the tanh potential into the WCA potential to represent the excluded volume effect;
        Therefore, the softcore potential will only be applied to the chromosomes (type 1, 2, 3, 4).
        Here, we use step function to select the atoms.
        When r < r_0SC, tanh potential; r_0SC < r < cutOffSC, WCA potential, r_0SC ~ 0.485
        """
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

    def addIdealPotential(self, ICFile, sigmatanhIC = 11.349, rctanhIC = 0.646, dendIC = 1000, cutOffIC = 4.0):
        '''
        Adds the Intra-chromosome Ideal Chromosome potential using custom values for interactions between beads separated by a genomic distance :math:`d`. 
        The parameters :sigmatanhIC and rctanhIC are part of the probability of crosslink function :math:`f(r_{i,j}) = \frac{1}{2}\left( 1 + tanh\left[sigmatanhIC(rctanhIC - r_{i,j}\right] \right)`, 
        where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.
        ICFile: the path to the ideal potential scaling factor file (txt)
        dendIC: cutoff of genomic separation, default: 1000 for 100KB model
        cutOffIC: cutoff of IC: default = 4.0
        '''
        energyIC = ("IClist(d) * step(dendIC - d) * step(cutOffIC - r) * (f * step(rctanhIC - r) + g * step(r-rctanhIC) - 0.5*(rctanhIC / cutOffIC) ^ 4);"
                        "f = 0.5*(1. + tanh(sigmatanhIC * (rctanhIC - r)));"
                        "g = 0.5*(rctanhIC / r) ^ 4;"
                        "d = abs(idx2-idx1)")

        IC = mm.CustomNonbondedForce(energyIC)
        IClist_listfromfile = np.loadtxt(ICFile)[:,1] 

        tabIClist = mm.Discrete1DFunction(IClist_listfromfile)
        IC.addTabulatedFunction('IClist', tabIClist) 

        IC.addGlobalParameter('dendIC', dendIC)

        IC.addGlobalParameter('sigmatanhIC', sigmatanhIC)  
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

    def addType2TypePotential(self, TypesTable, sigmatanhLP = 11.349, rctanhLP = 0.646, crossLPcutOff = 4.0):
        """
        Adds the type-to-type potential using custom values for interactions between the chromatin types. 
        The parameters :sigmatanhLP and rctanhLP are part of the probability of crosslink function :math:`f(r_{i,j}) = \frac{1}{2}\left( 1 + tanh\left[sigmatanhLP(rctanhLP - r_{i,j}\right] \right)`, where :math:`r_{i,j}` is the spatial distance between loci (beads) *i* and *j*.

        The function receives a txt/TSV/CSV file containing the upper triangular matrix of the type-to-type interactions. A file example can be found `here <https://www.ndb.rice.edu>`__.

        +---+-------+-------+-------+-------+
        |   |   A   |   B   |   C   |   D   |
        +---+-------+-------+-------+-------+
        | A | xxxxx | xxxxx | xxxxx |   0
        +---+-------+-------+-------+-------+
        | B |       | xxxxx | xxxxx |   0
        +---+-------+-------+-------+-------+
        | C |       |       | xxxxx |   0
        +---+-------+-------+-------+-------+
        | D |       |       |       |   0
        +---+-------+-------+-------+-------+

        sigmatanhLP (float, required):
          Parameter in the probability of crosslink function. (Default value = 11.349).
        rctanhLP (float, required):
          Parameter in the probability of crosslink function. (Default value = 0.646).
        TypesTable (file, required):
          A txt file containing the type-to-type interactions. (Default value: :code:`None`).
        """

        energy = ("mapType(t1,t2) * (f * step(rctanhLP - r) + g * step(r-rctanhLP) - 0.5*(rctanhLP/crossLPcutOff)^4 ) * step(crossLPcutOff - r);"
                        "f = 0.5*(1. + tanh(sigmatanhLP*(rctanhLP - r)));"
                        "g = 0.5*(rctanhLP / r) ^ 4;")

        crossLP = mm.CustomNonbondedForce(energy)

        crossLP.addGlobalParameter('sigmatanhLP', sigmatanhLP)
        crossLP.addGlobalParameter('rctanhLP', rctanhLP)
        crossLP.addGlobalParameter('crossLPcutOff', crossLPcutOff)
        crossLP.setCutoffDistance(crossLPcutOff)
        
        crossLP.setForceGroup(13)
        #generate the above shown 4*4 matrix, 4 is the number of chr bead types 
        tab = np.zeros((6, 6))
        with open(TypesTable, "r") as fin:
            for line in fin.readlines():
                line = line.split()
                i = int(line[1]) - 1
                j = int(line[2]) - 1
                tab[i, j] = tab[j, i] = float(line[5])

        tab = list(np.ravel(tab))
        diff_types_size = 6

        fTypes = mm.Discrete2DFunction(diff_types_size, diff_types_size, tab)
        crossLP.addTabulatedFunction('mapType', fTypes) 

        crossLP.addPerParticleParameter("t")

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            crossLP.addExclusion(idx1, idx2)

        crossLP.setNonbondedMethod(crossLP.CutoffNonPeriodic) #avoid periodic boundary

        for i in range(self.chrSystem.getNumParticles()):
            crossLP.addParticle([self.compartType[i]])

        self.chrSystem.addForce(crossLP)   

    def addhardwall(self, R, epshw = 0.5, sigmahw = 0.5, CutOffhw = 0.5*(2.**(1./6.))):
        '''
        If your system doesn't contain the lamina to prevent chromosomes going outside, then, add this force
        add external force (WCA repulsion) by adding a hard wall for chromosomes
        epshw: epsilon of lj/cut, default: 0.5
        sigmahw: sigma of lj/cut, default: 0.5
        CutOffhw: cutoff of lj/cut, default: 0.5*2. ** (1. / 6.)
        '''
        hw = mm.CustomExternalForce("4. * epshw * ((sigmahw/f)^12-(sigmahw/f)^6) * step(CutOffhw - f);"
                "f = R - sqrt(x^2+y^2+z^2)")
        hw.addGlobalParameter("epshw", epshw)
        hw.addGlobalParameter("sigmahw", sigmahw)
        hw.addGlobalParameter("CutOffhw", CutOffhw)
        hw.addGlobalParameter("R", R)

        hw.setForceGroup(14)
        for i in range(self.chrSystem.getNumParticles()):
            hw.addParticle(i,[])

        self.chrSystem.addForce(hw)

    def addparticle_hw(self, eps_p_hw = 0.5, sigma_p_hw = 0.75, cutoff_p_hw = 0.75*(2.**(1./6.))):
        """
        add nonbonded hard-wall potential between (chr, speckles, nucleoli) and lamina
        The size of lamina is 1.0 (diameter), so the sigma should be r_lamina + r_chromain = 0.75 for the 100KB model
        curOffLJ: cutoff, default: 0.75*1.12
        """
        #load the bead-specific rescaling factors

        particle_hw_energy = mm.CustomNonbondedForce("step(max(phw1,phw2)-N_hw) * step(N_hw-min(phw1,phw2)) * 4*((sigma_p_hw/r)^12-(sigma_p_hw/r)^6-(sigma_p_hw/cutoff_p_hw)^12+(sigma_p_hw/cutoff_p_hw)^6) * step(cutoff_p_hw-r)")

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
        """
        add nonbonded plain LJpotential between Nuc-Spe
        """
        CutOffLJplain = 2. ** (1 / 6) #1.12

        LJplain = mm.CustomNonbondedForce("(4.*epsLJplain * ((sigLJplain/r)^12-(sigLJplain/r)^6) + epsLJplain) * step(CutOffLJplain - r)")
        LJplain.setCutoffDistance(CutOffLJplain)
        LJplain.addGlobalParameter("epsLJplain", 1.0)
        LJplain.addGlobalParameter("sigLJplain", 1.0)
        LJplain.addGlobalParameter('CutOffLJplain', CutOffLJplain)

        LJplain.addInteractionGroup(self.beadGroups[4], self.beadGroups[5])
        LJplain.setNonbondedMethod(LJplain.CutoffNonPeriodic) #avoid periodic boundary

        LJplain.setForceGroup(15)
        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            LJplain.addExclusion(idx1, idx2)

        for i in range(self.chrSystem.getNumParticles()):
            LJplain.addParticle([])

        self.chrSystem.addForce(LJplain)

    def addLJNuc(self, epsNuc = 6.0, sigmaNuc = 0.5, cutOffNuc = 2.0):
        """
        add nonbonded rescaled LJpotential between Nuc-Nuc
        epsNuc: float, epsilon of LJ, default: 6.0
        sigNuc: float, sigma of LJ, default: 0.5
        cutOffNuc: float, cut off of LJ, default: 2.0
        """
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

    def addSpeLJ(self, epsSpe = 8.0, sigmaSpe = 0.5, cutOffSpe = 2.0, XiSpe = 0.95, cutOffYukawa = 6.0):
        """
        add nonbonded rescaled LJpotential between Spe-Spe
        epsSpe: float, epsilon of LJ, default: 8.0
        sigSPe: float, sigma of LJ, default: 0.5
        CutOffSpe: float, cut off of LJ, default: 2.0
        XiSpe: Xi for Yukawa potential
        cutOffYukawa: cutOff scaling for Yukawa potential
        """
        cutOffYukawa *= XiSpe

        LJSpe = mm.CustomNonbondedForce("LJspe + 2.5 * XiSpe * (exp(-r / XiSpe) / r - exp(-CutOffYukawa / XiSpe) / CutOffYukawa) * step(CutOffYukawa - r);"
                                          "LJspe = 4.*epsSpe*((sigmaSpe/r)^12-(sigmaSpe/r)^6 - (sigmaSpe/CutOffSpe)^12 + (sigmaSpe/CutOffSpe)^6) * step(CutOffSpe - r);")
        LJSpe.setCutoffDistance(cutOffYukawa)
        LJSpe.addGlobalParameter("epsSpe", epsSpe)
        LJSpe.addGlobalParameter("sigmaSpe", sigmaSpe)
        LJSpe.addGlobalParameter("XiSpe", XiSpe)
        LJSpe.addGlobalParameter('CutOffSpe', cutOffSpe)
        LJSpe.addGlobalParameter('CutOffYukawa', cutOffYukawa)

        LJSpe.addInteractionGroup(self.beadGroups[5], self.beadGroups[5])
        LJSpe.setNonbondedMethod(LJSpe.CutoffNonPeriodic) #avoid periodic boundary

        LJSpe.setForceGroup(17)

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            LJSpe.addExclusion(idx1, idx2)

        for i in range(self.chrSystem.getNumParticles()):
            LJSpe.addParticle([])
        self.chrSystem.addForce(LJSpe)
			
    def addNAD(self, rescalarFile, epsilon_nad = 0.5, sigma_nad = 0.5, cutOff_nad = 2.0):
        """
        add nonbonded rescaled LJpotential between chr and Nucleoli
        rescalarFile: path to the txt file of rescaling factors
        curOff_nad: cutoff, default: 2.0
        """
        #load the bead-specific rescaling factors

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

    def addTSA(self, tsaTable, sigmatanhtsa = 4., rctanhtsa = 0.75, cutOfftsa = 4.0):
        """
        add nonbonded rescaled LJpotential between chr and Speckles
        curOffLJ: cutoff, default: 4.0
        """
        #load the bead-specific rescaling factors

        tsa_energy = ("tsa(min(p1,p2))*step(max(p1,p2)-N_tsa_chr_nuc)*step(N_tsa_chr_nuc_spec-max(p1,p2))*(f*step(rctanhtsa-r)+g*step(r-rctanhtsa)-0.5*(rctanhtsa/cutOfftsa)^4)*step(cutOfftsa-r);"
                         "f = 0.5*(1.+tanh(sigmatanhtsa*(rctanhtsa-r)));"
                         "g = 0.5*(rctanhtsa/r)^4;")

        tsaMat = np.loadtxt(tsaTable)
        M = self.chrSystem.getNumParticles()

        tsa_energy = mm.CustomNonbondedForce(tsa_energy)

        tsa = mm.Discrete1DFunction(tsaMat)
        tsa_energy.addTabulatedFunction('tsa', tsa)
        tsa_energy.addGlobalParameter("sigmatanhtsa", sigmatanhtsa)
        tsa_energy.addGlobalParameter("rctanhtsa", rctanhtsa)
        tsa_energy.addGlobalParameter("cutOfftsa", cutOfftsa)
        tsa_energy.addGlobalParameter("N_tsa_chr_nuc", self.N_chr_nuc-0.5)
        tsa_energy.addGlobalParameter("N_tsa_chr_nuc_spec", self.N_chr_nuc_spec-0.5)
        tsa_energy.setCutoffDistance(cutOfftsa)

        tsa_energy.addPerParticleParameter("p")

        for idx1, idx2 in self.Allbonds: #exclude all nearest neighbor interactions
            tsa_energy.addExclusion(idx1, idx2)

        tsa_energy.setNonbondedMethod(tsa_energy.CutoffNonPeriodic) #avoid periodic boundary

        tsa_energy.setForceGroup(19)

        for i in range(self.chrSystem.getNumParticles()):
            tsa_energy.addParticle([i])

        self.chrSystem.addForce(tsa_energy)

    def add_p_DamID(self, DamIDTable, sigmatanh_p_damid = 4., rctanh_p_damid = 0.75, cutOff_p_damid = 2.0):
        """
        add nonbonded rescaled LJpotential between chr and Lamina
        curOffLJ: cutoff, default: 1.5
        """
        #load the bead-specific rescaling factors

        p_dam_energy = ("damid(min(pdam1,pdam2))*step(max(pdam1,pdam2)-N_damid)*(f*step(rctanh_p_damid-r)+g*step(r-rctanh_p_damid)-0.5*(rctanh_p_damid/cutOff_p_damid)^4)*step(cutOff_p_damid-r);"
                         "f = 0.5*(1.+tanh(sigmatanh_p_damid*(rctanh_p_damid-r)));"
                         "g = 0.5*(rctanh_p_damid/r)^4;")

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


    def addDamID(self, R, DamIDfile, sigmatanhdamid = 4.0, rctanhdamid = 0.75, cutoffdamid = 4.0):
        """
        If your system doesn't contact lamina beads, then, add this potential to do the optimization for DamID
        Add the potential of contact function (tanh potential) between chr and sphere confinement.
        DamIDfile: the txt file of alpha value
        sigmatanhdamid, rctanhdamid: results from imaging data,
        cutOffIC: cutoff of IC: default = 4.0
        """
        energyDamID = ("alpha * step(cutoffdamid - r) * (0.5 * (1. + tanh(sigmatanhdamid * (rctanhdamid - r))) * step(rctanhdamid - r) + 0.5 * (rctanhdamid / r) ^ 4 * step(r - rctanhdamid) - 0.5*(rctanhdamid / cutoffdamid) ^ 4);"
        "r = R - sqrt(x^2+y^2+z^2)")

        rescalar = np.loadtxt(DamIDfile) #alpha values for DamID

        DamID_energy = mm.CustomExternalForce(energyDamID)
        DamID_energy.addGlobalParameter("sigmatanhdamid", sigmatanhdamid)
        DamID_energy.addGlobalParameter("rctanhdamid", rctanhdamid)
        DamID_energy.addGlobalParameter("cutoffdamid", cutoffdamid)
        DamID_energy.addGlobalParameter("R", R)
        DamID_energy.setForceGroup(20)
        #add a new parameter (index) for each bead
        DamID_energy.addPerParticleParameter("alpha")

        #create the new parameter (index)
        for i in range(self.chrSystem.getNumParticles()):
            DamID_energy.addParticle(i, [rescalar[i]])

        self.chrSystem.addForce(DamID_energy)

    def addinter_chrom(self, InterFile, sigmatanh_inter = 11.349, rctanh_inter = 0.646, cutOff_inter = 4.0):
        inter_energy = ("mapid(m1,m2) * (f * step(rctanh_inter - r) + g * step(r - rctanh_inter) - 0.5*(rctanh_inter/cutOff_inter)^4) * step(cutOff_inter - r);"
                        "f = 0.5*(1. + tanh(sigmatanh_inter*(rctanh_inter - r)));"
                        "g = 0.5*(rctanh_inter / r) ^ 4;")

        cross_inter = mm.CustomNonbondedForce(inter_energy)
        
        cross_inter.addGlobalParameter('sigmatanh_inter', sigmatanh_inter)
        cross_inter.addGlobalParameter('rctanh_inter', rctanh_inter)
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

    def createSimulation(self, platform_type):
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
          
