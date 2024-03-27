## New feature

This tutorial explains how to introduce new features into the model. The following changes to the source code were made to include a new of particles that can binding with centromeres. 

- chromosome.py, nucleolus.py, speckle.py, lamina.py: source code in the original OpenNucleome software package and are left unchanged.

- newfeature.py: This new file defines the interaction potentials between the newly introduced particles and other nuclear components. Here, for simplicity, the code only defines potentials between the new particles and between them and centromeres.

- whole_nucleus_model.py: We made the following changes to this file to accommodate the new feature. 
  
    * 1\. Line 80: Change the number of bead types from 8 to 9 (Compartment A: 1, Compartment B: 2, Pericentromeres: 3, Centromeres: 4, Nucleoli: 5, dP Speckles: 6, P Speckles: 7, New Beads: 8, Lamina: 9), and always set lamina beads as the last type.

    * 2\. Line 121-139: In the old file, we calculate $N_{chr}$, $N_{chr}+N_{nuc}$, $N_{chr}+N_{nuc}+N_{spec}$. In the new file, in addition to these three numbers, we also calculate $N_{chr}+N_{nuc}+N_{spec}+N_{newbeads}$

    * 3\. Line 150: Similar to other components in the system, we create an object for the new bead types, and we can add corresponding potentials by calling the functions in the specific class for the new type beads.

    * 2.4\. Line 165: Add a temporary name for the new bead types.

    * 2.5\. Line 279-290: Create a function by adding the new potentials into the system, and as shown in line 306, 308, 310, 312, we also modify the indexes of potentials.

- configs_generation: To include the newly introduced particles in the system, scripts for generating initial configurations need to be updated as well. The revised scripts are included in this folder. 
