## New feature

As an example, this folder provides tutorials for including new features into the existing model.

- chromosome.py, nucleolus.py, speckle.py, lamina.py: source code.

- whole_nucleus_model.py: users should change some lines for adding new features, and a more comprehensive tutorials will be followed.

- newfeature.py: save the potentials between new beads and other components. Here, as an example, we add the potentials between 1. the new beads and new beads; 2. the new beads and centromeres.

- configs_generation: generate the configuration with new beads.

Here, we show how to add new features.

* 1\. Create the configurations with new beads, a comprehensive tutorial is included in the folder "./configs_generation".

* 2\. Modify the functions in the script "whole_nucleus_model.py".

    * 2.1\. Line 80: Change the number of beads from 8 to 9 (Compartment A: 1, Compartment B: 2, Pericentromeres: 3, Centromeres: 4, Nucleoli: 5, dP Speckles: 6, P Speckles: 7, New Beads: 8, Lamina: 9), and always set the lamina bead as the last type.

    * 2.2\. Line 92: Create 9 lists instead of 8 lists for saving the indexes of specific type beads.

    * 2.3\. Line 121-139: In the old file, we calculate $N_{chr}$, $N_{chr}+N_{nuc}$, $N_{chr}+N_{nuc}+N_{spec}$. In the new file, in addition to these three numbers, we also calculate $N_{chr}+N_{nuc}+N_{spec}+N_{newbeads}$

    * 2.4\. Line 150: Similar to other components in the system, we create an object for the new type beads, and we can add corresponding potentials by calling the functions in the specific class for the new type beads.

    * 2.5\. Line 165: Add a fake name for the new type beads.

    * 2.6\. Line 279-290: Create a function by adding the new potentials into the system, and as shown in line 306, 308, 310, 312, we also modify the indexes of potentials.

* 3\. Add a new file "newfeature.py" into the source code. As an example, we only include the potentials between 1. the new beads and the new beads; 2. the new beads and centromeres.
