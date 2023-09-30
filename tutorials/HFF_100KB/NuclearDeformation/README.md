## NuclearDeformation 

This folder provides tutorials for setting up and performing simulations of the HFF nucleus at the 100KB resolution with the presence of various nuclear landmarks, including nuclear lamina, speckles, and nucleoli. Different from the sphere situation, we applied a force to squeeze the nucleus.

- simulation.ipynb: The jupyter notebook showing how to create a system and run a simulation

- human.pdb: An initial configuration used for creating a simulation.

- coordinate_transformation.py: This file converts the trajectory with position in OpenMM default unit to reduced unit 

- final_frame.py: We output the final configuration according to the converted trajectory after running script "coordinate_transformation.py"

- bead_info.txt: This file contains the information of each chromatin bead, and is necessary when creating the final configuration. 1st column: the index of each bead; 2nd column: the index of belonging chromosome; 3rd column: the type of each chromatin bead.

- human_final.pdb: The final configuration after finishing the simulation.
