## Nuclear Deformation with Langevin Dynamics

With Langevin dynamics, this folder provides tutorials for setting up and performing Nuclear Deformation simulations of the HFF nucleus at the 100KB resolution with the presence of various nuclear landmarks, including nuclear lamina, speckles, and nucleoli. And here, we also show how to create the system with customized setup.

- Nuclear_Deformation.ipynb: The jupyter notebook showing how to create a system and run a simulation with the nucleus deformation. Make sure you put all the force field parameter files in the parameter folder.

- parameters: The folder saves all the parameters that the simulation needs

- input.csv: The file logs which potential that the users hope to add

- lamina_bond.txt: The file logs where the lamina connected bonds are

- step_100000.dcd: The trajectory from the Nuclear Deformation with Langevin Dynamics

- reduced_traj.dcd: The trajectory after converting the OpenMM default unit to the reduced unit

- compt_final_frame.txt: A file saves the type of each bead of the last frame

- human.pdb: An initial configuration used for creating a simulation

- human_final.pdb: The configuration of the last frame from the trajectory (reduced_traj.dcd)

- coordinate_transformation.py: This file converts the trajectory with position in OpenMM default unit to reduced unit 

- final_frame.py: We output the final configuration according to the converted trajectory after running script "coordinate_transformation.py"

- bead_info.txt: This file contains the information of each chromatin bead, and is necessary when creating the final configuration. 1st column: the index of each bead; 2nd column: the index of belonging chromosome; 3rd column: the type of each chromatin bead.

- strength_zero: This folder contains the initial structure and final structure with strength = 0.0 when squeezing the nucleus

- strength_half: This folder contains the initial structure and final structure with strength = 0.5 when squeezing the nucleus

- strength_one: This folder contains the initial structure and final structure with strength = 1.0 when squeezing the nucleus

The pipeline is as follows:

    1. Obtain the trajectories from Nuclear_Deformation.ipynb

    2. Convert the trajectory with position in OpenMM default unit to reduced unit

    python coordinate_transformation.py step_100000.dcd

    3. Create a configuration file with the last frame

    python coordinate_transformation.py reduced_traj.dcd
