## Configuration relaxation

This folder provides tutorials for relaxing the initial configurations from the folder "../configs_generation/init_config_pool"

- ConfigsRelax.ipynb: The jupyter notebook showing how to relax the configurations, same as the OpenNucleome/tutorials/HFF_100KB/SphericalNucleus/LangevinDynamics/Simulation_Langevin.ipynb. Make sure you put all the force field parameter files in the parameter folder.

- parameters: The folder saves all the parameters that the simulation needs

- step_100000.dcd: The trajectory from the Langevin Dynamics

- reduced_traj.dcd: The trajectory after converting the OpenMM default unit to the reduced unit

- compt_final_frame.txt: A file saves the type of each bead of the last frame

- human.pdb: An initial configuration used for creating a simulation, same as "../configs_generation/init_config_pool/human_1.pdb"

- human_final.pdb: The configuration of the last frame from the trajectory (reduced_traj.dcd)

- coordinate_transformation.py: This file converts the trajectory with position in OpenMM default unit to reduced unit 

- final_frame.py: We output the final configuration according to the converted trajectory after running script "coordinate_transformation.py"

- bead_info.txt: This file contains the information of each chromatin bead, and is necessary when creating the final configuration. 1st column: the index of each bead; 2nd column: the index of belonging chromosome; 3rd column: the type of each chromatin bead.

The pipeline is as follows:

    1. Obtain the trajectories from ConfigsRelax.ipynb

    2. Convert the trajectory with position in OpenMM default unit to reduced unit

    python coordinate_transformation.py step_100000.dcd

    3. Create a configuration file with the last frame

    python coordinate_transformation.py reduced_traj.dcd
