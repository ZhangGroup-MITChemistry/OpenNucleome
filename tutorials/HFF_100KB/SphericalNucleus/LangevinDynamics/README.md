## Spherical Nucleus with Langevin Dynamics

This folder provides tutorials for setting up and performing Langevin Dynamics simulations of the HFF nucleus at the 100KB resolution with the presence of various nuclear landmarks, including nuclear lamina, speckles, and nucleoli.

- Simulation_Langevin.ipynb: The jupyter notebook showing how to create a system and run a simulation. Make sure you put all the force field parameter files in the parameter folder.

- parameters: The folder saves all the parameters that the simulation needs

- step_100000.dcd: The trajectory from the Langevin Dynamics

- reduced_traj.dcd: The trajectory after converting the OpenMM default unit to the reduced unit

- compt_final_frame.txt: A file saves the type of each bead of the last frame

- human.pdb: An initial configuration used for creating a simulation

- human_final.pdb: The configuration of the last frame from the trajectory (reduced_traj.dcd)

- coordinate_transformation.py: This file converts the trajectory with position in OpenMM default unit to reduced unit 

- final_frame.py: We output the final configuration according to the converted trajectory after running script "coordinate_transformation.py"

- bead_info.txt: This file contains the information of each chromatin bead, and is necessary when creating the final configuration. 1st column: the index of each bead; 2nd column: the index of belonging chromosome; 3rd column: the type of each chromatin bead.

- initial_configs: This folder contains 50 initial configurations, which were chosen according to the "Section: Initial configurations for simulations"

- final_configs: This folder contains the corresponding 50 final configurations after running 10,000,000 steps for each initial configuration

The pipeline is as follows:

(i) Obtain the trajectories from Simulation_Langevin.ipynb


(ii) Convert the trajectory with position in OpenMM default unit to reduced unit

```
python coordinate_transformation.py step_100000.dcd
```

(iii) Create a configuration file with the last frame

```
python coordinate_transformation.py reduced_traj.dcd
```
