## Spherical Nucleus with Brownian Dynamics

This folder provides tutorials for setting up and performing Brownian Dynamics simulations of the HFF nucleus at the 100KB resolution with the presence of various nuclear landmarks, including nuclear lamina, speckles, and nucleoli.

- BrownianSimulation.ipynb: The jupyter notebook showing how to create a system and run a simulation. Make sure you put all the force field parameter files in the parameter folder.

- parameters: The folder saves all the parameters that the simulation needs

- step_100000.dcd: The trajectory from the Brownian Dynamics

- reduced_traj.dcd: The trajectory after converting the OpenMM default unit to the reduced unit

- compt_final_frame.txt: A file saves the type of each bead of the last frame

- human.pdb: An initial configuration used for creating a simulation

- human_final.pdb: The configuration of the last frame from the trajectory (reduced_traj.dcd)

- coordinate_transformation.py: This file converts the trajectory with position in OpenMM default unit to reduced unit 

- final_frame.py: We output the final configuration according to the converted trajectory after running script "coordinate_transformation.py"

- bead_info.txt: This file contains the information of each chromatin bead, and is necessary when creating the final configuration. 1st column: the index of each bead; 2nd column: the index of belonging chromosome; 3rd column: the type of each chromatin bead.

- initial_configs: This folder contains 8 initial configurations

- final_configs: This folder contains the corresponding 8 final configurations after running 30,000,000 steps for each initial configuration

The pipeline is as follows:

(i) Obtain the trajectories from BrownianSimulation.ipynb

(ii) Convert the trajectory with position in OpenMM default unit to reduced unit

```
python coordinate_transformation.py step_100000.dcd
```

(iii) Create a configuration file with the last frame

```
python coordinate_transformation.py reduced_traj.dcd
```
