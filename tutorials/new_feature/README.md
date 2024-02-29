## New feature

This folder provides tutorials for including new features into the existing model

- NewFeature.ipynb: The jupyter notebook shows how to run the simulations with the new features. Make sure you put all the force field parameter files in the parameter folder.

- parameters: The folder saves all the parameters that the simulation needs

- step_100000.dcd: The trajectory from the Langevin Dynamics

- compt_final_frame.txt: A file saves the type of each bead of the last frame

- human.pdb: An initial configuration used for creating a simulation, same as "../configs_generation/init_config_pool/human_1.pdb"

- whole_nucleus_model.py, chromosome.py, nucleolus.py, speckle.py, lamina.py: source code

- newfeature.py: save the potential between new beads and other components. Here, as an example, we add the potential between 1. new beads and new beads 2. new beads and centromeres

- configs_generation: generate the configuration with new beads
