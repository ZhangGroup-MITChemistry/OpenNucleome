## Nuclear Deformation with Langevin Dynamics

With Langevin dynamics, this folder provides tutorials for setting up and performing Nuclear Deformation simulations of the HFF nucleus at the 100KB resolution with the presence of various nuclear landmarks, including nuclear lamina, speckles, and nucleoli. And here, we also show how to create the system with customized setup.

- input_params: The folder saves all the customized parameters (input_params/input.csv: The file logs which potential that the users hope to add; input_params/lamina_bond.txt: The file logs where the lamina connected bonds are)

- NuclearDeformation.ipynb: The jupyter notebook showing how to create a system and run a simulation with the nucleus deformation.

- results: This folder contains the trajectory and new configuration files after running the notebook "NuclearDeformation.ipynb".

- strength_zero: This folder contains the initial structure and final structure with strength = 0.0 when squeezing the nucleus.

- strength_half: This folder contains the initial structure and final structure with strength = 0.5 when squeezing the nucleus.

- strength_one: This folder contains the initial structure and final structure with strength = 1.0 when squeezing the nucleus.
