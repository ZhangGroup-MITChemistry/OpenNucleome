## HFF model with nuclear landmarks at 100KB resolution

This example is the HFF nucleus with nuclear landmarks (lamina, speckle, nucleolus) at 100KB resolution.

The model has 60642 chromosome beads, 300 nucleolus beads, 600 speckle beads, and 8000 lamina beads.

We implemented the adaptive learning rate optimization (adam) to coarsely optimize the ideal (intra) potential, the compt-compt potential, and the inter potential with Hi-C data; chr-lamina and chr-speckle potential with sequencing data (DamID and TSA-Seq). All the optimization scripts are under the folder "/OpenNucleome/analysis_code".

To run the simulation:
```
chmod 777 run.sh

./run.sh
```

human.pdb: the initial structure

parameter.py: put all the potential parameter files into this folder

simulation.py: run the molecular dynamics simulation

openChrModel.py: the source code

coor_transformation.py: convert the position used in OpenMM to reduced unit

final_frame.py: extract the last frame of the simulation to restart the next iteration

/potential: folder that contains all the potential parameter files

## Files about the chr-nuclear landmarks
DamID_8900.txt: chr-lamina potential parameters

TSA_8900.txt: chr-speckle potential parameters

nuc_rescaling.txt: the probability for a chromatin bead to be in the nucleolus state as assigned by the SPIN algorithm

## Files about the chr-chr

ideal_chromosome.txt: intra-chromosome potential parameters

eij_compartment_uniform.txt: compt-compt potential parameters

inter_chromosome.txt: inter-chromosome potential parameters
