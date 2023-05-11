## HFF model with nuclear landmarks at 100KB resolution

This example is the HFF nucleus with nuclear landmarks (lamina, speckle, nucleolus) at 100KB resolution.

The model has 60642 chromosome beads, 300 nucleolus beads, 1600 speckle beads, and 8000 lamina beads.

We implemented the adaptive learning rate optimization (adam) to coarsely optimize the ideal (intra-chromosomal) potential, the compt-compt potential, and the inter-chromosomal potential with Hi-C data; chr-lamina and chr-speckle potential with sequencing data (DamID and TSA-Seq). All the optimization scripts are under the folder "/OpenNucleome/analysis_code".

To run the simulation:
```
python simulation.py
```

human.pdb: the initial structure

simulation.py: run the molecular dynamics simulation

coordinate_transformation.py: convert the position used in OpenMM to reduced unit

final_frame.py: extract the last frame of the simulation to restart the next iteration

## Files about the Chromosomes-Nuclear landmarks
chr_lam_spec.txt: Chr-lamina potential parameters

chr_spec_param.txt: Chr-speckle potential parameters

chr_nuc_param.txt: the probability for a chromatin bead to be in the nucleolus state as assigned by the SPIN algorithm

## Files about the Chromosomes-Chromosomes

ideal_param.txt: intra-chromosomal potential parameters

types_param.txt: compt-compt potential parameters

inter_param.txt: inter-chromosomal potential parameters
