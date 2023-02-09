### Optimization of potential between chromosomes and nuclear landmarks

calculate_contact_prob.f90: calculates the contact probabilities between chr-chr for each individual simulation

adam_training.py, update_alpha_chr_chr.py: update the potential parameters according to the comparison between simulations and experiments

compartment_genome_100KB_diploid_HFF.txt, mol_index_cell_type-HFF.txt, tad_index_genome_100KB_diploid_cell_type-HFF.txt: chromosome molecular information files

We used the Adam training to do the optimization, and the iter_num folder stores all the parameters used in the optimization.

To run the simulation:
```
chmod 777 run.sh

./run.sh
```
