### Optimization of potential between chromosomes and chromosomes

mol_info, compartment_genome_100KB_diploid_HFF.txt, mol_index_cell_type-HFF.txt, tad_index_genome_100KB_diploid_cell_type-HFF.txt: chromosome molecular information files

expt_constraints_HFF_100KB.txt: contact probabilities from the experimental diploid Hi-C map

We used the Adam optimizer to do the optimization, and the iter_num folder stores all the parameters used in the optimization.

To run the optimization (All the parameters in this file can be changed):
```
gfortran -o calculate_contact_prob calculate_contact_prob.f90

./calculate_contact_prob ../../examples/HFF_100KB/DUMP_FILE.dcd 501 -1 ./contact_prob/ 1 counter.txt

python adam.py 1 1
```

calculate_contact_prob.f90: Calculate the contact prob after getting the trajectory (DUMP_FILE.dcd); 501 represents the first frame chosen for analysis; -1 represents the last frame in the trajectory; ./contact_prob/ represents that the contact_prob results will be in the current folder; counter.txt will log numbers of contact for every CV.

adam.py: The first 1 represents the 1st iteration, and the second 1 represents the number of replica.

update_alpha_chr_chr.py: Update the potential parameters according to the comparison between simulations and experiments

After the optimization, this folder will generate the new chr-chr potential.
