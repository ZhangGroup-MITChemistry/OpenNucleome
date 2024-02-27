## chromosome-nuclear landmarks optimization code

chr_NL_optimization.ipynb logs the process of optimizing the chromosome-nuclear landmarks interactions. The users can use the notebook OpenNucleome/tutorials/compute_contact_prob/compute_contact_prob.ipynb to get the DamID and TSA-Seq, then, use this notebook to optimize the interactions with the Adam Optimizer. The loss function is similar to the chr-chr optimization $$L = \sum_i (\left< f_i\right> - f_i^\text{exp})^2$$

The pipline is the same as that in the folder "../chr_chr_optimization"

- DamID-OE.txt: experimental DamID sequencing data

- TSA-Seq-OE.txt: experimental TSASeq sequencing data

- mol_info: save all the molecular information files necessary for computing contact probabilities

- adam_chr_NL_param: save all the parameters used in the adam optimization for the chromosome-nuclear landmarks interactions

- potential: save all the potentials used in the chr-NL interactions

- frame_10.dcd: trajectory file example which contains 10 frames
