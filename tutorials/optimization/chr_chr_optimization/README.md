## chromosome-chromosome optimization code

chr_chr_optimization.ipynb logs the process of optimizing the chromosome-chromosome interactions. The users can use the notebook OpenNucleome/tutorials/compute_contact_prob/compute_contact_prob.ipynb to get the ideal, compt-compt, and inter-chromosomal contact probablities, then, use this notebook to optimize the interactions with the Adam Optimizer. The loss function is $$L = \sum_i (\left< f_i\right> - f_i^\text{exp})^2$$

- adam_chr_chr_param: save all the parameters used in the adam optimization for the chromosome-chromosome interactions

- potential: save all the potentials used in chr-chr interactions

- frame_10.dcd: trajectory file example which contains 10 frames

- contact_prob: store all the chromosome-chromosome contact probabilities, the number of contacts, and the number of analyzed frames

- expt_constraints_HFF_100KB.txt: experimental chromosome-chromosome contact probabilities, the first 2489 values are ideal contacts, the next 6 are compt-compt contacts, and the remaining 231 are interchromosomal contacts
