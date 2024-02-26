## Choose a set of proper initial configurations to run the simulations

Here, as an example, we only choose 20 initial configurations out of 200 relaxed configurations. 

- ConfigsSelection.ipynb: call and outputs the indexes of 20 initial configurations

- config_pool: save all the initial relaxed configurations

- DamID_haploid.txt: experimental haploid DamID sequencing data. ../analysis_code/mol_info/DAM-ID_HFF_100KB.txt is averaged over the homologous chromosomes,

- DamID_simulation.txt: simulated haploid DamID sequencing data. Each row with a specific row index represent the DamID result of the corresponding initial configuration with the same index.

- expt_constraints_HFF_100KB.txt: experimental chromosome-chromosome contact probabilities, the first 2489 values are ideal contacts, the next 6 are compt-compt contacts, and the remaining 231 are interchromosomal contacts

- inter_contact_simulation.txt: simulated haploid inter-chromosomal contact probabilities. Each row with a specific row index represent the DamID result of the corresponding initial configuration with the same index.

- choose_20.txt: the indexes of 20 proper initial configurations
