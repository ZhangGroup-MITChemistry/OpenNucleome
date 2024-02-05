## Analysis code

- compute_contact_prob.ipynb: call and computes the contact probabilities (DamID, TSASeq, and chr-chr)

- chr_chr_optimization.ipynb: call and optimizes the chromosome-chromosome interactions (ideal, compt-compt, interchr)

- chr_NL_optimization.ipynb: call and optimizes the chromosome-nuclear landmarks interactions (chr-spe, chr-lam)

- DamID-OE.txt: experimental DamID sequencing data

- TSA-Seq-OE.txt: experimental TSASeq sequencing data

- contact_calculation: execution file compiled from Fortran code "../../openNucleome/utility/chromosome_contact_calculation.f90"

- mol_info: save all the molecular information files necessary for computing contact probabilities

- adam_chr_chr_param: save all the parameters used in the adam optimization for the chromosome-chromosome interactions

- adam_chr_NL_param: save all the parameters used in the adam optimization for the chromosome-nuclear landmarks interactions

- potential: save all the potentials used in the simulations of different iterations

- frame_10.dcd: trajectory file example which contains 10 frames

- contact_prob: store all the chromosome-chromosome contact probabilities, the number of contacts, and the number of analyzed frames

- expt_constraints_HFF_100KB.txt: experimental chromosome-chromosome contact probabilities, the first 2489 values are ideal contacts, the next 6 are compt-compt contacts, and the remaining 231 are interchromosomal contacts
