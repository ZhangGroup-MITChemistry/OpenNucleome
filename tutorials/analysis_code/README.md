## Analysis code

- compute_contact_prob.ipynb: call and computes the contact probabilities (DamID, TSASeq, and chr-chr)

- contact_calculation: execution file compiled from Fortran code "../../openNucleome/utility/chromosome_contact_calculation.f90"

- mol_info: store all the molecular information files necessary for computing contact probabilities

- frame_10.dcd: trajectory file example which contains 10 frames

- contact_prob: store all the chromosome-chromosome contact probabilities, the number of contacts, and the number of analyzed frames

- expt_constraints_HFF_100KB.txt: experimental chromosome-chromosome contact probabilities, the first 2489 values are ideal contacts, the next 6 are compt-compt contacts, and the remaining 231 are interchromosomal contacts
