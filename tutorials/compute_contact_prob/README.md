## contact probability computation code

- ComputeContactProb.ipynb: this notebook includes the contact probablity calculation of chromosome-lamina (DamID) and chromosome-speckle (TSA-Seq), where the cluster analysis of speckles are also involved. The number of clusters with respect to the time are plotted. In addition, We show how to visualize the comparison between experimental contacts and simulated contacts, including DamID, TSA-Seq, ideal, compt-compt, and inter-chromosomal contact probs.

- contact_calculation: execution file compiled from Fortran code "../../openNucleome/utils/chromosome_contact_calculation.f90"

- mol_info: save all the molecular information files necessary for computing contact probabilities

- reduced_traj.dcd: trajectory file example which contains 50 frames from "OpenNucleome/tutorials/HFF_100KB/SphericalNucleus/LangevinDynamics/reduced_traj.dcd"

- contact_prob: store all the chromosome-chromosome contact probabilities, the number of contacts, and the number of analyzed frames

- expt_constraints_HFF_100KB.txt: experimental chromosome-chromosome contact probabilities, the first 2489 values are ideal contacts, the next 6 are compt-compt contacts, and the remaining 231 are interchromosomal contacts
