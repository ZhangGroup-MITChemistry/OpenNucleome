## chromosome-chromosome optimization code

chr_chr_optimization.ipynb logs the process of optimizing the chromosome-chromosome interactions. The users can use the notebook OpenNucleome/tutorials/compute_contact_prob/compute_contact_prob.ipynb to get the ideal, compt-compt, and inter-chromosomal contact probablities, then, use this notebook to optimize the interactions with the Adam Optimizer. The loss function is $$L = \sum_i (\left< f_i\right> - f_i^\text{exp})^2$$

The pipeline is as follows:

(i) Starting with a set of values for $\alpha_\mathrm{ideal}(|i-j|)$, $\alpha_\mathrm{compt}(T_i, T_j)$, and $\alpha_\mathrm{inter}(I,J)$,  we performed 50 independent 3-million-step long MD simulations to obtain an ensemble of nuclear configurations. The 500K steps of each trajectory are discarded as equilibration. We collected the configurations at every 2000 simulation steps from the rest of the simulation trajectories to compute the ensemble averages defined on the left-hand side of Eq. S13. 

(ii) Check the convergence of the optimization by calculating the percentage of error defined as $\sum_i(\left< f_i\right>-f_i^\text{exp})/\sum_i f_i^\text{exp}$. The summation over $i$ includes all the average contact probabilities defined in Eq. S28.

(iii) If the error is less than a tolerance value $e_\text{tol}$, the optimization has converged, and we stop the simulations. Otherwise, we update the parameters, $\alpha$, using the Adam optimizer. With the new parameter values, we return to step one and restart the iteration.

- adam_chr_chr_param: save all the parameters used in the adam optimization for the chromosome-chromosome interactions

- potential: save all the potentials used in chr-chr interactions

- frame_10.dcd: trajectory file example which contains 10 frames

- contact_prob: store all the chromosome-chromosome contact probabilities, the number of contacts, and the number of analyzed frames

- expt_constraints_HFF_100KB.txt: experimental chromosome-chromosome contact probabilities, the first 2489 values are ideal contacts, the next 6 are compt-compt contacts, and the remaining 231 are interchromosomal contacts
