### Optimization of potential between chromosomes and nuclear landmarks

calculate_contact_prob.py: calculates the contact probabilities between chr-specles and chr-lamina for each individual simulation

update_alpha.py: update the potential parameters according to the comparison between simulations and experiments

gLengthFile.txt, maternalIdxFile.txt, paternalIdxFile.txt: chromosome molecular information files

We used the Adam training to do the optimization, and the iter_num folder stores all the parameters used in the optimization.

To run the simulation:
```
chmod 777 run.sh

./run.sh
```
