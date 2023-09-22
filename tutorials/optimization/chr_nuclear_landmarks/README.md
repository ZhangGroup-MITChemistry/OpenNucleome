### Optimization of potential between chromosomes and nuclear landmarks

mol_info, gLengthFile.txt, maternalIdxFile.txt, paternalIdxFile.txt: chromosome molecular information files

DamID-OE.txt, TSA-Seq-OE.txt: experimental DamID and TSA-Seq results

We used the Adam training to do the optimization, and the iter_num folder stores all the parameters used in the optimization.

To run the simulation:
```
python calculate_contact_prob.py 1 1

python update_alpha 1 1
```

calculate_contact_prob.py: calculates the contact probabilities between chr-specles and chr-lamina for each individual simulation. The first 1 represents the 1st iteration, and the second 1 represents the 1st replica.

update_alpha.py: update the potential parameters according to the comparison between simulations and experiments. The first 1 represents the 1st iteration, and the second 1 represents the number of replica.
