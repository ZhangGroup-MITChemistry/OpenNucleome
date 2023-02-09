#!/bin/bash

gfortran -o calculate_contact_prob calculate_contact_prob.f90
calculate_contact_prob ../../examples/HFF_100KB/DUMP_FILE.dcd 500 -1 ./ 1 counter.txt

python adam_training.py 1 1 
