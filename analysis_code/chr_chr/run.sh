#!/bin/bash

# All the parameters in this file can be changed

gfortran -o calculate_contact_prob calculate_contact_prob.f90

# Calculate the contact prob after getting the trajectory (DUMP_FILE.dcd); 501 represents the first frame chosen for analysis; -1 represents the last frame in the trajectory; ./contact_prob represents that the contact_prob results will be in this folder; 1 represents the index of replica, and counter.txt will log numbers of contact for every CV.
./calculate_contact_prob ../../examples/HFF_100KB/DUMP_FILE.dcd 501 -1 ./contact_prob/ 1 counter.txt

# The first 1 represents 1st iteration, and the second 1 represents the number of replica.
python adam.py 1 1 
