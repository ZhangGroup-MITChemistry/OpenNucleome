#!/bin/bash

python parameter.py
python simulation.py
python coor_transformation.py
python final_frame.py

rm HFF_3e6_every2000.dcd
