import subprocess
import numpy as np
import sys

iter_num                  = int(sys.argv[1])

subprocess.call(["cp -p potential/%02d/ideal_chromosome_iter%02d.txt ideal_chromosome.txt"%(iter_num,iter_num)],shell=True,stdout=subprocess.PIPE)
subprocess.call(["cp -p potential/%02d/eij_compartment_uniform_iter%02d.txt eij_compartment_uniform.txt"%(iter_num,iter_num)],shell=True,stdout=subprocess.PIPE)
subprocess.call(["cp -p potential/%02d/inter_chromosome_iter%02d.txt inter_chromosome.txt"%(iter_num,iter_num)],shell=True,stdout=subprocess.PIPE)
subprocess.call(["cp -p potential/%02d/DamID_8900.txt ./"%(iter_num)],shell=True,stdout=subprocess.PIPE)
subprocess.call(["cp -p potential/%02d/TSA_8900.txt ./"%(iter_num)],shell=True,stdout=subprocess.PIPE)
