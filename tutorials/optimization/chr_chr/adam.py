import numpy as np 
import os 
import subprocess
import sys 
import copy

run_number                  = int(sys.argv[1]) 
N_replicas                  = int(sys.argv[2])

nAtom  = 60642
nIdeal = 2490-1
ncompt = 3
nChrom = 22
ninter = int(nChrom*(nChrom-1)/2)
ncv    = int(nIdeal + ncompt*(ncompt+1)/2 + ninter)

chop_ideals             = 1000                     

## Adam training parameters
m_dw                    = np.loadtxt('iter_num/%02d/mdw.txt'%(run_number-1))
v_dw                    = np.loadtxt('iter_num/%02d/vdw.txt'%(run_number-1))
m_db                    = np.loadtxt('iter_num/%02d/mdb.txt'%(run_number-1))
v_db                    = np.loadtxt('iter_num/%02d/vdb.txt'%(run_number-1))
beta1                   = 0.9
beta2                   = 0.999
epsilon                 = 1e-8
eta1                    = 0.01
eta2                    = 0.01
eta3                    = 0.01
t                       = int(np.loadtxt('iter_num/%02d/t.txt'%(run_number-1)))

start_cv                = 1     
end_cv                  = 2726

old_iter                = run_number-1

cvInd   = np.zeros((ncv, ), dtype=float)
irun    = 0
for replica in range(1,N_replicas+1,1):
    #If simulation didn't complete
    if os.path.exists("contact_prob/cvInd.txt_replica_%d"%replica):
        cvInd   += np.loadtxt("contact_prob/cvInd.txt_replica_%d"%replica)
        irun    += np.loadtxt("contact_prob/nframes.txt_replica_%d"%replica)
cvInd /= irun
np.savetxt('contact_prob/cvInd_iter%02d.txt'%(run_number), cvInd, fmt='%14.7e')

expt        = np.loadtxt("expt_constraints_HFF_100KB.txt")
expt        = cvInd[0]/expt[0]*expt

grad        = -cvInd + expt
## START TO DO THE ADAM TRAINING
## momentum beta 1
# *** weights *** #
m_dw        = beta1*m_dw + (1-beta1)*grad
# *** biases *** #
m_db        = beta1*m_db + (1-beta1)*grad
## rms beta 2
# *** weights *** #
v_dw        = beta2*v_dw + (1-beta2)*(grad**2)
# *** biases *** #
v_db        = beta2*v_db + (1-beta2)*grad

subprocess.call(["mkdir -p iter_num/%02d"%(run_number)],shell=True,stdout=subprocess.PIPE)
np.savetxt('iter_num/%02d/mdw.txt'%(run_number), m_dw.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/vdw.txt'%(run_number), v_dw.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/mdb.txt'%(run_number), m_db.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/vdb.txt'%(run_number), v_db.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/t.txt'%(run_number), np.array([t+1]).reshape((-1,1)), fmt='%d')

## bias correction
m_dw_corr   = m_dw/(1-beta1**t)
m_db_corr   = m_db/(1-beta1**t)
v_dw_corr   = v_dw/(1-beta2**t)
v_db_corr   = v_db/(1-beta2**t)

dalpha1     = m_dw_corr/(np.sqrt(v_dw_corr)+epsilon)
dalpha2     = m_db_corr/(np.sqrt(v_db_corr)+epsilon)

eigen_value_best        = 0 

np.savetxt("iter_num/%02d/dalpha.iter%02d.cutEig%d_noIdeal.txt"%(run_number-1,run_number,eigen_value_best),dalpha1.reshape((-1,1)),fmt='%15.12e')

subprocess.call(["python update_alpha_chr_chr.py ../../../ %d %d %.4f %.4f %.4f %d"%(run_number,eigen_value_best,eta1,eta2,eta3,chop_ideals)],shell=True,stdout=subprocess.PIPE)
