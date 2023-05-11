import numpy as np 
import os 
import subprocess
import sys 
import copy
import MDAnalysis as mda

run_number                  = int(sys.argv[1]) 
N_replicas                  = int(sys.argv[2])
update_chr_list             = np.array([i for i in range(1,23)])

update_chr_list -= 1

ncv                     = 30321             

start_cv                = 0
end_cv                  = 30321
first_frames            = 500

old_iter                = run_number-1

"""Info files"""
gLength             = np.loadtxt("mol_info/gLengthFile.txt",dtype=int)
maternalIdx         = np.loadtxt("mol_info/maternalIdxFile.txt",dtype=int)
paternalIdx         = np.loadtxt("mol_info/paternalIdxFile.txt",dtype=int)
damid_data_low_res  = np.loadtxt("DamID-OE.txt",usecols=[1])
tsa_data_low_res    = np.loadtxt("TSA-Seq-OE.txt",usecols=[1])
"""End Info Files"""

### Get number of frames data   ### 
n_frames        = np.zeros(N_replicas)
for i in range(N_replicas):
    traj_data   = mda.coordinates.LAMMPS.DCDReader("../../examples/HFF_100KB/DUMP_FILE.dcd")
    n_frames[i] = len(traj_data)-first_frames

cvInd       = np.zeros((ncv, ), dtype=float)
cvInd_tsa   = np.zeros((ncv, ), dtype=float)
irun    = 0
for replica in range(1,N_replicas+1,1):
    cvInd       += np.load("contact_prob/cvInd.txt_%d.npy"%replica)
    cvInd_tsa   += np.load("contact_prob/cvInd_tsa.txt_%d.npy"%replica)
    irun    += n_frames[replica-1]
cvInd       /= irun
cvInd_tsa   /= irun
np.savetxt('contact_prob/cvInd_iter%02d.txt'%(run_number), cvInd, fmt='%14.7e')
np.savetxt('contact_prob/cvInd_tsa_iter%02d.txt'%(run_number), cvInd_tsa, fmt='%14.7e')

### Experimental constraint portion 
# Speckles and Lamina
damid_data_low_res_haploid  = np.zeros(30321)
tsa_data_low_res_haploid    = np.zeros(30321)
for i in range(23):
    damid_data_low_res_haploid[gLength[i]:gLength[i+1]] = 0.5*(damid_data_low_res[maternalIdx[i][0]-1:maternalIdx[i][1]] + 
                                                       damid_data_low_res[paternalIdx[i][0]-1:paternalIdx[i][1]]
                                                       )
    tsa_data_low_res_haploid[gLength[i]:gLength[i+1]] = 0.5*(tsa_data_low_res[maternalIdx[i][0]-1:maternalIdx[i][1]] + 
                                                       tsa_data_low_res[paternalIdx[i][0]-1:paternalIdx[i][1]]
                                                       )

gw_lamina                   =   np.mean(cvInd)
gw_speckles                 =   np.mean(cvInd_tsa)

expt                        =   damid_data_low_res_haploid*gw_lamina
expt_tsa                    =   tsa_data_low_res_haploid*gw_speckles

m_dw_dam                    = np.loadtxt('iter_num/%02d/mdw_dam.txt'%(run_number-1))
v_dw_dam                    = np.loadtxt('iter_num/%02d/vdw_dam.txt'%(run_number-1))
m_db_dam                    = np.loadtxt('iter_num/%02d/mdb_dam.txt'%(run_number-1))
v_db_dam                    = np.loadtxt('iter_num/%02d/vdb_dam.txt'%(run_number-1))
beta1_dam                   = 0.9
beta2_dam                   = 0.999
epsilon_dam                 = 1e-8
eta_dam                     = 0.01
t_dam                       = int(np.loadtxt('iter_num/%02d/t_dam.txt'%(run_number-1)))

grad_dam        = -cvInd + expt
# START TO DO THE ADAM TRAINING
# momentum beta 1
# *** weights *** #
m_dw_dam        = beta1_dam*m_dw_dam + (1-beta1_dam)*grad_dam
# *** biases *** #
m_db_dam        = beta1_dam*m_db_dam + (1-beta1_dam)*grad_dam
# rms beta 2
# *** weights *** #
v_dw_dam        = beta2_dam*v_dw_dam + (1-beta2_dam)*(grad_dam**2)
# *** biases *** #
v_db_dam        = beta2_dam*v_db_dam + (1-beta2_dam)*grad_dam

subprocess.call(["mkdir -p iter_num/%02d"%run_number],shell=True,stdout=subprocess.PIPE)
np.savetxt('iter_num/%02d/mdw_dam.txt'%(run_number), m_dw_dam.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/vdw_dam.txt'%(run_number), v_dw_dam.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/mdb_dam.txt'%(run_number), m_db_dam.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/vdb_dam.txt'%(run_number), v_db_dam.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/t_dam.txt'%(run_number), np.array([t_dam+1]).reshape((-1,1)), fmt='%d')

## bias correction
m_dw_corr_dam   = m_dw_dam/(1-beta1_dam**t_dam)
m_db_corr_dam   = m_db_dam/(1-beta1_dam**t_dam)
v_dw_corr_dam   = v_dw_dam/(1-beta2_dam**t_dam)
v_db_corr_dam   = v_db_dam/(1-beta2_dam**t_dam)

dalpha1_dam     = m_dw_corr_dam/(np.sqrt(v_dw_corr_dam)+epsilon_dam)
dalpha2_dam     = m_db_corr_dam/(np.sqrt(v_db_corr_dam)+epsilon_dam)

damid = np.loadtxt("../../examples/HFF_100KB/chr_lam_param.txt")[:60642]

for i in update_chr_list:
    damid[maternalIdx[i][0]-1:maternalIdx[i][1]] -= eta_dam*dalpha1_dam[gLength[i]:gLength[i+1]]
    damid[paternalIdx[i][0]-1:paternalIdx[i][1]] -= eta_dam*dalpha1_dam[gLength[i]:gLength[i+1]]

m_dw_tsa                    = np.loadtxt('iter_num/%02d/mdw_tsa.txt'%(run_number-1))
v_dw_tsa                    = np.loadtxt('iter_num/%02d/vdw_tsa.txt'%(run_number-1))
m_db_tsa                    = np.loadtxt('iter_num/%02d/mdb_tsa.txt'%(run_number-1))
v_db_tsa                    = np.loadtxt('iter_num/%02d/vdb_tsa.txt'%(run_number-1))
beta1_tsa                   = 0.9
beta2_tsa                   = 0.999
epsilon_tsa                 = 1e-8
eta_tsa                     = 0.01
t_tsa                       = int(np.loadtxt('iter_num/%02d/t_tsa.txt'%(run_number-1)))

grad_tsa        = -cvInd_tsa + expt_tsa
# START TO DO THE ADAM TRAINING
# momentum beta 1
# *** weights *** #
m_dw_tsa        = beta1_tsa*m_dw_tsa + (1-beta1_tsa)*grad_tsa
# *** biases *** #
m_db_tsa        = beta1_tsa*m_db_tsa + (1-beta1_tsa)*grad_tsa
# rms beta 2
# *** weights *** #
v_dw_tsa        = beta2_tsa*v_dw_tsa + (1-beta2_tsa)*(grad_tsa**2)
# *** biases *** #
v_db_tsa        = beta2_tsa*v_db_tsa + (1-beta2_tsa)*grad_tsa

np.savetxt('iter_num/%02d/mdw_tsa.txt'%(run_number), m_dw_tsa.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/vdw_tsa.txt'%(run_number), v_dw_tsa.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/mdb_tsa.txt'%(run_number), m_db_tsa.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/vdb_tsa.txt'%(run_number), v_db_tsa.reshape((-1,1)), fmt='%15.12e')
np.savetxt('iter_num/%02d/t_tsa.txt'%(run_number), np.array([t_tsa+1]).reshape((-1,1)), fmt='%d')

## bias correction
m_dw_corr_tsa   = m_dw_tsa/(1-beta1_tsa**t_tsa)
m_db_corr_tsa   = m_db_tsa/(1-beta1_tsa**t_tsa)
v_dw_corr_tsa   = v_dw_tsa/(1-beta2_tsa**t_tsa)
v_db_corr_tsa   = v_db_tsa/(1-beta2_tsa**t_tsa)

dalpha1_tsa     = m_dw_corr_tsa/(np.sqrt(v_dw_corr_tsa)+epsilon_tsa)
dalpha2_tsa     = m_db_corr_tsa/(np.sqrt(v_db_corr_tsa)+epsilon_tsa)

tsaseq = np.loadtxt("../../examples/HFF_100KB/chr_spec_param.txt")[:60642]

for i in update_chr_list:
    tsaseq[maternalIdx[i][0]-1:maternalIdx[i][1]] -= eta_tsa*dalpha1_tsa[gLength[i]:gLength[i+1]]
    tsaseq[paternalIdx[i][0]-1:paternalIdx[i][1]] -= eta_tsa*dalpha1_tsa[gLength[i]:gLength[i+1]]


#Added portion to overide the parameters to 0.0 if no expt signal on segment
zero_signal_damid   = (damid_data_low_res[:]    == 0.0)
zero_signal_tsa     = (tsa_data_low_res[:]      == 0.0)
damid[zero_signal_damid]  = 0.0 
tsaseq[zero_signal_tsa]    = 0.0

np.savetxt("new_potential/chr_spec_param.txt",np.append(tsaseq,[0]*9900).reshape((-1,1)),fmt='%.6f')
np.savetxt("new_potential/chr_lam_param.txt",np.append(damid,[0]*9900).reshape((-1,1)),fmt='%.6f')
