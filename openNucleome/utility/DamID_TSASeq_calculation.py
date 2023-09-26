import numpy as np
import MDAnalysis as mda
import sklearn.cluster as sk
import matplotlib.pyplot as plt
import copy
import subprocess
from joblib import Parallel, delayed
import scipy.stats
import sys
import os

def process_traj_win(win_num):

    traj_data                   = mda.coordinates.LAMMPS.DCDReader("../HFF_100KB/SphericalNucleus/DUMP_FILE.dcd")
    N_frame                     = len(traj_data)-first_frame
    exp_tsa_seq                 = np.zeros(N_chr_beads)
    exp_damid                   = np.zeros(N_chr_beads)
    exp_tsa_seq_all_frames      = np.zeros((N_frame,N_chr_beads))
    exp_damid_all_frames        = np.zeros((N_frame,N_chr_beads))
    
    for frame_number in range(first_frame,len(traj_data),1):
        chr_I_data          = traj_data.trajectory[frame_number].positions[:N_chr_beads]
        speckles_data       = traj_data.trajectory[frame_number].positions[(N_chr_beads+N_nucleolus_particles):(N_chr_beads+N_nucleolus_particles+N_speckles_particles)]
        lamina_data         = traj_data.trajectory[frame_number].positions[(N_chr_beads+N_nucleolus_particles+N_speckles_particles):]

        #Following code identifies the speckles clusters in the system and calculates and saves center of mass positions
        """Code Snippet from DBSCAN Python Page"""
        """https://scikit-learn.org/stable/auto_examples/cluster/plot_dbscan.html#sphx-glr-auto-examples-cluster-plot-dbscan-py"""
        db = sk.DBSCAN(eps=1.0).fit(speckles_data)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)

        cluster_com_master = np.zeros((n_clusters_,3))
        radius_cluster_master = np.zeros(n_clusters_)
        i = 0
        while(i<=max(labels)):
            #Go cluster by cluster
            r_points_cluster            = copy.deepcopy(speckles_data[labels==i])
            cluster_com                 = np.mean(r_points_cluster,axis=0)
            cluster_com_master[i]       = copy.deepcopy(cluster_com)
            radius_cluster_master[i]    = (np.sum((r_points_cluster-cluster_com)**2,axis=None)/len(r_points_cluster))**0.5
            i +=1

        for i in range(len(chr_I_data)):
            ###Speckles
            distances_from_speckles = np.sum((cluster_com_master-chr_I_data[i])**2,axis=1)**0.5
            distances_from_speckles -= radius_cluster_master

            #Do this only for speckles
            exp_tsa_seq_all_frames[frame_number-first_frame,i] = np.sum(0.5*(1.0+np.tanh(sigma_tanh*(rc_tanh-distances_from_speckles))),axis=None)
            #exp_tsa_seq will be normalized by total (itself) while getting enrichment and so it does not matter if I divide by number of speckles
            #But have to divide exp_tsa_seq_all_frames by number of speckles so that the values can be interpretted as probabilities
            #when calculating correlations
            exp_tsa_seq_all_frames[frame_number-first_frame,i] /= float(n_clusters_)
            exp_tsa_seq[i] += copy.deepcopy(exp_tsa_seq_all_frames[frame_number-first_frame,i])        

            ###Lamina
            distances_from_lamina = np.sum((lamina_data-chr_I_data[i])**2,axis=1)**0.5

            #Do this only for lamina
            exp_damid_all_frames[frame_number-first_frame,i] = np.sum(0.5*(1.0+np.tanh(sigma_tanh*(rc_tanh-distances_from_lamina))),axis=None)
            exp_damid[i] += copy.deepcopy(exp_damid_all_frames[frame_number-first_frame,i])


        if frame_number%20==0:
            print("Window %d Frame %d done"%(win_num,frame_number))
            print(n_clusters_,n_noise_) #This will report speckles

    np.save("contact_prob/exp_tsa_seq_iter%d_win%d"%(iter_num,win_num),exp_tsa_seq)
    np.save("contact_prob/exp_damid_iter%d_win%d"%(iter_num,win_num),exp_damid)
    np.save("contact_prob/exp_tsa_seq_all_frames_iter%d_win%d"%(iter_num,win_num),exp_tsa_seq_all_frames)
    np.save("contact_prob/exp_damid_all_frames_iter%d_win%d"%(iter_num,win_num),exp_damid_all_frames)
    
    damid_all_frames_haploid = np.zeros((N_frame,int(N_chr_beads/2)))
    tsaseq_all_frames_haploid = np.zeros((N_frame,int(N_chr_beads/2)))
    for i in range(23):
        damid_all_frames_haploid[:,gLength[i]:gLength[i+1]] = 0.5*(exp_damid_all_frames[:,maternalIdx[i][0]-1:maternalIdx[i][1]]+exp_damid_all_frames[:,paternalIdx[i][0]-1:paternalIdx[i][1]])
        tsaseq_all_frames_haploid[:,gLength[i]:gLength[i+1]] = 0.5*(exp_tsa_seq_all_frames[:,maternalIdx[i][0]-1:maternalIdx[i][1]]+exp_tsa_seq_all_frames[:,paternalIdx[i][0]-1:paternalIdx[i][1]])

    np.save("contact_prob/cvInd.txt_%d"%win_num,np.sum(damid_all_frames_haploid,axis=0))
    np.save("contact_prob/cvInd_tsa.txt_%d"%win_num,np.sum(tsaseq_all_frames_haploid,axis=0))

