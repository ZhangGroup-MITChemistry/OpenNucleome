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

def DamID_TSASeq_calculation(traj_data, gLength, maternalIdx, paternalIdx, start_frame, end_frame):

    """
    The function computes the DamID and TSA-Seq signals. The returns of this function are DamID, TSA-Seq, and number of speckle clusters.

    Parameters
    ----------
    traj_data (string, required) :
        The trajectory file

    gLength (string, required) :
        The difference between each two neighboring values in the file represents the length of the chromosome

    maternalIdx (string, required) :
        Index of each maternal chromosome

    paternalIdx (string, required) :
        Index of each paternal chromosome

    start_frame (int, required) :
        The starting analyzed frame

    end_frame (int, required) :
        The ending analyzed frame

    """
    # The next two parameters are used when computing the DamID and TSASeq signals
    sigma_tanh          =  4
    rc_tanh             =  0.75

    N_chr_beads         = 60642            #Number of chromosome particles

    first_frame                 = start_frame      #The analysis starts from the first frame
    N_nucleolus_particles       = 300              #The number of nucleolus particles
    N_speckles_particles        = 1600             #The number of speckle particles
    N_lamina_particles          = 8000             #The number of lamina particles
    radius_nucleus              = 13.0             #The radius of cell nucleus, 13.0 (LJ unit) = 5.0 µm

    # Info files
    gLength             = np.loadtxt(gLength,dtype=int)      #The length of each chromosome
    maternalIdx         = np.loadtxt(maternalIdx,dtype=int)  #Index of each maternal chromosome
    paternalIdx         = np.loadtxt(paternalIdx,dtype=int)  #Index of each paternal chromosome

    traj_data                   = mda.coordinates.LAMMPS.DCDReader(traj_data)
    N_frame                     = end_frame-first_frame
    exp_tsa_seq                 = np.zeros(N_chr_beads)
    exp_damid                   = np.zeros(N_chr_beads)
    exp_tsa_seq_all_frames      = np.zeros((N_frame,N_chr_beads))
    exp_damid_all_frames        = np.zeros((N_frame,N_chr_beads))
    N_speckle                   = []
    
    for frame_number in range(first_frame,end_frame,1):
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

        N_speckle.append(n_clusters_)

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
            #Speckles
            distances_from_speckles = np.sum((cluster_com_master-chr_I_data[i])**2,axis=1)**0.5
            distances_from_speckles -= radius_cluster_master

            #Do this only for speckles
            exp_tsa_seq_all_frames[frame_number-first_frame,i] = np.sum(0.5*(1.0+np.tanh(sigma_tanh*(rc_tanh-distances_from_speckles))),axis=None)

            exp_tsa_seq_all_frames[frame_number-first_frame,i] /= float(n_clusters_)
            exp_tsa_seq[i] += copy.deepcopy(exp_tsa_seq_all_frames[frame_number-first_frame,i])        

            #Lamina
            distances_from_lamina = np.sum((lamina_data-chr_I_data[i])**2,axis=1)**0.5

            #Do this only for lamina
            exp_damid_all_frames[frame_number-first_frame,i] = np.sum(0.5*(1.0+np.tanh(sigma_tanh*(rc_tanh-distances_from_lamina))),axis=None)
            exp_damid[i] += copy.deepcopy(exp_damid_all_frames[frame_number-first_frame,i])
    
    damid_all_frames_haploid = np.zeros((N_frame,int(N_chr_beads/2)))
    tsaseq_all_frames_haploid = np.zeros((N_frame,int(N_chr_beads/2)))
    for i in range(23):
        damid_all_frames_haploid[:,gLength[i]:gLength[i+1]] = 0.5*(exp_damid_all_frames[:,maternalIdx[i][0]-1:maternalIdx[i][1]]+exp_damid_all_frames[:,paternalIdx[i][0]-1:paternalIdx[i][1]])
        tsaseq_all_frames_haploid[:,gLength[i]:gLength[i+1]] = 0.5*(exp_tsa_seq_all_frames[:,maternalIdx[i][0]-1:maternalIdx[i][1]]+exp_tsa_seq_all_frames[:,paternalIdx[i][0]-1:paternalIdx[i][1]])

    return (np.mean(damid_all_frames_haploid,axis=0), np.mean(tsaseq_all_frames_haploid,axis=0), N_speckle)
