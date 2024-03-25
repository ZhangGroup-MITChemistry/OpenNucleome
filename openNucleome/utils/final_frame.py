import numpy as np
import MDAnalysis as mda
import sys
import os

current_path = os.path.abspath(__file__)
father_path = os.path.abspath(os.path.dirname(current_path) + os.path.sep + '..')
config_path = os.path.join(father_path, 'parameters', 'HFF_100KB')

def final_frame(traj_data, type_final, output_config, n_chrom = 46, n_nuc = 300, n_spec = 1600, nuc_type = 5, lam_type = 8):
    """
    The function extracts the final frame of the previous simulation to restart a new simulation.

    Parameters
    ----------
    traj_data (string, required):
        The trajectory file from a simulation

    type_final (string, required):
        The file contains the type of each bead in the last frame

    output_config (string, required):
        The final frame of the trajectory

    n_chrom (int, required):
        The number of chromosomes (Default value = 46)

    n_nuc (int, required):
        The number of nucleolus (Default value = 300)

    n_spec (int, required):
        The number of speckles (Default value = 1600)

    nuc_type (int, required):
        The type index of nucleolus (Default value = 5)
        
    lam_type (int, required):
        The type index of lamina (Default value = 8)
    """

    traj = mda.coordinates.LAMMPS.DCDReader(traj_data)
    chr_info = np.loadtxt(os.path.join(config_path, 'chromatin_info.txt'), dtype=int)
    X, zero = 'X', 0
    type_final = np.loadtxt(type_final, dtype=int)
    coor = traj.trajectory[-1].positions
    n_all_beads = len(coor)
    n_chrom_beads = len(chr_info)

    file_write = open(output_config, 'w')
    for i in range(n_chrom_beads):
        if i == 0:
            file_write.write('CRYST1   26.400   26.400   26.400  90.00  90.00  90.00 P 1           1\n')
        file_write.write('ATOM%7s%3s%8s%4s%12.3f%8.3f%8.3f%6.2f%6.2f            \n'%(chr_info[i,0],chr_info[i,2],X,chr_info[i,1],coor[i,0],coor[i,1],coor[i,2],zero,zero))
    for k in range(1,n_all_beads-n_chrom_beads+1):
        idx_curr = n_chrom_beads+k
        mol_curr = n_chrom+k
        if k <= n_nuc:
            file_write.write('ATOM%7s%3s%8s%4s%12.3f%8.3f%8.3f%6.2f%6.2f            \n'%(idx_curr,nuc_type,X,mol_curr,coor[idx_curr-1,0],coor[idx_curr-1,1],coor[idx_curr-1,2],zero,zero))
        elif k <= n_nuc+n_spec:
            file_write.write('ATOM%7s%3s%8s%4s%12.3f%8.3f%8.3f%6.2f%6.2f            \n'%(idx_curr,type_final[idx_curr-1],X,mol_curr,coor[idx_curr-1,0],coor[idx_curr-1,1],coor[idx_curr-1,2],zero,zero))
        else:
            file_write.write('ATOM%7s%3s%8s%4s%12.3f%8.3f%8.3f%6.2f%6.2f            \n'%(idx_curr,lam_type,X,mol_curr,coor[idx_curr-1,0],coor[idx_curr-1,1],coor[idx_curr-1,2],zero,zero))
    file_write.write("END")
    file_write.close()
