import numpy as np
import MDAnalysis as mda

traj = mda.coordinates.LAMMPS.DCDReader("DUMP_FILE.dcd")

index  = np.loadtxt('bead_info.txt', dtype='int') # 1st column: bead index; 2nd column: chain index; 3rd column: type index

X = np.array(['X']*len(index)).reshape((-1,1))

last = np.array([0.00]*len(index))

coor  = traj.trajectory[-1].positions
n_chrom = 46
n_chrom_beads = len(index)
spe_type = 6
nuc_type = 5
lam_type = 7

file_write = open('human_final.pdb', 'w')
for i in range(len(index)):
    if i == 0:
        file_write.write('CRYST1   26.400   26.400   26.400  90.00  90.00  90.00 P 1           1\n')
    file_write.write('ATOM%7s%3s%8s%4s%12.3f%8.3f%8.3f%6.2f%6.2f            \n'%(index[i,0],index[i,2],X[i,0],index[i,1],coor[i,0],coor[i,1],coor[i,2],last[i],last[i]))
for k in range(1,len(coor)-len(index)+1):
    if k<=300:
        file_write.write('ATOM%7s%3s%8s%4s%12.3f%8.3f%8.3f%6.2f%6.2f            \n'%(n_chrom_beads+k,nuc_type,X[i,0],n_chrom+k,coor[len(index)+k-1,0],coor[len(index)+k-1,1],coor[len(index)+k-1,2],last[i],last[i]))
    elif k<=900:
        file_write.write('ATOM%7s%3s%8s%4s%12.3f%8.3f%8.3f%6.2f%6.2f            \n'%(n_chrom_beads+k,spe_type,X[i,0],n_chrom+k,coor[len(index)+k-1,0],coor[len(index)+k-1,1],coor[len(index)+k-1,2],last[i],last[i]))
    else:
        file_write.write('ATOM%7s%3s%8s%4s%12.3f%8.3f%8.3f%6.2f%6.2f            \n'%(n_chrom_beads+k,lam_type,X[i,0],n_chrom+k,coor[len(index)+k-1,0],coor[len(index)+k-1,1],coor[len(index)+k-1,2],last[i],last[i]))
file_write.write("END")
file_write.close()
