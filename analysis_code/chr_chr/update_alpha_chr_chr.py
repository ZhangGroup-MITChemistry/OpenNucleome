from numpy import *
from subprocess import *
import subprocess
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import copy
import sys
import os

nAtom       = 2490-1
nIdeal      = 2489
ncompt      = 3
nChrom      = 22
inter_ncv   = int(nChrom*(nChrom-1)/2)
ncv_interTAD_intraChrom = int(ncompt*(ncompt+1)/2)

r_cut =  0.702
sigma =  6.218

def update_alpha(iterId,egcut):
    cutoff  = egcut

    old_iter    = iterId-1
    afile_old   = '%s/%02d/alpha_iter%02d.txt'%(write_potential_path,old_iter,old_iter)
    afile_new   = '%s/%02d/alpha_iter%02d.txt'%(write_potential_path,iterId,iterId)

    alpha_old   = loadtxt(afile_old)
    alpha       = copy.copy(alpha_old)
    
    separation_cut = 1000
    
    dalpha                      = loadtxt('../../../../analysis_code/chr_chr/iter_num/%02d/dalpha.iter%02d.cutEig%d_noIdeal.txt'%(old_iter,iterId,cutoff))
    
    alpha[1:separation_cut]                         -= dalpha[1:separation_cut]*eta1
    alpha[nIdeal:nIdeal+ncv_interTAD_intraChrom]    -= dalpha[nIdeal:nIdeal+ncv_interTAD_intraChrom]*eta2
    alpha[nIdeal+ncv_interTAD_intraChrom:]          -= dalpha[nIdeal+ncv_interTAD_intraChrom:]*eta3

    savetxt(afile_new, alpha, fmt='%15.12e')

    ## Update the intra-chrom potential
    fo = open('ideal_chromosome_iter%02d.txt'%(iterId), 'w')
    fo.write('%8d '%(0))
    fo.write('%12.6f '%(0))
    fo.write('\n')
    for ib in range(nIdeal):
        fo.write('%8d '%(ib+1))
        if ib<chop_ideals:
            fo.write('%12.6f '%(alpha[ib]))
        else:
            fo.write('%12.6f '%(0.0))
        fo.write('\n')
    fo.close()

    ## Update the compt-compt potential
    eij_intraTAD_intraChrom = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    eij_interTAD_intraChrom = alpha[nIdeal:nIdeal+ncv_interTAD_intraChrom]
    eij_interTAD_interChrom = copy.deepcopy(eij_interTAD_intraChrom)
    fo = open('eij_compartment_uniform_iter%02d.txt'%(iterId), 'w')
    icv = 0
    for ic in range(ncompt+1):
        for jc in range(ic, ncompt+1):
            if ic < ncompt and jc < ncompt:
                if ic == ncompt-1 or jc == ncompt-1:
                    fo.write('pair_coeff      %4d %4d  tanhlr/cut/domainab    %12.6f %12.6f %12.6f ${rctanh}  ${sigmatanh}\n'\
                        %(ic+1, jc+1, eij_intraTAD_intraChrom[icv],\
                            eij_interTAD_intraChrom[icv], eij_interTAD_interChrom[icv]))
                else:
                    fo.write('pair_coeff      %4d %4d  tanhlr/cut/domainab    %12.6f %12.6f %12.6f ${rctanh}  ${sigmatanh}\n'\
                        %(ic+1,jc+1,eij_intraTAD_intraChrom[icv],eij_interTAD_intraChrom[icv],eij_interTAD_interChrom[icv]))
                icv += 1
            else:
                fo.write('pair_coeff      %4d %4d  tanhlr/cut/domainab    %12.6f %12.6f %12.6f ${rctanh}  ${sigmatanh}\n'\
                    %(ic+1,jc+1,0.0,0.0,0.0))
    fo.close()

    ## Update the inter-chrom potential
    fo = open('inter_chromosome_iter%02d.txt'%(iterId), 'w')
    icv = 2495
    for ic in range(1,22,1):
        for jc in range(ic+1,23,1):
            fo.write('%d %d %12.6f\n'%(ic,jc,alpha[icv]))
            icv += 1

if __name__ == "__main__":

    run_path    = sys.argv[1]
    iterId      = int(sys.argv[2])
    eigcutoff   = int(sys.argv[3])
    write_potential_path    = "%s/examples/HFF_100KB/potential"%run_path
    os.chdir("%s/%02d"%(write_potential_path,iterId))

    eta1        = float(sys.argv[4])
    eta2        = float(sys.argv[5])
    eta3        = float(sys.argv[6])
    chop_ideals = int(sys.argv[7])

    update_alpha(iterId,eigcutoff)
