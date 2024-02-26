Program main

!
! part of the code for reading DCD file is copied from
! http://www.ks.uiuc.edu/Research/namd/wiki/index.cgi?ReadingDCDinFortran
!

    ! =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~ !
    ! variable definitions
    ! =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~ !

    implicit none
    character*500      :: dcdFileName, instr, molFileName, comptFileName, domainFileName, counterFileName

    ! variables to read in dcd file
    double precision  :: d
    real              :: t,dummyr
    real,allocatable  :: x(:,:),y(:,:),z(:,:)
    integer           :: nset,natom,dummyi,i,j,nframes
    character*4       :: dummyc

    ! parameters for contact map
    double precision, parameter :: r_cut=0.540

    integer           :: nIdeal, nCV, ndomain, intra, maxLoop, nAtomSum,    &
                        icm, jcm, tadI, tadJ, chrSta, chrEnd, natom_chr, imol, jmol, nChrom, ceil1, ceil2
    integer,allocatable:: boundary(:,:), domainIndx(:,:)

    ! optimization variables
    double precision :: g_a, g_b, g_norm, icount
    double precision :: dx, dy, dz, rij, r2, pij, r_cut4, p_cut
    double precision :: nideal_10kb_boundary, nideal_50kb_boundary, i50kb, i10kb,&
                        nideal_10kb_start, nideal_50kb_start
    integer,allocatable:: active(:), ideal_index(:), domain(:), mol_index(:), &
                          ctcfList(:,:), counter(:), counter1(:), compt(:), &
                          intra_TAD_index_map(:,:), cv_index_map(:,:), &
                          interTAD_intra_chrom_index_map(:,:), &
                          interTAD_inter_chrom_index_map(:,:), &
                          inter_chrom_index_map(:,:)
    double precision, allocatable :: cvRaw(:,:),      &   ! cmap from HiC
                                     cvRaw1(:,:),     &
                                     cvCorr(:,:),     &
                                     cvCorr1(:,:),    &
                                     cvRaw_avg(:),    &
                                     cvRaw1_avg(:)

    ! general variables
    integer           :: ib, jb, kb, ifr, iCS,jCS, gpos, iCV, jmi, lb, ub, scaleId, id, jpd, kpd
    integer           :: reclength, polynomial_degree, pdPlus1, ncompt, ipd, cvInd
    double precision  :: polyscale(5)
    
    ! blas variable
    CHARACTER*1       :: TRANSA, TRANSB
    INTEGER           :: M, N, K, LDA, LDB, LDC, first_frame, last_frame
    REAL*8            :: ALPHA, BETA

    ! for files
    CHARACTER*4       :: replica_number
    CHARACTER*500      :: write_file_path,write_file_1,write_file_2,write_file_3,write_file_4,write_file_5,write_file_6
    
    ! =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~ !
    ! function starts here
    ! =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~ !
    
    ! getarg retrieves command line arguments, and sets first to "filename"
    call getarg(1,dcdFileName)
    call getarg(2,instr)
    read (instr, '(i10)') first_frame
    call getarg(3,instr)
    read (instr, '(i10)') last_frame
    call getarg(4,write_file_path)
    !
    ! read dcd file
    open(10,file=trim(dcdFileName),status='old',form='unformatted')
    read(10) dummyc, nframes, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
    read(10) dummyi, dummyr
    read(10) natom

    allocate(x(natom,nframes))
    allocate(y(natom,nframes))
    allocate(z(natom,nframes))

    do i = 1, nframes
       read(10) (d, j=1, 6)
       read(10) (x(j,i),j=1,natom)
       read(10) (y(j,i),j=1,natom)
       read(10) (z(j,i),j=1,natom)
    end do

    print*
    print *, 'The dcd contains ',nframes,' frames and ',natom,' atoms'
    print*
    close(10)

    natom     = 57520 !exclude the sex chromosome
    natom_chr = 57520 !exclude the sex chromosome
    nIdeal = 2489
    ncompt = 3
    nChrom = 22 ! 22 is 23-1, exclude the sex chromosome in the mouse cell nuclei
    
    nCV = nIdeal+ncompt*(ncompt+1)/2  + nChrom*(nChrom-1)/2 

    allocate(interTAD_intra_chrom_index_map(ncompt,ncompt))
    cvInd = nIdeal
    do ib = 1, ncompt
        do jb = ib, ncompt
            cvInd = cvInd +1
            interTAD_intra_chrom_index_map(ib,jb) = cvInd
            interTAD_intra_chrom_index_map(jb,ib) = cvInd
        enddo
    enddo

    allocate(inter_chrom_index_map(nChrom,nChrom))
    do ib = 1, nChrom-1
        do jb = ib+1, nChrom
            cvInd = cvInd +1
            inter_chrom_index_map(ib,jb) = cvInd
            inter_chrom_index_map(jb,ib) = cvInd
        enddo
    enddo

    !
    ! read chromatin state
    allocate(compt(natom_chr))
    !open(unit=10, file=trim(comptFileName))
    open(unit=10, & 
    file='mol_info/compartment_genome_100KB_diploid_HFF.txt')
    do ib = 1, natom_chr
        read(10,*) compt(ib)
    enddo
    close(10)

    !
    ! read molecule index
    allocate(mol_index(natom))
    !open(unit=10, file=trim(molFileName))
    open(unit=10, &
    file='mol_info/mol_index_cell_type-HFF.txt')
    do ib = 1, natom
        read(10,*) mol_index(ib)
    enddo
    close(10)

    !
    ! read domain list
    allocate(domain(natom_chr))
    !open(unit=10, file=trim(domainFileName))
    open(unit=10, &
    file='mol_info/tad_index_genome_100KB_diploid_cell_type-HFF.txt')
    do ib = 1, natom_chr
        read(10,*) gpos, domain(ib)
    enddo
    close(10)
    
    r_cut4 = r_cut**4
    if (last_frame .ne. -1) then
        nframes = last_frame
    endif

    allocate(cvRaw(nCV, nframes))
    allocate(cvRaw_avg(nCV))
    allocate(counter(nCV))

    do ifr = first_frame, nframes
        !
        ! initialize the contact map
        cvraw(:,ifr) = 0.0
        counter      = 0.0

        do ib = 1, natom -1
            do jb = ib+1, natom

                tadI = domain(ib)
                tadJ = domain(jb)

                if ( tadI .eq. tadJ ) then
                    cycle
                endif

                icm = compt(ib)
                jcm = compt(jb)

                if ( (icm .eq. ncompt+1) .or. (jcm .eq. ncompt+1) ) then
                   cycle
                endif

                dx = x(ib,ifr)-x(jb,ifr)
                dy = y(ib,ifr)-y(jb,ifr)
                dz = z(ib,ifr)-z(jb,ifr)
                r2 = dx**2 + dy**2 + dz**2
                rij = sqrt(r2)

                if (rij <= r_cut) then
                    pij = 1.0
                else
                    pij = r_cut4 / r2 / r2
                endif

                imol = mol_index(ib)
                jmol = mol_index(jb)


                ceil1 = ceiling(imol/2.0)
                ceil2 = ceiling(jmol/2.0)

                if (imol .eq. jmol) then
                    jmi = jb - ib
                    !Changed for male cell
                    if (imol < 47) then
                        cvRaw(jmi, ifr) = cvRaw(jmi, ifr) + pij
                        counter(jmi) = counter(jmi) + 1
                    endif
                    cvRaw(interTAD_intra_chrom_index_map(icm,jcm), ifr)  = cvRaw(interTAD_intra_chrom_index_map(icm,jcm), ifr) + pij
                    counter(interTAD_intra_chrom_index_map(icm,jcm))  = counter(interTAD_intra_chrom_index_map(icm,jcm)) + 1
                elseif (ceil1 .ne. ceil2) then
                    cvRaw(inter_chrom_index_map(ceil1,ceil2),ifr) = cvRaw(inter_chrom_index_map(ceil1,ceil2),ifr) + pij
                    counter(inter_chrom_index_map(ceil1,ceil2)) = counter(inter_chrom_index_map(ceil1,ceil2)) + 1
                endif

            enddo
        enddo                


        do icv = 1, ncv
            if (counter(icv) > 0) then
                cvRaw(icv,ifr) = cvRaw(icv,ifr) / counter(icv)
            endif
        enddo

    enddo
    print *, "finished fine"

    cvRaw_avg = sum(cvRaw(:,first_frame:nframes),2) 

    write_file_1 = trim(write_file_path)//'contact_prob.txt'
    write_file_2 = trim(write_file_path)//'nframes.txt'
    write_file_3 = trim(write_file_path)//'counter.txt'

    Open (unit=11, File=trim(write_file_1))
    do ib = 1, ncv
        write(11,'(E15.8)') cvRaw_avg(ib)/(nframes - first_frame + 1)
    enddo
    close(11)

    open(11, file=trim(write_file_2))
    write(11, *) nframes - first_frame + 1
    close(11)

    open(11, file=trim(write_file_3))
    do icv = 1, ncv
        write(11, '(I10)') counter(icv)
    enddo
    close(11)

    deallocate(x)
    deallocate(y)
    deallocate(z)

    deallocate(cvRaw)
    deallocate(cvRaw_avg)
    deallocate(counter)

    deallocate(interTAD_intra_chrom_index_map)
end program main
