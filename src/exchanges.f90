! Copyright (C) Dmitry Korotin dmitry@korotin.name

program exchange_parameters

  use general
  use parameters
  use green_mod
  use iomodule
  use meminfo
  use omp_lib

  implicit none

  integer :: i, j, time_start, time_end, count_rate, nz, ia, ja, idim, jdim, iz, istart, jstart, iend, jend
  character(len=3) :: fmt='   ' ! is used for pretty output only
  real(dp) :: pos_delta
  real(dp) :: est_gb, est_mb, rss_gb
  real(dp) :: elapsed
  character(len=20) :: exec_ts, exec_start_ts
  character(len=32) :: s_mev, s_k, s_dist
  complex(dp), allocatable :: z(:), Gz(:,:,:,:,:), delta(:,:), tmp1(:,:), hksum(:,:,:)
  real(dp), allocatable :: occ(:,:,:)
  integer, allocatable :: istart_idx(:), idim_idx(:), iend_idx(:)
  integer :: ii, jj
  complex(dp) :: zstep, tmp2
  complex(dp), allocatable :: Jexc(:,:), Jorb(:,:,:,:)
  integer(kind=8) :: rss_kb, peak_kb, total_elems, est_bytes
  integer, parameter :: bytes_per_complex = 16
  logical :: diag_pass
  logical :: dbg_print
  integer :: dbg_ia, dbg_ja, dbg_i, dbg_j, dbg_ispin
  complex(dp), allocatable :: Gtest(:,:,:,:,:,:)
  complex(dp), allocatable :: ztest(:)
  integer :: ia2,ja2,i2,j2,ispin2
  real(dp) :: max_abs, max_rel, a, b
  real(dp) :: max_sum
  integer :: diag_iz_max
  real(dp) :: sumJ_baseline, sumJ_stream
  complex(dp), allocatable :: tmp_loc(:,:)
  

  ! Input parameters are below:
  
  ! Intergation mesh:
  integer :: nz1, nz2, nz3
  real(dp) :: emin, emax, height
  
  integer :: atom_of_interest
  character(len=1) :: l
  
  real(dp) :: distance ! distance for the nearest neighbours search
                       ! could be used as integer do define a coordination sphere number
                       ! if mode='csphere'
  character(len=10) :: mode ! mode for the nearest neighbours search

  namelist /exchanges/ nz1, nz2, nz3, height, emin, emax, distance, iverbosity, mode, l, atom_of_interest

  ! default values
  nz1 = 150
  nz2 = 3000
  nz3 = 150
  height = 0.01
  emin = -30.0
  emax = 0.0
  distance = 8.d-1
  mode = 'distance'
  l = 'd'
  atom_of_interest = -100
  
  call system_clock(time_start,count_rate)

  write(stdout,'(/,5x,a66)') '------------------------------------------------------------------'
  write(stdout,'(5x,a66)')   '                      Program EXCHANGES                           '
  write(stdout,'(5x,a66)')   ' for calculation of exchange parameters of the Heisenberg model.  '
  write(stdout,'(/,5x,a66)') 'Please cite "D. M. Korotin et al., Phys. Rev. B 91, 224405 (2015)"'
  write(stdout,'(5x,a66)')   '    in publications or presentations arising from this work.      '
  write(stdout,'(5x,a66,/)') '------------------------------------------------------------------'
  write(stdout,'(5x,a66,/)') '   We are using the model with the exchange term defined as:      '
  write(stdout,'(5x,a66,/)') '                   H = \sum_ij J_{ij} e_i e_j,                    '
  write(stdout,'(5x,a66,/)') ' where e_i,j are unit vectors and sum runs once over ions pairs   '
  write(stdout,'(5x,a66,/)') '------------------------------------------------------------------'

  write(stdout,'(/,5x,a11,i3,a8,/)') 'Running in ', OMP_get_max_threads(), ' threads'
  call get_timestamp(exec_start_ts)
  write(stdout,'(5x,A)') trim(exec_start_ts)//'  Start execution'
  write(stdout,'(5x,A)') ''

  read(stdin, exchanges, iostat=ios)
  if( ios .ne. 0 ) stop "Can't read input"

  write(stdout,'(5x,a34)') 'Parameters of integration contour:'
  write(stdout,'(5x,a7,f7.3,a9,f6.3,a11,f5.3,a5)') 'emin = ', emin, '  emax = ', emax, '  height = ', height,' (eV)'
  write(stdout,'(5x,a6,i4,a8,i4,a8,i4)') 'nz1 = ', nz1, '  nz2 = ', nz2, '  nz3 = ', nz3
  
  call read_hamilt()
  call read_crystal()

  ! report memory after reading input H and crystal
  call print_mem_status('After reading inputs')

  if(atom_of_interest == -100) atom_of_interest = block_atom(1)

  ! debug
  if( iverbosity .ge. 3) then
    ! output H(0)
    allocate( hksum(hdim,hdim,nspin) )
    hksum = cmplx(0.0,0.0,dp)

    do i=1, nkp
      hksum(:,:,:) = hksum(:,:,:) + h(:,:,i,:)*wk(i)
    end do

    do i=1, nspin
      write(stdout,'(/,5x,a13,i2,a1)') 'H(0) for spin', i, ':'
      call output_matrix_by_blocks(hdim,dreal(hksum(:,:,i)))
    end do

    deallocate( hksum ) 
  end if
  ! end of debug
  
  ! Output of the readed values
  write(stdout,'(/,5x,a20,f13.9,a5)') 'Cell constant (alat):', alat, 'Bohr'
  write(stdout,'(5x,20a)') 'Cell vectors (rows):'
  do i = 1, 3
    write(stdout,'(7x,3f9.5)') cell(:,i)
  end do
  write(stdout,'(5x,a26)') 'Atoms of the initial cell:'
  do i = 1, natoms
    write(stdout,'(7x,a3,x,3f9.5)') atomlabel(i), tau(:,i)
  end do
  
  write(stdout,'(/,5x,a71)') 'We have the following basis for the hamiltonian and the green function:'
  write(stdout,'(5x,i3,a20,i3,a8)') hdim, ' orbitals grouped in', nblocks, ' blocks:'
  do i=1, nblocks
    write(fmt,'(i1)') block_dim(i)
    write(stdout,'(5x,a3,a7,i2,a2,3x,i1,x,a,a11,'//adjustl(fmt)//'a11,a1)') &
        adjustr(atomlabel(block_atom(i))), ' (atom ', block_atom(i), '):', &
        block_dim(i),block_l(i), '-orbitals (', (orbitals(block_orbitals(i,j)), j=1, block_dim(i)), ')'
  end do
  
  call atoms_list(mode,distance,atom_of_interest,l)

  ! Pretty output
  write(stdout,'(/5x,a17,i4,a8)') 'We will consider ', nnnbrs,' atoms: '

  do i = 1, nnnbrs
    pos_delta = sqrt( (taunew(1,i)-tau(1,atom_of_interest))**2 + &
                  (taunew(2,i)-tau(2,atom_of_interest))**2 + &
                  (taunew(3,i)-tau(3,atom_of_interest))**2 )
    write(stdout,'(5x,i2,a2,a3,x,3f9.5,3x,a12,f9.5,a7)') i, ': ', &
      atomlabel( block_atom( parent(i) ) ), taunew(:,i), '( distance =', pos_delta, ' alat )'
  end do


  ! define the z mesh
  ! i will generate nz+1 mesh points as a trick for zstep calculation below
  ! only nz points will be used in fact.

  nz = nz1 + nz2 + nz3
  allocate(z(nz+1))
  z = cmplx(0.0,0.0,dp)

  z(1) = cmplx(emin,0.d0,dp)

  zstep = cmplx(0.d0,(height/nz1),dp)
  do i = 2, nz1+1
    z(i) = z(i-1) + zstep
  end do

  zstep = cmplx((emax-emin)/nz2,0.d0,dp)
  do i = nz1+2, nz1+nz2
    z(i) = z(i-1) + zstep
  end do

  zstep = -1.d0*cmplx(0.d0,(height/nz3),dp)
  do i = nz1+nz2+1, nz+1
    z(i) = z(i-1) + zstep
  end do

  !Now compute full Green function for all atoms
  ! estimate memory required for G
  total_elems = int(nz,8) * int(nnnbrs,8) * int(nnnbrs,8) * int(MAXVAL(block_dim),8) * int(MAXVAL(block_dim),8) * int(nspin,8)
  est_bytes = total_elems * bytes_per_complex
  est_mb = real(est_bytes,dp) / 1024.0_dp / 1024.0_dp
  est_gb = real(est_bytes,dp) / 1024.0_dp / 1024.0_dp / 1024.0_dp
  call print_estimated_G(est_bytes)

  ! allocate temporary storage for single-z Green function and occupations (streaming)
  allocate(Gz(nnnbrs,nnnbrs,MAXVAL(block_dim),MAXVAL(block_dim),nspin))
  allocate(occ(nnnbrs,nspin,MAXVAL(block_dim)))
  occ = 0.0_dp
  ! memory status after allocating per-z buffers omitted to match reference format

  ! Diagnostics removed for clean output (kept timestamps only)

  allocate(delta(hdim,hdim))
  call compute_delta(h,delta)

  ! debug
  if( iverbosity .ge. 3) then
    write(stdout,'(/,5x,a15)') 'Computed delta:'
    call output_matrix_by_blocks(hdim,dreal(delta))
  end if
  ! end of debug

  allocate(Jexc(nnnbrs,nnnbrs))
  Jexc = cmplx(0.0,0.0,dp)
  allocate( Jorb(nnnbrs,nnnbrs,MAXVAL(block_dim),MAXVAL(block_dim)))
  Jorb = cmplx(0.0,0.0,dp)

  allocate( tmp1(MAXVAL(block_dim),MAXVAL(block_dim)) )

  ! prepare per-atom indices
  allocate(istart_idx(nnnbrs), idim_idx(nnnbrs), iend_idx(nnnbrs))
  do ia = 1, nnnbrs
    istart_idx(ia) = block_start(parent(ia))
    idim_idx(ia) = block_dim(parent(ia))
    iend_idx(ia) = istart_idx(ia) + idim_idx(ia) - 1
  end do

  ! --- Per-iz Jorb checksum diagnostics (recompute now that delta is available)
  if (diag_iz_max .gt. 0) then
    allocate(Gtest(1,nnnbrs,nnnbrs,MAXVAL(block_dim),MAXVAL(block_dim),nspin))
    do iz = 1, diag_iz_max
      zstep = z(iz+1) - z(iz)
      allocate(tmp_loc(MAXVAL(block_dim),MAXVAL(block_dim)))
      call compute_g(1,nnnbrs,nblocks,MAXVAL(block_dim),Gtest,H,(/z(iz)/),parent,taunew,block_start,block_dim)
      call compute_g_onez(nnnbrs,nblocks,MAXVAL(block_dim),Gz,H,z(iz),parent,taunew,block_start,block_dim)
      sumJ_baseline = 0.0_dp
      sumJ_stream = 0.0_dp
      do ia2 = 1, nnnbrs
        istart = block_start(parent(ia2))
        idim = block_dim(parent(ia2))
        iend = istart + idim - 1
        do ja2 = ia2+1, nnnbrs
          jstart = block_start(parent(ja2))
          jdim = block_dim(parent(ja2))
          jend = jstart + jdim - 1
          if (idim .ne. jdim) cycle
          tmp_loc = cmplx(0.0,0.0,dp)
          tmp_loc(1:idim,1:idim) = MATMUL( MATMUL(delta(istart:iend,istart:iend), Gtest(1,ia2,ja2,1:idim,1:jdim,2)), MATMUL(delta(jstart:jend,jstart:jend), Gtest(1,ja2,ia2,1:jdim,1:idim,1)) )
          sumJ_baseline = sumJ_baseline + sum( DIMAG(tmp_loc(1:idim,1:idim)* zstep ) )
          tmp_loc(1:idim,1:idim) = MATMUL( MATMUL(delta(istart:iend,istart:iend), Gz(ia2,ja2,1:idim,1:jdim,2)), MATMUL(delta(jstart:jend,jstart:jend), Gz(ja2,ia2,1:jdim,1:idim,1)) )
          sumJ_stream = sumJ_stream + sum( DIMAG(tmp_loc(1:idim,1:idim)* zstep ) )
        end do
      end do
      write(stdout,'(5x,a,i4,1x,a,2(1x,f12.6))') 'DIAG_iz: iz=', iz, ' sumJ_baseline, sumJ_stream =', sumJ_baseline, sumJ_stream
      deallocate(tmp_loc)
    end do
    deallocate(Gtest)
  end if

  ! streaming over z: compute G for one z, accumulate Jorb and occupations
  do iz = 1, nz
    zstep = z(iz+1) - z(iz)
    call compute_g_onez(nnnbrs,nblocks,MAXVAL(block_dim),Gz,H,z(iz),parent,taunew,block_start,block_dim)

    do ia = 1, nnnbrs
      istart = istart_idx(ia)
      idim = idim_idx(ia)
      iend = iend_idx(ia)
      do ja = ia+1, nnnbrs
        jstart = istart_idx(ja)
        jdim = idim_idx(ja)
        jend = iend_idx(ja)

        if (idim .ne. jdim) then
          write(stdout,*) ia,ja,idim,jdim
          stop 'Not equal subblocks size'
        end if

        tmp1 = cmplx(0.0,0.0,dp)
        tmp1(1:idim,1:idim) = MATMUL( &
                          MATMUL(delta(istart:iend,istart:iend),Gz(ia,ja,1:idim,1:jdim,2)), &
                          MATMUL(delta(jstart:jend,jstart:jend),Gz(ja,ia,1:jdim,1:idim,1)) &
                          )

            if (ia==dbg_ia .and. ja==dbg_ja) then
              write(stdout,'(5x,a,i4,a,i4)') 'DIAG_accum: ia=',ia,' ja=',ja
              write(stdout,'(5x,a)') ' DIAG_accum: delta block (re,im):'
              do ii = 1, idim
                do jj = 1, idim
                  write(stdout,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: delta(',ii,',',jj,')=', real(delta(istart+ii-1,istart+jj-1)), aimag(delta(istart+ii-1,istart+jj-1))
                end do
              end do

              write(stdout,'(5x,a)') ' DIAG_accum: Gz(ia,ja,*,*,spin=2) (re,im):'
              do ii = 1, idim
                do jj = 1, jdim
                  write(stdout,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: Gz(ia,ja)(',ii,',',jj,')=', real(Gz(ia,ja,ii,jj,2)), aimag(Gz(ia,ja,ii,jj,2))
                end do
              end do

              write(stdout,'(5x,a)') ' DIAG_accum: Gz(ja,ia,*,*,spin=1) (re,im):'
              do ii = 1, jdim
                do jj = 1, idim
                  write(stdout,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: Gz(ja,ia)(',ii,',',jj,')=', real(Gz(ja,ia,ii,jj,1)), aimag(Gz(ja,ia,ii,jj,1))
                end do
              end do

              write(stdout,'(5x,a)') ' DIAG_accum: tmp1 (re,im):'
              do ii = 1, idim
                do jj = 1, idim
                  write(stdout,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: tmp1(',ii,',',jj,')=', real(tmp1(ii,jj)), aimag(tmp1(ii,jj))
                end do
              end do

              write(stdout,'(5x,a,1x,2(f12.6))') ' DIAG_accum: DIMAG(tmp1*zstep) sample =', DIMAG(tmp1(1,1)*zstep), DIMAG(tmp1(idim,idim)*zstep)
            end if

            ! accumulate orbital-resolved contribution
            Jorb(ia,ja,1:idim,1:idim) = Jorb(ia,ja,1:idim,1:idim) + DIMAG(tmp1(1:idim,1:idim)*zstep)
            ! accumulate scalar exchange (sum over orbital block)
            Jexc(ia,ja) = Jexc(ia,ja) + sum( DIMAG(tmp1(1:idim,1:idim)*zstep) )

      end do
    end do

    ! accumulate orbital occupations from Gz diagonal
    do ia = 1, nnnbrs
      do j = 1, nspin
        do i = 1, idim_idx(ia)
          occ(ia,j,i) = occ(ia,j,i) + ((-1.d0/pi) * DIMAG( Gz(ia,ia,i,i,j) * zstep ))
        end do
      end do
    end do

  end do

  deallocate(istart_idx, idim_idx, iend_idx)

  deallocate(tmp1)

  Jexc = (-1.d0/tpi)*Jexc
  Jorb = (-1.d0/tpi)*Jorb

  DO ia = 1, nnnbrs
    DO ja = ia+1, nnnbrs
        ! Compute distance between atoms for pretty output
        pos_delta= SQRT((taunew(1,ia) - taunew(1,ja))**2+ &
                        (taunew(2,ia) - taunew(2,ja))**2+ &
                        (taunew(3,ia) - taunew(3,ja))**2 )
        !
        write(s_mev,'(F12.6)') real(Jexc(ia,ja))*1.0d3
        write(s_k,'(F7.2)') real(Jexc(ia,ja))/kb_ev
        write(s_dist,'(F7.3)') pos_delta
        write(stdout,'(5x,A,3x,I3,2x,A,2x,I3)') 'Exchange interaction between atoms', ia, 'and', ja
        write(stdout,'(7x,A,1x,A,1x,A)') trim(s_mev)//' meV =', trim(s_k)//' K', '(distance: '//trim(s_dist)//')'

        write(stdout,*) 'Orbital exchange interaction matrix J_{i,j,m,n} (in K and meV)'
        do i=1,block_dim(parent(ia))
          write(stdout,'(7x,5(1x,F11.2),5x,5(1x,F12.6))') (real(Jorb(ia,ja,i,j)/kb_ev), j=1,block_dim(parent(ia))), (real(Jorb(ia,ja,i,j)*1.0d3), j=1,block_dim(parent(ia)))
        end do

    END DO
  END DO

  write(stdout,*)
  write(stdout,*) '    Computed orbitals occupations should coincide with your DFT results'
  write(stdout,*) '    If they differ significantly - check your integration contour.'
  DO ia=1,nnnbrs
    DO j=1,nspin
      write(stdout,'(/5x,a8,i3,a5,i2)') 'For atom', ia, 'spin', j
      DO i = 1, block_dim(parent(ia))
        write(stdout,'(7x,a8,i2,a13,f6.3)') 'Orbital', i, ' occupation: ', occ(ia,j,i)
      END DO
    END DO
  END DO


  if( allocated(z) ) deallocate(z)
  if( allocated(Gz) ) deallocate(Gz)
  if( allocated(occ) ) deallocate(occ)
  if( allocated(delta) ) deallocate(delta)
  if( allocated(Jorb) ) deallocate(Jorb)
  if( allocated(Jexc) ) deallocate(Jexc)
  call clear()

  call system_clock(time_end,count_rate)

  elapsed = real(time_end - time_start, dp) / real(count_rate, dp)
  call get_timestamp(exec_ts)
  write(stdout,'(5x,A)') ''
  write(stdout,'(5x,A,F8.3,A)') trim(exec_ts)//'  Execution time: ', elapsed, ' seconds'

end program exchange_parameters

subroutine compute_delta(h,delta)
  use parameters, only : dp
  use general, only : hdim, nspin, nkp, wk
  
  implicit none
  complex(dp) :: h(hdim,hdim,nkp,nspin), delta(hdim,hdim)
  integer :: ik

  delta = cmplx(0.0,0.0,dp)
  
  ! delta = h(spin_up) - h(spin_down)
  do ik=1, nkp
    delta(:,:) = delta(:,:) + wk(ik)*( h(:,:,ik,1) - h(:,:,ik,2) )
  end do

  delta = dreal(delta) 

END SUBROUTINE compute_delta

subroutine atoms_list(mode,distance,atom_of_interest,l_of_interest)

  use iomodule
  use parameters, only : dp, maxnnbrs
  use general
  
  implicit none
  character(len=10), intent(in) :: mode
  real(dp), intent(in) :: distance
  character(len=1), intent(in) :: l_of_interest
  integer, intent(in) :: atom_of_interest

  ! temp storage
  real(dp) :: taunew_(3,maxnnbrs), vect(3), d
  integer :: parent_(maxnnbrs), nnnbrs_

  integer :: i,j, iblock, nvect, index
  logical :: have_atom_already

  character(len=40) :: currentline

  nnnbrs = 0
  parent = -1
  taunew = -1000

  !call find_nnbrs(natoms,tau,cell,block_atom(1),distance,nnnbrs_,taunew_,parent_)
  call find_nnbrs(natoms,tau,cell,atom_of_interest,distance,nnnbrs_,taunew_,parent_)

  select case ( trim(mode) )
    case ('list') ! list mode
      ! a primitive way to find ATOMS_LIST string in stdin
      currentline = ''
      do while( ios .ge. 0 .AND. trim(currentline) .ne. 'ATOMS_LIST' )
        read(stdin, fmt="(a40)", iostat=ios ) currentline
        if ( ios .gt. 0 ) stop "Can not find ATOMS_LIST in stdin"
      end do
      !call find_section(stdin,'ATOMS_LIST')
      read(stdin,*) nvect

      do i = 1, nvect
        read(stdin,*) vect(1:3)
        do j = 1, nnnbrs_
          call haa(taunew_, vect, have_atom_already, index)
          if( .not. have_atom_already ) then
            write(stdout,'(/,"Error: Can not find atom in position", 3f9.5, " within sphere of ",f9.5," with ",a2," orbitals"/)') &
              vect,distance,l_of_interest
            stop ' '
          else
            nnnbrs = nnnbrs + 1
            parent(nnnbrs) = parent_(index)
            taunew(:,nnnbrs) = vect
            exit
          end if
        end do
      end do
      ! reset arrays
      nnnbrs_ = nnnbrs
      parent_ = parent
      taunew_ = taunew

    case default ! 'distance mode'

      ! We want to consider all hamilt atoms anyway
      do i = 1, natoms
        call haa(taunew_, tau(:,i), have_atom_already, index)
        if( .not. have_atom_already ) then
          nnnbrs_ = nnnbrs_+1
          taunew_(:,nnnbrs_) = tau(:,i)
          parent_(nnnbrs_) = i
        end if
      end do

  end select

  nnnbrs = 0
  parent = -1
  taunew = -1000

  !filter atoms by l
  do i = 1, nnnbrs_
    do iblock = 1, nblocks
      if( block_atom(iblock) .eq. parent_(i) .and. block_l(iblock) .eq. l_of_interest ) then
        nnnbrs = nnnbrs + 1
        parent(nnnbrs) = iblock
        taunew(:,nnnbrs) = taunew_(:,i)
      end if
    end do 
  end do

end subroutine atoms_list

subroutine read_crystal()

  use iomodule
  use parameters, only : dp
  use general
  
  implicit none
  integer :: i,j
  character(len=3) :: dummy
  character :: l_symbol

  call open_input_file(iunsystem,'system.am')

  call find_section(iunsystem,'&cell')
  read(iunsystem,*) alat
  
  do i = 1,3
    read(iunsystem,*) cell(:,i)
  end do
  
  call find_section(iunsystem,'&atoms')
  read(iunsystem,*) natoms
  allocate( tau(3,natoms) )
  allocate( atomlabel(natoms) )

  do i=1, natoms
    read(iunsystem,*) atomlabel(i), tau(:,i)
  end do

  call find_section(iunsystem,'&basis')
  read(iunsystem,*) i, nblocks
  if(i .ne. hdim) stop 'Basis dimention in system.am and hamilt.am are unequal'

  allocate( block_atom(nblocks) )
  allocate( block_l(nblocks) )
  allocate( block_dim(nblocks) )
  allocate( block_orbitals(nblocks,7) )
  allocate( block_start(nblocks) )
  block_orbitals = 0

  do i = 1, nblocks
    read(iunsystem,*) dummy, block_atom(i), block_l(i), block_dim(i), block_start(i), block_orbitals(i,1:block_dim(i))
  end do

  call find_section(iunsystem,'&efermi')
  read(iunsystem,*) efermi

  close(iunsystem)

end subroutine read_crystal

subroutine read_hamilt()

  use iomodule
  use parameters, only : dp
  use general
  
  implicit none
  integer :: i,j,ispin, ik
  real(dp) :: Hre, Him

  call open_input_file(iunhamilt,'hamilt.am')

  call find_section(iunhamilt,'&nspin')
  read(iunhamilt,*) nspin
  if (nspin .ne. 2) stop 'A spin-polarized hamiltonian is neccessary'

  call find_section(iunhamilt,'&nkp')
  read(iunhamilt,*) nkp

  call find_section(iunhamilt,'&dim')
  read(iunhamilt,*) hdim

  allocate( h(hdim,hdim,nkp,nspin) )
  allocate( wk(nkp) )
  allocate( xk(3,nkp) )

  h = cmplx(0.0,0.0)
  wk = 0.0
  xk = 0.0

  call find_section(iunhamilt,'&kpoints')
  do ik=1, nkp
    read(iunhamilt,*) wk(ik), xk(:,ik)
  end do
      
  call find_section(iunhamilt,'&hamiltonian')
  
  do ispin = 1, nspin
    do ik = 1, nkp
      do i=1, hdim
        do j=i, hdim
          read(iunhamilt,*) Hre, Him
          h(i,j,ik,ispin) = cmplx(Hre,Him)
          h(j,i,ik,ispin) = dconjg( h(i,j,ik,ispin) )
        end do
      end do
    end do
  end do
  
  close(iunhamilt)

end subroutine read_hamilt

subroutine clear()

  use general

  if( allocated(tau) ) deallocate(tau)
  if( allocated(atomlabel) ) deallocate(atomlabel)
  if( allocated(h) ) deallocate(h)
  if( allocated(wk) ) deallocate(wk)
  if( allocated(xk) ) deallocate(xk)
  if( allocated(block_atom) ) deallocate( block_atom )
  if( allocated(block_l) ) deallocate( block_l )
  if( allocated(block_dim) ) deallocate( block_dim )
  if( allocated(block_orbitals) ) deallocate( block_orbitals )
  if( allocated(block_start) ) deallocate(block_start)

end subroutine clear

subroutine inverse_complex_matrix(dim,a)
  
  use parameters, only : dp
  implicit none

  integer :: dim, info
  integer, allocatable :: ipiv(:)
  complex(dp), intent(inout) :: a(dim,dim)
  complex(dp), allocatable :: work(:)

  allocate(ipiv(dim))
  allocate(work(dim))

  call ZGETRF(dim,dim,a,dim,ipiv,info)
  if(info /= 0) stop "inverse_complex_matrix Error in ZGETRF"
  call ZGETRI(dim,a,dim,ipiv,work,dim,info)
  if(info /= 0) stop "inverse_complex_matrix Error in ZGETRI"

  if (allocated(ipiv)) deallocate(ipiv)
  if (allocated(work)) deallocate(work)

end subroutine
