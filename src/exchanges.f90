! Copyright (C) Dmitry Korotin dmitry@korotin.name

program exchange_parameters

  use general
  use parameters
  use green_mod
  use iomodule
  use hamiltonian_mod
  use meminfo
  use mesh_mod
  use exchange_utils
  use io_mod
  use diag_mod
  use streaming_mod
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
  ! read_hamilt moved to hamiltonian_mod

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



  ! define the z mesh (moved to mesh_mod)
  call build_zmesh(nz1,nz2,nz3,emin,emax,height,z,nz)

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
      call print_diag_iz(iz, sumJ_baseline, sumJ_stream)
      deallocate(tmp_loc)
    end do
    deallocate(Gtest)
  end if

  ! streaming over z: compute G for one z, accumulate Jorb and occupations
  call compute_streaming(nz, z, nnnbrs, nblocks, MAXVAL(block_dim), H, parent, taunew, block_start, block_dim, delta, occ, Jorb, Jexc, dbg_ia, dbg_ja)

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
        ! delegate printing to io_mod
        call print_exchange_pair(ia, ja, Jexc(ia,ja), Jorb(ia,ja,1:block_dim(parent(ia)),1:block_dim(parent(ia))), kb_ev, pos_delta)

    END DO
  END DO

  write(stdout,*)
  write(stdout,*) '    Computed orbitals occupations should coincide with your DFT results'
  write(stdout,*) '    If they differ significantly - check your integration contour.'
  call print_occupations(occ,parent,block_dim,nnnbrs,nspin)


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
