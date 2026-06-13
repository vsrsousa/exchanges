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
  use setup_mod
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
  logical :: diag_verbose
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
  diag_iz_max = 0
  diag_verbose = .false.

  call system_clock(time_start,count_rate)
  call get_timestamp(exec_start_ts)
  call print_program_header(exec_start_ts, OMP_get_max_threads())

  read(stdin, exchanges, iostat=ios)
  if( ios .ne. 0 ) stop "Can't read input"

  write(stdout,'(5x,a34)') 'Parameters of integration contour:'
  write(stdout,'(5x,a7,f7.3,a9,f6.3,a11,f5.3,a5)') 'emin = ', emin, '  emax = ', emax, '  height = ', height,' (eV)'
  write(stdout,'(5x,a6,i4,a8,i4,a8,i4)') 'nz1 = ', nz1, '  nz2 = ', nz2, '  nz3 = ', nz3

  call setup_system(mode, distance, l, atom_of_interest)



  ! define the z mesh (moved to mesh_mod)
  call build_zmesh(nz1,nz2,nz3,emin,emax,height,z,nz)

  !Now compute full Green function for all atoms
  ! estimate memory required for G
  total_elems = int(nz,8) * int(nnnbrs,8) * int(nnnbrs,8) * int(MAXVAL(block_dim),8) * int(MAXVAL(block_dim),8) * int(nspin,8)
  est_bytes = total_elems * bytes_per_complex
  est_mb = real(est_bytes,dp) / 1024.0_dp / 1024.0_dp
  est_gb = real(est_bytes,dp) / 1024.0_dp / 1024.0_dp / 1024.0_dp
  call print_estimated_G(est_bytes)

  ! allocate exchange/occupation buffers
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

  call init_exchange_buffers(nnnbrs, MAXVAL(block_dim), occ, Jorb, Jexc, tmp1)

  ! prepare per-atom indices
  allocate(istart_idx(nnnbrs), idim_idx(nnnbrs), iend_idx(nnnbrs))
  call prepare_indices(nnnbrs, parent, block_start, block_dim, istart_idx, idim_idx, iend_idx)

  ! --- Per-iz Jorb checksum diagnostics (recompute now that delta is available)
  if (diag_iz_max .gt. 0) then
    call run_diag_iz(diag_iz_max, z, nz, nnnbrs, nblocks, MAXVAL(block_dim), H, parent, taunew, block_start, block_dim, delta)
  end if

  ! streaming over z: compute G for one z, accumulate Jorb and occupations
  call compute_streaming(nz, z, nnnbrs, nblocks, MAXVAL(block_dim), H, parent, taunew, block_start, block_dim, delta, occ, Jorb, Jexc, dbg_ia, dbg_ja)

  deallocate(istart_idx, idim_idx, iend_idx)

  deallocate(tmp1)

  Jexc = (-1.d0/tpi)*Jexc
  Jorb = (-1.d0/tpi)*Jorb

  DO ia = 1, nnnbrs
    DO ja = ia+1, nnnbrs
        ! delegate printing to io_mod
    END DO
  END DO

  call print_all_exchanges(nnnbrs, parent, block_dim, Jexc, Jorb, taunew)

  write(stdout,*)
  write(stdout,*) '    Computed orbitals occupations should coincide with your DFT results'
  write(stdout,*) '    If they differ significantly - check your integration contour.'
  call print_occupations(occ,parent,block_dim,nnnbrs,nspin)


  call finalize_exchanges(z, occ, delta, Jorb, Jexc, tmp1, istart_idx, idim_idx, iend_idx)

  call system_clock(time_end,count_rate)

  elapsed = real(time_end - time_start, dp) / real(count_rate, dp)
  call get_timestamp(exec_ts)
  write(stdout,'(5x,A)') ''
  write(stdout,'(5x,A,F8.3,A)') trim(exec_ts)//'  Execution time: ', elapsed, ' seconds'

end program exchange_parameters
