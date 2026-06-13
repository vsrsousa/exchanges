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
  use input_mod
  use calc_mod
  use omp_lib

  implicit none

  integer :: time_start, time_end, count_rate, nz, ia, ja
  real(dp) :: est_gb, est_mb, rss_gb, elapsed
  character(len=20) :: exec_ts, exec_start_ts
  complex(dp), allocatable :: z(:)
  real(dp), allocatable :: occ(:,:,:)
  integer, allocatable :: istart_idx(:), idim_idx(:), iend_idx(:)
  complex(dp), allocatable :: Jexc(:,:), Jorb(:,:,:,:)
  complex(dp), allocatable :: delta(:,:), tmp1(:,:)
  integer :: dbg_ia, dbg_ja
  integer(kind=8) :: rss_kb, peak_kb, total_elems, est_bytes
  integer, parameter :: bytes_per_complex = 16

  ! Input parameters are below:

  ! namelist and defaults are in `input_mod`


  call system_clock(time_start,count_rate)
  call get_timestamp(exec_start_ts)
  call print_program_header(exec_start_ts, OMP_get_max_threads())

  call read_input()

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

  call calculate_exchanges(z, nz, H)

  call system_clock(time_end,count_rate)

  elapsed = real(time_end - time_start, dp) / real(count_rate, dp)
  call get_timestamp(exec_ts)
  write(stdout,'(5x,A)') ''
  write(stdout,'(5x,A,F8.3,A)') trim(exec_ts)//'  Execution time: ', elapsed, ' seconds'

end program exchange_parameters
