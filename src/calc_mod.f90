module calc_mod
  use parameters, only: dp, tpi
  use general, only: hdim, nspin, nnnbrs, nblocks, block_dim, block_start, parent, taunew, blocks
  use exchange_utils, only: compute_delta, init_exchange_buffers, prepare_indices
  use streaming_mod, only: compute_streaming
  use io_mod, only: print_all_exchanges, print_occupations
  use iomodule, only: stdout
  use diag_mod, only: run_diag_iz
  use meminfo, only: finalize_exchanges, print_estimated_G
  use input_mod, only: diag_iz_max
  implicit none

contains

  subroutine calculate_exchanges(z, nz, h)
    implicit none
    complex(dp), allocatable, intent(inout) :: z(:)
    integer, intent(in) :: nz
    complex(dp), intent(in) :: h(:,:,:,:)

    complex(dp), allocatable :: delta(:,:), tmp1(:,:)
    real(dp), allocatable :: occ(:,:,:)
    complex(dp), allocatable :: Jexc(:,:), Jorb(:,:,:,:)
    integer, allocatable :: istart_idx(:), idim_idx(:), iend_idx(:)
    integer(kind=8) :: total_elems, est_bytes
    integer :: maxbd, i
    integer, parameter :: bytes_per_complex = 16

    ! estimate memory for G (informational)
    maxbd = 0
    if (nblocks .gt. 0) then
      do i = 1, nblocks
        if (blocks(i)%dim > maxbd) maxbd = blocks(i)%dim
      end do
    end if
    total_elems = int(nz,8) * int(nnnbrs,8) * int(nnnbrs,8) * int(maxbd,8) * int(maxbd,8) * int(nspin,8)
    est_bytes = total_elems * bytes_per_complex
    call print_estimated_G(est_bytes)

    allocate(delta(hdim,hdim))
    call compute_delta(h,delta)

    call init_exchange_buffers(nnnbrs, maxbd, occ, Jorb, Jexc, tmp1)

    allocate(istart_idx(nnnbrs), idim_idx(nnnbrs), iend_idx(nnnbrs))
    call prepare_indices(nnnbrs, parent, block_start, block_dim, istart_idx, idim_idx, iend_idx)

    if (diag_iz_max .gt. 0) then
      call run_diag_iz(diag_iz_max, z, nz, nnnbrs, nblocks, maxbd, h, parent, taunew, block_start, block_dim, delta)
    end if

    call compute_streaming(nz, z, nnnbrs, nblocks, maxbd, h, parent, taunew, block_start, block_dim, delta, occ, Jorb, Jexc, -1, -1)

    deallocate(istart_idx, idim_idx, iend_idx)
    deallocate(tmp1)

    Jexc = (1.d0/tpi) * Jexc
    Jorb = (1.d0/tpi) * Jorb

    call print_all_exchanges(nnnbrs, parent, block_dim, Jexc, Jorb, taunew)
    write(stdout,*)
    write(stdout,*) '    Computed orbitals occupations should coincide with your DFT results'
    write(stdout,*) '    If they differ significantly - check your integration contour.'
    call print_occupations(occ,parent,block_dim,nnnbrs,nspin)

    call finalize_exchanges(z, occ, delta, Jorb, Jexc, tmp1, istart_idx, idim_idx, iend_idx)

  end subroutine calculate_exchanges

end module calc_mod
