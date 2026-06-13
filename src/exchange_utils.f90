module exchange_utils
  use parameters, only: dp, pi
  use general, only: hdim, nspin, nkp, wk
  implicit none

contains

  subroutine compute_delta(h,delta)
    use parameters, only: dp
    use general, only: hdim, nspin, nkp, wk
    implicit none
    complex(dp), intent(in) :: h(hdim,hdim,nkp,nspin)
    complex(dp), intent(out) :: delta(hdim,hdim)
    integer :: ik

    delta = cmplx(0.0,0.0,dp)

    ! delta = h(spin_up) - h(spin_down)
    do ik=1, nkp
      delta(:,:) = delta(:,:) + wk(ik)*( h(:,:,ik,1) - h(:,:,ik,2) )
    end do

    delta = dreal(delta)

  end subroutine compute_delta

  subroutine accumulate_pair(ia, ja, idim, jdim, delta_a, delta_b, Gz12, Gz21, zstep, Jorb_block, Jexc_scalar, dbg_ia, dbg_ja)
    use diag_mod
    implicit none
    integer, intent(in) :: ia, ja, idim, jdim
    complex(dp), intent(in) :: delta_a(idim,idim)
    complex(dp), intent(in) :: delta_b(jdim,jdim)
    complex(dp), intent(in) :: Gz12(idim,jdim)
    complex(dp), intent(in) :: Gz21(jdim,idim)
    complex(dp), intent(in) :: zstep
    complex(dp), intent(inout) :: Jorb_block(idim,idim)
    complex(dp), intent(inout) :: Jexc_scalar
    integer, intent(in) :: dbg_ia, dbg_ja
    complex(dp), allocatable :: tmp1(:,:)

    allocate(tmp1(idim,idim))
    tmp1 = cmplx(0.0,0.0,dp)

    tmp1(1:idim,1:idim) = MATMUL( MATMUL(delta_a(1:idim,1:idim), Gz12(1:idim,1:jdim)), MATMUL(delta_b(1:jdim,1:jdim), Gz21(1:jdim,1:idim)) )

    if (ia == dbg_ia .and. ja == dbg_ja) then
      call print_diag_accum(ia,ja,idim,jdim, delta_a, Gz12, Gz21, tmp1, zstep)
    end if

    Jorb_block(1:idim,1:idim) = Jorb_block(1:idim,1:idim) + DIMAG(tmp1(1:idim,1:idim)*zstep)
    Jexc_scalar = Jexc_scalar + sum( DIMAG(tmp1(1:idim,1:idim)*zstep) )

    deallocate(tmp1)

  end subroutine accumulate_pair

  subroutine accumulate_occupations(nnnbrs, idim_idx, nspin, Gz, zstep, occ)
    implicit none
    integer, intent(in) :: nnnbrs, nspin
    integer, intent(in) :: idim_idx(:)
    complex(dp), intent(in) :: Gz(:,:,:,:,:)
    complex(dp), intent(in) :: zstep
    real(dp), intent(inout) :: occ(:,:,:)
    integer :: ia, j, i, idim

    do ia = 1, nnnbrs
      idim = idim_idx(ia)
      do j = 1, nspin
        do i = 1, idim
          occ(ia,j,i) = occ(ia,j,i) + ((-1.d0/pi) * DIMAG( Gz(ia,ia,i,i,j) * zstep ))
        end do
      end do
    end do

  end subroutine accumulate_occupations

  subroutine prepare_indices(nnnbrs, parent, block_start, block_dim, istart_idx, idim_idx, iend_idx)
    implicit none
    integer, intent(in) :: nnnbrs
    integer, intent(in) :: parent(:)
    integer, intent(in) :: block_start(:)
    integer, intent(in) :: block_dim(:)
    integer, intent(out) :: istart_idx(:), idim_idx(:), iend_idx(:)
    integer :: ia

    do ia = 1, nnnbrs
      istart_idx(ia) = block_start(parent(ia))
      idim_idx(ia) = block_dim(parent(ia))
      iend_idx(ia) = istart_idx(ia) + idim_idx(ia) - 1
    end do

  end subroutine prepare_indices

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

  end subroutine inverse_complex_matrix

end module exchange_utils
