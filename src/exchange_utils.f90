module exchange_utils
  use parameters, only: dp
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

end module exchange_utils
