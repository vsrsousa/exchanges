module mesh_mod
  use parameters, only: dp
  implicit none

contains

  subroutine build_zmesh(nz1,nz2,nz3,emin,emax,height,z,nz)
    use parameters, only: dp
    implicit none
    integer, intent(in) :: nz1, nz2, nz3
    real(dp), intent(in) :: emin, emax, height
    complex(dp), allocatable, intent(out) :: z(:)
    integer, intent(out) :: nz
    complex(dp) :: zstep
    integer :: i

    nz = nz1 + nz2 + nz3
    allocate(z(nz+1))
    z = cmplx(0.0,0.0,dp)

    z(1) = cmplx(emin,0.d0,dp)

    zstep = cmplx(0.d0,(height/dble(nz1)),dp)
    do i = 2, nz1+1
      z(i) = z(i-1) + zstep
    end do

    zstep = cmplx((emax-emin)/dble(nz2),0.d0,dp)
    do i = nz1+2, nz1+nz2
      z(i) = z(i-1) + zstep
    end do

    zstep = -1.d0*cmplx(0.d0,(height/dble(nz3)),dp)
    do i = nz1+nz2+1, nz+1
      z(i) = z(i-1) + zstep
    end do

  end subroutine build_zmesh

end module mesh_mod
