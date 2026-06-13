module diag_mod
  use parameters, only: dp
  implicit none

contains

  subroutine print_diag_iz(iz, sumJ_baseline, sumJ_stream)
    integer, intent(in) :: iz
    real(dp), intent(in) :: sumJ_baseline, sumJ_stream
    write(stdout,'(5x,a,i4,1x,a,2(1x,f12.6))') 'DIAG_iz: iz=', iz, ' sumJ_baseline, sumJ_stream =', sumJ_baseline, sumJ_stream
  end subroutine print_diag_iz

end module diag_mod
