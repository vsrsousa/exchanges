module diag_mod
  use parameters, only: dp
  implicit none

contains

  subroutine print_diag_iz(iz, sumJ_baseline, sumJ_stream)
    integer, intent(in) :: iz
    real(dp), intent(in) :: sumJ_baseline, sumJ_stream
    write(*,'(5x,a,i4,1x,a,2(1x,f12.6))') 'DIAG_iz: iz=', iz, ' sumJ_baseline, sumJ_stream =', sumJ_baseline, sumJ_stream
  end subroutine print_diag_iz

  subroutine print_diag_accum(ia,ja, idim, jdim, delta_block, Gz12, Gz21, tmp1, zstep)
    use parameters, only: dp
    implicit none
    integer, intent(in) :: ia, ja, idim, jdim
    complex(dp), intent(in) :: delta_block(idim,idim)
    complex(dp), intent(in) :: Gz12(idim,jdim)
    complex(dp), intent(in) :: Gz21(jdim,idim)
    complex(dp), intent(in) :: tmp1(idim,idim)
    complex(dp), intent(in) :: zstep
    integer :: ii, jj

    write(*,'(5x,a,i4,a,i4)') 'DIAG_accum: ia=',ia,' ja=',ja
    write(*,'(5x,a)') ' DIAG_accum: delta block (re,im):'
    do ii = 1, idim
      do jj = 1, idim
        write(*,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: delta(',ii,',',jj,')=', real(delta_block(ii,jj)), aimag(delta_block(ii,jj))
      end do
    end do

    write(*,'(5x,a)') ' DIAG_accum: Gz(ia,ja,*,*,spin=2) (re,im):'
    do ii = 1, idim
      do jj = 1, jdim
        write(*,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: Gz(ia,ja)(',ii,',',jj,')=', real(Gz12(ii,jj)), aimag(Gz12(ii,jj))
      end do
    end do

    write(*,'(5x,a)') ' DIAG_accum: Gz(ja,ia,*,*,spin=1) (re,im):'
    do ii = 1, jdim
      do jj = 1, idim
        write(*,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: Gz(ja,ia)(',ii,',',jj,')=', real(Gz21(ii,jj)), aimag(Gz21(ii,jj))
      end do
    end do

    write(*,'(5x,a)') ' DIAG_accum: tmp1 (re,im):'
    do ii = 1, idim
      do jj = 1, idim
        write(*,'(5x,a,i3,a,i3,a,1x,f12.6,1x,f12.6)') ' DIAG_accum: tmp1(',ii,',',jj,')=', real(tmp1(ii,jj)), aimag(tmp1(ii,jj))
      end do
    end do

    write(*,'(5x,a,1x,2(f12.6))') ' DIAG_accum: DIMAG(tmp1*zstep) sample =', DIMAG(tmp1(1,1)*zstep), DIMAG(tmp1(idim,idim)*zstep)

  end subroutine print_diag_accum

end module diag_mod
