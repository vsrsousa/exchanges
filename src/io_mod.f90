module io_mod
  use parameters, only: dp
  implicit none

contains

  subroutine print_exchange_pair(ia,ja,Jexc_scalar,Jorb_block,kb_ev,pos_delta)
    use parameters, only: dp
    implicit none
    integer, intent(in) :: ia, ja
    complex(dp), intent(in) :: Jexc_scalar
    complex(dp), intent(in) :: Jorb_block(:,:)
    real(dp), intent(in) :: kb_ev, pos_delta
    character(len=32) :: s_mev, s_k, s_dist
    integer :: i, j, idim

    idim = size(Jorb_block,1)
    write(s_mev,'(F12.6)') real(Jexc_scalar)*1.0d3
    write(s_k,'(F7.2)') real(Jexc_scalar)/kb_ev
    write(s_dist,'(F7.3)') pos_delta
    write(*,*)
    write(*,'(5x,A,3x,I3,2x,A,2x,I3)') 'Exchange interaction between atoms', ia, 'and', ja
    write(*,'(5x,A,1x,A,1x,A)') trim(s_mev)//' meV =', trim(s_k)//' K', '(distance: '//trim(s_dist)//')'

    write(*,'(5x,A)') 'Orbital exchange interaction matrix J_{i,j,m,n} (in K and meV)'
    do i = 1, idim
      write(*,'(5x,5(1x,F11.2),5x,5(1x,F12.6))') (real(Jorb_block(i,j)/kb_ev), j=1,idim), (real(Jorb_block(i,j)*1.0d3), j=1,idim)
    end do

  end subroutine print_exchange_pair

  subroutine print_occupations(occ,parent,block_dim,nnnbrs,nspin)
    use parameters, only: dp
    implicit none
    real(dp), intent(in) :: occ(:,:,:)
    integer, intent(in) :: parent(:)
    integer, intent(in) :: block_dim(:)
    integer, intent(in) :: nnnbrs, nspin
    integer :: ia, j, i, id

    do ia=1, nnnbrs
      do j=1, nspin
        write(*,'(/5x,a8,i3,a5,i2)') 'For atom', ia, 'spin', j
        id = block_dim(parent(ia))
        do i = 1, id
            write(*,'(5x,a8,i2,a13,f6.3)') 'Orbital', i, ' occupation: ', occ(ia,j,i)
        end do
      end do
    end do

  end subroutine print_occupations

  subroutine print_program_header(exec_start_ts, nthreads)
    character(len=*), intent(in) :: exec_start_ts
    integer, intent(in) :: nthreads
    write(*,'(/,5x,a66)') '------------------------------------------------------------------'
    write(*,'(5x,a66)')   '                      Program EXCHANGES                           '
    write(*,'(5x,a66)')   ' for calculation of exchange parameters of the Heisenberg model.  '
    write(*,'(/,5x,a66)') 'Please cite "D. M. Korotin et al., Phys. Rev. B 91, 224405 (2015)"'
    write(*,'(5x,a66)')   '    in publications or presentations arising from this work.      '
    write(*,'(5x,a66,/)') '------------------------------------------------------------------'
    write(*,'(/,5x,a66)') '   We are using the model with the exchange term defined as:      '
    write(*,'(5x,a66,/)') '                   H = \sum_ij J_{ij} e_i e_j,                    '
    write(*,'(5x,a66,/)') ' where e_i,j are unit vectors and sum runs once over ions pairs   '
    write(*,'(5x,a66,/)') '------------------------------------------------------------------'
    write(*,'(/,5x,a11,i3,a8,/)') 'Running in ', nthreads, ' threads'
    write(*,'(5x,A)') trim(exec_start_ts)//'  Start execution'
    write(*,'(5x,A)') ''
  end subroutine print_program_header

end module io_mod
