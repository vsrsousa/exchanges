module io_mod
  use parameters, only: dp
  use iomodule, only: out_unit, stdout
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
    write(s_dist,'(F7.3)') pos_delta
    write(stdout,*)
    write(stdout,'(5x,A,3x,I3,2x,A,2x,I3)') 'Exchange interaction between atoms', ia, 'and', ja

    if (trim(adjustl(out_unit)) == 'K') then
      write(s_mev,'(F12.6)') real(Jexc_scalar)/kb_ev
      write(stdout,'(5x,A,1x,A)') trim(s_mev)//' K', '(distance: '//trim(s_dist)//')'
      write(stdout,'(5x,A)') 'Orbital exchange interaction matrix J_{i,j,m,n} (in K)'
      do i = 1, idim
        write(stdout,'(5x)',advance='no')
        do j = 1, idim
          write(stdout,'(1x,F11.2)',advance='no') real(Jorb_block(i,j)/kb_ev)
        end do
        write(stdout,*)
      end do
    else
      ! default: meV
      write(s_mev,'(F12.6)') real(Jexc_scalar)*1.0d3
      write(stdout,'(5x,A,1x,A)') trim(s_mev)//' meV', '(distance: '//trim(s_dist)//')'
      write(stdout,'(5x,A)') 'Orbital exchange interaction matrix J_{i,j,m,n} (in meV)'
      do i = 1, idim
        write(stdout,'(5x)',advance='no')
        do j = 1, idim
          write(stdout,'(1x,F11.2)',advance='no') real(Jorb_block(i,j)*1.0d3)
        end do
        write(stdout,*)
      end do
    end if

  end subroutine print_exchange_pair

  subroutine print_occupations(occ,parent,block_dim,nnnbrs,nspin)
    use parameters, only: dp
    use general, only: blocks
    implicit none
    real(dp), intent(in) :: occ(:,:,:)
    integer, intent(in) :: parent(:)
    integer, intent(in) :: block_dim(:)
    integer, intent(in) :: nnnbrs, nspin
    integer :: ia, j, i, id

    do ia=1, nnnbrs
      do j=1, nspin
        write(stdout,'(/5x,a8,i3,a5,i2)') 'For atom', ia, 'spin', j
        id = blocks(parent(ia))%dim
        do i = 1, id
            write(stdout,'(5x,a8,i2,a13,f6.3)') 'Orbital', i, ' occupation: ', occ(ia,j,i)
        end do
      end do
    end do

  end subroutine print_occupations

  subroutine print_program_header(exec_start_ts, nthreads)
    character(len=*), intent(in) :: exec_start_ts
    integer, intent(in) :: nthreads
    write(stdout,'(/,5x,a66)') '------------------------------------------------------------------'
    write(stdout,'(5x,a66)')   '                      Program EXCHANGES                           '
    write(stdout,'(5x,a66)')   ' for calculation of exchange parameters of the Heisenberg model.  '
    write(stdout,'(/,5x,a66)') 'Please cite "D. M. Korotin et al., Phys. Rev. B 91, 224405 (2015)"'
    write(stdout,'(5x,a66)')   '    in publications or presentations arising from this work.      '
    write(stdout,'(5x,a66,/)') '------------------------------------------------------------------'
    write(stdout,'(/,5x,a66)') '   We are using the model with the exchange term defined as:      '
    write(stdout,'(5x,a66,/)') '                  H = -\sum_ij J_{ij} e_i e_j,                    '
    write(stdout,'(5x,a66,/)') ' where e_i,j are unit vectors and sum runs once over ions pairs   '
    write(stdout,'(5x,a66,/)') '------------------------------------------------------------------'
    write(stdout,'(/,5x,a11,i3,a8,/)') 'Running in ', nthreads, ' threads'
    write(stdout,'(5x,A)') trim(exec_start_ts)//'  Start execution'
    write(stdout,'(5x,A)') ''
  end subroutine print_program_header

  subroutine print_all_exchanges(nnnbrs, parent, block_dim, Jexc, Jorb, taunew)
    use parameters, only: dp, kb_ev
    use general, only: blocks
    implicit none
    integer, intent(in) :: nnnbrs
    integer, intent(in) :: parent(:)
    integer, intent(in) :: block_dim(:)
    complex(dp), intent(in) :: Jexc(:,:)
    complex(dp), intent(in) :: Jorb(:,:,:,:)
    real(dp), intent(in) :: taunew(3,*)
    integer :: ia, ja
    real(dp) :: pos_delta

    do ia = 1, nnnbrs
      do ja = ia+1, nnnbrs
        pos_delta = sqrt( (taunew(1,ia) - taunew(1,ja))**2 + (taunew(2,ia) - taunew(2,ja))**2 + (taunew(3,ia) - taunew(3,ja))**2 )
        call print_exchange_pair(ia, ja, Jexc(ia,ja), Jorb(ia,ja,1:blocks(parent(ia))%dim,1:blocks(parent(ia))%dim), kb_ev, pos_delta)
      end do
    end do

  end subroutine print_all_exchanges

  subroutine sync_out_unit_from_input()
    ! (removed) sync handled by iomodule:set_out_unit to avoid naming conflicts
  end subroutine sync_out_unit_from_input

end module io_mod
