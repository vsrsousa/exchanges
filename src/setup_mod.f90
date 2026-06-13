module setup_mod
  use parameters, only: dp
  use general
  use iomodule
  use hamiltonian_mod
  use meminfo
  use io_mod
  implicit none
contains

  subroutine setup_system(mode, distance, l, atom_of_interest)
    character(len=*), intent(in) :: mode
    real(dp), intent(in) :: distance
    character(len=1), intent(in) :: l
    integer, intent(inout) :: atom_of_interest
    implicit none
    complex(dp), allocatable :: hksum(:,:,:)
    integer :: i

    call read_hamilt()
    call read_crystal()

    call print_mem_status('After reading inputs')

    if(atom_of_interest == -100) atom_of_interest = block_atom(1)

    if( iverbosity .ge. 3) then
      allocate( hksum(hdim,hdim,nspin) )
      hksum = cmplx(0.0,0.0,dp)
      do i=1, nkp
        hksum(:,:,:) = hksum(:,:,:) + h(:,:,i,:)*wk(i)
      end do

      do i=1, nspin
        write(stdout,'(/,5x,a13,i2,a1)') 'H(0) for spin', i, ':'
        call output_matrix_by_blocks(hdim,dreal(hksum(:,:,i)))
      end do

      deallocate( hksum )
    end if

    write(stdout,'(/,5x,a20,f13.9,a5)') 'Cell constant (alat):', alat, 'Bohr'
    write(stdout,'(5x,20a)') 'Cell vectors (rows):'
    do i = 1, 3
      write(stdout,'(7x,3f9.5)') cell(:,i)
    end do
    write(stdout,'(5x,a26)') 'Atoms of the initial cell:'
    do i = 1, natoms
      write(stdout,'(7x,a3,x,3f9.5)') atomlabel(i), tau(:,i)
    end do

    write(stdout,'(/,5x,a71)') 'We have the following basis for the hamiltonian and the green function:'
    write(stdout,'(5x,i3,a20,i3,a8)') hdim, ' orbitals grouped in', nblocks, ' blocks:'
    do i=1, nblocks
      write(fmt,'(i1)') block_dim(i)
      write(stdout,'(5x,a3,a7,i2,a2,3x,i1,x,a,a11,'//adjustl(fmt)//'a11,a1)') &
          adjustr(atomlabel(block_atom(i))), ' (atom ', block_atom(i), '):', &
          block_dim(i),block_l(i), '-orbitals (', (orbitals(block_orbitals(i,j)), j=1, block_dim(i)), ')'
    end do

    call atoms_list(mode,distance,atom_of_interest,l)

    write(stdout,'(/5x,a17,i4,a8)') 'We will consider ', nnnbrs,' atoms: '

    do i = 1, nnnbrs
      write(stdout,'(5x,i2,a2,a3,x,3f9.5,3x,a12,f9.5,a7)') i, ': ', &
        atomlabel( block_atom( parent(i) ) ), taunew(:,i), '( distance =', &
        sqrt( (taunew(1,i)-tau(1,atom_of_interest))**2 + &
              (taunew(2,i)-tau(2,atom_of_interest))**2 + &
              (taunew(3,i)-tau(3,atom_of_interest))**2 ), ' alat )'
    end do

  end subroutine setup_system

end module setup_mod
