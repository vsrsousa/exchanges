module setup_mod
  use parameters, only: dp, bohr_to_ang
  use general
  use input_mod
  use iomodule
  use hamiltonian_mod
  use meminfo
  use io_mod
  implicit none
contains

  subroutine setup_system(mode, distance, l, atom_of_interest, dist_unit)
    implicit none
    character(len=*), intent(in) :: mode
    real(dp), intent(in) :: distance
    character(len=1), intent(in) :: l
    integer, intent(inout) :: atom_of_interest
    character(len=6), intent(in) :: dist_unit
    complex(dp), allocatable :: hksum(:,:,:)
    integer :: i, j
    character(len=8) :: fmt
    real(dp) :: dist_alat, dist_bohr, dist_print
    character(len=16) :: s_num

    call read_hamilt()
    call read_crystal()

    call print_mem_status('After reading inputs')

    if(atom_of_interest == -100) atom_of_interest = blocks(1)%atom

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

    if (trim(cell_info%unit) == 'ang') then
      write(stdout,'(/,5x,a20,f13.9,a5)') 'Cell constant (alat):', alat*bohr_to_ang, 'Ang'
    else
      write(stdout,'(/,5x,a20,f13.9,a5)') 'Cell constant (alat):', alat, 'Bohr'
    end if
    write(stdout,'(5x,20a)') 'Cell vectors (rows):'
    do i = 1, 3
      write(stdout,'(7x,3f9.5)') cell_info%vec(:,i)
    end do
    write(stdout,'(5x,a26)') 'Atoms of the initial cell:'
    do i = 1, natoms
      write(stdout,'(7x,a3,x,3f9.5)') trim(adjustl(atoms(i)%label)), atoms(i)%pos
    end do

    write(stdout,'(/,5x,a71)') 'We have the following basis for the hamiltonian and the green function:'
    write(stdout,'(5x,i3,a20,i3,a8)') hdim, ' orbitals grouped in', nblocks, ' blocks:'
    do i=1, nblocks
      write(fmt,'(i1)') blocks(i)%dim
        write(stdout,'(5x,a3,a7,i2,a2,3x,i1,x,a,a11,'//adjustl(fmt)//'a11,a1)') &
          trim(adjustl(atoms(blocks(i)%atom)%label)), ' (atom ', blocks(i)%atom, '):', &
          blocks(i)%dim,blocks(i)%l, '-orbitals (', (orbitals(blocks(i)%orbitals(j)), j=1, blocks(i)%dim), ')'
    end do

    call atoms_list(mode,distance,atom_of_interest,l)

    write(stdout,'(/5x,a17,i4,a8)') 'We will consider ', nnnbrs,' atoms: '

    do i = 1, nnnbrs
      dist_alat = sqrt( (taunew(1,i)-atoms(atom_of_interest)%pos(1))**2 + &
                (taunew(2,i)-atoms(atom_of_interest)%pos(2))**2 + &
                (taunew(3,i)-atoms(atom_of_interest)%pos(3))**2 )
      dist_bohr = dist_alat * alat
      dist_print = dist_bohr * bohr_to_ang  ! always print in Angstrom
      write(s_num,'(F9.5)') dist_print
      write(stdout,'(5x,i2,a2,a3,x,3f9.5,3x,a)') i, ': ', &
        trim(adjustl(atoms(blocks(parent(i))%atom)%label)), taunew(:,i), '( distance = '//trim(s_num)//' Ang)'
    end do

  end subroutine setup_system

end module setup_mod
