module hamiltonian_mod
  use iomodule
  use parameters, only : dp, ang_to_bohr
  use general
  implicit none

contains

  subroutine read_hamilt()
    implicit none
    integer :: i,j,ispin, ik
    real(dp) :: Hre, Him

    call open_input_file(iunhamilt,'hamilt.am')

    call find_section(iunhamilt,'&nspin')
    read(iunhamilt,*) nspin
    if (nspin .ne. 2) stop 'A spin-polarized hamiltonian is neccessary'

    call find_section(iunhamilt,'&nkp')
    read(iunhamilt,*) nkp

    call find_section(iunhamilt,'&dim')
    read(iunhamilt,*) hdim

    allocate( h(hdim,hdim,nkp,nspin) )
    allocate( wk(nkp) )
    allocate( xk(3,nkp) )

    h = cmplx(0.0,0.0)
    wk = 0.0
    xk = 0.0

    call find_section(iunhamilt,'&kpoints')
    do ik=1, nkp
      read(iunhamilt,*) wk(ik), xk(:,ik)
    end do
        
    call find_section(iunhamilt,'&hamiltonian')
    
    do ispin = 1, nspin
      do ik = 1, nkp
        do i=1, hdim
          do j=i, hdim
            read(iunhamilt,*) Hre, Him
            h(i,j,ik,ispin) = cmplx(Hre,Him)
            h(j,i,ik,ispin) = dconjg( h(i,j,ik,ispin) )
          end do
        end do
      end do
    end do
    
    close(iunhamilt)

  end subroutine read_hamilt

  subroutine read_crystal()

    implicit none
    integer :: i,j
    character(len=3) :: dummy
    character :: l_symbol
    character(len=32) :: tmp_label
    character(len=128) :: cell_line
    real(dp) :: tmp_alat
    character(len=8) :: tmp_unit
    integer :: ios

    call open_input_file(iunsystem,'system.am')

    call find_section(iunsystem,'&cell')
    cell_line = ''
    read(iunsystem,'(A)', iostat=ios) cell_line
    if (ios /= 0) stop 'Can''t read &cell line in system.am'
    tmp_unit = 'bohr'
    tmp_alat = 0.0_dp
    read(cell_line,*, iostat=ios) tmp_alat, tmp_unit
    if (ios /= 0) then
      tmp_unit = 'bohr'
      read(cell_line,*, iostat=ios) tmp_alat
      if (ios /= 0) stop 'Bad &cell specification'
    end if
    ! normalize tmp_unit
    tmp_unit = adjustl(trim(tmp_unit))
    if ( tmp_unit(1:3) == 'ang' .or. tmp_unit == 'angstrom' ) then
      alat = tmp_alat * ang_to_bohr
      cell_info%unit = 'ang'
    else
      alat = tmp_alat
      cell_info%unit = 'bohr'
    end if

    do i = 1,3
      read(iunsystem,*) cell_info%vec(:,i)
    end do
    cell_info%alat = alat
    
    call find_section(iunsystem,'&atoms')
    read(iunsystem,*) natoms

    if (allocated(atoms)) deallocate(atoms)
    allocate(atoms(natoms))

    do i=1, natoms
      read(iunsystem,*) tmp_label, atoms(i)%pos
      atoms(i)%label = adjustl(trim(tmp_label))
    end do

    call find_section(iunsystem,'&basis')
    read(iunsystem,*) i, nblocks
    if(i .ne. hdim) stop 'Basis dimention in system.am and hamilt.am are unequal'

    allocate( block_atom(nblocks) )
    allocate( block_l(nblocks) )
    allocate( block_dim(nblocks) )
    allocate( block_orbitals(nblocks,7) )
    allocate( block_start(nblocks) )
    block_orbitals = 0

    do i = 1, nblocks
      read(iunsystem,*) dummy, block_atom(i), block_l(i), block_dim(i), block_start(i), block_orbitals(i,1:block_dim(i))
    end do

    ! Populate derived-type blocks for easier handling (and keep legacy arrays for compatibility)
    allocate(blocks(nblocks))
    do i = 1, nblocks
      blocks(i)%atom = block_atom(i)
      blocks(i)%l = adjustl(block_l(i))
      blocks(i)%dim = block_dim(i)
      blocks(i)%start = block_start(i)
      if (blocks(i)%dim > 0) then
        allocate(blocks(i)%orbitals(blocks(i)%dim))
        blocks(i)%orbitals = block_orbitals(i,1:blocks(i)%dim)
      else
        allocate(blocks(i)%orbitals(0))
      end if
    end do

    call find_section(iunsystem,'&efermi')
    read(iunsystem,*) efermi

    close(iunsystem)

  end subroutine read_crystal


  subroutine atoms_list(mode,distance,atom_of_interest,l_of_interest)

    implicit none
    character(len=10), intent(in) :: mode
    real(dp), intent(in) :: distance
    character(len=1), intent(in) :: l_of_interest
    integer, intent(in) :: atom_of_interest

    ! temp storage
    real(dp) :: taunew_(3,maxnnbrs), vect(3), d
    integer :: parent_(maxnnbrs), nnnbrs_

    integer :: i,j, iblock, nvect, idx
    logical :: have_atom_already

    character(len=40) :: currentline

    nnnbrs = 0
    parent = -1
    taunew = -1000

    call find_nnbrs(natoms,atoms,atom_of_interest,distance,nnnbrs_,taunew_,parent_)

    select case ( trim(mode) )
      case ('list') ! list mode
        currentline = ''
        do while( ios .ge. 0 .AND. trim(currentline) .ne. 'ATOMS_LIST' )
          read(stdin, fmt="(a40)", iostat=ios ) currentline
          if ( ios .gt. 0 ) stop "Can not find ATOMS_LIST in stdin"
        end do
        read(stdin,*) nvect

        do i = 1, nvect
          read(stdin,*) vect(1:3)
          do j = 1, nnnbrs_
            call haa(taunew_, vect, have_atom_already, idx)
            if( .not. have_atom_already ) then
              write(stdout,'(/,"Error: Can not find atom in position", 3f9.5, " within sphere of ",f9.5," with ",a2," orbitals"/)') &
                vect,distance,l_of_interest
              stop ' '
            else
              nnnbrs = nnnbrs + 1
              parent(nnnbrs) = parent_(idx)
              taunew(:,nnnbrs) = vect
              exit
            end if
          end do
        end do
        ! reset arrays
        nnnbrs_ = nnnbrs
        parent_ = parent
        taunew_ = taunew

      case default ! 'distance mode'

        ! We want to consider all hamilt atoms anyway
        do i = 1, natoms
          call haa(taunew_, atoms(i)%pos, have_atom_already, idx)
          if( .not. have_atom_already ) then
            nnnbrs_ = nnnbrs_+1
            taunew_(:,nnnbrs_) = atoms(i)%pos
            parent_(nnnbrs_) = i
          end if
        end do

    end select

    nnnbrs = 0
    parent = -1
    taunew = -1000

    !filter atoms by l (support mixed labels like 'sp','spd')
    do i = 1, nnnbrs_
      do iblock = 1, nblocks
        if( blocks(iblock)%atom .eq. parent_(i) ) then
          if ( trim(l_of_interest) == 'all' .or. index(adjustl(blocks(iblock)%l), trim(l_of_interest)) > 0 ) then
            nnnbrs = nnnbrs + 1
            parent(nnnbrs) = iblock
            taunew(:,nnnbrs) = taunew_(:,i)
          end if
        end if
      end do
    end do

  end subroutine atoms_list

end module hamiltonian_mod
