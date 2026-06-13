module hamiltonian_mod
  use iomodule
  use parameters, only : dp
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

    call open_input_file(iunsystem,'system.am')

    call find_section(iunsystem,'&cell')
    read(iunsystem,*) alat
    
    do i = 1,3
      read(iunsystem,*) cell(:,i)
    end do
    
    call find_section(iunsystem,'&atoms')
    read(iunsystem,*) natoms
    allocate( tau(3,natoms) )
    allocate( atomlabel(natoms) )

    do i=1, natoms
      read(iunsystem,*) atomlabel(i), tau(:,i)
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

    integer :: i,j, iblock, nvect, index
    logical :: have_atom_already

    character(len=40) :: currentline

    nnnbrs = 0
    parent = -1
    taunew = -1000

    call find_nnbrs(natoms,tau,cell,atom_of_interest,distance,nnnbrs_,taunew_,parent_)

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
            call haa(taunew_, vect, have_atom_already, index)
            if( .not. have_atom_already ) then
              write(stdout,'(/,"Error: Can not find atom in position", 3f9.5, " within sphere of ",f9.5," with ",a2," orbitals"/)') &
                vect,distance,l_of_interest
              stop ' '
            else
              nnnbrs = nnnbrs + 1
              parent(nnnbrs) = parent_(index)
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
          call haa(taunew_, tau(:,i), have_atom_already, index)
          if( .not. have_atom_already ) then
            nnnbrs_ = nnnbrs_+1
            taunew_(:,nnnbrs_) = tau(:,i)
            parent_(nnnbrs_) = i
          end if
        end do

    end select

    nnnbrs = 0
    parent = -1
    taunew = -1000

    !filter atoms by l
    do i = 1, nnnbrs_
      do iblock = 1, nblocks
        if( block_atom(iblock) .eq. parent_(i) .and. block_l(iblock) .eq. l_of_interest ) then
          nnnbrs = nnnbrs + 1
          parent(nnnbrs) = iblock
          taunew(:,nnnbrs) = taunew_(:,i)
        end if
      end do 
    end do

  end subroutine atoms_list

end module hamiltonian_mod
