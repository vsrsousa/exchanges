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

end module hamiltonian_mod
