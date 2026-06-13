! Copyright (C) Dmitry Korotin dmitry@korotin.name

module general
 
  use parameters, only : dp, maxnnbrs 	 
  implicit none
  save

  !Crystal cell
  real(dp) :: alat, cell(3,3)
  integer :: natoms

  type :: cell_t
    real(dp) :: alat
    real(dp) :: vec(3,3)
    character(len=8) :: unit  ! 'bohr' or 'ang'
  end type cell_t

  type(cell_t) :: cell_info

  type :: atom_t
    character(len=:), allocatable :: label
    real(dp) :: pos(3)
  end type atom_t

  type(atom_t), allocatable :: atoms(:)

  !Hamiltonian
  integer :: &
            nkp, &  ! number of k-points
            hdim, & ! hamiltonian dimention
            nspin	!

  real(dp), allocatable :: &
            wk(:), & ! k-points weight (nkp)
            xk(:,:) ! k-points coordinates in 2pi/a units (3,nkp)
  complex(dp), allocatable :: h(:,:,:,:) ! Hamiltonian matrix (hdim,hdim,nkp,nspin)

  
  !Basis
  integer :: nblocks
  type :: block_t
    integer :: atom
    character(len=:), allocatable :: l
    integer :: dim
    integer, allocatable :: orbitals(:)
    integer :: start
  end type block_t

  type(block_t), allocatable :: blocks(:)

  integer, allocatable :: block_atom(:), block_dim(:), block_orbitals(:,:)
  character, allocatable :: block_l(:)

  integer, allocatable :: block_start(:)

  !Other
  real(dp) :: efermi

  ! Nearest neighbours
  real(dp) :: &
             taunew(3,maxnnbrs) ! positions of nearest neighbours
  integer :: nnnbrs, & ! total number of nearest neighbours
             parent(maxnnbrs) ! parent atom for every neighbour


  ! That is how the atomic orbitals are encoded in the system.am file
  character(len=11) :: orbitals(16)
  data orbitals/'s', 'y', 'z', 'x', 'xy', 'yz', '3z^2-1', 'xz', 'x^2-y^2', 'y(3x^2-y^2)', 'xyz', &
                'y(5z^2-1)', 'z(5z^2-3)', 'x(5z^2-1)', 'z(x^2-y^2)', 'x(3y^2-x^2)'/


end module general
