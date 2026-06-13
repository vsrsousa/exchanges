module input_mod
  use parameters, only: dp
  use iomodule, only: iverbosity, stdin  ! brings `iverbosity` and `stdin` into scope for namelist
  implicit none
  save

  integer :: nz1 = 150, nz2 = 3000, nz3 = 150
  real(dp) :: emin = -30.0_dp, emax = 0.0_dp, height = 0.01_dp
  real(dp) :: distance = 8.d-1
  character(len=6) :: distance_unit = 'bohr'  ! 'bohr' (default) or 'ang'
  character(len=10) :: mode = 'distance'
  character(len=1) :: l = 'd'
  character(len=3) :: out_unit = 'meV'  ! 'meV' (default) or 'K'
  integer :: atom_of_interest = -100

  integer :: diag_iz_max = 0
  logical :: diag_verbose = .false.

  namelist /exchanges/ nz1, nz2, nz3, height, emin, emax, distance, distance_unit, iverbosity, mode, l, atom_of_interest, diag_iz_max, diag_verbose, out_unit

contains

  subroutine read_input()
    implicit none
    integer :: ios
    read(stdin, exchanges, iostat=ios)
    if (ios .ne. 0) stop 'Can''t read input'
  end subroutine read_input

  function get_out_unit() result(u)
    character(len=3) :: u
    u = out_unit
  end function get_out_unit

end module input_mod
