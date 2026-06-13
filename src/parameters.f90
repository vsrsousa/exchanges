! Copyright (C) Dmitry Korotin dmitry@korotin.name

module parameters
  
  implicit none
  save

  integer, parameter :: &
    dp = selected_real_kind(14,200), &
    maxnnbrs = 300 ! maximum number of nearest neighbours
  
  real(dp), parameter :: pi               = 3.14159265358979323846_dp
  real(dp), parameter :: tpi              = 2.0_dp * pi
  real(dp), parameter :: kb_ev = 8.6173324d-5
  real(dp), parameter :: bohr_to_ang = 0.529177210903_dp
  real(dp), parameter :: ang_to_bohr = 1.0_dp / bohr_to_ang
  

end module parameters
