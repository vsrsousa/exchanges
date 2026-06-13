module meminfo
  implicit none

contains

  subroutine get_mem_kb(rss_kb, vmpeak_kb)
    implicit none
    integer(kind=8), intent(out) :: rss_kb, vmpeak_kb
    character(len=256) :: line
    integer :: ios

    rss_kb = -1_8
    vmpeak_kb = -1_8

    open(unit=99, file='/proc/self/status', status='old', action='read', iostat=ios)
    if (ios /= 0) return

    do
      read(99,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (index(line,'VmRSS:') /= 0) then
        read(line(7:),*) rss_kb
      else if (index(line,'VmPeak:') /= 0) then
        read(line(8:),*) vmpeak_kb
      end if
    end do

    close(99)
  end subroutine get_mem_kb

end module meminfo
