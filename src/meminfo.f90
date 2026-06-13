module meminfo
  use general
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

  subroutine print_mem_status(stage)
    implicit none
    character(len=*), intent(in) :: stage
    integer(kind=8) :: rss_kb, vmpeak_kb
    real(8) :: rss_gb
    character(len=128) :: s_rss_kb, s_vmpeak_kb, s_rss_gb, line
      ! character(len=20) :: ts

    rss_kb = -1_8
    vmpeak_kb = -1_8

    call get_mem_kb(rss_kb, vmpeak_kb)
    rss_gb = 0.0_8
    if (rss_kb .gt. 0_8) rss_gb = real(rss_kb,8) / (1024.0_8*1024.0_8)

    write(s_rss_kb,'(I0)') rss_kb
    write(s_vmpeak_kb,'(I0)') vmpeak_kb
    write(s_rss_gb,'(F8.3)') rss_gb

        write(line,'(A)') trim(adjustl(stage)) // ': RSS=' // trim(adjustl(s_rss_kb)) // &
          ' kB  Peak=' // trim(adjustl(s_vmpeak_kb)) // &
          ' kB  (' // trim(adjustl(s_rss_gb)) // ' GB)'
    write(*,'(5x,A)') ''
      write(*,'(5x,A)') trim(line)
  end subroutine print_mem_status

  subroutine print_estimated_G(est_bytes)
    implicit none
    integer(kind=8), intent(in) :: est_bytes
    real(8) :: est_mb, est_gb
    character(len=128) :: s_est_bytes, s_est_mb, s_est_gb, line
    character(len=20) :: ts

    est_mb = real(est_bytes,8) / 1024.0_8 / 1024.0_8
    est_gb = est_mb / 1024.0_8

    write(s_est_bytes,'(I0)') est_bytes
    write(s_est_mb,'(F8.2)') est_mb
    write(s_est_gb,'(F8.3)') est_gb

        write(line,'(A)') 'Estimated memory for G: ' // trim(adjustl(s_est_bytes)) // &
          ' bytes (' // trim(adjustl(s_est_mb)) // ' MB, ' // trim(adjustl(s_est_gb)) // ' GB)'
    write(*,'(5x,A)') ''
      write(*,'(5x,A)') trim(line)
  end subroutine print_estimated_G
  pure function int4_to_str(val) result(str)
    integer, intent(in) :: val
    character(len=4) :: str
    integer :: tmp, i, d
    tmp = val
    do i = 4, 1, -1
      d = mod(tmp,10)
      str(i:i) = achar(48 + d)
      tmp = tmp / 10
    end do
  end function int4_to_str

  pure function int2_to_str(val) result(str)
    integer, intent(in) :: val
    character(len=2) :: str
    integer :: tmp, i, d
    tmp = val
    do i = 2, 1, -1
      d = mod(tmp,10)
      str(i:i) = achar(48 + d)
      tmp = tmp / 10
    end do
  end function int2_to_str

  subroutine get_timestamp(ts)
    implicit none
    character(len=20), intent(out) :: ts
    integer :: v(8)
    character(len=4) :: y
    character(len=2) :: mo, da, hh, mm, ss

    call date_and_time(values=v)
    y = int4_to_str(v(1))
    mo = int2_to_str(v(2))
    da = int2_to_str(v(3))
    hh = int2_to_str(v(5))
    mm = int2_to_str(v(6))
    ss = int2_to_str(v(7))

    ts = y//'-'//mo//'-'//da//' '//hh//':'//mm//':'//ss

  end subroutine get_timestamp

  subroutine finalize_exchanges(z, occ, delta, Jorb, Jexc, tmp1, istart_idx, idim_idx, iend_idx)
    use parameters, only: dp
    implicit none
    complex(dp), allocatable, intent(inout) :: z(:)
    real(dp), allocatable, intent(inout) :: occ(:,:,:)
    complex(dp), allocatable, intent(inout) :: delta(:,:)
    complex(dp), allocatable, intent(inout) :: Jorb(:,:,:,:)
    complex(dp), allocatable, intent(inout) :: Jexc(:,:)
    complex(dp), allocatable, intent(inout) :: tmp1(:,:)
    integer, allocatable, intent(inout) :: istart_idx(:), idim_idx(:), iend_idx(:)

    if (allocated(z)) deallocate(z)
    if (allocated(occ)) deallocate(occ)
    if (allocated(delta)) deallocate(delta)
    if (allocated(Jorb)) deallocate(Jorb)
    if (allocated(Jexc)) deallocate(Jexc)
    if (allocated(tmp1)) deallocate(tmp1)
    if (allocated(istart_idx)) deallocate(istart_idx)
    if (allocated(idim_idx)) deallocate(idim_idx)
    if (allocated(iend_idx)) deallocate(iend_idx)

    call clear()

  end subroutine finalize_exchanges

end module meminfo

subroutine clear()

  use general
  if( allocated(atoms) ) deallocate(atoms)
  if( allocated(h) ) deallocate(h)
  if( allocated(wk) ) deallocate(wk)
  if( allocated(xk) ) deallocate(xk)
  if( allocated(block_atom) ) deallocate( block_atom )
  if( allocated(block_l) ) deallocate( block_l )
  if( allocated(block_dim) ) deallocate( block_dim )
  if( allocated(block_orbitals) ) deallocate( block_orbitals )
  if( allocated(block_start) ) deallocate(block_start)
  if( allocated(blocks) ) deallocate(blocks)

end subroutine clear
