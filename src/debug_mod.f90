module debug_mod
  implicit none
  save

  integer :: diag_unit = -1
  logical :: diag_use_stdout = .true.

contains

  subroutine diag_print(msg)
    character(len=*), intent(in) :: msg
    if (diag_use_stdout) then
      write(*,'(A)') trim(msg)
    else
      write(diag_unit,'(A)') trim(msg)
    end if
  end subroutine diag_print

  subroutine diag_open_file(filename, iostat)
    character(len=*), intent(in) :: filename
    integer, intent(out), optional :: iostat
    integer :: ios

    if (.not. diag_use_stdout .and. diag_unit > 0) then
      close(diag_unit, iostat=ios)
    end if

    open(newunit=diag_unit, file=filename, status='replace', action='write', iostat=ios)
    if (present(iostat)) iostat = ios
    if (ios == 0) then
      diag_use_stdout = .false.
    else
      diag_use_stdout = .true.
    end if
  end subroutine diag_open_file

  subroutine diag_close()
    integer :: ios
    if (.not. diag_use_stdout .and. diag_unit > 0) then
      close(diag_unit, iostat=ios)
      diag_unit = -1
    end if
    diag_use_stdout = .true.
  end subroutine diag_close

end module debug_mod
