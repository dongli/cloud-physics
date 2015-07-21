module ui

    implicit none

contains

    subroutine report_error(message, file_path, line_number)

        character(*), intent(in) :: message
        character(*), intent(in), optional :: file_path
        integer, intent(in), optional :: line_number

        write(*, "('[Error]: ', A, ': ', I4, ': ', A)") file_path, line_number, message
        stop 1

    end subroutine report_error

end module ui
