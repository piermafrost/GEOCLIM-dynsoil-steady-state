module ascii_io_functions
    implicit none

    contains


        !------------------------------------------------------------------------------------------------------------------!


        subroutine read_comment( funit, ierr )
            ! read file until it finds uncommented line
            ! IMPORTANT: commented lines are BLANK LINES or LINES BEGINNING BY #
            integer, intent(in):: funit
            integer, intent(out), optional:: ierr
            character(len=1):: line

            line = '#'

            if (present(ierr)) then

                ! read until finding uncommented line or error is raised
                ierr = 0
                do while (line=='#' .and. ierr==0)
                read(unit=funit, fmt=*, iostat=ierr) line
                end do

                ! if no error raised, get back to the previous record (in that case, the previous line)
                if (ierr==0) then
                backspace(unit=funit)
                end if

            else

                ! read until finding uncommented line
                do while (line=='#')
                read(unit=funit, fmt=*) line
                end do

                ! get back to the previous record (in that case, the previous line)
                backspace(unit=funit)

            end if

        end subroutine


        !------------------------------------------------------------------------------------------------------------------!


        function file_length(filename, fileunit)
            ! return the length (ie, number of lines) of a text file. Input argument: file path (char)
            character(len=*), intent(in), optional:: filename
            integer, intent(in), optional:: fileunit
            integer:: ierr, k, funit, file_length

            ! Open file (if file name was given):
            if (present(filename)) then
                funit = 1
                open(unit=funit, file=filename, status='old', action='read')
            elseif (present(fileunit)) then
                    funit = fileunit
            else
                print *, 'Error while calling "file_length" function.'
                print *, 'Not enough input arguments (one at least is expected)'
                stop
            end if

            ! Read file until end-of-file:
            file_length = 0
            read(unit=funit, fmt=*, iostat=ierr)
            do while (ierr==0)
                file_length = file_length + 1
                read(unit=funit, fmt=*, iostat=ierr)
            end do

            ! If error NOT type end-of-file is reported (ierr>0), rewind and read the file again without iostat argument
            ! to force printing the error message:
            if (ierr>0) then
                rewind(unit=funit)
                do k = 0,file_length
                    read(unit=funit, fmt=*)
                end do
            end if

            ! If file name was given, close file. If not, rewind file:
            if (present(filename)) then
                close(unit=funit)
            elseif (present(fileunit)) then
                rewind(unit=funit)
            end if

        end function


        !------------------------------------------------------------------------------------------------------------------!


end module
