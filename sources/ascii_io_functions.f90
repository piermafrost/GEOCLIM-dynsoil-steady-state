module ascii_io_functions
        implicit none

        contains


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
