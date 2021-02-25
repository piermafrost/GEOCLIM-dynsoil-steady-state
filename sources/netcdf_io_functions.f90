module netcdf_io_functions
implicit none

contains


  !-------------------------------------------------------------------------------------------------------------------------------!

  function netcdf_get_size( filename, dimname )
    use netcdf
    character(len=*), intent(in):: filename, dimname
    integer:: netcdf_get_size
    integer:: ierr, fid, did

      ! open intput file
      ierr = nf90_open( filename, NF90_NOWRITE, fid )
      call nf90_check(ierr, 'Error while openning file '//filename)

      ! get dimension ID
      ierr = nf90_inq_dimid( fid, dimname, did )
      call nf90_check(ierr, 'Error while inquiring ID of dimension '//dimname//' in file '//filename)

      ! get dimension length
      ierr = nf90_inquire_dimension( fid, did, len=netcdf_get_size )
      call nf90_check(ierr, 'Error while getting length of dimension '//dimname//'in file '//filename)

      ! close file
      ierr = nf90_close(fid)
      call nf90_check(ierr, 'Error while closing file '//filename)

  end function


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine load_netcdf_0D( filename, varname, vardata, units )
    use netcdf
    character(len=*), intent(in):: filename, varname
    double precision, intent(out), dimension(:):: vardata
    character(len=*), intent(out), optional:: units
    integer:: fid, vid, ierr

      ! open intput file
      ierr = nf90_open( filename, NF90_NOWRITE, fid )
      call nf90_check(ierr, 'Error while openning file '//filename)

      ! get variable ID
      ierr = nf90_inq_varid( fid, varname, vid )
      call nf90_check(ierr, 'Error while inquiring ID of variable '//varname//' in file '//filename)

      ! get variable units
      if (present(units)) then
        ierr = nf90_get_att( fid, vid, 'units', units )
        call nf90_check( ierr, 'Warning: not able to get attribute "units" of variable '//varname//' in file '//filename, &
                         kill=.false. )
        if (ierr/=NF90_NOERR) units = 'unkown' ! default units
      end if

      ! get variable
      ierr = nf90_get_var( fid, vid, vardata )
      call nf90_check(ierr, 'Error while getting variable '//varname//'in file '//filename)

      ! close file
      ierr = nf90_close(fid)
      call nf90_check(ierr, 'Error while closing file '//filename)

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine load_netcdf_1D( filename, varname, vardata, units, fillvalue, starts, counts )
    use netcdf
    character(len=*), intent(in):: filename, varname
    double precision, intent(out), dimension(:):: vardata
    character(len=*), intent(out), optional:: units
    double precision, intent(out), optional:: fillvalue
    integer, dimension(:), intent(in), optional:: starts, counts
    integer:: fid, vid, ierr

      ! open intput file
      ierr = nf90_open( filename, NF90_NOWRITE, fid )
      call nf90_check(ierr, 'Error while openning file '//filename)

      ! get variable ID
      ierr = nf90_inq_varid( fid, varname, vid )
      call nf90_check(ierr, 'Error while inquiring ID of variable '//varname//' in file '//filename)

      ! get variable units
      if (present(units)) then
        ierr = nf90_get_att( fid, vid, 'units', units )
        call nf90_check( ierr, 'Warning: not able to get attribute "units" of variable '//varname//' in file '//filename, &
                         kill=.false. )
        if (ierr/=NF90_NOERR) units = 'unkown' ! default units
      end if

      ! get variable fillvalue
      if (present(fillvalue)) then
        ierr = nf90_get_att( fid, vid, '_FillValue', fillvalue )
        call nf90_check(ierr, 'Warning: not able to get attribute "_FillValue" of variable '//varname//' in file '//filename, &
                        kill=.false.)
        if (ierr/=NF90_NOERR) fillvalue = 0 ! default fillvalue
      end if

      ! get variable
      if (present(starts) .and. present(counts)) then
        ierr = nf90_get_var( fid, vid, vardata, start=starts, count=counts )
      else
        ierr = nf90_get_var( fid, vid, vardata )
      end if
      call nf90_check(ierr, 'Error while getting variable '//varname//'in file '//filename)

      ! close file
      ierr = nf90_close(fid)
      call nf90_check(ierr, 'Error while closing file '//filename)

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine load_netcdf_2D( filename, varname, vardata, units, fillvalue, starts, counts )
    use netcdf
    character(len=*), intent(in):: filename, varname
    double precision, intent(out), dimension(:,:):: vardata
    character(len=*), intent(out), optional:: units
    double precision, intent(out), optional:: fillvalue
    integer, dimension(:), intent(in), optional:: starts, counts
    integer:: fid, vid, ierr

      ! open intput file
      ierr = nf90_open( filename, NF90_NOWRITE, fid )
      call nf90_check(ierr, 'Error while openning file '//filename)

      ! get variable ID
      ierr = nf90_inq_varid( fid, varname, vid )
      call nf90_check(ierr, 'Error while inquiring ID of variable '//varname//' in file '//filename)

      ! get variable units
      if (present(units)) then
        ierr = nf90_get_att( fid, vid, 'units', units )
        call nf90_check( ierr, 'Warning: not able to get attribute "units" of variable '//varname//' in file '//filename, &
                         kill=.false. )
        if (ierr/=NF90_NOERR) units = 'unkown' ! default units
      end if

      ! get variable fillvalue
      if (present(fillvalue)) then
        ierr = nf90_get_att( fid, vid, '_FillValue', fillvalue )
        call nf90_check(ierr, 'Warning: not able to get attribute "_FillValue" of variable '//varname//' in file '//filename, &
                        kill=.false.)
        if (ierr/=NF90_NOERR) fillvalue = 0 ! default fillvalue
      end if

      ! get variable
      if (present(starts) .and. present(counts)) then
        ierr = nf90_get_var( fid, vid, vardata, start=starts, count=counts )
      else
        ierr = nf90_get_var( fid, vid, vardata )
      end if
      call nf90_check(ierr, 'Error while getting variable '//varname//'in file '//filename)

      ! close file
      ierr = nf90_close(fid)
      call nf90_check(ierr, 'Error while closing file '//filename)

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine load_netcdf_3D( filename, varname, vardata, units, fillvalue, starts, counts )
    use netcdf
    character(len=*), intent(in):: filename, varname
    double precision, intent(out), dimension(:,:,:):: vardata
    character(len=*), intent(out), optional:: units
    double precision, intent(out), optional:: fillvalue
    integer, dimension(:), intent(in), optional:: starts, counts
    integer:: fid, vid, ierr

      ! open intput file
      ierr = nf90_open( filename, NF90_NOWRITE, fid )
      call nf90_check(ierr, 'Error while openning file '//filename)

      ! get variable ID
      ierr = nf90_inq_varid( fid, varname, vid )
      call nf90_check(ierr, 'Error while inquiring ID of variable '//varname//' in file '//filename)

      ! get variable units
      if (present(units)) then
        ierr = nf90_get_att( fid, vid, 'units', units )
        call nf90_check( ierr, 'Warning: not able to get attribute "units" of variable '//varname//' in file '//filename, &
                         kill=.false. )
        if (ierr/=NF90_NOERR) units = 'unkown' ! default units
      end if

      ! get variable fillvalue
      if (present(fillvalue)) then
        ierr = nf90_get_att( fid, vid, '_FillValue', fillvalue )
        call nf90_check(ierr, 'Warning: not able to get attribute "_FillValue" of variable '//varname//' in file '//filename, &
                        kill=.false.)
        if (ierr/=NF90_NOERR) fillvalue = 0 ! default fillvalue
      end if

      ! get variable
      if (present(starts) .and. present(counts)) then
        ierr = nf90_get_var( fid, vid, vardata, start=starts, count=counts )
      else
        ierr = nf90_get_var( fid, vid, vardata )
      end if
      call nf90_check(ierr, 'Error while getting variable '//varname//'in file '//filename)

      ! close file
      ierr = nf90_close(fid)
      call nf90_check(ierr, 'Error while closing file '//filename)

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine load_netcdf_4D( filename, varname, vardata, units, fillvalue, starts, counts )
    use netcdf
    character(len=*), intent(in):: filename, varname
    double precision, intent(out), dimension(:,:,:,:):: vardata
    character(len=*), intent(out), optional:: units
    double precision, intent(out), optional:: fillvalue
    integer, dimension(:), intent(in), optional:: starts, counts
    integer:: fid, vid, ierr

      ! open intput file
      ierr = nf90_open( filename, NF90_NOWRITE, fid )
      call nf90_check(ierr, 'Error while openning file '//filename)

      ! get variable ID
      ierr = nf90_inq_varid( fid, varname, vid )
      call nf90_check(ierr, 'Error while inquiring ID of variable '//varname//' in file '//filename)

      ! get variable units
      if (present(units)) then
        ierr = nf90_get_att( fid, vid, 'units', units )
        call nf90_check( ierr, 'Warning: not able to get attribute "units" of variable '//varname//' in file '//filename, &
                         kill=.false. )
        if (ierr/=NF90_NOERR) units = 'unkown' ! default units
      end if

      ! get variable fillvalue
      if (present(fillvalue)) then
        ierr = nf90_get_att( fid, vid, '_FillValue', fillvalue )
        call nf90_check(ierr, 'Warning: not able to get attribute "_FillValue" of variable '//varname//' in file '//filename, &
                        kill=.false.)
        if (ierr/=NF90_NOERR) fillvalue = 0 ! default fillvalue
      end if

      ! get variable
      if (present(starts) .and. present(counts)) then
        ierr = nf90_get_var( fid, vid, vardata, start=starts, count=counts )
      else
        ierr = nf90_get_var( fid, vid, vardata )
      end if
      call nf90_check(ierr, 'Error while getting variable '//varname//'in file '//filename)

      ! close file
      ierr = nf90_close(fid)
      call nf90_check(ierr, 'Error while closing file '//filename)

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine put_var_real0D( ID, varid, VAR, begin, length )
    use netcdf

    integer, intent(in):: ID, varid
    real, intent(in):: VAR
    integer, dimension(:), optional:: begin, length
    integer:: ierr
    character(len=3):: num

    !==============================================================!
    if ( present(begin) .and. present(length) ) then
      ierr=nf90_put_var( ID, varid ,(/VAR/), start=begin, count=length )
    else
      ierr=nf90_put_var( ID, varid ,(/VAR/) )
    end if
    write(num,fmt="(I3.3)") varid
    call nf90_check(ierr, 'Error while putting variable #'//num)
    !==============================================================!

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine put_var_real1D( ID, varid, VAR, begin, length )
    use netcdf

    integer, intent(in):: ID, varid
    real, dimension(:), intent(in):: VAR
    integer, dimension(:), optional:: begin, length
    integer:: ierr
    character(len=3):: num

    !==============================================================!
    if ( present(begin) .and. present(length) ) then
      ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
    else
      ierr=nf90_put_var( ID, varid ,VAR )
    end if
    write(num,fmt="(I3.3)") varid
    call nf90_check(ierr, 'Error while putting variable #'//num)
    !==============================================================!

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine put_var_real2D( ID, varid, VAR, begin, length )
    use netcdf

    integer, intent(in):: ID, varid
    real, dimension(:,:), intent(in):: VAR
    integer, dimension(:), optional:: begin, length
    integer:: ierr
    character(len=3):: num

    !==============================================================!
    if ( present(begin) .and. present(length) ) then
      ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
    else
      ierr=nf90_put_var( ID, varid ,VAR )
    end if
    write(num,fmt="(I3.3)") varid
    call nf90_check(ierr, 'Error while putting variable #'//num)
    !==============================================================!

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine put_var_real3D( ID, varid, VAR, begin, length )
    use netcdf

    integer, intent(in):: ID, varid
    real, dimension(:,:,:), intent(in):: VAR
    integer, dimension(:), optional:: begin, length
    integer:: ierr
    character(len=3):: num

    !==============================================================!
    if ( present(begin) .and. present(length) ) then
      ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
    else
      ierr=nf90_put_var( ID, varid ,VAR )
    end if
    write(num,fmt="(I3.3)") varid
    call nf90_check(ierr, 'Error while putting variable #'//num)
    !==============================================================!

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine put_var_real4D( ID, varid, VAR, begin, length )
    use netcdf

    integer, intent(in):: ID, varid
    real, dimension(:,:,:,:), intent(in):: VAR
    integer, dimension(:), optional:: begin, length
    integer:: ierr
    character(len=3):: num

    !==============================================================!
    if ( present(begin) .and. present(length) ) then
      ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
    else
      ierr=nf90_put_var( ID, varid ,VAR )
    end if
    write(num,fmt="(I3.3)") varid
    call nf90_check(ierr, 'Error while putting variable #'//num)
    !==============================================================!

  end subroutine

  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine put_var_int1D( ID, varid, VAR, begin, length )
    use netcdf

    integer, intent(in):: ID, varid
    integer, dimension(:), intent(in):: VAR
    integer, dimension(:), optional:: begin, length
    integer:: ierr
    character(len=3):: num

    !==============================================================!
    if ( present(begin) .and. present(length) ) then
      ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
    else
      ierr=nf90_put_var( ID, varid ,VAR )
    end if
    write(num,fmt="(I3.3)") varid
    call nf90_check(ierr, 'Error while putting variable #'//num)
    !==============================================================!

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine put_var_int2D( ID, varid, VAR, begin, length )
    use netcdf

    integer, intent(in):: ID, varid
    integer, dimension(:,:), intent(in):: VAR
    integer, dimension(:), optional:: begin, length
    integer:: ierr
    character(len=3):: num

    !==============================================================!
    if ( present(begin) .and. present(length) ) then
      ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
    else
      ierr=nf90_put_var( ID, varid ,VAR )
    end if
    write(num,fmt="(I3.3)") varid
    call nf90_check(ierr, 'Error while putting variable #'//num)
    !==============================================================!

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine put_var_int3D( ID, varid, VAR, begin, length )
    use netcdf

    integer, intent(in):: ID, varid
    integer, dimension(:,:,:), intent(in):: VAR
    integer, dimension(:), optional:: begin, length
    integer:: ierr
    character(len=3):: num

    !==============================================================!
    if ( present(begin) .and. present(length) ) then
      ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
    else
      ierr=nf90_put_var( ID, varid ,VAR )
    end if
    write(num,fmt="(I3.3)") varid
    call nf90_check(ierr, 'Error while putting variable #'//num)
    !==============================================================!

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine put_var_int4D( ID, varid, VAR, begin, length )
    use netcdf

    integer, intent(in):: ID, varid
    integer, dimension(:,:,:,:), intent(in):: VAR
    integer, dimension(:), optional:: begin, length
    integer:: ierr
    character(len=3):: num

    !==============================================================!
    if ( present(begin) .and. present(length) ) then
      ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
    else
      ierr=nf90_put_var( ID, varid ,VAR )
    end if
    write(num,fmt="(I3.3)") varid
    call nf90_check(ierr, 'Error while putting variable #'//num)
    !==============================================================!

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine close_file(fid)
    use netcdf
    integer, intent(in):: fid
    integer:: ierr
    ierr = nf90_close(fid)
    call nf90_check(ierr, 'Error while closing output file')
  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine create_output_variable(fid, varname, dimnames, defdim, units, fillvalue, varid)
    ! Define a variable type NF90_FLOAT in an output file
    ! Expected input arguments:
    !   1) ID of netCDF output file
    !   2) name of variable
    !   3) complete list of output file dimensions names
    !   4) list of dimension on which the variable is defined (INT, 1/0)
    !   5) OPTIONAL: units of the variable
    !   6) OPTIONAL: variable fillvalue
    ! Output argument: variable ID
    ! WARNING: the output file has to be in definition mode already

    use netcdf

    integer, intent(in):: fid
    character(len=*), intent(in):: varname
    character(len=*), dimension(:):: dimnames
    integer, dimension(:), intent(in):: defdim
    character(len=*), optional:: units
    double precision, intent(in), optional:: fillvalue
    integer, intent(out):: varid
    integer, dimension(:), allocatable:: listdim, dimids
    integer:: ndim, ndefdim, k, ierr

    ndim = size(dimnames)

    allocate(listdim(ndim))

    ! find variable-defined dimensions:
    ndefdim = 0
    do k = 1,ndim
      if (defdim(k)==1) then
        ndefdim = ndefdim+1
        listdim(ndefdim)=k
      end if
    end do

    allocate(dimids(ndefdim))

    ! Get dimensions ID
    do k = 1,ndefdim
      ierr = nf90_inq_dimid(fid, dimnames(listdim(k)), dimid=dimids(k))
      call nf90_check(ierr, 'Error while inquiring ID of output file dimension '//dimnames(listdim(k)))
    end do

    ! Define variable
    ierr = nf90_def_var(fid, varname, NF90_FLOAT, dimids, varid)
    call nf90_check(ierr, 'Error in output file while defining variable '//varname)

    ! Attributes
    if (present(units)) then
      ierr = nf90_put_att(fid, varid, 'units', units)
      call nf90_check(ierr, 'Error in output file while putting attribute "units" of variable '//varname)
    end if
    if (present(fillvalue)) then
      ierr = nf90_put_att(fid, varid, '_FillValue', real(fillvalue))
      call nf90_check(ierr, 'Error in output file while putting attribute "_FillValue" of variable '//varname)
    end if


  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


  function empty_string(string)
    character(len=*), intent(in):: string
    logical:: empty_string
    integer:: n
    empty_string = .true.
    do n = 1,len(string)
      if (string(n:n)/=' ') then
        empty_string = .false.
      end if
    end do
  end function


  subroutine define_dimension(fid, dimname, dlen, dtype, varid)
    use netcdf
    integer, intent(in):: fid, dlen
    character(len=*), intent(in):: dimname
    integer, intent(in), optional:: dtype
    integer, intent(out), optional:: varid
    integer:: ierr, dimid

    ierr = nf90_def_dim(fid, dimname, dlen, dimid)
    call nf90_check(ierr, 'Error while defining dimension '//dimname)

    if (present(dtype) .and. present(varid)) then
      ierr = nf90_def_var(fid, dimname, dtype, dimid, varid)
      call nf90_check(ierr, 'Error while defining variable '//dimname)
    end if

  end subroutine


  subroutine create_output_file(filename, fid, dimnames, dimunits, &
                                x1,nx1, x2,nx2, x3,nx3, x4,nx4, x5,nx5, x6,nx6, x7,nx7, x8,nx8, x9,nx9, x10,nx10)
    ! Create output netCDF file, define dimensions and put dimension variables
    ! Expected input arguments:
    !   1) file path (char)
    !   2) vector of dimensions names (char) => determine the number of dimensions
    !   3) OPTIONAL: vector of dimensions units (char). Must be the same size than dimension names vector. The units attribute is
    !   not defined if the corresponding string is blank
    !   4) for each dimension *:
    !           'x*' double precision 1D array. In that case, the array will be dimension variable and the dimension type will be FLOAT)
    !       OR  'nx*' integer scalar. In that case, this integer is interpreted as the dimension length. dimension type wil be INTEGER
    !       OR  nothing. In that case, the dimension will be unlimited, and type FLOAT
    ! 10 dimesions at most can be created

    use netcdf
    character(len=*), intent(in):: filename
    integer, intent(out):: fid
    character(len=*), dimension(:):: dimnames
    character(len=*), dimension(:), optional:: dimunits
    double precision, dimension(:), optional:: x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
    integer, optional:: nx1, nx2, nx3, nx4, nx5, nx6, nx7, nx8, nx9, nx10
    integer:: k, ndim, ierr, varid(10)


      !=============!
      ! Create file !
      !=============!
      ierr = nf90_create( filename, NF90_CLOBBER, fid )
      call nf90_check(ierr, 'Error while creating output file '//filename)


      !===================!
      ! Define dimensions !
      !===================!
      ndim = size(dimnames)

      ! 1st dimension:
      !-----------------
      if (ndim>=1) then
        if (present(x1)) then
          call define_dimension(fid, dimnames(1), size(x1), NF90_FLOAT, varid(1))
        elseif (present(nx1)) then
          call define_dimension(fid, dimnames(1), nx1, NF90_INT, varid(1))
        else
          call define_dimension(fid, dimnames(1), NF90_UNLIMITED, NF90_FLOAT, varid(1))
        end if
        ! Axis attribute
        ierr = nf90_put_att(fid, varid(1), 'axis', 'X')
        call nf90_check(ierr, 'Error while putting attributes "axis" "X" in variable '//dimnames(1))
      end if

      ! 2nd dimension:
      !-----------------
      if (ndim>=2) then
        if (present(x2)) then
          call define_dimension(fid, dimnames(2), size(x2), NF90_FLOAT, varid(2))
        elseif (present(nx2)) then
          call define_dimension(fid, dimnames(2), nx2, NF90_INT, varid(2))
        else
          call define_dimension(fid, dimnames(2), NF90_UNLIMITED, NF90_FLOAT, varid(2))
        end if
        ! Axis attribute
        ierr = nf90_put_att(fid, varid(2), 'axis', 'Y')
        call nf90_check(ierr, 'Error while putting attributes "axis" "Y" in variable '//dimnames(2))
      end if

      ! 3rd dimension:
      !-----------------
      if (ndim>=3) then
        if (present(x3)) then
          call define_dimension(fid, dimnames(3), size(x3), NF90_FLOAT, varid(3))
        elseif (present(nx3)) then
          call define_dimension(fid, dimnames(3), nx3, NF90_INT, varid(3))
        else
          call define_dimension(fid, dimnames(3), NF90_UNLIMITED, NF90_FLOAT, varid(3))
        end if
        ! Axis attribute
        ierr = nf90_put_att(fid, varid(3), 'axis', 'Z')
        call nf90_check(ierr, 'Error while putting attributes "axis" "Z" in variable '//dimnames(3))
      end if

      ! 4th dimension:
      !-----------------
      if (ndim>=4) then
        if (present(x4)) then
          call define_dimension(fid, dimnames(4), size(x4), NF90_FLOAT, varid(4))
        elseif (present(nx4)) then
          call define_dimension(fid, dimnames(4), nx4, NF90_INT, varid(4))
        else
          call define_dimension(fid, dimnames(4), NF90_UNLIMITED, NF90_FLOAT, varid(4))
        end if
      end if

      ! 5th dimension:
      !-----------------
      if (ndim>=5) then
        if (present(x5)) then
          call define_dimension(fid, dimnames(5), size(x5), NF90_FLOAT, varid(5))
        elseif (present(nx5)) then
          call define_dimension(fid, dimnames(5), nx5, NF90_INT, varid(5))
        else
          call define_dimension(fid, dimnames(5), NF90_UNLIMITED, NF90_FLOAT, varid(5))
        end if
      end if

      ! 6th dimension:
      !-----------------
      if (ndim>=6) then
        if (present(x6)) then
          call define_dimension(fid, dimnames(6), size(x6), NF90_FLOAT, varid(6))
        elseif (present(nx6)) then
          call define_dimension(fid, dimnames(6), nx6, NF90_INT, varid(6))
        else
          call define_dimension(fid, dimnames(6), NF90_UNLIMITED, NF90_FLOAT, varid(6))
        end if
      end if

      ! 7th dimension:
      !-----------------
      if (ndim>=7) then
        if (present(x7)) then
          call define_dimension(fid, dimnames(7), size(x7), NF90_FLOAT, varid(7))
        elseif (present(nx7)) then
          call define_dimension(fid, dimnames(7), nx7, NF90_INT, varid(7))
        else
          call define_dimension(fid, dimnames(7), NF90_UNLIMITED, NF90_FLOAT, varid(7))
        end if
      end if

      ! 8th dimension:
      !-----------------
      if (ndim>=8) then
        if (present(x8)) then
          call define_dimension(fid, dimnames(8), size(x8), NF90_FLOAT, varid(8))
        elseif (present(nx8)) then
          call define_dimension(fid, dimnames(8), nx8, NF90_INT, varid(8))
        else
          call define_dimension(fid, dimnames(8), NF90_UNLIMITED, NF90_FLOAT, varid(8))
        end if
      end if

      ! 9th dimension:
      !-----------------
      if (ndim>=9) then
        if (present(x9)) then
          call define_dimension(fid, dimnames(9), size(x9), NF90_FLOAT, varid(9))
        elseif (present(nx9)) then
          call define_dimension(fid, dimnames(9), nx9, NF90_INT, varid(9))
        else
          call define_dimension(fid, dimnames(9), NF90_UNLIMITED, NF90_FLOAT, varid(9))
        end if
      end if

      ! 10th dimension:
      !-----------------
      if (ndim>=10) then
        if (present(x10)) then
          call define_dimension(fid, dimnames(10), size(x10), NF90_FLOAT, varid(10))
        elseif (present(nx10)) then
          call define_dimension(fid, dimnames(10), nx10, NF90_INT, varid(10))
        else
          call define_dimension(fid, dimnames(10), NF90_UNLIMITED, NF90_FLOAT, varid(10))
        end if
      end if

      !=================!
      ! Units attribute !
      !=================!
      if (present(dimunits)) then

        do k = 1,ndim
          if (.not. empty_string(dimunits(k))) then
            ierr = nf90_put_att(fid, varid(k), 'units', dimunits(k))
            call nf90_check(ierr, 'Error while putting attributes "units" '//dimunits(1)//' in variable '//dimnames(k))
          end if
        end do

      end if

      !===============!
      ! Put variables !
      !===============!
      ierr = nf90_enddef(fid)
      call nf90_check(ierr, 'Error while end of definition of output file '//filename)

      ! 1st dimension:
      !---------------
      if (ndim>=1) then
        if (present(x1)) then
          call put_var_real1D(fid, varid(1), real(x1))
        elseif (present(nx1)) then
          call put_var_int1D(fid, varid(1), (/(k,k=1,nx1)/))
        end if
      end if

      ! 2nd dimension:
      !---------------
      if (ndim>=2) then
        if (present(x2)) then
          call put_var_real1D(fid, varid(2), real(x2))
        elseif (present(nx2)) then
          call put_var_int1D(fid, varid(2), (/(k,k=1,nx2)/))
        end if
      end if

      ! 3rd dimension:
      !---------------
      if (ndim>=3) then
        if (present(x3)) then
          call put_var_real1D(fid, varid(3), real(x3))
        elseif (present(nx3)) then
          call put_var_int1D(fid, varid(3), (/(k,k=1,nx3)/))
        end if
      end if

      ! 4th dimension:
      !---------------
      if (ndim>=4) then
        if (present(x4)) then
          call put_var_real1D(fid, varid(4), real(x4))
        elseif (present(nx4)) then
          call put_var_int1D(fid, varid(4), (/(k,k=1,nx4)/))
        end if
      end if

      ! 5th dimension:
      !---------------
      if (ndim>=5) then
        if (present(x5)) then
          call put_var_real1D(fid, varid(5), real(x5))
        elseif (present(nx5)) then
          call put_var_int1D(fid, varid(5), (/(k,k=1,nx5)/))
        end if
      end if

      ! 6th dimension:
      !---------------
      if (ndim>=6) then
        if (present(x6)) then
          call put_var_real1D(fid, varid(6), real(x6))
        elseif (present(nx6)) then
          call put_var_int1D(fid, varid(6), (/(k,k=1,nx6)/))
        end if
      end if

      ! 7th dimension:
      !---------------
      if (ndim>=7) then
        if (present(x7)) then
          call put_var_real1D(fid, varid(7), real(x7))
        elseif (present(nx7)) then
          call put_var_int1D(fid, varid(7), (/(k,k=1,nx7)/))
        end if
      end if

      ! 8th dimension:
      !---------------
      if (ndim>=8) then
        if (present(x8)) then
          call put_var_real1D(fid, varid(8), real(x8))
        elseif (present(nx8)) then
          call put_var_int1D(fid, varid(8), (/(k,k=1,nx8)/))
        end if
      end if

      ! 9th dimension:
      !---------------
      if (ndim>=9) then
        if (present(x9)) then
          call put_var_real1D(fid, varid(9), real(x9))
        elseif (present(nx9)) then
          call put_var_int1D(fid, varid(9), (/(k,k=1,nx9)/))
        end if
      end if

      ! 10th dimension:
      !---------------
      if (ndim>=10) then
        if (present(x10)) then
          call put_var_real1D(fid, varid(10), real(x10))
        elseif (present(nx10)) then
          call put_var_int1D(fid, varid(10), (/(k,k=1,nx10)/))
        end if
      end if


      !=========================!
      ! Back to definition mode !
      !=========================!
      ierr = nf90_redef(fid)
      call nf90_check(ierr, 'Error while getting back to definition mode in output file '//filename)


  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!
  !-------------------------------------------------------------------------------------------------------------------------------!


  subroutine nf90_check(ierr, message, kill)

    use netcdf
    integer, intent(in):: ierr
    character(len=*), intent(in):: message
    logical, optional, intent(in):: kill
    logical:: loc_kill

    if (present(kill)) then
      loc_kill = kill
    else
      loc_kill = .true.
    end if

    if (ierr/=NF90_NOERR) then
      print *
      print *, message
      print *, nf90_strerror(ierr)
      if (loc_kill) stop
    end if

  end subroutine


  !-------------------------------------------------------------------------------------------------------------------------------!


end module
