module io_module
  implicit none

  contains

    subroutine make_input_output( io_fname, cell_area,land_area,temp,runoff,slope,lith_frac, CaMg_rock, CO2_levels, &
                                  curr_temp, curr_runf, h, E, xs, W, WCaMg, &
                                  nlon,nlat,nlith,nrun, ForwBckw, ncout_ID )
      ! Main subroutine. Read input-output file, allocate variables, load variables and create output file

      use netcdf
      use netcdf_io_functions, only: netcdf_get_size, create_output_file, load_netcdf_1D, load_netcdf_2D, &
                                     load_netcdf_3D, nf90_check
      use ascii_io_functions, only: file_length

      ! internal file units (parameter file, forcing file and output variables ID scratch file) and fillvalue:
      ! define: IPARAM, IFORC, IOUT and DEFFILLVAL
      include 'common_parameters.inc'
      
      ! in/out variables
      character(len=*), intent(in):: io_fname
      double precision, intent(inout), dimension(:,:), allocatable:: cell_area, land_area, slope, curr_temp, curr_runf, h, E, xs, &
                                                                     W, WCaMg
      double precision, intent(inout), dimension(:,:,:), allocatable:: temp, runoff, lith_frac
      double precision, intent(inout), dimension(:), allocatable:: CaMg_rock, CO2_levels
      integer, intent(out):: nlon, nlat, nlith, nrun, ForwBckw, ncout_ID

      ! local variables
      double precision, dimension(:), allocatable:: lon, lat
      double precision:: T_fillval, R_fillval, S_fillval
      character(len=50):: varname, expect_units, units, dimname(10)
      character(len=50):: axunits(3), nounits
      character(len=500):: path
      integer:: nCO2

      ! dynsoil parameter variables
      include 'dynsoil_physical_parameters.inc'

      ! technical variables
      double precision:: forcing
      logical:: multirun
      integer:: ndim, ierr, k, length

      do k=1,len(nounits)
        nounits(k:k) = ' '
      end do


      !!=====================================!!
      !!  open input-output conditions file  !!
      !!=====================================!!

      open(unit=1, file=io_fname, status='old', action='read')




      !!===========================================================!!
      !!  read "input fields" variables names and load variables:  !!
      !!===========================================================!!


      !-----------------!
      ! Total cell area !
      !-----------------!

      ! read names
      call read_comment(1)
      read(unit=1, fmt=*) varname, path, expect_units, dimname(1:2)

      ! Sizing
      !**************************************!
      nlon = netcdf_get_size(path, dimname(1))
      nlat = netcdf_get_size(path, dimname(2))
      allocate( lon(nlon) )
      allocate( lat(nlat) )
      allocate( cell_area(nlon,nlat) )
      !**************************************!

      ! get variable
      call load_netcdf_1D(path, dimname(1), lon, units=axunits(1))
      call load_netcdf_1D(path, dimname(2), lat, units=axunits(2))
      !
      call load_netcdf_2D(path, varname, cell_area, units=units)
      call check_units(varname, expect_units, units)


      !----------------!
      ! Continent area !
      !----------------!

      ! read names
      call read_comment(1)
      read(unit=1, fmt=*) varname, path, expect_units, dimname(1:2)

      ! Sizing
      call check_sizing( path, varname, dimname(1:2), (/nlon,nlat/) )
      !******************************!
      allocate( land_area(nlon,nlat) )
      !******************************!

      ! get variable
      call check_coordinates(path, varname, dimname(1:2), lon, lat)
      call load_netcdf_2D(path, varname, land_area, units=units)
      call check_units(varname, expect_units, units)



      !-------------!
      ! Temperature !
      !-------------!

      ! read names
      call read_comment(1)
      read(unit=1, fmt=*) varname, path, expect_units, dimname(1:3)

      ! Sizing
      call check_sizing( path, varname, dimname(1:2), (/nlon,nlat/) )
      nCO2 = netcdf_get_size(path, dimname(3))
      !******************************!
      allocate( CO2_levels(nCO2) )
      allocate( temp(nlon,nlat,nCO2) )
      !******************************!

      ! get variable
      call load_netcdf_1D(path, dimname(3), CO2_levels, units=axunits(3))
      !
      call check_coordinates(path, varname, dimname(1:2), lon, lat)
      call load_netcdf_3D(path, varname, temp, units=units, fillvalue=T_fillval)
      call check_units(varname, expect_units, units)



      !--------!
      ! Runoff !
      !--------!

      ! read names
      call read_comment(1)
      read(unit=1, fmt=*) varname, path, expect_units, dimname(1:3)

      ! Sizing
      call check_sizing( path, varname, dimname(1:3), (/nlon,nlat,nCO2/) )
      !********************************!
      allocate( runoff(nlon,nlat,nCO2) )
      !********************************!

      ! get variable
      call check_coordinates(path, varname, dimname(1:3), lon, lat, CO2_levels)
      call load_netcdf_3D(path, varname, runoff, units=units, fillvalue=R_fillval)
      call check_units(varname, expect_units, units)



      !-------!
      ! Slope !
      !-------!

      ! read names
      call read_comment(1)
      read(unit=1, fmt=*) varname, path, expect_units, dimname(1:2)

      ! Sizing
      call check_sizing( path, varname, dimname(1:2), (/nlon,nlat/) )
      !**************************!
      allocate( slope(nlon,nlat) )
      !**************************!

      ! get variable
      call check_coordinates(path, varname, dimname(1:2), lon, lat)
      call load_netcdf_2D(path, varname, slope, units=units, fillvalue=S_fillval)
      call check_units(varname, expect_units, units)



      !--------------------!
      ! lithology fraction !
      !--------------------!

      ! read names
      call read_comment(1)
      read(unit=1, fmt=*) varname, path, expect_units, dimname(1:3)

      ! Sizing
      call check_sizing( path, varname, dimname(1:2), (/nlon,nlat/) )
      !*******************************************!
      nlith = netcdf_get_size(path, dimname(3)) - 1 ! Ignore the first 'lithology' class: water/ice
      allocate( lith_frac(nlon,nlat,nlith) )
      !*******************************************!

      ! get variable
      call check_coordinates(path, varname, dimname(1:2), lon, lat)
      call load_netcdf_3D(path, varname, lith_frac, units=units, starts=(/1,1,2/), counts=(/nlon,nlat,nlith/))
      call check_units(varname, expect_units, units)



      !!===============================!!
      !!  Allocate internal variables  !!
      !!===============================!!

      allocate( curr_temp(nlon,nlat) )
      allocate( curr_runf(nlon,nlat) )
      allocate(         h(nlon,nlat) )
      allocate(         E(nlon,nlat) )
      allocate(        xs(nlon,nlat) )
      allocate(         W(nlon,nlat) )
      allocate(     WCaMg(nlon,nlat) )
      allocate( CaMg_rock(nlith)     )

      ! fillvalues
      curr_temp = DEFFILLVAL
      curr_runf = DEFFILLVAL
      h         = DEFFILLVAL
      E         = DEFFILLVAL
      xs        = DEFFILLVAL
      W         = DEFFILLVAL
      WCaMg     = DEFFILLVAL
      where(   temp==T_fillval )  temp   = DEFFILLVAL
      where( runoff==R_fillval )  runoff = DEFFILLVAL
      where(  slope==S_fillval )  slope  = DEFFILLVAL



      !!===========================================================================!!
      !!  Check consistency of land area, missing-value and lithological fraction  !!
      !!===========================================================================!!

      call check_continental_cells( cell_area, land_area, temp, runoff, slope, lith_frac )



      !!==================================================!!
      !!  read other variables names and load variables:  !!
      !!==================================================!!


      !---------------------------!
      ! Single run/Multirun flag: !
      !---------------------------!

      call read_comment(1)
      read(unit=1, fmt=*) multirun


      !-----------------------!
      ! Forward/Backward run: !
      !-----------------------!

      call read_comment(1)
      read(unit=1, fmt=*) ForwBckw


      !-------------!
      ! Parameters: !
      !-------------!

      call read_comment(1)
      read(unit=1, fmt=*) path
      open(unit=IPARAM, file=path, status='old', action='read')

      if (multirun) then
              
        call read_comment(IPARAM)
        nrun = file_length(fileunit=IPARAM)
        
        ! Check: read all the parameterizations:
        call read_comment(IPARAM)
        do k = 1,nrun
          read(unit=IPARAM, fmt=*)  ke, a, b, krp, Ea_rp, T0_rp, h0, kd, kw, Ea, T0, sigma, CaMg_rock
        end do

        rewind(unit=IPARAM)
        call read_comment(IPARAM)
        ! parameters will be read step by step in the main file

      else

        nrun = 0

        ! Get parameters values:
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) ke
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) a
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) b
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) krp
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) Ea_rp
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) T0_rp
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) h0
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) kd
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) kw
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) Ea
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) T0
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) sigma
        call read_comment(IPARAM)
        read(unit=IPARAM, fmt=*) CaMg_rock

        close(unit=IPARAM)

      end if


      !---------------!
      ! Forcing file: !
      !---------------!

      call read_comment(1)
      read(unit=1, fmt=*) path

      open(unit=IFORC, file=path, status='old', action='read')
      call read_comment(IFORC)

      if (multirun) then

        length = file_length(fileunit=IFORC)
        call read_comment(IFORC)

        ! Compare the length of forcing file and number of runs (from parameters file)
        if (length/=nrun) then

          ! In case of mismatch, if the length of forcing file is 1, create a "fake" scratch forcing file by putting the single
          ! forcing for each parameterization (nrun)
          if (length==1) then
            read(unit=IFORC, fmt=*) forcing
            close(unit=IFORC)
            open(unit=IFORC, status='scratch', action='readwrite')
            do k = 1,nrun
              write(unit=IFORC, fmt=*) forcing
            end do
            rewind(unit=IFORC)

            ! Otherwise, raise error
          else
            print *, 'ERROR: inconsistent nubmer of parameterizations and forcings. File lengtha mismatch'
            stop
          end if

        end if
      end if



      !!======================================================================!!
      !!  read output file and variables names and create output netCDF file  !!
      !!======================================================================!!


      !----------------!
      ! read file name !
      !----------------!
      call read_comment(1)
      read(unit=1, fmt=*) path

      !---------------------------------------------------!
      ! read dimensions name. 4 dimensions area expected: !
      !---------------------------------------------------!
      ! longitude
      call read_comment(1)
      read(unit=1, fmt=*) dimname(1)
      ! latitude
      call read_comment(1)
      read(unit=1, fmt=*) dimname(2)
      ! lithology
      call read_comment(1)
      read(unit=1, fmt=*) dimname(3)
      ! model runs (this line MUST be present, but the value is ignore if multirun=.false.)
      call read_comment(1)
      read(unit=1, fmt=*) dimname(4)


      !--------------------!
      ! create output file !
      !--------------------!
      if (multirun) then
        call create_output_file(path, ncout_ID, dimname(1:4), (/axunits(1),axunits(2),nounits,nounits/), x1=lon, x2=lat, nx3=nlith)
      else
        call create_output_file(path, ncout_ID, dimname(1:3), (/axunits(1),axunits(2),nounits/), x1=lon, x2=lat, nx3=nlith)
      end if


      !-------------------------!
      ! Define output variables !
      !-------------------------!

      ! Scratch file recording the status of output variables (ID, write or not...)
      open(unit=IOUT, status='scratch', action='readwrite')

      if (multirun) then
        ndim = 4
      else
        ndim = 3
      end if

      ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !
      call define_variables( ncout_ID, dimname(1:ndim), DEFFILLVAL, axunits, multirun )
      ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !


      !--------------------!
      ! End of define mode !
      !--------------------!
      ierr = nf90_enddef(ncout_ID)
      call nf90_check(ierr, 'Error while end of definition mode in output file '//path)



      !!======================================!!
      !!  close input-output conditions file  !!
      !!======================================!!

      close(unit=1)



    end subroutine



    !-------------------------------------------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------------------------------------------!



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


    !-------------------------------------------------------------------------------------------------------------------------------!


    subroutine define_variables( fid, dimname, fillval, axunits, multirun )

      use netcdf_io_functions, only: create_output_variable

      include 'common_parameters.inc' ! define IOUT, IPARAM and IFORC
      integer, intent(in):: fid
      character(len=*), dimension(:), intent(in):: dimname, axunits
      double precision, intent(in):: fillval
      logical, intent(in):: multirun
      character(len=100):: varname
      character(len=50):: units
      integer, dimension(:), allocatable:: defdim
      logical:: write_var
      integer:: varid, k, n, ierr


      n = size(dimname)
      allocate(defdim(n))


      call read_comment(1, ierr=ierr)

      do while(ierr==0) ! read until the end of the interface file!

        read(unit=1, fmt=*) varname, units, write_var, defdim

        ! A given units $* means: use the units of the axis number *
        if (units(1:1)=='$') then
          read(units(2:50), fmt=*) k
          units = axunits(k)
        end if

        if (write_var) then
          call create_output_variable(fid, varname, dimname, defdim, units, fillval, varid)
        else
          varid = 0
          defdim(n) = 0 ! make sure the two booleans in IOUT file will be .false. and the program will not try to write the variable
        end if

        if (multirun) then
          write(unit=IOUT, fmt=*) write_var, (defdim(n)==1), varid
        else
          write(unit=IOUT, fmt=*) write_var, varid
        end if

        call read_comment(1, ierr=ierr)

      end do


      if (ierr>0) then ! Error other than end-of-file has been raised
        backspace(1)
        call read_comment(1)
      end if


      deallocate(defdim)


    end subroutine


    !-------------------------------------------------------------------------------------------------------------------------------!


    subroutine check_sizing( path, varname, dimnames, sizes )

      use netcdf_io_functions, only: netcdf_get_size

      character(len=*), intent(in):: path, varname
      character(len=*), dimension(:), intent(in):: dimnames
      integer, dimension(:):: sizes
      integer:: n, l

      do n = 1,size(dimnames)

        l = netcdf_get_size(path, dimnames(n))

        if (l/=sizes(n)) then

          print *, 'ERROR: size mismatch of dimension '//dimnames(n)//' of variable '//varname
          print *, '  Expected: ',sizes(n)
          print *, '       Got: ',l
          stop

        end if

      end do

    end subroutine


    !-------------------------------------------------------------------------------------------------------------------------------!


    subroutine check_units( varname, expected_units, units )

      use netcdf_io_functions, only: empty_string

      character(len=*), intent(in):: varname, expected_units, units
      integer:: l0, l
      logical:: success

      if (expected_units(1:12)/='DO_NOT_CHECK') then

        success = .false.
        l0 = 1

        ! search for "OR" separator, that are "|" character
        do l = 1,len(expected_units)

          if (expected_units(l:l)=='|') then
            if ( (units(1:l-l0) == expected_units(l0:l-1)) .and. empty_string(units(l-l0+1:)) ) then
              success = .true. ! units matches one of the expected units
            end if
          l0 = l+1
          end if

        end do

        ! if no match was found
        if (.not. success) then

          ! try to match the rest of the string (or the entire string if no | was found) 
          if ( units(1:l-l0+1) /= expected_units(l0:) ) then ! if not, raise error
            print *, 'ERROR: incorrect units of input variable '//varname
            print *, '  Expected: '//expected_units
            print *, '       Got: '//units
            print*, 'Note: | characters are interpreted as "or" separators'
            stop
          end if

        end if

      end if

    end subroutine


    !-------------------------------------------------------------------------------------------------------------------------------!


    subroutine check_coordinates_singlevar( path, varname, dimname, x0 )

      use netcdf_io_functions, only: load_netcdf_1D

      include 'common_parameters.inc' !to get the value of MAX_ALLOWED_INACC

      character(len=*), intent(in):: path, varname, dimname
      double precision, dimension(:), intent(in):: x0
      double precision, dimension(:), allocatable:: x

      allocate(x(size(x0)))
      call load_netcdf_1D(path, dimname, x)
      if ( maxval(abs((x-x0)/x0)) > MAX_ALLOWED_INACC ) then
        print *, 'ERROR: coordinates mismatch found for variable '//varname
        stop
      end if
      deallocate(x)

    end subroutine

    !-----------------------!

    subroutine check_coordinates( path, varname, dimnames, x1, x2, x3, x4, x5, x6, x7 )

      use netcdf_io_functions, only: load_netcdf_1D

      character(len=*), intent(in):: path, varname
      character(len=*), dimension(:), intent(in):: dimnames
      double precision, dimension(:), intent(in), optional:: x1, x2, x3, x4, x5, x6, x7

      if (present(x1)) call check_coordinates_singlevar( path, varname, dimnames(1), x1 )
      if (present(x2)) call check_coordinates_singlevar( path, varname, dimnames(2), x2 )
      if (present(x3)) call check_coordinates_singlevar( path, varname, dimnames(3), x3 )
      if (present(x4)) call check_coordinates_singlevar( path, varname, dimnames(4), x4 )
      if (present(x5)) call check_coordinates_singlevar( path, varname, dimnames(5), x5 )
      if (present(x6)) call check_coordinates_singlevar( path, varname, dimnames(6), x6 )
      if (present(x7)) call check_coordinates_singlevar( path, varname, dimnames(7), x7 )

    end subroutine


    !-------------------------------------------------------------------------------------------------------------------------------!


    function check_continental_cells_single2Dvar( varname, land_area, var )

      include 'common_parameters.inc' ! to get the value of DEFFILLVAL

      logical:: check_continental_cells_single2Dvar
      character(len=*), intent(in):: varname
      double precision, dimension(:,:), intent(in):: land_area, var
      double precision:: area_err, tot_land_area
      integer:: nerr

      nerr = count( (land_area>0 .and. var==DEFFILLVAL) )
      area_err = sum( land_area, mask=(land_area>0 .and. var==DEFFILLVAL) )
      tot_land_area = sum(land_area)

      if (nerr > 0) then
        print *, 'WARNING: found missing values on continental cells of variable '//varname
        print *, 'Number of continent cells with missing values:     ', nerr
        print *, 'Fraction of continental cells with missing values: ', dble(nerr)/dble(count((land_area>0)))
        print *, 'Total area of those cells (m2):                    ', area_err
        print *, 'Which is a fraction of total land area:            ', area_err/tot_land_area
        print *
        check_continental_cells_single2Dvar = .false.
      else
        check_continental_cells_single2Dvar = .true.
      end if

    end function

    !-----------------------!

    function check_continental_cells_single3Dvar( varname, land_area, var )

      include 'common_parameters.inc' ! to get the value of DEFFILLVAL

      logical:: check_continental_cells_single3Dvar
      character(len=*), intent(in):: varname
      double precision, dimension(:,:), intent(in):: land_area
      double precision, dimension(:,:,:), intent(in):: var
      double precision:: area_err, tot_land_area
      integer:: k, nerr

      check_continental_cells_single3Dvar = .true.
      tot_land_area = sum(land_area)

      do k = 1,size( var, 3 )

        nerr = count( (land_area>0 .and. var(:,:,k)==DEFFILLVAL) )
        area_err = sum( land_area, mask=(land_area>0 .and. var(:,:,k)==DEFFILLVAL) )

        if (nerr > 0) then
          print *, 'WARNING: found missing values on continental cells of variable '//varname
          print *, '"Vertical" (3rd axis) level number: ', k
          print *, 'Number of continent cells with missing values:     ', nerr
          print *, 'Fraction of continental cells with missing values: ', dble(nerr)/dble(count((land_area>0)))
          print *, 'Total area of those cells (m2):                    ', area_err
          print *, 'Which is a fraction of total land area:            ', area_err/tot_land_area
          print *
          check_continental_cells_single3Dvar = .false.
        end if

      end do

    end function

    !-----------------------!

    function check_continental_lithology( cell_area, land_area, lith_frac )

      include 'common_parameters.inc' !to get the value of MAX_ALLOWED_INACC

      logical:: check_continental_lithology
      double precision, dimension(:,:), intent(in):: cell_area, land_area
      double precision, dimension(:,:,:), intent(in):: lith_frac
      double precision:: area_err, tot_land_area, areadiff, max_area_diff
      integer:: nerr,i,j

      nerr = 0
      max_area_diff = 0d0
      area_err = 0d0
      tot_land_area = sum(land_area)

      do j=1,size(cell_area,2)
        do i=1,size(cell_area,1)

          if (land_area(i,j)>0) then

            areadiff = abs(  sum(lith_frac(i,j,:)) * cell_area(i,j)  -  land_area(i,j)  )   /   land_area(i,j)

            if (areadiff > MAX_ALLOWED_INACC) then
              nerr = nerr+1
              area_err = area_err + land_area(i,j)
              if (areadiff > max_area_diff) max_area_diff = areadiff
            end if

          end if

        end do
      end do

      if (nerr > 0) then
        print *, 'WARNING: inconsistency of lithogolical class fraction and continental area.'
        print *, 'Found cells where the sum of all lithology fractions times the cell area differs from the cell land area.'
        print *, 'Number of cells affected:               ', nerr
        print *, 'Fraction of continental cells affected: ', dble(nerr)/dble(count((land_area>0)))
        print *, 'Total area of those cells (m2):         ', area_err
        print *, 'Which is a fraction of total land area: ', area_err/tot_land_area
        print *, 'Maximum relative difference found:      ', max_area_diff
        print *
        check_continental_lithology = .false.
      else
        check_continental_lithology = .true.
      end if

    end function

    !-----------------------!

    subroutine check_continental_cells( cell_area, land_area, temp, runoff, slope, lith_frac )

      double precision, dimension(:,:), intent(in):: cell_area, land_area, slope
      double precision, dimension(:,:,:), intent(in):: temp, runoff, lith_frac
      logical, dimension(4):: checkpoints
      integer:: kill, ierr
      character(len=1):: killchar

      checkpoints(1) = check_continental_cells_single3Dvar( 'temperature', land_area, temp   )
      checkpoints(2) = check_continental_cells_single3Dvar( 'runoff',      land_area, runoff )
      checkpoints(3) = check_continental_cells_single2Dvar( 'slope',       land_area, slope  )
      checkpoints(4) = check_continental_lithology( cell_area, land_area, lith_frac )

      ! Kill the program if the checkpoints are not validated
      !if (.not. all(checkpoints)) stop

      ! Other option: ask for killing (with the possibility to call the program with external argument 0|1):
      if (.not. all(checkpoints)) then

        print *, 'Do you want to continue running the program? (1:yes, 0:no)'

        ! read potential argument
        call get_command_argument(1, killchar, status=ierr)
        read(killchar, fmt=*, iostat=ierr) kill
        if (ierr/=0) kill=-1

        ! reading loop
        do while (kill/=0 .and. kill/=1)
          read(unit=*, fmt=*) kill
        end do

        if (kill==0) stop

      end if

    end subroutine


    !-------------------------------------------------------------------------------------------------------------------------------!



end module
