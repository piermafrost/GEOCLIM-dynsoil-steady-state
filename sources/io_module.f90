module io_module
  implicit none

  contains

    subroutine make_input_output( io_fname, cell_area,land_area,temp,runoff,slope,lith_frac, CaMg_rock, CO2_levels, &
                                  curr_temp, curr_runf, h, E, xs, W, WCaMg, &
                                  nlon,nlat,nlith,nrun, ForwBckw, ncout_ID )
      ! Main subroutine. Read input-output file, allocate variables, load variables and create output file

      use netcdf
      use netcdf_io_functions, only: netcdf_get_size, create_output_file, nf90_check
      use ascii_io_functions, only: file_length, read_comment

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
      double precision, dimension(:), allocatable:: lon, lat, uniform_lithfrac
      double precision:: T_fillval, R_fillval, S_fillval, L_fillval, xi
      character(len=1000):: line, varname, varname2
      character(len=50):: dimname(10)
      character(len=50):: axunits(2), nounits
      character(len=1000):: path
      character(len=500), dimension(:), allocatable:: multipath
      integer:: nCO2

      ! dynsoil parameter variables
      include 'dynsoil_physical_parameters.inc'

      ! technical variables
      double precision:: forcing
      logical:: multirun, fixed_CO2
      integer:: ndim, ierr, k, length

      nounits = ''


      print *
      print *
      print *, 'Load input data'
      print *, '==============='
      print *


      !!=====================================!!
      !!  open input-output conditions file  !!
      !!=====================================!!

      open(unit=1, file=io_fname, status='old', action='read')




      !!===========================================================!!
      !!  read "input fields" variables names and load variables:  !!
      !!===========================================================!!


      !----------!
      ! CO2 axis !
      !----------!

      print *
      print *, 'Read CO2 axis'

      call read_comment(1)
      read(unit=1, fmt=*) path
      
      ! Try open taken file name
      call add_rootpath(path)
      open(unit=2, file=path, status='old', action='read', iostat=ierr)

      if (ierr==0) then

        ! Read CO2 axis from file
        print *, '--- read from separate file: '//trim(path)
        call read_comment(2)
        nCO2 = file_length(fileunit=2)
        !****************************!
        allocate( CO2_levels(nCO2+2) )
        allocate( multipath(nCO2) )
        !****************************!
        call read_comment(2)
        print *, nCO2
        read(unit=2, fmt=*) CO2_levels(2:nCO2+1)
        close(unit=2)


      else ! Try read CO2 axis directly from main file

        print *, '--- direct reading from main IO file'
        backspace(unit=1)
        read(unit=1, fmt='(A)') line
        nCO2 = 1
        do k = 1,len(line)
          if (line(k:k)==',') nCO2 = nCO2 + 1
        end do
        !****************************!
        allocate( CO2_levels(nCO2+2) )
        allocate( multipath(nCO2) )
        !****************************!
        read(line, fmt=*) CO2_levels(2:nCO2+1)

      end if

      fixed_CO2 = (nCO2==1)



      !-----------------!
      ! Total cell area !
      !-----------------!

      print *
      print *, 'load total cell area'

      ! read names
      call read_comment(1)
      read(unit=1, fmt=*) path, dimname(1:2), varname
      call add_rootpath(path)

      ! Sizing
      !**************************************!
      nlon = netcdf_get_size(path, dimname(1))
      nlat = netcdf_get_size(path, dimname(2))
      allocate( lon(nlon) )
      allocate( lat(nlat) )
      allocate( cell_area(nlon,nlat) )
      !**************************************!

      ! get variable (and variable axis) + check units
      call load_variable('area', varname, dimname(1), dimname(2), single_input_file=path, &
                         varout2D=cell_area, x=lon, y=lat, fill_missing=.true.            )



      !----------------!
      ! Continent area !
      !----------------!

      print *
      print *, 'load cell land fraction/area'

      ! read names
      call read_comment(1)
      read(unit=1, fmt=*) path, dimname(1:2), varname
      call add_rootpath(path)

      ! Sizing
      call check_sizing( path, dimname(1:2), (/nlon,nlat/) )
      !******************************!
      allocate( land_area(nlon,nlat) )
      !******************************!

      ! get variable + check coordinates and units
      call load_variable('landarea', varname, dimname(1), dimname(2), single_input_file=path,           &
                         varout2D=land_area, xref=lon, yref=lat, totarea=cell_area, fill_missing=.true. )



      !-------------------------!
      ! GCM outputs information !
      !-------------------------!

      print *
      print *, 'Read GCM output files:'
      print *, '    * variables and dimension names'

      ! read dimension names
      call read_comment(1)
      read(unit=1, fmt=*) dimname(1)
      call read_comment(1)
      read(unit=1, fmt=*) dimname(2)

      ! read temperature and runoff variable names
      call read_comment(1)
      read(unit=1, fmt='(A)') varname
      call read_comment(1)
      read(unit=1, fmt='(A)') varname2

      print *, '    * files names for each CO2 level'

      ! read multiple GCM input file names
      !
      !   read START tag
      call read_comment(1)
      read(unit=1, fmt=*) line
      if (line /= '<<--START-->>') then
        print *
        print *, 'ERROR: Expected "START" tag not found while reading GCM output file names'
        stop
      end if
      !
      !   pre-read files names
      open(unit=333, status='scratch', action='readwrite')
      k = 0
      call read_comment(1)
      read(unit=1, fmt='(A)') line
      do while (line /= '<<--STOP-->>')
        k = k + 1
        write(unit=333, fmt='(A)') line
        call read_comment(1)
        read(unit=1, fmt='(A)') line
      end do
      !
      ! Check consistency between number of files and CO2 axis lenth
      if (k/=nCO2) then
        print *
        print *, 'ERROR: number of input GCM files differs from size of CO2 axis'
        stop
      end if
      !
      !   read file names
      rewind(unit=333)
      do k = 1,nCO2
        read(unit=333, fmt=*) multipath(k)
        call add_rootpath(multipath(k))
      end do
      close(unit=333)


      print *, '    * load data'


      ! Sizing
      do k = 1,nCO2
        call check_sizing(multipath(k), dimname(1:2), (/nlon,nlat/) )
      end do


      ! Temperature
      !------------

      print *
      print *, '          - temperature'

      !********************************!
      allocate( temp(nlon,nlat,nCO2+2) )
      !********************************!

      ! get variable + check coordinates and units
      call load_variable('temperature', varname, dimname(1), dimname(2), multiple_input_file=multipath, &
                         varout3D=temp(:,:,2:nCO2+1), xref=lon, yref=lat, fillvalue=T_fillval           )


      ! Runoff
      !-------

      print *
      print *, '          - runoff'

      !**********************************!
      allocate( runoff(nlon,nlat,nCO2+2) )
      !**********************************!

      ! get variable + check coordinates and units
      call load_variable('runoff', varname2, dimname(1), dimname(2), multiple_input_file=multipath, &
                         varout3D=runoff(:,:,2:nCO2+1), xref=lon, yref=lat, fillvalue=R_fillval    )



      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Expand CO2 axis => extrapolate data outside the given range !
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      if (.not. fixed_CO2) then

        ! 0 PAL level:
        CO2_levels(1) = 0
        !
        ! linear extrapolation coefficient (for temperature only)
        xi = - CO2_levels(2) / (CO2_levels(3) - CO2_levels(2))
        !
        ! temperature
        where (temp(:,:,2) == T_fillval .or. temp(:,:,3) == T_fillval)
          temp(:,:,1) = T_fillval
        else where
          temp(:,:,1) = (1-xi)*temp(:,:,2) + xi*temp(:,:,3)
        end where
        !
        ! runoff => zero
        runoff(:,:,1) = 0

        ! Upper bound for CO2 levels: 128 times the highest level
        CO2_levels(nCO2+2) = 128*CO2_levels(nCO2+1)
        !
        ! linear extrapolation coefficient
        xi = (CO2_levels(nCO2+2) - CO2_levels(nCO2)) / (CO2_levels(nCO2+1) - CO2_levels(nCO2))
        !
        ! temperature
        where (temp(:,:,nCO2) == T_fillval .or. temp(:,:,nCO2+1) == T_fillval)
          temp(:,:,nCO2+2) = T_fillval
        else where
          temp(:,:,nCO2+2) = (1-xi)*temp(:,:,nCO2) + xi*temp(:,:,nCO2+1)
        end where
        !
        ! runoff
        where (runoff(:,:,nCO2) == R_fillval .or. runoff(:,:,nCO2+1) == R_fillval)
          runoff(:,:,nCO2+2) = R_fillval
        else where
          runoff(:,:,nCO2+2) = (1-xi)*runoff(:,:,nCO2) + xi*runoff(:,:,nCO2+1)
          ! Avoid negative runoff:
          where (runoff(:,:,nCO2+2) < 0)
            runoff(:,:,nCO2+2) = 0
          end where
        end where

      end if



      !-------!
      ! Slope !
      !-------!

      print *
      print *, 'Load slope'

      ! read names
      call read_comment(1)
      read(unit=1, fmt=*) path, dimname(1:2), varname
      call add_rootpath(path)

      ! Sizing
      call check_sizing( path, dimname(1:2), (/nlon,nlat/) )
      !**************************!
      allocate( slope(nlon,nlat) )
      !**************************!

      ! get variable + check coordinates and units
      call load_variable('slope', varname, dimname(1), dimname(2), single_input_file=path, &
                         varout2D=slope, xref=lon, yref=lat, fillvalue=S_fillval           )



      !--------------------!
      ! lithology fraction !
      !--------------------!

      print *
      print *, 'Load lithology fraction'

      call read_comment(1)

      ! Try direct reading of lithology fraction (assuming uniform lithology)

      read(unit=1, fmt='(A)') line
      nlith = 1
      do k = 1,len(line)
        if (line(k:k)==',') nlith = nlith + 1
      end do
      allocate( uniform_lithfrac(nlith) )
      read(line, fmt=*, iostat=ierr) uniform_lithfrac

      if (ierr==0) then

        print *, '--- direct reading from main IO file'

        !************************************!
        allocate( lith_frac(nlon,nlat,nlith) )
        !************************************!

        do k = 1,nlith
          lith_frac(:,:,k) = uniform_lithfrac(k)
        end do

        L_fillval = -1d99


      else ! Expect netCDF file name. Load lith_frac from file

        print *, '--- load from separate netCDF file'

        ! read names
        backspace(unit=1)
        read(unit=1, fmt=*) path, dimname(1:3), varname
        call add_rootpath(path)

        ! Sizing
        call check_sizing( path, dimname(1:2), (/nlon,nlat/) )
        !***************************************!
        nlith = netcdf_get_size(path, dimname(3))
        allocate( lith_frac(nlon,nlat,nlith) )
        !***************************************!

        ! get variable + check coordinates and units
        call load_variable('lithology', varname, dimname(1), dimname(2), z_varname=dimname(3), single_input_file=path, &
                           varout3D=lith_frac, xref=lon, yref=lat, fillvalue=L_fillval                                 )

      end if



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
      where (      temp==T_fillval )  temp      = DEFFILLVAL
      where (    runoff==R_fillval )  runoff    = DEFFILLVAL
      where (     slope==S_fillval )  slope     = DEFFILLVAL
      where ( lith_frac==L_fillval )  lith_frac = DEFFILLVAL



      !!===========================================================================!!
      !!  Check consistency of land area, missing-value and lithological fraction  !!
      !!===========================================================================!!

      print *
      print *
      print *, 'Check data consistency'
      print *, '======================'

      call check_continental_cells( land_area, temp, runoff, slope, lith_frac )



      !!==================================================!!
      !!  read other variables names and load variables:  !!
      !!==================================================!!

      print *
      print *
      print *, 'Read run parameters and forcings'
      print *, '================================'


      !---------------------------!
      ! Single run/Multirun flag: !
      !---------------------------!

      print *
      print *, 'Parameterization'

      call read_comment(1)
      read(unit=1, fmt=*) multirun

      if (multirun) then
        print *, '--- single'
      else
        print *, '--- multiple'
      end if


      !-----------------------!
      ! Forward/Backward run: !
      !-----------------------!

      print *
      print *, 'Computation mode'

      call read_comment(1)
      read(unit=1, fmt=*) ForwBckw

      if (ForwBckw == 1) then
        print *, '--- forward (forced by CO2)'
      elseif (ForwBckw == -1) then
        print *, '--- backward (forced by degassing)'
      else
        print *
        print *, 'ERROR: illegal mode', ForwBckw
        print *, 'expect "+1" (forward mode) or "-1" (inverse mode)'
        stop
      end if


      !-------------!
      ! Parameters: !
      !-------------!

      print *
      print *, 'Read model parameters'

      call read_comment(1)
      read(unit=1, fmt=*) path
      call add_rootpath(path)
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

      print *
      if (ForwBckw == 1) then
        print *, 'Read model forcing (CO2)'
      else
        print *, 'Read model forcing (degassing)'
      end if

      call read_comment(1)
      read(unit=1, fmt=*) path
      call add_rootpath(path)

      print *, '  from file: '//trim(path)


      if (fixed_CO2) then ! No CO2 dimension => fixed CO2 climate => write "fake" CO2 forcing to skip climate interpolation

        if (ForwBckw==-1) then
          print *
          print *, 'Error: cannot do a fixed-CO2 run (climate fields not defined on a CO2 axis) in inverse mode.'
          stop
        end if

        open(unit=IFORC, status='scratch', action='readwrite')

        if (nrun==0) then
          write(unit=IFORC, fmt=*) CO2_levels(1)
        else
          do k = 1,nrun
            write(unit=IFORC, fmt=*) CO2_levels(1)
          end do
        end if

        rewind(unit=IFORC)


      else

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
              print *
              print *, 'ERROR: inconsistent number of parameterizations and forcings. File length mismatch'
              stop
            end if

          end if
        end if

      end if



      !!======================================================================!!
      !!  read output file and variables names and create output netCDF file  !!
      !!======================================================================!!

      print *
      print *
      print *, 'Read information for creating output file'
      print *, '========================================='


      !----------------!
      ! read file name !
      !----------------!
      call read_comment(1)
      read(unit=1, fmt=*) path
      call add_rootpath(path)

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


      print *
      print *, 'Create output file'

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




    !==============================================================================================================================!
    !------------------------------------------------------------------------------------------------------------------------------!
    !==============================================================================================================================!




    subroutine add_rootpath(filename)
      character(len=*), intent(inout):: filename
      include 'path.inc' ! => 'gdss_root_path' variable, path of root directory
      !
      ! test if absolute path given (file name starts with a '/')
      ! if not, add current root path to the name:
      if (filename(1:1) /= '/') then
        filename = gdss_root_path//filename
      end if
    end subroutine



    !==============================================================================================================================!



    subroutine load_variable(internal_varname, varname, x_varname, y_varname, z_varname, &
                             single_input_file, multiple_input_file,                     &
                             varout2D, varout3D, x, y, xvec, yvec, xref, yref, totarea,  &
                             fill_missing, fillvalue                                     )
        use netcdf

        character(len=*), intent(in):: internal_varname
        character(len=*),  intent(in):: varname, x_varname, y_varname
        character(len=*),  intent(in), optional:: z_varname
        character(len=*), intent(in), optional:: single_input_file
        character(len=*), dimension(:),  intent(in),  optional:: multiple_input_file
        double precision, dimension(:,:),   intent(in),  optional:: totarea
        double precision, dimension(:,:),   intent(out), optional:: varout2D
        double precision, dimension(:,:,:), intent(out), optional:: varout3D
        double precision, intent(out), optional:: x(:), y(:), xvec(:,:), yvec(:,:)
        double precision, intent(in),  optional:: xref(:), yref(:)
        double precision, intent(out), optional:: fillvalue
        logical, intent(in), optional:: fill_missing

        double precision, allocatable, dimension(:):: loc_x, loc_y
        double precision, allocatable, dimension(:,:):: dummyvar2D
        double precision, allocatable, dimension(:,:,:):: dummyvar3D
        double precision, allocatable, dimension(:):: fillvalue_vec1D
        double precision, allocatable, dimension(:,:):: fillvalue_vec2D
        character(len=500):: var_units, vname
        character(len=1):: oper
        logical:: loc_fill_missing
        integer:: nx, ny, nz
        integer:: ierr, check
        integer:: k, n, noperations


        ! How to handle missing-value
        if (present(fill_missing)) then
            loc_fill_missing = fill_missing
        else
            loc_fill_missing = .false. ! default behaviour
        end if


        ! Check consistency between 2D variable or 3D (ie, multiple CO2 levels) variable
        ! => the right "input_file" and "var" arguments must be given.
        check = -11
        if (present(varout2D)) then
            check = check + 1
        end if
        if (present(single_input_file)) then
            check = check + 10
        end if
        if (present(z_varname)) then
            check = check - 2
        end if
        if (present(varout3D)) then
            check = check + 3
        end if
        if (present(multiple_input_file)) then
            check = check + 8
        end if
        if (check /= 0) then
            print *
            print *, 'INTERNAL ERROR: inconsistent set of optional variable passed to the subroutine'
            print *, '"load_variable", in module "GCM_io_module"'
            stop
        end if


        ! Instead of reading just 1 variable name, they may be several variables,
        ! with arithmetic operation between them (addition or substraction)
        noperations = get_arithmetic_operations(varname)
        ! => write the scratch files:
        !        unit=334: list of individual variables
        !        unit=335: list of mathematics operators (+ or -)


        if (present(varout2D)) then

            ! Shaping
            nx = size(varout2D, 1)
            ny = size(varout2D, 2)
            allocate( dummyvar2D(nx,ny) )
            allocate( fillvalue_vec1D(noperations) )
            allocate( loc_x(nx) )
            allocate( loc_y(ny) )

            ! Initialization
            varout2D = 0d0

            ! Loop to get variables and perform arithmetic operations:
            do n = 1,noperations

                ! load current variable
                read(unit=334, fmt=*) vname
                if (loc_fill_missing) then
                    ! load variable
                    call load_netcdf_dble2D(single_input_file, x_varname, y_varname, vname, dummyvar2D, &
                                            x=loc_x, y=loc_y, varunits=var_units, fillval=fillvalue_vec1D(n), fillval_iostat=ierr)
                    ! set var=0 on "missing" cells
                    if (ierr==NF90_NOERR) then
                        where (dummyvar2D==fillvalue_vec1D(n)) dummyvar2D = 0d0
                    else
                        fillvalue_vec1D(n) = -1d99 ! dummy fillvalue
                    end if
                else
                    ! load variable
                    call load_netcdf_dble2D(single_input_file, x_varname, y_varname, vname, dummyvar2D, &
                                            x=loc_x, y=loc_y, varunits=var_units, fillval=fillvalue_vec1D(n))
                end if
                ! check axis
                if (present(xref) .and. present(yref)) then
                    call check_coordinates(vname, loc_x, xref, loc_y, yref)
                end if
                ! check variable units
                if (present(totarea)) then
                    call check_units(internal_varname, vname, dummyvar2D, var_units, fillvalue=fillvalue_vec1D(n), totarea=totarea)
                else
                    call check_units(internal_varname, vname, dummyvar2D, var_units, fillvalue=fillvalue_vec1D(n))
                end if

                ! perform arithmetic operation
                read(unit=335, fmt=*) oper
                select case (oper)
                    case ("+")
                        where (varout2D/=fillvalue_vec1D(1) .and. dummyvar2D/=fillvalue_vec1D(n))
                            varout2D = varout2D + dummyvar2D
                        else where
                            varout2D = fillvalue_vec1D(1)
                        end where
                    case ("-")
                        where (varout2D/=fillvalue_vec1D(1) .and. dummyvar2D/=fillvalue_vec1D(n))
                            varout2D = varout2D - dummyvar2D
                        else where
                            varout2D = fillvalue_vec1D(1)
                        end where
                end select

            end do

            if (present(x))  x = loc_x
            if (present(y))  y = loc_y
            if (present(fillvalue))  fillvalue = fillvalue_vec1D(1)
            deallocate(dummyvar2D)
            deallocate(fillvalue_vec1D)
            deallocate(loc_x)
            deallocate(loc_y)


        else ! varout3D case


            ! Shaping
            nx = size(varout3D, 1)
            ny = size(varout3D, 2)
            nz = size(varout3D, 3)
            allocate( dummyvar3D(nx,ny,nz) )
            allocate( loc_x(nx) )
            allocate( loc_y(ny) )

            ! Initialization
            varout3D = 0d0


            if (present(z_varname)) then ! direct (one-shot) load of 3D variable

                
                allocate( fillvalue_vec1D(noperations) )

                ! Loop to get variables and perform arithmetic operations:
                do n = 1,noperations

                    ! load current variable
                    read(unit=334, fmt=*) vname
                    if (loc_fill_missing) then
                        ! load variable
                        call load_netcdf_dble3D(single_input_file, x_varname, y_varname, z_varname, vname, dummyvar3D, &
                                                x=loc_x, y=loc_y, varunits=var_units, fillval=fillvalue_vec1D(n), &
                                                fillval_iostat=ierr)
                        ! set var=0 on "missing" cells
                        if (ierr==NF90_NOERR) then
                            where (dummyvar3D==fillvalue_vec1D(n)) dummyvar3D = 0d0
                        else
                            fillvalue_vec1D(n) = -1d99 ! dummy fillvalue
                        end if
                    else
                        ! load variable
                        call load_netcdf_dble3D(single_input_file, x_varname, y_varname, z_varname, vname, dummyvar3D, &
                                                x=loc_x, y=loc_y, varunits=var_units, fillval=fillvalue_vec1D(n))
                    end if
                    ! check axis
                    if (present(xref) .and. present(yref)) then
                        call check_coordinates(vname, loc_x, xref, loc_y, yref)
                    end if
                    ! check variable units
                    if (present(totarea)) then
                        call check_units_3Dvar(internal_varname, vname, dummyvar3D, var_units, fillvalue=fillvalue_vec1D(n), &
                                               totarea=totarea)
                    else
                        call check_units_3Dvar(internal_varname, vname, dummyvar3D, var_units, fillvalue=fillvalue_vec1D(n))
                    end if

                    ! perform arithmetic operation
                    read(unit=335, fmt=*) oper
                    select case (oper)
                        case ("+")
                            where (varout3D/=fillvalue_vec1D(1) .and. dummyvar3D/=fillvalue_vec1D(n))
                                varout3D = varout3D + dummyvar3D
                            else where
                                varout3D = fillvalue_vec1D(1)
                            end where
                        case ("-")
                            where (varout3D/=fillvalue_vec1D(1) .and. dummyvar3D/=fillvalue_vec1D(n))
                                varout3D = varout3D - dummyvar3D
                            else where
                                varout3D = fillvalue_vec1D(1)
                            end where
                    end select

                end do

                if (present(fillvalue))  fillvalue = fillvalue_vec1D(1)
                deallocate(fillvalue_vec1D)


            else ! 3rd dimension = CO2 axis => load var for each CO2 level, one level per file


                allocate( fillvalue_vec2D(nz,noperations) )

                ! Loop to get variables and perform arithmetic operations:
                do n = 1,noperations

                    ! load current variable
                    read(unit=334, fmt=*) vname
                    do k = 1,nz
                        if (loc_fill_missing) then
                            ! load variable
                            call load_netcdf_dble2D(multiple_input_file(k), x_varname, y_varname, vname, dummyvar3D(:,:,k), &
                                            x=loc_x, y=loc_y, varunits=var_units, fillval=fillvalue_vec2D(k,n), fillval_iostat=ierr)
                            ! set var=0 on "missing" cells
                            if (ierr==NF90_NOERR) then
                                where (dummyvar3D(:,:,k)==fillvalue_vec2D(k,n)) dummyvar3D(:,:,k) = 0d0
                            else
                                fillvalue_vec2D(k,n) = -1d99 ! dummy fillvalue
                            end if
                        else
                            ! load variable
                            call load_netcdf_dble2D(multiple_input_file(k), x_varname, y_varname, vname, dummyvar3D(:,:,k), &
                                                    x=loc_x, y=loc_y, varunits=var_units, fillval=fillvalue_vec2D(k,n))
                        end if
                        ! check axis
                        if (present(xref) .and. present(yref)) then
                            call check_coordinates(vname, loc_x, xref, loc_y, yref)
                        end if
                        ! check variable units
                        if (present(totarea)) then
                            call check_units(internal_varname, vname, dummyvar3D(:,:,k), var_units, fillvalue=fillvalue_vec2D(k,n),&
                                            totarea=totarea)
                        else
                            call check_units(internal_varname, vname, dummyvar3D(:,:,k), var_units, fillvalue=fillvalue_vec2D(k,n))
                        end if
                        !
                        if (present(xvec))  xvec(:,k) = loc_x
                        if (present(yvec))  yvec(:,k) = loc_y

                    end do

                    ! perform arithmetic operation
                    read(unit=335, fmt=*) oper
                    select case (oper)
                        case ("+")
                            do k = 1,nz
                                where (varout3D(:,:,k)/=fillvalue_vec2D(1,1) .and. dummyvar3D(:,:,k)/=fillvalue_vec2D(k,n))
                                    varout3D(:,:,k) = varout3D(:,:,k) + dummyvar3D(:,:,k)
                                else where
                                    varout3D(:,:,k) = fillvalue_vec2D(1,1)
                                end where
                            end do
                        case ("-")
                            do k = 1,nz
                                where (varout3D(:,:,k)/=fillvalue_vec2D(1,1) .and. dummyvar3D(:,:,k)/=fillvalue_vec2D(k,n))
                                    varout3D(:,:,k) = varout3D(:,:,k) - dummyvar3D(:,:,k)
                                else where
                                    varout3D(:,:,k) = fillvalue_vec2D(1,1)
                                end where
                            end do
                    end select

                end do

                if (present(fillvalue))  fillvalue = fillvalue_vec2D(1,1)
                deallocate(fillvalue_vec2D)


            end if


            deallocate(dummyvar3D)
            deallocate(loc_x)
            deallocate(loc_y)


        end if


        ! close scratch files
        close(unit=334)
        close(unit=335)


    end subroutine


    !=======================!


    function get_arithmetic_operations(varstring)
    ! => create scratch files 334 and 335 containing (respectively) the lists of variables and operators
    ! Note: the first operators (corresponding to the first variable) is automatically '+'
    ! return the number of variables/operations
        character(len=*), intent(in):: varstring
        integer:: get_arithmetic_operations
        character(len=1), dimension(2), parameter:: operators = (/'+', '-'/) ! Legal arithmetic operators
        integer:: i0, i1
        open(unit=334, status='scratch', action='readwrite') ! individual variables
        open(unit=335, status='scratch', action='readwrite') ! operators
        i0 = 1
        i1 = 1
        get_arithmetic_operations = 1
        write(unit=335, fmt='(A)') '+'
        do i1 = 1,len(varstring)
            if (any(operators == varstring(i1:i1))) then
                write(unit=334, fmt='(A)') trim(adjustl(varstring(i0:i1-1)))
                write(unit=335, fmt='(A)') varstring(i1:i1)
                i0 = i1+1
                get_arithmetic_operations = get_arithmetic_operations + 1
            end if
        end do
        write(unit=334, fmt='(A)') trim(adjustl(varstring(i0:)))
        rewind(unit=334)
        rewind(unit=335)
    end function


    !=======================!


    subroutine load_netcdf_dble2D(fname, x_varname, y_varname, varname, var, x, y, varunits, fillval, fillval_iostat)

        use netcdf
        use netcdf_io_functions, only: nf90_check

        character(len=*), intent(in):: fname, x_varname, y_varname, varname
        double precision, dimension(:,:), intent(out):: var
        double precision, dimension(:), intent(out), optional:: x, y
        character(len=*), intent(out), optional:: varunits
        double precision, intent(out), optional:: fillval
        integer, intent(out), optional:: fillval_iostat
        integer:: nx, ny
        integer:: ierr, fid, varid, truedimids(2), shp_idx(2)
        integer, dimension(:), allocatable:: dimids
        logical:: transp
        double precision, dimension(:,:), allocatable:: loc_var

        nx = size(var, 1)
        ny = size(var, 2)
        ! must also match the size of "x" and "y"


        ! variable shape-independent operations (open file, get variable ID, ...)
        call load_netcdf_generic(fname, varname, fid, varid, dimids, truedimids, shp_idx, &
                                 x_varname=x_varname, y_varname=y_varname, x=x, y=y,      &
                                 varunits=varunits, fillval=fillval, fillval_iostat=fillval_iostat)

        ! Check that variable is defined on the given dimensions
        if (dimids(shp_idx(1))==truedimids(1) .and. dimids(shp_idx(2))==truedimids(2)) then
            transp = .false.
            allocate(loc_var(nx,ny))
        elseif (dimids(shp_idx(1))==truedimids(2) .and. dimids(shp_idx(2))==truedimids(1)) then
            transp = .true.
            allocate(loc_var(ny,nx))
        else
            print *
            print *
            print *, 'Error: variable "'//trim(varname)//'" of file "'//trim(fname)// &
            '" is not defined on the given dimensions "'//trim(x_varname)//'" and "'//trim(y_varname)//'"'
            stop
        end if

        ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !
        ierr = nf90_get_var(fid, varid, loc_var)
        call nf90_check(ierr, 'Error while getting variable "'//trim(varname)//'" of file "'//trim(fname)//'"')
        ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !
        if (transp) then
            var = transpose(loc_var)
        else
            var = loc_var
        end if

        ! close netCDF file
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//trim(fname), kill=.false.)

        deallocate(dimids)

    end subroutine


    !======================================================================!


    subroutine load_netcdf_dble3D(fname, x_varname, y_varname, z_varname, varname, var, x, y, z, varunits, fillval, fillval_iostat)

        use netcdf
        use netcdf_io_functions, only: nf90_check

        character(len=*), intent(in):: fname, x_varname, y_varname, z_varname, varname
        double precision, dimension(:,:,:), intent(out):: var
        double precision, dimension(:), intent(out), optional:: x, y, z
        character(len=*), intent(out), optional:: varunits
        double precision, intent(out), optional:: fillval
        integer, intent(out), optional:: fillval_iostat
        integer:: nx, ny, nz
        integer:: ierr, fid, varid, truedimids(3), shp_idx(3)
        integer, dimension(:), allocatable:: dimids

        nx = size(var, 1)
        ny = size(var, 2)
        nz = size(var, 3)
        ! must also match the size of "x", "y" and "z"


        ! variable shape-independent operations (open file, get variable ID, ...)
        call load_netcdf_generic(fname, varname, fid, varid, dimids, truedimids, shp_idx,                      &
                                 x_varname=x_varname, y_varname=y_varname, z_varname=z_varname, x=x, y=y, z=z, &
                                 varunits=varunits, fillval=fillval, fillval_iostat=fillval_iostat             )

        ! Check that variable is defined on the given dimensions
        ! NOTE: load_netcdf_dble3D does not allow transposition (ie, permuted dimensions)
        if (.not. (dimids(shp_idx(1))==truedimids(1) .and. &
                   dimids(shp_idx(2))==truedimids(2) .and. &
                   dimids(shp_idx(3))==truedimids(3))) then
            print *
            print *
            print *, 'Error: variable "'//trim(varname)//'" of file "'//trim(fname)// &
            '" is not defined on the given dimensions {"'//trim(x_varname)//'", "'//trim(y_varname)//'", "'//trim(z_varname)//'"}'
            stop
        end if

        ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !
        ierr = nf90_get_var(fid, varid, var)
        call nf90_check(ierr, 'Error while getting variable "'//trim(varname)//'" of file "'//trim(fname)//'"')
        ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !

        ! close netCDF file
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//trim(fname), kill=.false.)

        deallocate(dimids)

    end subroutine


    !======================================================================!


    subroutine load_netcdf_generic(fname, varname, fid, varid, dimids, truedimids, shp_idx, &
                                   x_varname, y_varname, z_varname, x, y, z, varunits, fillval, fillval_iostat)

        use netcdf
        use netcdf_io_functions, only: nf90_check

        character(len=*), intent(in):: fname, varname
        integer, intent(out):: fid, varid
        integer, dimension(:), intent(out):: truedimids ! must be dimension(2) or dimension(3)
        integer, dimension(:), intent(out), allocatable:: dimids
        integer, dimension(:), intent(out):: shp_idx
        character(len=*), intent(in), optional:: x_varname, y_varname, z_varname
        double precision, dimension(:), intent(out), optional:: x, y, z
        character(len=*), intent(out), optional:: varunits
        double precision, intent(out), optional:: fillval
        integer, intent(out), optional:: fillval_iostat
        integer, dimension(:), allocatable:: shp
        integer:: ndim_true, ndim, k, n
        integer:: ierr

        ! Get number of dimension of the variable "in the code" (ie, without degenerated dimensions),
        ! defined by the length of output argument "truedimids"
        ndim_true = size(truedimids)
        if (ndim_true/=2 .and. ndim_true/=3) then
            print *
            print *, 'INTERNAL ERROR in subroutine "load_netcdf_generic" of module "io_module":'
            print *, 'argument "truedimids" must be length 2 or 3 (=> shape of the variable to load)'
            stop
        end if
        ! "true" number of dimension must also be the length of "truedimids"
        if (ndim_true /= size(shp_idx)) then
            print *
            print *, 'INTERNAL ERROR in subroutine "load_netcdf_generic" of module "io_module":'
            print *, 'argument "shp_idx" shape does not match argument "truedimids" one'
        end if

        ! open netCDF file
        ierr = nf90_open(fname, NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//trim(fname))

        ! Get axis data
        if (present(x_varname)) then
            if (present(x)) then
                call load_axis(fname, fid, x_varname, truedimids(1), x)
            else
                call load_axis(fname, fid, x_varname, truedimids(1))
            end if
        end if
        if (present(y_varname)) then
            if (present(y)) then
                call load_axis(fname, fid, y_varname, truedimids(2), y)
            else
                call load_axis(fname, fid, y_varname, truedimids(2))
            end if
        end if
        if (present(z_varname) .and. ndim_true==3) then
            if (present(z)) then
                call load_axis(fname, fid, z_varname, truedimids(3), z)
            else
                call load_axis(fname, fid, z_varname, truedimids(3))
            end if
        end if

        ! Get variable ID
        ierr = nf90_inq_varid(fid, varname, varid)
        call nf90_check(ierr, 'Error while inquiring variable "'//trim(varname)//'" ID in file "'//trim(fname)//'"')

        ! Get variable number of dimension
        ierr = nf90_inquire_variable(fid, varid, ndims=ndim)
        call nf90_check(ierr, 'Error while inquiring number of dim. of variable "'//trim(varname)//'" in file "'//trim(fname)//'"')

        allocate(dimids(ndim))
        allocate(shp(ndim))

        ! Get variable dimension IDs
        ierr = nf90_inquire_variable(fid, varid, dimids=dimids)
        call nf90_check(ierr, 'Error while inquiring dimensions IDs of variable "'//trim(varname)//'" in file "'//trim(fname)//'"')

        ! Get variable shape (inquire length of all dimensions)
        do k = 1,ndim
            ierr = nf90_inquire_dimension(fid, dimids(k), len=shp(k))
            call nf90_check(ierr, 'Error while inquiring length of dimension in file "'//trim(fname)//'"')
        end do

        ! Check variable shape (must have exactly 3 non degenerated dimensions)
        if (count(shp > 1) /= ndim_true) then
            print *
            write(*, fmt='(A,I2,A)') ' Error: variable "'//trim(varname)//'" of file "'//trim(fname)//'" must have exactly ', &
                                     ndim_true, ' non-degenerated (ie, size-1) dimensions'
            stop
        else
            k = 0
            do n = 1,ndim_true
                k = k + 1
                do while (shp(k)==1)
                    k = k + 1
                end do
                shp_idx(n) = k
            end do
        end if

        ! Notify for degenerated dimensions
        if (ndim>ndim_true) then
            print *
            print *,                 'Variable "'//trim(varname)//'" of file "'//trim(fname)//'"'
            write(*, fmt='(A,I2,A)') '    Ignore ', ndim-ndim_true, ' degenerated (size-1) dimension(s).'
        end if

        ! Get variable units (if asked)
        if (present(varunits)) then
            varunits = ''
            ierr = nf90_get_att(fid, varid, 'units', varunits)
            call nf90_check(ierr, 'Warning: unable to get attribute "units" of variable "'//trim(varname)//'" in file "' &
                                   //trim(fname)//'". Variable assumed to be dimensionless.', kill=.false.)
        end if

        ! Get variale fill-value (if asked)
        if (present(fillval)) then
            ierr = nf90_get_att(fid, varid, '_FillValue', fillval)
            if (present(fillval_iostat)) then
                fillval_iostat = ierr
            else
                call nf90_check(ierr, 'Error while getting attribute "_FillValue" of variable "'//trim(varname)//'" in file "' &
                                                                                                          //trim(fname)//'"')
            end if
        end if

        deallocate(shp)

    end subroutine


    !=======================!


    subroutine load_axis(fname, fid, axname, axdimid, ax)
        use netcdf
        use netcdf_io_functions, only: nf90_check
        character(len=*), intent(in):: fname, axname
        integer, intent(in):: fid
        integer, intent(out):: axdimid
        double precision, dimension(:), intent(out), optional:: ax
        integer:: varid, ierr

        ! Inquire dim ID
        ierr = nf90_inq_dimid(fid, axname, axdimid)
        call nf90_check(ierr, 'Error while inquiring dimension "'//trim(axname)//'" ID in file "'//trim(fname)//'"')
        ! Inquire var ID
        ierr = nf90_inq_varid(fid, axname, varid)
        call nf90_check(ierr, 'Error while inquiring variable "'//trim(axname)//'" ID in file "'//trim(fname)//'"')
        ! Get variable
        if (present(ax)) then
            ierr = nf90_get_var(fid, varid, ax)
            call nf90_check(ierr, 'Error while getting variable "'//trim(axname)//'" in file "'//trim(fname)//'"')
        end if

    end subroutine



    !==============================================================================================================================!



    subroutine define_variables( fid, dimname, fillval, axunits, multirun )

      use netcdf_io_functions, only: create_output_variable
      use ascii_io_functions, only: read_comment

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



    !==============================================================================================================================!



    subroutine check_sizing( path, dimnames, sizes )

      use netcdf_io_functions, only: netcdf_get_size

      character(len=*), intent(in):: path
      character(len=*), dimension(:), intent(in):: dimnames
      integer, dimension(:):: sizes
      integer:: n, l

      do n = 1,size(dimnames)

        l = netcdf_get_size(path, dimnames(n))

        if (l/=sizes(n)) then

          print *
          print *, 'ERROR: size mismatch of dimension '//trim(dimnames(n))//' in file '//trim(path)
          print *, '  Expected: ',sizes(n)
          print *, '       Got: ',l
          stop

        end if

      end do

    end subroutine



    !==============================================================================================================================!



    subroutine compare_units(units_string, known_units, conversion_factor, conversion_offset, passed)

        use physical_units, only: units

        character(len=*), intent(in):: units_string
        type(units), intent(in):: known_units
        double precision, intent(out):: conversion_factor, conversion_offset
        logical, intent(out):: passed
        integer:: k

        conversion_factor = 1d0
        conversion_offset = 0d0

        passed = .false.

        if (units_string == known_units%reference) then
            passed = .true.
        else
            do k = 1,known_units%naccepted
                if (units_string == known_units%accepted(k)%string) then
                    passed = .true.
                    conversion_factor = known_units%accepted(k)%conversion(1)
                    conversion_offset = known_units%accepted(k)%conversion(2)
                end if
            end do
        end if

    end subroutine


    !=======================!


    subroutine check_units(which_variable, varname, var, varunits, fillvalue, totarea)

        use physical_units, only: units

        character(len=*), intent(in):: which_variable, varname, varunits
        double precision, dimension(:,:), intent(inout):: var
        double precision, dimension(:,:), intent(in), optional:: totarea
        double precision, intent(in), optional:: fillvalue
        double precision:: factor, offset
        logical:: multiply_by_area

        ! Get univt conversion parameters:
        call check_units_generic(which_variable, varname, varunits, factor, offset, multiply_by_area)

        ! Conversion
        if (present(fillvalue)) then
            where (var/=fillvalue)  var = factor*(var + offset)
        else
            var = factor*(var + offset)
        end if

        ! Multiply value by cell area
        if (multiply_by_area) then
            if (present(totarea)) then
                print *, 'Variable "'//trim(varname)//'" multiplied by total area'
                if (present(fillvalue)) then
                    where (var/=fillvalue) var = var * totarea
                else
                    var = var * totarea
                end if
            else
                print *
                print *, 'INTERNAL ERROR: total area not passed to the function "check_units" in module "GCM_io_module".'
                print *, 'You have the right to be mad at the developer.'
                stop
            end if
        end if

    end subroutine


    !=======================!


    subroutine check_units_3Dvar(which_variable, varname, var, varunits, fillvalue, totarea)

        use physical_units, only: units

        character(len=*), intent(in):: which_variable, varname, varunits
        double precision, dimension(:,:,:), intent(inout):: var
        double precision, dimension(:,:), intent(in), optional:: totarea
        double precision, intent(in), optional:: fillvalue
        double precision:: factor, offset
        logical:: multiply_by_area
        integer:: k

        ! Get univt conversion parameters:
        call check_units_generic(which_variable, varname, varunits, factor, offset, multiply_by_area)

        ! Conversion
        if (present(fillvalue)) then
            where (var/=fillvalue)  var = factor*(var + offset)
        else
            var = factor*(var + offset)
        end if

        ! Multiply value by cell area
        if (multiply_by_area) then
            if (present(totarea)) then
                print *, 'Variable "'//trim(varname)//'" multiplied by total area'
                if (present(fillvalue)) then
                    do k = 1,size(var, 3)
                        where (var(:,:,k)/=fillvalue)  var(:,:,k) = var(:,:,k) * totarea
                    end do
                else
                    do k = 1,size(var, 3)
                        var(:,:,k) = var(:,:,k) * totarea
                    end do
                end if
            else
                print *
                print *, 'INTERNAL ERROR: total area not passed to the function "check_units_3Dvar" in module "GCM_io_module".'
                print *, 'You have the right to be mad at the developer.'
                stop
            end if
        end if

    end subroutine


    !=======================!


    subroutine check_units_generic(which_variable, varname, varunits, factor, offset, multiply_by_area)

        use physical_units, only: units, area_units, fraction_units, slope_units, temperature_units, runoff_units

        character(len=*), intent(in):: which_variable, varname, varunits
        double precision, intent(out):: factor, offset
        logical, intent(out):: multiply_by_area
        type(units):: known_units
        logical:: passed
        integer:: go_on


        multiply_by_area = .false.

        ! case-dependent statement
        select case (which_variable)

            case ("area")
                known_units = area_units()

            case ("landarea")
                known_units = area_units()

            case ("temperature")
                known_units = temperature_units()

            case ("runoff")
                known_units = runoff_units()

            case ("slope")
                known_units = slope_units()

            case ("lithology")
                known_units = fraction_units()

            case default
                print *
                print *, 'INTERNAL ERROR in function "check_units" of module "GCM_io_module": unkown variable case "' &
                                                                                            //trim(which_variable)//'"'
                stop

        end select


        ! Compare given units to reference units, and get conversion values if not equals
        call compare_units(varunits, known_units, factor, offset, passed)
        ! in the case "landarea", also try "fraction" units
        if (which_variable == "landarea" .and. (.not. passed)) then
            call compare_units(varunits, fraction_units(), factor, offset, passed)
            multiply_by_area = passed
        end if

        ! Print message
        if (passed) then
            if (factor/=1 .or. offset/=0 .or. multiply_by_area) then
                print *
                print *, 'Automatic conversion of variable "'//trim(varname)//'":'
                print *, '    "'//trim(varunits)//'" => "'//trim(known_units%reference)//'"'
            end if
        else
            print *
            print *, 'WARNING: unkown units of variable "'//trim(varname)//'".'
            print *, '    got units:      "'//trim(varunits)//'"'
            print *, '    expected units: "'//trim(known_units%reference)//'"'
            go_on = -1
            do while (go_on/=1)
                print *, 'Enter 1 to continue the execution without modification, 0 to abort it:'
                read *, go_on
                if (go_on==0) stop
            end do
        end if

    end subroutine



    !==============================================================================================================================!



    subroutine check_coordinates_singlevar(varname, x, xref)
        include 'common_parameters.inc' !to get the value of MAX_ALLOWED_INACC
        character(len=*), intent(in):: varname
        double precision, dimension(:), intent(in):: x, xref
        double precision:: dx
        integer:: nx
        !
        nx = size(xref)
        dx = minval(xref(2:nx) - xref(1:nx-1))
        if ( maxval(abs((x-xref)/dx)) > MAX_ALLOWED_INACC ) then
            print *
            print *, 'ERROR: coordinates mismatch found for variable '//trim(varname)
            stop
        end if
    end subroutine


    !=======================!


    subroutine check_coordinates(varname, x1, xref1, x2, xref2, x3, xref3, x4, xref4, x5, xref5, x6, xref6, x7, xref7)
        character(len=*), intent(in):: varname
        double precision, dimension(:), intent(in), optional:: x1, x2, x3, x4, x5, x6, x7
        double precision, dimension(:), intent(in), optional:: xref1, xref2, xref3, xref4, xref5, xref6, xref7
        !
        if (present(x1) .and. present(xref1)) call check_coordinates_singlevar(varname, x1, xref1)
        if (present(x2) .and. present(xref2)) call check_coordinates_singlevar(varname, x2, xref2)
        if (present(x3) .and. present(xref3)) call check_coordinates_singlevar(varname, x3, xref3)
        if (present(x4) .and. present(xref4)) call check_coordinates_singlevar(varname, x4, xref4)
        if (present(x5) .and. present(xref5)) call check_coordinates_singlevar(varname, x5, xref5)
        if (present(x6) .and. present(xref6)) call check_coordinates_singlevar(varname, x6, xref6)
        if (present(x7) .and. present(xref7)) call check_coordinates_singlevar(varname, x7, xref7)
    end subroutine



    !==============================================================================================================================!



    subroutine check_continental_cells_single2Dvar( varname, land_area, var, check, new_land_area )

      include 'common_parameters.inc' ! to get the value of DEFFILLVAL

      character(len=*), intent(in):: varname
      double precision, dimension(:,:), intent(in):: land_area, var
      logical, intent(out):: check
      double precision, dimension(:,:), intent(inout), optional:: new_land_area
      double precision:: area_err, tot_land_area
      integer:: nerr

      nerr = count( (land_area>0 .and. var==DEFFILLVAL) )
      area_err = sum( land_area, mask=(land_area>0 .and. var==DEFFILLVAL) )
      tot_land_area = sum(land_area)
      if (present(new_land_area)) then
        where (land_area>0 .and. var==DEFFILLVAL) new_land_area = 0
      end if

      if (nerr > 0) then
        print *
        print *, 'WARNING: found missing values on continental cells of variable '//trim(varname)
        print *, 'Number of continent cells with missing values:     ', nerr
        print *, 'Fraction of continental cells with missing values: ', dble(nerr)/dble(count((land_area>0)))
        print *, 'Total area of those cells (m2):                    ', area_err
        print *, 'Which is a fraction of total land area:            ', area_err/tot_land_area
        print *
        check = .false.
      else
        check = .true.
      end if

    end subroutine


    !=======================!


    subroutine check_continental_cells_single3Dvar( varname, land_area, var, check, new_land_area )

      include 'common_parameters.inc' ! to get the value of DEFFILLVAL

      character(len=*), intent(in):: varname
      double precision, dimension(:,:), intent(in):: land_area
      double precision, dimension(:,:,:), intent(in):: var
      logical, intent(out):: check
      double precision, dimension(:,:), intent(inout), optional:: new_land_area
      double precision:: area_err, tot_land_area
      integer:: k, nerr

      check = .true.
      tot_land_area = sum(land_area)

      do k = 1,size( var, 3 )

        nerr = count( (land_area>0 .and. var(:,:,k)==DEFFILLVAL) )
        area_err = sum( land_area, mask=(land_area>0 .and. var(:,:,k)==DEFFILLVAL) )
        if (present(new_land_area)) then
          where (land_area>0 .and. var(:,:,k)==DEFFILLVAL) new_land_area = 0
        end if

        if (nerr > 0) then
          print *
          print *, 'WARNING: found missing values on continental cells of variable '//trim(varname)
          print *, '"Vertical" (3rd axis) level number: ', k
          print *, 'Number of continent cells with missing values:     ', nerr
          print *, 'Fraction of continental cells with missing values: ', dble(nerr)/dble(count((land_area>0)))
          print *, 'Total area of those cells (m2):                    ', area_err
          print *, 'Which is a fraction of total land area:            ', area_err/tot_land_area
          print *
          check = .false.
        end if

      end do

    end subroutine


    !=======================!


    subroutine check_continental_lithology(land_area, lith_frac, check, new_land_area )

      include 'common_parameters.inc' !to get the value of MAX_ALLOWED_INACC

      double precision, dimension(:,:), intent(in):: land_area
      double precision, dimension(:,:,:), intent(in):: lith_frac
      logical, intent(out):: check
      double precision, dimension(:,:), intent(inout), optional:: new_land_area
      double precision:: area_err, tot_land_area, sumdiff, max_sumdiff
      integer:: nerr,i,j

      nerr = 0
      max_sumdiff = 0d0
      area_err = 0d0
      tot_land_area = sum(land_area)

      if (present(new_land_area)) then

        do j=1,size(lith_frac,2)
          do i=1,size(lith_frac,1)

            if (land_area(i,j)>0) then

              sumdiff = abs( sum(lith_frac(i,j,:))  -  1 )

              if (sumdiff > MAX_ALLOWED_INACC) then
                nerr = nerr+1
                new_land_area(i,j) = 0
                area_err = area_err + land_area(i,j)
                if (sumdiff > max_sumdiff) max_sumdiff = sumdiff
              end if

            end if

          end do
        end do

      else

        do j=1,size(land_area,2)
          do i=1,size(land_area,1)

            if (land_area(i,j)>0) then

              sumdiff = abs( sum(lith_frac(i,j,:))  -  1 )

              if (sumdiff > MAX_ALLOWED_INACC) then
                nerr = nerr+1
                area_err = area_err + land_area(i,j)
                if (sumdiff > max_sumdiff) max_sumdiff = sumdiff
              end if

            end if

          end do
        end do

      end if


      if (nerr > 0) then
        print *
        print *, 'WARNING: inconsistent lithology class fraction'
        print *, 'Found cells where the sum of all lithology fractions differs from 1'
        print *, 'Number of cells affected:               ', nerr
        print *, 'Fraction of continental cells affected: ', dble(nerr)/dble(count((land_area>0)))
        print *, 'Total area of those cells (m2):         ', area_err
        print *, 'Which is a fraction of total land area: ', area_err/tot_land_area
        print *, 'Maximum absolute difference:            ', max_sumdiff
        print *
        check = .false.
      else
        check = .true.
      end if

    end subroutine


    !=======================!


    subroutine check_continental_cells( land_area, temp, runoff, slope, lith_frac )

      double precision, dimension(:,:), intent(inout):: land_area
      double precision, dimension(:,:), intent(in):: slope
      double precision, dimension(:,:,:), intent(in):: temp, runoff, lith_frac
      double precision, dimension(:,:), allocatable:: land_area_corr, land_area_corr_2
      logical, dimension(5):: checkpoints
      integer:: answer, ierr
      character(len=1):: answerchar

      allocate(land_area_corr(   size(land_area,1), size(land_area,2) ))
      allocate(land_area_corr_2( size(land_area,1), size(land_area,2) ))

      land_area_corr = land_area
      call check_continental_cells_single3Dvar( 'temperature', land_area, temp,      checkpoints(1), land_area_corr )
      call check_continental_cells_single3Dvar( 'runoff',      land_area, runoff,    checkpoints(2), land_area_corr )
      call check_continental_cells_single2Dvar( 'slope',       land_area, slope,     checkpoints(3), land_area_corr )
      call check_continental_cells_single3Dvar( 'lithology',   land_area, lith_frac, checkpoints(4), land_area_corr )
      land_area_corr_2 = land_area_corr
      call check_continental_lithology( land_area_corr, lith_frac, checkpoints(5), land_area_corr_2 )

      if (.not. all(checkpoints)) then

        print *, 'Hit one of the following options:'
        print *, '    0: abort the program'
        print *, '    1: remove all the erratic points (set land_area=0)'
        print *, '    2: remove only the missing points (and ignore the lithology issue)'
        print *, '    3: ignore all issues and run the program as it is'

        ! read potential argument
        call get_command_argument(1, answerchar, status=ierr)
        read(answerchar, fmt=*, iostat=ierr) answer
        if (ierr==0) then
          print *, answer
        else
          answer=-1
        end if

        ! reading loop (in case nor valid argument found)
        do while (answer/=0 .and. answer/=1 .and. answer/=2 .and. answer/=3)
          read(unit=*, fmt=*) answer
        end do

        ! According to the user answer:
        select case (answer)
          case(0)
            stop
          case(1)
            land_area = land_area_corr_2
          case(2)
            land_area = land_area_corr
        end select

      end if

      deallocate(land_area_corr)
      deallocate(land_area_corr_2)

    end subroutine



    !==============================================================================================================================!




end module

