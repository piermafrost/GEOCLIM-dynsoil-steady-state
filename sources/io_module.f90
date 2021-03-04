module io_module
  implicit none


  !<><><><><><><><><><><><><><><><><><>!
  ! input-output parameters and types: !
  !<><><><><><><><><><><><><><><><><><>!

  integer, parameter:: N_TOT_DIM = 4 ! (lon, lat, litho, runs)


  ! Technical parameters:
  integer, parameter, private:: LINE_CHARLEN=1000, FILE_CHARLEN=500, OTHR_CHARLEN=50
  !character(len=*), parameter, private:: LINE_CHARFMT='A1000', FILE_CHARFMT='A500', OTHR_CHARFMT='A50'
  double precision, parameter, private:: DEFFILLVAL = 9.96921d+36
  double precision, parameter, private:: MAX_ALLOWED_INACC = 1.d-6


  ! Type of object storing the information about output file
  !    sub-type:
  type varinfo
    logical:: unlimited, write_var
    integer:: varid
  end type
  !    main type:
  type outinfo
    character(len=FILE_CHARLEN):: file_name
    integer:: file_id
    type(varinfo), dimension(:), allocatable:: variables
  end type



  !<><><><><><><><><><><><><><>!
  ! input-output subroutines:  !
  !<><><><><><><><><><><><><><>!

  contains



    !# Main subroutine: read IO file, load input data, create output file
    !####################################################################

    subroutine make_input_output( io_fname, cell_area,land_area,temp,runoff,slope,lith_frac, CO2_levels, GMST, &
                                  curr_temp, curr_runf, h, E, xs, W, W_all, &
                                  nlon,nlat,nlith,nrun, ForwBckw, forcing, params, output_info )
      ! Main subroutine. Read input-output file, allocate variables, load variables and create output file

      use netcdf
      use netcdf_io_functions, only: netcdf_get_size, create_output_file, nf90_check
      use ascii_io_functions, only: file_length, read_comment
      use climate_module, only: INTERPOLATION_MODE, interp_coeff

      ! in/out variables
      character(len=*), intent(in):: io_fname
      double precision, intent(out), dimension(:,:),   allocatable:: cell_area, land_area, slope, curr_temp, curr_runf
      double precision, intent(out), dimension(:,:,:), allocatable:: temp, runoff, lith_frac, h, E, xs, W
      double precision, intent(out), dimension(:,:),   allocatable:: W_all
      double precision, intent(out), dimension(:),     allocatable:: CO2_levels, GMST, forcing
      integer, intent(out):: nlon, nlat, nlith, nrun, ForwBckw
      double precision, dimension(:,:,:), allocatable, intent(out):: params ! nparams x nlith x nrun
      type(outinfo), intent(out):: output_info

      ! local variables
      double precision, dimension(:),     allocatable:: lon, lat, uniform_lithfrac
      double precision, dimension(:,:,:), allocatable:: glob_temp
      double precision:: T_fillval, R_fillval, S_fillval, L_fillval, GT_fillval, xi
      character(len=LINE_CHARLEN):: line, varname, varname2
      character(len=OTHR_CHARLEN):: dimname(10)
      character(len=OTHR_CHARLEN):: axunits(2), nounits
      character(len=FILE_CHARLEN):: path
      character(len=FILE_CHARLEN), dimension(:), allocatable:: multipath
      integer:: nCO2, nCO2bis, lev0, lev1

      ! technical variables
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
        fixed_CO2 = (nCO2==1)
        if (fixed_CO2) then
          nCO2bis = 1
          lev0 = 1
          lev1 = 1
        else
          nCO2bis = nCO2+2 ! CO2 axis will be extrapolated
          lev0 = 2
          lev1 = nCO2+1
        end if
        !*****************************!
        allocate( CO2_levels(nCO2bis) )
        allocate( multipath(nCO2) )
        !*****************************!
        call read_comment(2)
        read(unit=2, fmt=*) CO2_levels(lev0:lev1)
        close(unit=2)


      else ! Try read CO2 axis directly from main file

        print *, '--- direct reading from main IO file'
        backspace(unit=1)
        read(unit=1, fmt='(A)') line
        nCO2 = 1
        do k = 1,len(line)
          if (line(k:k)==',') nCO2 = nCO2 + 1
        end do
        fixed_CO2 = (nCO2==1)
        if (fixed_CO2) then
          nCO2bis = 1
          lev0 = 1
          lev1 = 1
        else
          nCO2bis = nCO2+2 ! CO2 axis will be extrapolated
          lev0 = 2
          lev1 = nCO2+1
        end if
        !*****************************!
        allocate( CO2_levels(nCO2bis) )
        allocate( multipath(nCO2) )
        !*****************************!
        read(line, fmt=*) CO2_levels(lev0:lev1)

      end if




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
      call load_variable('area', varname, dimname(1), dimname(2), single_input_file=path,                           &
                         varout2D=cell_area, x=lon, y=lat, xunits=axunits(1), yunits=axunits(2), fill_missing=.true.)



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

      !*********************************!
      allocate( temp(nlon,nlat,nCO2bis) )
      !*********************************!

      ! get variable + check coordinates and units
      call load_variable('temperature', varname, dimname(1), dimname(2), multiple_input_file=multipath, &
                         varout3D=temp(:,:,lev0:lev1), xref=lon, yref=lat, fillvalue=T_fillval          )


      ! Runoff
      !-------

      print *
      print *, '          - runoff'

      !***********************************!
      allocate( runoff(nlon,nlat,nCO2bis) )
      !***********************************!

      ! get variable + check coordinates and units
      call load_variable('runoff', varname2, dimname(1), dimname(2), multiple_input_file=multipath, &
                         varout3D=runoff(:,:,lev0:lev1), xref=lon, yref=lat, fillvalue=R_fillval    )



      ! Global temperature
      !-------------------

      print *
      print *, 'Read GCM output file (for global temperature):'
      print *, '    * variable and dimension names'

      ! read dimension names
      call read_comment(1)
      read(unit=1, fmt=*) dimname(1)
      call read_comment(1)
      read(unit=1, fmt=*) dimname(2)

      ! read global temperature variable name
      call read_comment(1)
      read(unit=1, fmt='(A)') varname

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


      if (varname == 'none') then

        ! if no data specified, put fill-value on GMST
        print *, 'NO DATA'
        allocate( GMST(nCO2bis) )
        GMST = DEFFILLVAL


      else

        if (k==0) then
          print *, '      => keep previous files'
        else
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
        end if


        print *, '    * load global temperature'


        ! Sizing
        do k = 1,nCO2
            call check_sizing(multipath(k), dimname(1:2), (/nlon,nlat/) )
        end do

        !***********************************!
        allocate( glob_temp(nlon,nlat,nCO2) )
        allocate(      GMST(nCO2bis)    )
        !***********************************!

        ! get variable + check coordinates and units
        call load_variable('temperature', varname, dimname(1), dimname(2), multiple_input_file=multipath, &
                            varout3D=glob_temp, xref=lon, yref=lat, fillvalue=GT_fillval                  )

        do k = 1,nCO2
          GMST(lev0+k-1) = sum(glob_temp(:,:,k)*cell_area, mask=(glob_temp(:,:,k)/=GT_fillval))  / &
                           sum(cell_area, mask=(glob_temp(:,:,k)/=GT_fillval))
        end do

      end if
      
      close(unit=333)



      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Expand CO2 axis => extrapolate data outside the given range !
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      if (.not. fixed_CO2) then

        ! Extraploate upper and lower bounds for CO2 levels
        select case (INTERPOLATION_MODE)
          case ('linear')
            CO2_levels(1) = 0
            CO2_levels(nCO2+2) = 4*CO2_levels(nCO2+1)
          case ('log')
            CO2_levels(1) = CO2_levels(2)/16
            CO2_levels(nCO2+2) = 32*CO2_levels(nCO2+1)
        end select

        ! linear extrapolation coefficient
        xi = interp_coeff(CO2_levels(1), CO2_levels(2), CO2_levels(3))
        !
        ! temperature
        where (temp(:,:,2) == T_fillval .or. temp(:,:,3) == T_fillval)
          temp(:,:,1) = T_fillval
        else where
          temp(:,:,1) = (1-xi)*temp(:,:,2) + xi*temp(:,:,3)
        end where
        !
        ! runoff
        where (runoff(:,:,2) == R_fillval .or. runoff(:,:,3) == R_fillval)
          runoff(:,:,1) = R_fillval
        else where
          runoff(:,:,1) = (1-xi)*runoff(:,:,2) + xi*runoff(:,:,3)
          ! Avoid negative runoff:
          where (runoff(:,:,1) < 0)
            runoff(:,:,1) = 0
          end where
        end where
        !
        ! global mean surface temperature
        GMST(1) = (1-xi)*GMST(2) + xi*GMST(3)

        ! linear extrapolation coefficient
        xi = interp_coeff(CO2_levels(nCO2+2), CO2_levels(nCO2), CO2_levels(nCO2+1))
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
        !
        ! global mean surface temperature
        GMST(nCO2+2) = (1-xi)*GMST(nCO2) + xi*GMST(nCO2+1)

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

      allocate( curr_temp(nlon,nlat)       )
      allocate( curr_runf(nlon,nlat)       )
      allocate(         h(nlith,nlon,nlat) )
      allocate(         E(nlith,nlon,nlat) )
      allocate(        xs(nlith,nlon,nlat) )
      allocate(         W(nlith,nlon,nlat) )
      allocate(     W_all(nlon,nlat)       )

      ! fillvalues
      curr_temp = DEFFILLVAL
      curr_runf = DEFFILLVAL
      h         = DEFFILLVAL
      E         = DEFFILLVAL
      xs        = DEFFILLVAL
      W         = DEFFILLVAL
      W_all     = DEFFILLVAL
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
        print *, '--- multiple'
      else
        print *, '--- single'
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

      ! params(:, i_litho, i_run) =
      ! #1:  ke
      ! #2:  a
      ! #3:  b
      ! #4:  krp
      ! #5:  EA_rp
      ! #6:  T0_rp
      ! #7:  h0
      ! #8:  kd
      ! #9:  kw
      ! #10: EA
      ! #11: T0
      ! #12: sigma
      ! #13: CaMg

      print *
      print *, 'Read model parameters'

      call read_comment(1)
      read(unit=1, fmt=*) path
      call add_rootpath(path)
      open(unit=2, file=path, status='old', action='read')

      if (multirun) then
              
        call read_comment(2)
        nrun = file_length(fileunit=2)

        !*********************************!
        allocate( params(13, nlith, nrun) )
        !*********************************!
        
        ! read all the parameterizations:
        call read_comment(2)
        do k = 1,nrun
          read(unit=2, fmt=*) params(1,:,k),  params(2,:,k),  params(3,:,k),  params(4,:,k),   params(5,:,k),   params(6,:,k),  &
                              params(7,:,k),  params(8,:,k),  params(9,:,k),  params(10,:,k),  params(11,:,k),  params(12,:,k), &
                              params(13,:,k)
        end do

        close(unit=2)

      else

        nrun = 1

        !*********************************!
        allocate( params(13, nlith, nrun) )
        !*********************************!

        ! Get parameters values:
        call read_comment(2)
        read(unit=2, fmt=*) params(1,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(2,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(3,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(4,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(5,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(6,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(7,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(8,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(9,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(10,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(11,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(12,:,1)
        call read_comment(2)
        read(unit=2, fmt=*) params(13,:,1)

        close(unit=2)

      end if


      !---------------!
      ! Forcing file: !
      !---------------!

      !*********************!
      allocate(forcing(nrun))
      !*********************!

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

        forcing = CO2_levels(1)


      else

        open(unit=2, file=path, status='old', action='read')
        call read_comment(2)
        length = file_length(fileunit=2)
        call read_comment(2)

        ! Compare the length of forcing file and number of runs (from parameters file)
        if (length==nrun) then

          do k = 1,nrun
            read(unit=2, fmt=*) forcing(k)
          end do
          close(unit=2)

        else

          ! In case of length mismatch, if the length of forcing file is 1, consider constant forcing
          if (length==1) then
            read(unit=2, fmt=*) forcing(1)
            close(unit=2)
            forcing = forcing(1)

          else ! Otherwise, raise error
            print *
            print *, 'ERROR: inconsistent number of parameterizations and forcings. File length mismatch'
            stop
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
      read(unit=1, fmt=*) output_info%file_name
      call add_rootpath(output_info%file_name)

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
        call create_output_file(output_info%file_name, output_info%file_id, &
                                dimname(1:4), (/axunits(1),axunits(2),nounits,nounits/), x1=lon, x2=lat, nx3=nlith)
      else
        call create_output_file(output_info%file_name, output_info%file_id, &
                                dimname(1:3), (/axunits(1),axunits(2),nounits/), x1=lon, x2=lat, nx3=nlith)
      end if

      ierr = nf90_put_att(output_info%file_id, NF90_GLOBAL, 'CO2_interpolation', INTERPOLATION_MODE)
      call nf90_check(ierr, 'Error while putting global attribute "CO2 interpolation"')


      !-------------------------!
      ! Define output variables !
      !-------------------------!

      if (multirun) then
        ndim = 4
      else
        ndim = 3
      end if

      ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !
      call define_variables( output_info, dimname(1:ndim), DEFFILLVAL, axunits )
      ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !


      !--------------------!
      ! End of define mode !
      !--------------------!
      ierr = nf90_enddef(output_info%file_id)
      call nf90_check(ierr, 'Error while end of definition mode in output file '//path)



      !!======================================!!
      !!  close input-output conditions file  !!
      !!======================================!!

      close(unit=1)



    end subroutine




    !==============================================================================================================================!
    !------------------------------------------------------------------------------------------------------------------------------!
    !==============================================================================================================================!




    !# Secondary subroutines:
    !########################


    subroutine add_rootpath(filename)
      character(len=FILE_CHARLEN), intent(inout):: filename
      include 'path.inc' ! => 'gdss_root_path' variable, path of root directory
      !
      ! test if absolute path given (file name starts with a '/')
      ! if not, add current root path to the name:
      if (filename(1:1) /= '/') then
        filename = gdss_root_path//trim(filename)
      end if
    end subroutine



    !==============================================================================================================================!



    function read_answer(valid_range, command_argnum)

      integer:: read_answer
      integer, intent(in), dimension(2):: valid_range
      integer, intent(in), optional:: command_argnum
      character(len=1):: answerchar
      integer:: ierr

      ! initialize with invalid answer:
      read_answer = valid_range(1) - 1

      ! read potential command argument
      if (present(command_argnum)) then
        call get_command_argument(command_argnum, answerchar, status=ierr)
        read(answerchar, fmt=*, iostat=ierr) read_answer
        if (ierr==0)  print *, read_answer
      end if

      ! reading loop (in case nor valid argument found)
      do while (read_answer<valid_range(1) .or. read_answer>valid_range(2))
        read(unit=*, fmt=*) read_answer
      end do

    end function



    !==============================================================================================================================!



    subroutine load_variable(internal_varname, varname, x_varname, y_varname, z_varname,       &
                             single_input_file, multiple_input_file,                           &
                             varout2D, varout3D, x, y, xvec, yvec, xref, yref, xunits, yunits, &
                             totarea, fill_missing, fillvalue                                  )
        use netcdf

        character(len=*), intent(in):: internal_varname
        character(len=LINE_CHARLEN), intent(in):: varname
        character(len=OTHR_CHARLEN), intent(in):: x_varname, y_varname
        character(len=OTHR_CHARLEN), intent(in), optional:: z_varname
        character(len=FILE_CHARLEN), intent(in), optional:: single_input_file
        character(len=FILE_CHARLEN), dimension(:),  intent(in),  optional:: multiple_input_file
        character(len=OTHR_CHARLEN), intent(out), optional:: xunits, yunits
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
        character(len=OTHR_CHARLEN):: var_units, loc_xunits, loc_yunits
        character(len=LINE_CHARLEN):: vname
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
                    call load_netcdf_dble2D(single_input_file, x_varname, y_varname, vname, dummyvar2D,                 &
                                            x=loc_x, y=loc_y, xunits=loc_xunits, yunits=loc_yunits, varunits=var_units, &
                                            fillval=fillvalue_vec1D(n), fillval_iostat=ierr                             )
                    ! set var=0 on "missing" cells
                    if (ierr==NF90_NOERR) then
                        where (dummyvar2D==fillvalue_vec1D(n)) dummyvar2D = 0d0
                    else
                        fillvalue_vec1D(n) = -1d99 ! dummy fillvalue
                    end if
                else
                    ! load variable
                    call load_netcdf_dble2D(single_input_file, x_varname, y_varname, vname, dummyvar2D,                 &
                                            x=loc_x, y=loc_y, xunits=loc_xunits, yunits=loc_yunits, varunits=var_units, &
                                            fillval=fillvalue_vec1D(n)                                                  )
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
                        call load_netcdf_dble3D(single_input_file, x_varname, y_varname, z_varname, vname, dummyvar3D,      &
                                                x=loc_x, y=loc_y, xunits=loc_xunits, yunits=loc_yunits, varunits=var_units, &
                                                fillval=fillvalue_vec1D(n), &
                                                fillval_iostat=ierr)
                        ! set var=0 on "missing" cells
                        if (ierr==NF90_NOERR) then
                            where (dummyvar3D==fillvalue_vec1D(n)) dummyvar3D = 0d0
                        else
                            fillvalue_vec1D(n) = -1d99 ! dummy fillvalue
                        end if
                    else
                        ! load variable
                        call load_netcdf_dble3D(single_input_file, x_varname, y_varname, z_varname, vname, dummyvar3D,      &
                                                x=loc_x, y=loc_y, xunits=loc_xunits, yunits=loc_yunits, varunits=var_units, &
                                                fillval=fillvalue_vec1D(n))
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
                                            x=loc_x, y=loc_y, xunits=loc_xunits, yunits=loc_yunits, varunits=var_units,     &
                                            fillval=fillvalue_vec2D(k,n), fillval_iostat=ierr                               )
                            ! set var=0 on "missing" cells
                            if (ierr==NF90_NOERR) then
                                where (dummyvar3D(:,:,k)==fillvalue_vec2D(k,n)) dummyvar3D(:,:,k) = 0d0
                            else
                                fillvalue_vec2D(k,n) = -1d99 ! dummy fillvalue
                            end if
                        else
                            ! load variable
                            call load_netcdf_dble2D(multiple_input_file(k), x_varname, y_varname, vname, dummyvar3D(:,:,k),     &
                                                    x=loc_x, y=loc_y, xunits=loc_xunits, yunits=loc_yunits, varunits=var_units, &
                                                    fillval=fillvalue_vec2D(k,n)                                                )
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

        if (present(xunits))  xunits = loc_xunits
        if (present(yunits))  yunits = loc_yunits


        ! close scratch files
        close(unit=334)
        close(unit=335)


    end subroutine


    !=======================!


    function get_arithmetic_operations(varstring)
    ! => create scratch files 334 and 335 containing (respectively) the lists of variables and operators
    ! Note: the first operators (corresponding to the first variable) is automatically '+'
    ! return the number of variables/operations
        character(len=LINE_CHARLEN), intent(in):: varstring
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


    subroutine load_netcdf_dble2D(fname, x_varname, y_varname, varname, var, x, y, &
                                  xunits, yunits, varunits, fillval, fillval_iostat)

        use netcdf
        use netcdf_io_functions, only: nf90_check

        character(len=FILE_CHARLEN), intent(in):: fname
        character(len=OTHR_CHARLEN), intent(in):: x_varname, y_varname
        character(len=LINE_CHARLEN), intent(in):: varname
        double precision, dimension(:,:), intent(out):: var
        double precision, dimension(:), intent(out), optional:: x, y
        character(len=OTHR_CHARLEN), intent(out), optional:: xunits, yunits, varunits
        double precision, intent(out), optional:: fillval
        integer, intent(out), optional:: fillval_iostat
        integer:: nx, ny
        integer:: ierr, fid, varid, truedimids(2), shp_idx(2)
        integer, dimension(:), allocatable:: dimids
        logical:: transp
        double precision, dimension(:,:), allocatable:: loc_var
        character(len=OTHR_CHARLEN):: loc_xunits, loc_yunits

        nx = size(var, 1)
        ny = size(var, 2)
        ! must also match the size of "x" and "y"


        ! variable shape-independent operations (open file, get variable ID, ...)
        call load_netcdf_generic(fname, varname, fid, varid, dimids, truedimids, shp_idx, &
                                 x_varname=x_varname, y_varname=y_varname, x=x, y=y,      &
                                 xunits=loc_xunits, yunits=loc_yunits, varunits=varunits,  &
                                 fillval=fillval, fillval_iostat=fillval_iostat)

        if (present(xunits))  xunits = loc_xunits
        if (present(yunits))  yunits = loc_yunits

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


    subroutine load_netcdf_dble3D(fname, x_varname, y_varname, z_varname, varname, var, x, y, z, &
                                  xunits, yunits, varunits, fillval, fillval_iostat)

        use netcdf
        use netcdf_io_functions, only: nf90_check

        character(len=FILE_CHARLEN), intent(in):: fname
        character(len=OTHR_CHARLEN), intent(in):: x_varname, y_varname, z_varname
        character(len=LINE_CHARLEN), intent(in):: varname
        double precision, dimension(:,:,:), intent(out):: var
        double precision, dimension(:), intent(out), optional:: x, y, z
        character(len=OTHR_CHARLEN), intent(out), optional:: xunits, yunits, varunits
        double precision, intent(out), optional:: fillval
        integer, intent(out), optional:: fillval_iostat
        integer:: nx, ny, nz
        integer:: ierr, fid, varid, truedimids(3), shp_idx(3)
        integer, dimension(:), allocatable:: dimids
        character(len=OTHR_CHARLEN):: loc_xunits, loc_yunits

        nx = size(var, 1)
        ny = size(var, 2)
        nz = size(var, 3)
        ! must also match the size of "x", "y" and "z"


        ! variable shape-independent operations (open file, get variable ID, ...)
        call load_netcdf_generic(fname, varname, fid, varid, dimids, truedimids, shp_idx,                      &
                                 x_varname=x_varname, y_varname=y_varname, z_varname=z_varname, x=x, y=y, z=z, &
                                 xunits=loc_xunits, yunits=loc_yunits, varunits=varunits,                      &
                                 fillval=fillval, fillval_iostat=fillval_iostat                                )

        if (present(xunits))  xunits = loc_xunits
        if (present(yunits))  yunits = loc_yunits

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
                                   x_varname, y_varname, z_varname, x, y, z, xunits, yunits, varunits, fillval, fillval_iostat)

        use netcdf
        use netcdf_io_functions, only: nf90_check

        character(len=FILE_CHARLEN), intent(in):: fname
        character(len=LINE_CHARLEN), intent(in):: varname
        integer, intent(out):: fid, varid
        integer, dimension(:), intent(out):: truedimids ! must be dimension(2) or dimension(3)
        integer, dimension(:), intent(out), allocatable:: dimids
        integer, dimension(:), intent(out):: shp_idx
        character(len=OTHR_CHARLEN), intent(in), optional:: x_varname, y_varname, z_varname
        double precision, dimension(:), intent(out), optional:: x, y, z
        character(len=OTHR_CHARLEN), intent(out), optional:: xunits, yunits, varunits
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
            call load_axis(fname, fid, x_varname, truedimids(1))
            if (present(xunits)) call load_axis(fname, fid, x_varname, truedimids(1), axunits=xunits)
            if (present(x))      call load_axis(fname, fid, x_varname, truedimids(1), ax=x)
        end if
        if (present(y_varname)) then
            call load_axis(fname, fid, y_varname, truedimids(2))
            if (present(yunits)) call load_axis(fname, fid, y_varname, truedimids(2), axunits=yunits)
            if (present(y))      call load_axis(fname, fid, y_varname, truedimids(2), ax=y)
        end if
        if (present(z_varname) .and. ndim_true==3) then
            call load_axis(fname, fid, z_varname, truedimids(3))
            if (present(z)) call load_axis(fname, fid, z_varname, truedimids(3), ax=z)
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


    subroutine load_axis(fname, fid, axname, axdimid, axunits, ax)
        use netcdf
        use netcdf_io_functions, only: nf90_check
        character(len=FILE_CHARLEN), intent(in):: fname
        character(len=OTHR_CHARLEN), intent(in):: axname
        integer, intent(in):: fid
        integer, intent(out):: axdimid
        character(len=OTHR_CHARLEN), intent(out), optional:: axunits
        double precision, dimension(:), intent(out), optional:: ax
        integer:: varid, ierr

        ! Inquire dim ID
        ierr = nf90_inq_dimid(fid, axname, axdimid)
        call nf90_check(ierr, 'Error while inquiring dimension "'//trim(axname)//'" ID in file "'//trim(fname)//'"')
        ! Inquire var ID
        ierr = nf90_inq_varid(fid, axname, varid)
        call nf90_check(ierr, 'Error while inquiring variable "'//trim(axname)//'" ID in file "'//trim(fname)//'"')
        ! Get variable units
        if (present(axunits)) then
            ierr = nf90_get_att(fid, varid, 'units', axunits)
            call nf90_check(ierr, 'Warning: unable to get attribute "units" of variable "'//trim(axname)//'" in file "' &
                            //trim(fname)//'"', kill=.false.)
        end if
        ! Get variable
        if (present(ax)) then
            ierr = nf90_get_var(fid, varid, ax)
            call nf90_check(ierr, 'Error while getting variable "'//trim(axname)//'" in file "'//trim(fname)//'"')
        end if

    end subroutine



    !==============================================================================================================================!



    subroutine define_variables( output_info, dimname, fillval, axunits )

      use netcdf_io_functions, only: create_output_variable
      use ascii_io_functions, only: read_comment

      type(outinfo), intent(inout):: output_info
      character(len=OTHR_CHARLEN), dimension(:), intent(in):: dimname, axunits
      double precision, intent(in):: fillval
      character(len=100):: varname
      character(len=500):: varlongname
      character(len=50):: units
      integer, dimension(N_TOT_DIM):: defdim
      logical:: write_var
      integer:: fid, varid, k, ndim, nvar, ierr

      character(len=*), parameter:: loc_fmt = '(A100, A500, A50, L2, 4I2)'
      !       IMPORTANT: update this if N_TOT_DIM is not 4 ----------^^^


      fid = output_info%file_id

      ndim = size(dimname)


      call read_comment(1, ierr=ierr)

      ! Scratch file storing variables info before number of variable is known
      open(unit=333, status='scratch', action='readwrite')

      nvar = 0
      do while(ierr==0) ! read until the end of the interface file!

        nvar = nvar + 1

        read(unit=1, fmt=*) varname, units, write_var, defdim, varlongname

        ! A given units $* means: use the units of the axis number *
        if (units(1:1)=='$') then
          read(units(2:50), fmt=*) k
          units = axunits(k)
        end if

        write(unit=333, fmt=loc_fmt) varname, varlongname, units, write_var, defdim

        call read_comment(1, ierr=ierr)

      end do

      if (ierr>0) then ! Error other than end-of-file has been raised
        backspace(1)
        call read_comment(1) ! => trigger error printing
      end if


      !+++++++++++++++++++++++++++++++++++++!
      allocate( output_info%variables(nvar) )
      !+++++++++++++++++++++++++++++++++++++!

      ! Loop again to create variables
      rewind(unit=333)
      do k = 1,nvar

        read(unit=333, fmt=loc_fmt) varname, varlongname, units, write_var, defdim

        if (write_var) then
          call create_output_variable(fid, varname, dimname, defdim(1:ndim), units, fillval, long_name=varlongname, varid=varid)
        else
          varid = 0
          ! make sure the two booleans in output_info file will be .false. and the program will not try to write the variable
          defdim(N_TOT_DIM) = 0
        end if

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        output_info%variables(k)%write_var = write_var
        output_info%variables(k)%unlimited = (defdim(N_TOT_DIM)==1)
        output_info%variables(k)%varid     = varid
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      end do

      close(unit=333)


    end subroutine



    !==============================================================================================================================!



    subroutine check_sizing( path, dimnames, sizes )

      use netcdf_io_functions, only: netcdf_get_size

      character(len=FILE_CHARLEN), intent(in):: path
      character(len=OTHR_CHARLEN), dimension(:), intent(in):: dimnames
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

        character(len=OTHR_CHARLEN), intent(in):: units_string
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

        character(len=*), intent(in):: which_variable
        character(len=LINE_CHARLEN), intent(in):: varname
        character(len=OTHR_CHARLEN), intent(in):: varunits
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

        character(len=*), intent(in):: which_variable
        character(len=LINE_CHARLEN), intent(in):: varname
        character(len=OTHR_CHARLEN), intent(in):: varunits
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

        character(len=*), intent(in):: which_variable
        character(len=LINE_CHARLEN), intent(in):: varname
        character(len=OTHR_CHARLEN), intent(in):: varunits
        double precision, intent(out):: factor, offset
        logical, intent(out):: multiply_by_area
        type(units):: known_units
        logical:: passed
        integer:: answer


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
            print *
            print *, 'Hit one of the following options:'
            print *, '    0: abort the program'
            print *, '    1: ignore issue and continue execution'

            answer = read_answer((/0, 2/), command_argnum=1)
            if (answer==0) stop

        end if

    end subroutine



    !==============================================================================================================================!



    subroutine check_coordinates_singlevar(varname, x, xref)
        character(len=LINE_CHARLEN), intent(in):: varname
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
        character(len=LINE_CHARLEN), intent(in):: varname
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



    subroutine check_continental_cells_single2Dvar( varname, landarea, var, check, new_landarea )

      character(len=*), intent(in):: varname
      double precision, dimension(:,:), intent(in):: landarea, var
      logical, intent(out):: check
      double precision, dimension(:,:), intent(inout), optional:: new_landarea
      double precision:: area_err, tot_landarea
      integer:: nerr

      nerr = count( (landarea>0 .and. var==DEFFILLVAL) )
      area_err = sum( landarea, mask=(landarea>0 .and. var==DEFFILLVAL) )
      tot_landarea = sum(landarea)
      if (present(new_landarea)) then
        where (landarea>0 .and. var==DEFFILLVAL) new_landarea = 0
      end if

      if (nerr > 0) then
        print *
        print *, 'WARNING: found missing values on continental cells of variable '//trim(varname)
        print *, 'Number of continent cells with missing values:     ', nerr
        print *, 'Fraction of continental cells with missing values: ', dble(nerr)/dble(count((landarea>0)))
        print *, 'Total area of those cells (m2):                    ', area_err
        print *, 'Which is a fraction of total land area:            ', area_err/tot_landarea
        print *
        check = .false.
      else
        check = .true.
      end if

    end subroutine


    !=======================!


    subroutine check_continental_cells_single3Dvar( varname, landarea, var, check, new_landarea )

      character(len=*), intent(in):: varname
      double precision, dimension(:,:), intent(in):: landarea
      double precision, dimension(:,:,:), intent(in):: var
      logical, intent(out):: check
      double precision, dimension(:,:), intent(inout), optional:: new_landarea
      double precision:: area_err, tot_landarea
      integer:: k, nerr

      check = .true.
      tot_landarea = sum(landarea)

      do k = 1,size( var, 3 )

        nerr = count( (landarea>0 .and. var(:,:,k)==DEFFILLVAL) )
        area_err = sum( landarea, mask=(landarea>0 .and. var(:,:,k)==DEFFILLVAL) )
        if (present(new_landarea)) then
          where (landarea>0 .and. var(:,:,k)==DEFFILLVAL) new_landarea = 0
        end if

        if (nerr > 0) then
          print *
          print *, 'WARNING: found missing values on continental cells of variable '//trim(varname)
          print *, '"Vertical" (3rd axis) level number: ', k
          print *, 'Number of continent cells with missing values:     ', nerr
          print *, 'Fraction of continental cells with missing values: ', dble(nerr)/dble(count((landarea>0)))
          print *, 'Total area of those cells (m2):                    ', area_err
          print *, 'Which is a fraction of total land area:            ', area_err/tot_landarea
          print *
          check = .false.
        end if

      end do

    end subroutine


    !=======================!


    subroutine check_invalid(landarea, runoff, slope, lith_frac)

      double precision, dimension(:,:), intent(inout):: landarea, slope
      double precision, dimension(:,:,:), intent(inout):: runoff, lith_frac
      logical, dimension(:,:,:), allocatable:: errormask
      logical, dimension(:,:), allocatable:: errormask2D
      double precision:: area_err, tot_landarea, ex_invalid
      integer:: nerr, answer, k, nx, ny, nCO2, nlith


      nx    = size(runoff, 1)
      ny    = size(runoff, 2)
      nCO2  = size(runoff, 3)
      nlith = size(lith_frac, 3)
      allocate( errormask(nx, ny, nCO2) )
      allocate( errormask2D(nx, ny) )

      tot_landarea = sum(landarea)


      ! Check if negative runoff
      ! ------------------------

      do k = 1,nCO2
        errormask(:,:,k) = (landarea>0 .and. runoff(:,:,k)<0)
      end do
      errormask2D = any(errormask, dim=3)
      nerr        = count(errormask2D)
      area_err    = sum(landarea, mask=errormask2D)
      ex_invalid  = minval(runoff, mask=errormask)

      if (nerr > 0) then
        print *
        print *, 'WARNING: found negative runoff'
        print *, 'Number of cells affected:               ', nerr
        print *, 'Fraction of continental cells affected: ', dble(nerr)/dble(count((landarea>0)))
        print *, 'Total area of those cells (m2):         ', area_err
        print *, 'Which is a fraction of total land area: ', area_err/tot_landarea
        print *, 'Minimum runoff found (m/y):             ', ex_invalid
        print *
        print *, 'Hit one of the following options:'
        print *, '    0: abort the program'
        print *, '    1: remove all the erratic points (set land_area=0)'
        print *, '    2: replace negative runoff by 0'

        answer = read_answer((/0, 2/), command_argnum=3)

        ! According to the user answer:
        select case (answer)
          case(0)
            stop
          case(1)
            where (errormask2D)  landarea = 0
          case(2)
            where (errormask)  runoff = 0
        end select

      end if


      ! Check if negative or null slope
      ! -------------------------------

      errormask2D = (landarea>0 .and. slope<=0)
      nerr       = count(errormask2D)
      area_err   = sum(landarea, mask=errormask2D) 
      ex_invalid = minval(slope, mask=errormask2D)

      if (nerr > 0) then
        print *
        print *, 'WARNING: found negative or null slope'
        print *, 'Number of cells affected:               ', nerr
        print *, 'Fraction of continental cells affected: ', dble(nerr)/dble(count((landarea>0)))
        print *, 'Total area of those cells (m2):         ', area_err
        print *, 'Which is a fraction of total land area: ', area_err/tot_landarea
        print *, 'Minimum slope found (m/m):              ', ex_invalid
        print *
        print *, 'Hit one of the following options:'
        print *, '    0: abort the program'
        print *, '    1: remove all the erratic points (set land_area=0)'
        print *, '    2: replace by minimum non-null slope'

        answer = read_answer((/0, 2/), command_argnum=3)

        ! According to the user answer:
        select case (answer)
          case(0)
            stop
          case(1)
            where (errormask2D)  landarea = 0
          case(2)
            where (errormask2D)  slope = minval(slope, mask=(landarea>0 .and. slope>0))
        end select

      end if


      ! Check lithology fractions
      ! -------------------------------

      errormask2D = (landarea > 0  .and.  abs(sum(lith_frac, dim=3) - 1) > MAX_ALLOWED_INACC)
      nerr       = count(errormask2D)
      area_err   = sum(landarea, mask=errormask2D) 
      ex_invalid = maxval(abs(sum(lith_frac, dim=3) - 1), mask=errormask2D)

      if (nerr > 0) then
        print *
        print *, 'WARNING: inconsistent lithology class fraction'
        print *, 'Found cells where the sum of all lithology fractions differs from 1'
        print *, 'Number of cells affected:               ', nerr
        print *, 'Fraction of continental cells affected: ', dble(nerr)/dble(count((landarea>0)))
        print *, 'Total area of those cells (m2):         ', area_err
        print *, 'Which is a fraction of total land area: ', area_err/tot_landarea
        print *, 'Maximum absolute difference:            ', ex_invalid
        print *
        print *, 'Hit one of the following options:'
        print *, '    0: abort the program'
        print *, '    1: remove all the erratic points (set land_area=0)'
        print *, '    2: ignore issue'

        answer = read_answer((/0, 2/), command_argnum=3)

        ! According to the user answer:
        select case (answer)
          case(0)
            stop
          case(1)
            where (errormask2D)  landarea = 0
        end select

      end if


    deallocate(errormask)
    deallocate(errormask2D)

    end subroutine


    !=======================!


    subroutine check_continental_cells( landarea, temp, runoff, slope, lith_frac )

      double precision, dimension(:,:),   intent(inout):: landarea, slope
      double precision, dimension(:,:,:), intent(inout):: temp, runoff, lith_frac
      double precision, dimension(:,:), allocatable:: landarea_corr
      logical, dimension(4):: checkpoints
      integer:: answer

      allocate( landarea_corr(size(landarea,1), size(landarea,2)) )

      landarea_corr = landarea
      call check_continental_cells_single3Dvar( 'temperature', landarea, temp,      checkpoints(1), landarea_corr )
      call check_continental_cells_single3Dvar( 'runoff',      landarea, runoff,    checkpoints(2), landarea_corr )
      call check_continental_cells_single2Dvar( 'slope',       landarea, slope,     checkpoints(3), landarea_corr )
      call check_continental_cells_single3Dvar( 'lithology',   landarea, lith_frac, checkpoints(4), landarea_corr )

      if (.not. all(checkpoints)) then

        print *, 'Hit one of the following options:'
        print *, '    0: abort the program'
        print *, '    1: remove all the erratic points (set land_area=0)'
        print *, '    2: ignore issues'

        answer = read_answer((/0, 2/), command_argnum=2)

        ! According to the user answer:
        select case (answer)
          case(0)
            stop
          case(1)
            landarea = landarea_corr
        end select

      end if

      deallocate(landarea_corr)

      call check_invalid(landarea, runoff, slope, lith_frac)

    end subroutine



    !==============================================================================================================================!




end module

