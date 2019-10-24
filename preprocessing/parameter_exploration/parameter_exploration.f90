program parameter_exploration

        use dynsoil, only: dynsoil_steady_state
        use dynsoil_empirical_laws, only: erosion, reg_prod_opt, eq_reg_thick
        use netcdf_io_functions, only: nf90_check
        use netcdf
        implicit none



        !================================!
        !=========  DECLARATION =========!
        !================================!

        include 'dynsoil_physical_parameters.inc'

        integer:: n, i, j, k, l, i1,i2,i3,i4,i5

        ! Activate it to make a test with only one parameterization
        logical, parameter:: TEST_SINGLE_PARAMETERIZATION = .true.

        ! Geographic fields, input files and variables names

        integer, parameter:: nvar=8

        logical, parameter:: extend_output_file=.false.  ! .false. to create a new output file  -  .true. to extend an existing output file
        character(len=*), parameter:: output_file_name='output/parameter_exploration_test.nc'

        character(len=*), dimension(nvar), parameter:: fname = &
                (/ &
                '/Users/hematite/0000_Github_Projects/GEOCLIM_Modern/Code/Output/cell_area_360_720.nc',       &
                '/Users/hematite/0000_Github_Projects/GEOCLIM_Modern/Code/Output/cell_area_360_720.nc',       &
                '/Users/hematite/0000_Github_Projects/GEOCLIM_Modern/Code/Output/modern_temp_interp.nc',      &
                '/Users/hematite/0000_Github_Projects/GEOCLIM_Modern/Code/Output/modern_runoff_v2_interp.nc', &
                '/Users/hematite/0000_Github_Projects/GEOCLIM_Modern/Code/Output/modern_slope.nc',            &
                '/Users/hematite/0000_Github_Projects/GEOCLIM_Modern/Code/Output/cell_area_360_720.nc',       &
                '/Users/hematite/0000_Github_Projects/GEOCLIM_Modern/Code/Output/lith_mask_360_720.nc',       &
                '/Users/hematite/0000_Github_Projects/GEOCLIM_Modern/Code/Output/land_area_360_720.nc'        &
                /)
        character(len=*), dimension(nvar), parameter:: varname = &
                (/ 'lon', 'lat', 'tmp', 'rnf', 'slope', 'area', 'frac', 'area' /)
        character(len=*), parameter:: litho_dimname = 'lith'
        character(len=*), dimension(nvar), parameter:: missvalname = (/ ('_FillValue', i=1,nvar) /)

        character(len=*), parameter:: Kwest_fname    = 'parameter_list/Kwest.txt' 
        character(len=*), parameter:: kw_fname       = 'parameter_list/kw.txt'
        character(len=*), parameter:: sigma_fname    = 'parameter_list/sigma.txt'
        character(len=*), parameter:: CaMgBulk_fname = 'parameter_list/CaMg_BulkContCrust.txt'
        character(len=*), parameter:: CaMgSed_fname  = 'parameter_list/CaMg_Sediment.txt'
        

        ! Missing-value
        double precision, parameter:: deffillval = 9.96921d+36
        double precision, dimension(nvar):: fillval

        ! Variables:
        double precision, dimension(:), allocatable:: lon, lat
        double precision, dimension(:,:), allocatable:: temperature, runoff, slope, h, E, cell_area, land_area, W
        double precision, dimension(:,:,:), allocatable:: litho_frac, litho_frac2
        double precision:: RPopt, Esum, WSI, WCaMg
        logical, dimension(:,:), allocatable:: missingpoints
        double precision, parameter:: veget_factor=1d0, veget_eros_factor=1d0
        integer:: nlon0, nlon, nlat0, nlat, nlith, ncontpxl

        integer:: NKwest, Nkw, Nsigma, NCaMgBulk, NCaMgSed
        double precision, dimension(:), allocatable:: list_Kwest, list_kw, list_sigma, list_CaMgBulk, list_CaMgSed
        integer, dimension(:), allocatable:: list_cont_i, list_cont_j

        integer:: fid, dimid(4), varid(23), ierr





        !==================================================================================================!
        !===================================   LOAD PARAMETERS LISTS   ====================================!
        !==================================================================================================!


        !*******************************************************!
        !**************  Fixed-value parameters:  **************!
        !*******************************************************!
        Ea_rp = 42000.d0
        T0_rp = 286.d0
        Ea = 42000.d0
        T0 = 286.d0
        krp = 1.d-2 !1.d-2 !5.d-3 !3.d-3 !2.d-3 !1.2d-3
        h0 = 2.73
        ke = 0.0029110 ! for 0.5deg x 0.5deg
        a = 0.5
        b = 1.0
        CaMg_rock = (/ -1.0, 1521.0, 4759.0, 10317.0, 0.0, -1.0 /)
        !*******************************************************!


        !**** read lists of parameters to test in files ****!

        open(unit=1,file=Kwest_fname,status='old',action='read')
        !
        NKwest = -1
        ierr = 0
        do while (ierr==0)
          read(unit=1,fmt=*,iostat=ierr)
          NKwest = NKwest + 1
        end do
        if (ierr>0) then
          print *, 'Error while reading file '//Kwest_fname
          stop
        end if
        !
        allocate( list_Kwest(NKwest) )
        !
        rewind(unit=1)
        do k = 1,NKwest
          read(unit=1,fmt=*) list_Kwest(k)
        end do
        !
        close(unit=1)

        open(unit=1,file=kw_fname,status='old',action='read')
        !
        Nkw = -1
        ierr = 0
        do while (ierr==0)
          read(unit=1,fmt=*,iostat=ierr)
          Nkw = Nkw + 1
        end do
        if (ierr>0) then
          print *, 'Error while reading file '//kw_fname
          stop
        end if
        !
        allocate( list_kw(Nkw) )
        !
        rewind(unit=1)
        do k = 1,Nkw
          read(unit=1,fmt=*) list_kw(k)
        end do
        !
        close(unit=1)

        open(unit=1,file=sigma_fname,status='old',action='read')
        !
        Nsigma = -1
        ierr = 0
        do while (ierr==0)
          read(unit=1,fmt=*,iostat=ierr)
          Nsigma = Nsigma + 1
        end do
        if (ierr>0) then
          print *, 'Error while reading file '//sigma_fname
          stop
        end if
        !
        allocate( list_sigma(Nsigma) )
        !
        rewind(unit=1)
        do k = 1,Nsigma
          read(unit=1,fmt=*) list_sigma(k)
        end do
        !
        close(unit=1)

        open(unit=1,file=CaMgBulk_fname,status='old',action='read')
        !
        NCaMgBulk = -1
        ierr = 0
        do while (ierr==0)
          read(unit=1,fmt=*,iostat=ierr)
          NCaMgBulk = NCaMgBulk + 1
        end do
        if (ierr>0) then
          print *, 'Error while reading file '//CaMgBulk_fname
          stop
        end if
        !
        allocate( list_CaMgBulk(NCaMgBulk) )
        !
        rewind(unit=1)
        do k = 1,NCaMgBulk
          read(unit=1,fmt=*) list_CaMgBulk(k)
        end do
        !
        close(unit=1)

        open(unit=1,file=CaMgSed_fname,status='old',action='read')
        !
        NCaMgSed = -1
        ierr = 0
        do while (ierr==0)
          read(unit=1,fmt=*,iostat=ierr)
          NCaMgSed = NCaMgSed + 1
        end do
        if (ierr>0) then
          print *, 'Error while reading file '//CaMgSed_fname
          stop
        end if
        !
        allocate( list_CaMgSed(NCaMgSed) )
        !
        rewind(unit=1)
        do k = 1,NCaMgSed
          read(unit=1,fmt=*) list_CaMgSed(k)
        end do
        !
        close(unit=1)

        !***************************************************!





        !==================================================================================================!
        !==========================   LOAD CLIMATIC AND GEOGRAPHIC CONDITIONS   ===========================!
        !==================================================================================================!


        ! Get and check dimension sizes:

        ierr = nf90_open(fname(1), NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//fname(1))
        ierr = nf90_inq_dimid( fid, varname(1), dimid(1) )
        call nf90_check(ierr, 'Error while inquiring ID of dimension '//varname(1)//' in file '//fname(1))
        ierr = nf90_inquire_dimension( fid, dimid(1), len=nlon0 )
        call nf90_check(ierr, 'Error while inquiring length of dimension '//varname(1)//' in file '//fname(1))
        ierr = nf90_inq_dimid( fid, varname(2), dimid(1) )
        call nf90_check(ierr, 'Error while inquiring ID of dimension '//varname(2)//' in file '//fname(1))
        ierr = nf90_inquire_dimension( fid, dimid(1), len=nlat0 )
        call nf90_check(ierr, 'Error while inquiring length of dimension '//varname(2)//' in file '//fname(1))
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//fname(1))

        do k = 2,nvar

          ierr = nf90_open(fname(k), NF90_NOWRITE, fid)
          call nf90_check(ierr, 'Error while openning file '//fname(k))
          ierr = nf90_inq_dimid( fid, varname(1), dimid(1) )
          call nf90_check(ierr, 'Error while inquiring ID of dimension '//varname(1)//' in file '//fname(k))
          ierr = nf90_inquire_dimension( fid, dimid(1), len=nlon )
          call nf90_check(ierr, 'Error while inquiring length of dimension '//varname(1)//' in file '//fname(k))
          ierr = nf90_inq_dimid( fid, varname(2), dimid(1) )
          call nf90_check(ierr, 'Error while inquiring ID of dimension '//varname(2)//' in file '//fname(k))
          ierr = nf90_inquire_dimension( fid, dimid(1), len=nlat )
          call nf90_check(ierr, 'Error while inquiring length of dimension '//varname(2)//' in file '//fname(k))
          ierr = nf90_close(fid)
          call nf90_check(ierr, 'Error while closing file '//fname(k))

          if ( nlon/=nlon0 .or. nlat/=nlat0 ) then
            print *, 'Error: inconsistent file dimension'
            stop
          end if

        end do

        ierr = nf90_open(fname(7), NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//fname(k))
        ierr = nf90_inq_dimid( fid, litho_dimname, dimid(1) )
        call nf90_check(ierr, 'Error while inquiring ID of dimension '//litho_dimname//' in file '//fname(7))
        ierr = nf90_inquire_dimension( fid, dimid(1), len=nlith )
        call nf90_check(ierr, 'Error while inquiring length of dimension '//litho_dimname//' in file '//fname(7))
        nlith = nlith-1
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//fname(7))

        if ( nlith/=size(CaMg_rock) ) then
          print '(A,I2,A,I2,A)', 'Error: inconsistent number of lithology classes. Expected ', &
                                  size(CaMg_rock),', got ',nlith,' (water and ice removed)'
          stop
        end if


        ! Allocate variables
        allocate( lon(nlon) )
        allocate( lat(nlat) )
        allocate( temperature(nlon,nlat) )
        allocate( runoff(nlon,nlat) )
        allocate( slope(nlon,nlat) )
        allocate( h(nlon,nlat) )
        allocate( E(nlon,nlat) )
        allocate( cell_area(nlon,nlat) )
        allocate( land_area(nlon,nlat) )
        allocate( W(nlon,nlat) )
        allocate( litho_frac(nlon,nlat,nlith) )
        allocate( litho_frac2(nlith,nlon,nlat) )
        allocate( missingpoints(nlon,nlat) )
        allocate( list_cont_i(nlon*nlat) )
        allocate( list_cont_j(nlon*nlat) )


        ! Get longitude and latitude:
        ierr = nf90_open(fname(1), NF90_NOWRITE, fid)
        ierr = nf90_inq_varid( fid, varname(1), varid(1) )
        call nf90_check(ierr, 'Error while inquiring ID of variable '//varname(1)//' in file '//fname(1))
        ierr = nf90_get_var( fid, varid(1), lon )
        call nf90_check(ierr, 'Error while getting '//varname(1)//' in file '//fname(1))
        ierr = nf90_inq_varid( fid, varname(2), varid(2) )
        call nf90_check(ierr, 'Error while inquiring ID of variable '//varname(2)//' in file '//fname(1))
        ierr = nf90_get_var( fid, varid(2), lat )
        call nf90_check(ierr, 'Error while getting '//varname(2)//' in file '//fname(1))
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//fname(1))
        

        ! Load main inputs (geographic and climatic)
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        ! temperature:
        ierr = nf90_open(fname(3), NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//fname(3))
        ierr = nf90_inq_varid( fid, varname(3), varid(1) )
        call nf90_check(ierr, 'Error while inquiring ID of variable '//varname(3)//' in file '//fname(3))
        ierr = nf90_get_att( fid, varid(1), missvalname(3), fillval(3) )
        call nf90_check(ierr, 'Error while getting attribute '//missvalname(3)//' of variable '//varname(3)//' in file '//fname(3))
        ierr = nf90_get_var( fid, varid(1), temperature )
        call nf90_check(ierr, 'Error while getting variable '//varname(3)//' in file '//fname(3))
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//fname(3))

        ! runoff:
        ierr = nf90_open(fname(4), NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//fname(4))
        ierr = nf90_inq_varid( fid, varname(4), varid(1) )
        call nf90_check(ierr, 'Error while inquiring ID of variable '//varname(4)//' in file '//fname(4))
        ierr = nf90_get_att( fid, varid(1), missvalname(4), fillval(4) )
        call nf90_check(ierr, 'Error while getting attribute '//missvalname(4)//' of variable '//varname(4)//' in file '//fname(4))
        ierr = nf90_get_var( fid, varid(1), runoff )
        call nf90_check(ierr, 'Error while getting variable '//varname(4)//' in file '//fname(4))
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//fname(4))

        ! slope:
        ierr = nf90_open(fname(5), NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//fname(5))
        ierr = nf90_inq_varid( fid, varname(5), varid(1) )
        call nf90_check(ierr, 'Error while inquiring ID of variable '//varname(5)//' in file '//fname(5))
        ierr = nf90_get_att( fid, varid(1), missvalname(5), fillval(5) )
        call nf90_check(ierr, 'Error while getting attribute '//missvalname(5)//' of variable '//varname(5)//' in file '//fname(5))
        ierr = nf90_get_var( fid, varid(1), slope )
        call nf90_check(ierr, 'Error while getting variable '//varname(5)//' in file '//fname(5))
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//fname(5))

        ! cell area:
        ierr = nf90_open(fname(6), NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//fname(6))
        ierr = nf90_inq_varid( fid, varname(6), varid(1) )
        call nf90_check(ierr, 'Error while inquiring ID of variable '//varname(6)//' in file '//fname(6))
        !ierr = nf90_get_att( fid, varid(1), missvalname(6), fillval(6) )
        !call nf90_check(ierr, 'Error while getting attribute '//missvalname(6)//' of variable '//varname(6)//' in file '//fname(6))
        ierr = nf90_get_var( fid, varid(1), cell_area )
        call nf90_check(ierr, 'Error while getting variable '//varname(6)//' in file '//fname(6))
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//fname(6))

        ! lithology:
        ierr = nf90_open(fname(7), NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//fname(7))
        ierr = nf90_inq_varid( fid, varname(7), varid(1) )
        call nf90_check(ierr, 'Error while inquiring ID of variable '//varname(7)//' in file '//fname(7))
        !ierr = nf90_get_att( fid, varid(1), missvalname(7), fillval(7) )
        !call nf90_check(ierr, 'Error while getting attribute '//missvalname(7)//' of variable '//varname(7)//' in file '//fname(7))
        ierr = nf90_get_var( fid, varid(1), litho_frac, start=(/1,1,2/), count=(/nlon,nlat,nlith/) )
        call nf90_check(ierr, 'Error while getting variable '//varname(7)//' in file '//fname(7))
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//fname(7))

        ! land area:
        ierr = nf90_open(fname(8), NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//fname(8))
        ierr = nf90_inq_varid( fid, varname(8), varid(1) )
        call nf90_check(ierr, 'Error while inquiring ID of variable '//varname(8)//' in file '//fname(8))
        !ierr = nf90_get_att( fid, varid(1), missvalname(8), fillval(8) )
        !call nf90_check(ierr, 'Error while getting attribute '//missvalname(8)//' of variable '//varname(8)//' in file '//fname(8))
        ierr = nf90_get_var( fid, varid(1), land_area )
        call nf90_check(ierr, 'Error while getting variable '//varname(8)//' in file '//fname(8))
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//fname(8))
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        ! missingpoints and continental pixels:
        k = 0
        do j = 1,nlat
          do i = 1,nlon
            if ( temperature(i,j) == fillval(3) &
                .or. runoff(i,j)  == fillval(4) &
                .or. slope(i,j)   == fillval(5) &
                .or. land_area(i,j) == 0        &
                .or. all(litho_frac(i,j,:)==0)   ) then
              missingpoints(i,j) = .true.
            else
              k = k+1
              list_cont_i(k) = i
              list_cont_j(k) = j
            end if
          end do
        end do
        ncontpxl = k

        ! Null-slope security:
        where (slope==0) slope = minval( slope, (slope>0) )





        !==================================================================================================!
        !====================================   CREATE NETCDF OUTPUT   ====================================!
        !==================================================================================================!
        

        if (extend_output_file) then ! open existing output file and get data ID



                ierr = nf90_open( output_file_name, NF90_WRITE, fid )
                call nf90_check( ierr, 'Error while openning existing output file' )

                ierr = nf90_inq_dimid( fid, 'lon', dimid(1) )
                call nf90_check(ierr, 'Error while inquiring ID of dimension "lon" in output file')
                ierr = nf90_inquire_dimension( fid, dimid(1), len=nlon0 )
                call nf90_check(ierr, 'Error while inquiring length of dimension "lon" in output file')
                ierr = nf90_inq_dimid( fid, 'lat', dimid(2) )
                call nf90_check(ierr, 'Error while inquiring ID of dimension "lat" in output file')
                ierr = nf90_inquire_dimension( fid, dimid(2), len=nlat0 )
                call nf90_check(ierr, 'Error while inquiring length of dimension "lat" in output file')

                if ( nlon/=nlon0 .or. nlat/=nlat0 ) then
                        print *, 'Error: inconsistent output file dimension'
                        stop
                end if

                ierr = nf90_inq_dimid( fid, 'parameterization', dimid(4) )
                call nf90_check(ierr, 'Error while inquiring ID of dimension "parameterization" in output file')
                ierr = nf90_inquire_dimension( fid, dimid(4), len=n )
                call nf90_check(ierr, 'Error while inquiring length of dimension "parameterization" in output file')
                print *, n
                stop

                ierr = nf90_inq_varid( fid, 'parameterization', varid(4) )
                call nf90_check( ierr, 'Error while inquiring variable "parameterization" ID in output file' )
                !
                ierr = nf90_inq_varid( fid, 'weathering', varid(9) )
                call nf90_check( ierr, 'Error while inquiring variable "weathering" ID in output file' )
                !
                ierr = nf90_inq_varid( fid, 'Ea', varid(12) )
                call nf90_check( ierr, 'Error while inquiring variable "Ea" ID in output file' )
                !
                ierr = nf90_inq_varid( fid, 'T0', varid(13) )
                call nf90_check( ierr, 'Error while inquiring variable "T0" ID in output file' )
                !
                ierr = nf90_inq_varid( fid, 'Kwest', varid(19) )
                call nf90_check( ierr, 'Error while inquiring variable "Kwest" ID in output file' )
                !
                ierr = nf90_inq_varid( fid, 'kw', varid(20) )
                call nf90_check( ierr, 'Error while inquiring variable "kw" ID in output file' )
                !
                ierr = nf90_inq_varid( fid, 'sigma', varid(21) )
                call nf90_check( ierr, 'Error while inquiring variable "sigma" ID in output file' )
                !
                ierr = nf90_inq_varid( fid, 'CaMg', varid(22) )
                call nf90_check( ierr, 'Error while inquiring variable "CaMg" ID in output file' )
                !
                ierr = nf90_inq_varid( fid, 'WSI', varid(23) )
                call nf90_check( ierr, 'Error while inquiring variable "WSI" ID in output file' )




        else ! create a new output file (erase old one if it exists)



                ierr = nf90_create( output_file_name, NF90_CLOBBER, fid )
                call nf90_check( ierr, 'Error while creating output file' )
                
                ierr = nf90_def_dim( fid, 'lon', nlon, dimid(1) )
                call nf90_check( ierr, 'Error while defining dimension "lon" in output file' )
                !
                ierr = nf90_def_dim( fid, 'lat', nlat, dimid(2) )
                call nf90_check( ierr, 'Error while defining dimension "lat" in output file' )
                !
                ierr = nf90_def_dim( fid, 'lithology', 6, dimid(3) )
                call nf90_check( ierr, 'Error while defining dimension "lithology" in output file' )
                !
                ierr = nf90_def_dim( fid, 'parameterization', NF90_UNLIMITED, dimid(4) )
                call nf90_check( ierr, 'Error while defining dimension "parameterization" in output file' )

                ierr = nf90_def_var( fid, 'lon', NF90_FLOAT, dimid(1), varid(1) )
                call nf90_check( ierr, 'Error while defining variable "lon" in output file' )
                ierr = nf90_def_var( fid, 'lat', NF90_FLOAT, dimid(2), varid(2) )
                call nf90_check( ierr, 'Error while defining variable "lat" in output file' )
                ierr = nf90_def_var( fid, 'lithology', NF90_INT, dimid(3), varid(3) )
                call nf90_check( ierr, 'Error while defining variable "lithology" in output file' )
                ierr = nf90_def_var( fid, 'parameterization', NF90_INT, dimid(4), varid(4) )
                call nf90_check( ierr, 'Error while defining variable "parameterization" in output file' )
                !
                ierr = nf90_put_att( fid, varid(1), 'axis', 'X' )
                call nf90_check( ierr, 'Error while putting attribute "axis" "X" in output file variable "lon"' )
                ierr = nf90_put_att( fid, varid(2), 'axis', 'Y' )
                call nf90_check( ierr, 'Error while putting attribute "axis" "Y" in output file variable "lat"' )
                !
                ierr = nf90_def_var( fid, 'cell_area', NF90_FLOAT, dimid(1:2), varid(5) )
                call nf90_check( ierr, 'Error while defining variable "cell_area" in output file' )
                ierr = nf90_put_att( fid, varid(5), 'units', 'm2' )
                call nf90_check( ierr, 'Error while putting attribute "units" "m2" in output file variable "cell_area"' )
                ierr = nf90_put_att( fid, varid(5), '_FillValue', real(fillval(6)) )
                call nf90_check( ierr, 'Error while putting attribute "_FillValue" in output file variable "cell_area"' )
                !
                ierr = nf90_def_var( fid, 'reg_thickness', NF90_FLOAT, dimid(1:2), varid(7) )
                call nf90_check( ierr, 'Error while defining variable "reg_thickness" in output file' )
                ierr = nf90_put_att( fid, varid(7), 'units', 'm' )
                call nf90_check( ierr, 'Error while putting attribute "units" "m" in output file variable "reg_thickness"' )
                ierr = nf90_put_att( fid, varid(7), '_FillValue', real(deffillval) )
                call nf90_check( ierr, 'Error while putting attribute "_FillValue" in output file variable "reg_thickness"' )
                !
                ierr = nf90_def_var( fid, 'erosion', NF90_FLOAT, dimid(1:2), varid(8) )
                call nf90_check( ierr, 'Error while defining variable "erosion" in output file' )
                ierr = nf90_put_att( fid, varid(8), 'units', 'm/yr' )
                call nf90_check( ierr, 'Error while putting attribute "units" "m/yr" in output file variable "erosion"' )
                ierr = nf90_put_att( fid, varid(8), '_FillValue', real(deffillval) )
                call nf90_check( ierr, 'Error while putting attribute "_FillValue" in output file variable "erosion"' )
                !
                ierr = nf90_def_var( fid, 'weathering', NF90_FLOAT, (/dimid(1),dimid(2),dimid(4)/), varid(9) )
                call nf90_check( ierr, 'Error while defining variable "weathering" in output file' )
                ierr = nf90_put_att( fid, varid(9), 'units', 'mol(CaMg)/m2/yr' )
                call nf90_check( ierr, 'Error while putting attribute "units" "mol(CaMg)/m2/yr" in output file variable "weathering"' )
                ierr = nf90_put_att( fid, varid(9), '_FillValue', real(deffillval) )
                call nf90_check( ierr, 'Error while putting attribute "_FillValue" in output file variable "weathering"' )
                !
                ierr = nf90_def_var( fid, 'Ea_rp', NF90_FLOAT, dimid(3:2), varid(10) )
                call nf90_check( ierr, 'Error while defining variable "Ea_rp" in output file' )
                ierr = nf90_put_att( fid, varid(10), 'long_name', 'Activation energy for regolith production' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "Ea_rp"' )
                ierr = nf90_put_att( fid, varid(10), 'units', 'J/mol' )
                call nf90_check( ierr, 'Error while putting attribute "units" "J/mol" in output file variable "Ea_rp"' )
                !
                ierr = nf90_def_var( fid, 'T0_rp', NF90_FLOAT, dimid(3:2), varid(11) )
                call nf90_check( ierr, 'Error while defining variable "T0_rp" in output file' )
                ierr = nf90_put_att( fid, varid(11), 'long_name', 'Reference temperature of Arrhemius law (reg. prod.)' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "T0_rp"' )
                ierr = nf90_put_att( fid, varid(11), 'units', 'K' )
                call nf90_check( ierr, 'Error while putting attribute "units" "K" in output file variable "T0_rp"' )
                !
                ierr = nf90_def_var( fid, 'Ea', NF90_FLOAT, dimid(4), varid(12) )
                call nf90_check( ierr, 'Error while defining variable "Ea" in output file' )
                ierr = nf90_put_att( fid, varid(12), 'long_name', 'Activation energy for primary phases dissolution' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "Ea"' )
                ierr = nf90_put_att( fid, varid(12), 'units', 'J/mol' )
                call nf90_check( ierr, 'Error while putting attribute "units" "J/mol" in output file variable "Ea"' )
                !
                ierr = nf90_def_var( fid, 'T0', NF90_FLOAT, dimid(4), varid(13) )
                call nf90_check( ierr, 'Error while defining variable "T0" in output file' )
                ierr = nf90_put_att( fid, varid(13), 'long_name', 'Reference temperature of Arrhemius law (prim. ph. dissolution)' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "T0"' )
                ierr = nf90_put_att( fid, varid(13), 'units', 'K' )
                call nf90_check( ierr, 'Error while putting attribute "units" "K" in output file variable "T0"' )
                !
                ierr = nf90_def_var( fid, 'krp', NF90_FLOAT, dimid(3:2), varid(14) )
                call nf90_check( ierr, 'Error while defining variable "krp" in output file' )
                ierr = nf90_put_att( fid, varid(14), 'long_name', 'Main constant for regolith production' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "krp"' )
                ierr = nf90_put_att( fid, varid(14), 'units', '-' )
                call nf90_check( ierr, 'Error while putting attribute "units" "-" in output file variable "krp"' )
                !
                ierr = nf90_def_var( fid, 'h0', NF90_FLOAT, dimid(3:2), varid(15) )
                call nf90_check( ierr, 'Error while defining variable "h0" in output file' )
                ierr = nf90_put_att( fid, varid(15), 'long_name', 'Characteristic depth of Soil Production Function' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "h0"' )
                ierr = nf90_put_att( fid, varid(15), 'units', 'm' )
                call nf90_check( ierr, 'Error while putting attribute "units" "m" in output file variable "h0"' )
                !
                ierr = nf90_def_var( fid, 'ke', NF90_FLOAT, dimid(3:2), varid(16) )
                call nf90_check( ierr, 'Error while defining variable "ke" in output file' )
                ierr = nf90_put_att( fid, varid(16), 'long_name', 'Erodibility constant' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "ke"' )
                ierr = nf90_put_att( fid, varid(16), 'units', '-' )
                call nf90_check( ierr, 'Error while putting attribute "units" "-" in output file variable "ke"' )
                !
                ierr = nf90_def_var( fid, 'a', NF90_FLOAT, dimid(3:2), varid(17) )
                call nf90_check( ierr, 'Error while defining variable "a" in output file' )
                ierr = nf90_put_att( fid, varid(17), 'long_name', 'Runoff exponant for erosion' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "a"' )
                ierr = nf90_put_att( fid, varid(17), 'units', '-' )
                call nf90_check( ierr, 'Error while putting attribute "units" "-" in output file variable "a"' )
                !
                ierr = nf90_def_var( fid, 'b', NF90_FLOAT, dimid(3:2), varid(18) )
                call nf90_check( ierr, 'Error while defining variable "b" in output file' )
                ierr = nf90_put_att( fid, varid(18), 'long_name', 'Slope exponant for erosion' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "b"' )
                ierr = nf90_put_att( fid, varid(18), 'units', '-' )
                call nf90_check( ierr, 'Error while putting attribute "units" "-" in output file variable "b"' )
                !
                ierr = nf90_def_var( fid, 'Kwest', NF90_FLOAT, dimid(4), varid(19) )
                call nf90_check( ierr, 'Error while defining variable "Kwest" in output file' )
                ierr = nf90_put_att( fid, varid(19), 'long_name', 'Main constant for primary phases dissolution' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "Kwest"' )
                ierr = nf90_put_att( fid, varid(19), 'units', 'yr-1' )
                call nf90_check( ierr, 'Error while putting attribute "units" "yr-1" in output file variable "Kwest"' )
                !
                ierr = nf90_def_var( fid, 'kw', NF90_FLOAT, dimid(4), varid(20) )
                call nf90_check( ierr, 'Error while defining variable "kw" in output file' )
                ierr = nf90_put_att( fid, varid(20), 'long_name', 'Runoff sensitivity constant for primary phases dissolution' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "kw"' )
                ierr = nf90_put_att( fid, varid(20), 'units', 'yr/m' )
                call nf90_check( ierr, 'Error while putting attribute "units" "yr/m" in output file variable "kw"' )
                !
                ierr = nf90_def_var( fid, 'sigma', NF90_FLOAT, dimid(4), varid(21) )
                call nf90_check( ierr, 'Error while defining variable "sigma" in output file' )
                ierr = nf90_put_att( fid, varid(21), 'long_name', 'Age exponant for primary phases dissolution' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "sigma"' )
                ierr = nf90_put_att( fid, varid(21), 'units', '-' )
                call nf90_check( ierr, 'Error while putting attribute "units" "-" in output file variable "sigma"' )
                !
                ierr = nf90_def_var( fid, 'CaMg', NF90_FLOAT, dimid(3:4), varid(22) )
                call nf90_check( ierr, 'Error while defining variable "CaMg" in output file' )
                ierr = nf90_put_att( fid, varid(22), 'long_name', 'Molar abundancy of Ca and Mg in primary phases' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "CaMg"' )
                ierr = nf90_put_att( fid, varid(22), 'units', 'mol/m3' )
                call nf90_check( ierr, 'Error while putting attribute "units" "mol/m3" in output file variable "CaMg"' )
                !
                ierr = nf90_def_var( fid, 'WSI', NF90_FLOAT, dimid(4), varid(23) )
                call nf90_check( ierr, 'Error while defining variable "WSI" in output file' )
                ierr = nf90_put_att( fid, varid(23), 'long_name', 'Weathering Saturation Index' )
                call nf90_check( ierr, 'Error while putting attribute "long name" in output file variable "WSI"' )
                ierr = nf90_put_att( fid, varid(23), 'units', '-' )
                call nf90_check( ierr, 'Error while putting attribute "units" "mol/m3" in output file variable "WSI"' )
                !
                ierr = nf90_def_var( fid, 'land_area', NF90_FLOAT, dimid(1:2), varid(24) )
                call nf90_check( ierr, 'Error while defining variable "land_area" in output file' )
                ierr = nf90_put_att( fid, varid(24), 'units', 'm2' )
                call nf90_check( ierr, 'Error while putting attribute "units" "m2" in output file variable "land_area"' )
                ierr = nf90_put_att( fid, varid(24), '_FillValue', real(fillval(8)) )
                call nf90_check( ierr, 'Error while putting attribute "_FillValue" in output file variable "land_area"' )


                ierr = nf90_enddef( fid )
                call nf90_check( ierr, 'Error while end of definition of output file' )


                ierr = nf90_put_var( fid, varid(1), real(lon) )
                call nf90_check( ierr, 'Error while putting variable "lon" in output file' )
                !
                ierr = nf90_put_var( fid, varid(2), real(lat) )
                call nf90_check( ierr, 'Error while putting variable "lat" in output file' )
                !
                ierr = nf90_put_var( fid, varid(3), (/ (i,i=1,6) /)  )
                call nf90_check( ierr, 'Error while putting variable "lithology" in output file' )


                ierr = nf90_put_var( fid, varid(10), real(Ea_rp) )
                call nf90_check( ierr, 'Error while putting variable "Ea_rp" in output file' )
                !
                ierr = nf90_put_var( fid, varid(11), real(T0_rp) )
                call nf90_check( ierr, 'Error while putting variable "T0_rp" in output file' )
                !
                ierr = nf90_put_var( fid, varid(14), real(krp) )
                call nf90_check( ierr, 'Error while putting variable "krp" in output file' )
                !
                ierr = nf90_put_var( fid, varid(15), real(h0) )
                call nf90_check( ierr, 'Error while putting variable "h0" in output file' )
                !
                ierr = nf90_put_var( fid, varid(16), real(ke) )
                call nf90_check( ierr, 'Error while putting variable "ke" in output file' )
                !
                ierr = nf90_put_var( fid, varid(17), real(a) )
                call nf90_check( ierr, 'Error while putting variable "a" in output file' )
                !
                ierr = nf90_put_var( fid, varid(18), real(b) )
                call nf90_check( ierr, 'Error while putting variable "b" in output file' )



        end if





        !==================================================================================================!
        !==================================    LOOP TO TEST PARAMETERS   ==================================!
        !==================================================================================================!


        ! redefine litho to get rid of the oceanic part (it is define as the fraction of the CELL covered echo litho class, and we
        ! want the fraction of the LAND part of the cell covered by each litho class
        ! + reshape litho matrix
        do k = 1,ncontpxl

          i = list_cont_i(k)
          j = list_cont_j(k)

          litho_frac(i,j,:) = litho_frac(i,j,:)/sum(litho_frac(i,j,:))
          litho_frac2(:,i,j) = litho_frac(i,j,:)

        end do


        ! Put fillvalues:
        h = deffillval
        E = deffillval


        ! Erosion and regolith thickness
        Esum = 0
        do k = 1,ncontpxl
          i = list_cont_i(k)
          j = list_cont_j(k)

          ! optimal regolith production rate:
          RPopt = reg_prod_opt(temperature(i,j),runoff(i,j))
          ! regolith erosion rate:
          E(i,j) = erosion(temperature(i,j),runoff(i,j),slope(i,j))
          ! regolith thickness
          h(i,j) = eq_reg_thick(RPopt,E(i,j))

          Esum = Esum + land_area(i,j)*E(i,j)

        end do

        print *, 'Total erosion (1e9 m3/yr):'
        print *, 1e-9 * Esum



        if (.not. extend_output_file) then


          ! Record litho, h, E and areas

          ierr = nf90_put_var( fid, varid(7), real(h) )
          call nf90_check( ierr, 'Error while putting variable "h" in output file' ) 

          ierr = nf90_put_var( fid, varid(8), real(E) )
          call nf90_check( ierr, 'Error while putting variable "E" in output file' ) 

          ierr = nf90_put_var( fid, varid(5), real(cell_area) )
          call nf90_check( ierr, 'Error while putting variable "cell_area" in output file' ) 

          ierr = nf90_put_var( fid, varid(24), real(land_area) )
          call nf90_check( ierr, 'Error while putting variable "land_area" in output file' ) 


          n = 0


        end if


        !*******************************************************************************************************!
        !*******************************************************************************************************!


        ! Flag "TEST_SINGLE_PARAMETERIZATION" activated if you want to test the code for only 1 parameterization
        if (TEST_SINGLE_PARAMETERIZATION) then
          NKwest = 1
          Nkw = 1
          Nsigma = 1
          NCaMgBulk = 1
          NCaMgSed = 1
        end if


        ! Loops on Dynsoil parameters:
        print *, '0% done'

        do i1 = 1,NKwest
          Kwest = list_Kwest(i1)

          do i2 = 1,Nkw
            kw = list_kw(i2)

            do i3 = 1,Nsigma
              sigma = list_sigma(i3)

              ! Geographical loop and dynsoil model run
              WSI = 0
              do k = 1,ncontpxl
                i = list_cont_i(k)
                j = list_cont_j(k)
                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                call dynsoil_steady_state( temperature(i,j), runoff(i,j), slope(i,j), h(i,j), E(i,j), W(i,j) )
                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                WSI = WSI + W(i,j)*land_area(i,j)
              end do
              WSI = WSI / Esum

              ! Loops on Ca-Mg abundancy:
              do i4 = 1,NCaMgBulk
                CaMg_rock(1) = list_CaMgBulk(i4)

                do i5 = 1,NCaMgSed
                  CaMg_rock(6) = list_CaMgSed(i5)

                  if (CaMg_rock(6) < CaMg_rock(1)) then

                    n = n+1

                    ! loop on continental points
                    do k = 1,ncontpxl

                      i = list_cont_i(k)
                      j = list_cont_j(k)

                      ! CaMg weathering
                      WCaMg = 0.
                      do l = 1,nlith
                        WCaMg = WCaMg + litho_frac2(l,i,j)*CaMg_rock(l)*W(i,j)
                      end do

                      ! Storage:
                      !---------

                      ! Weathering:
                      ierr = nf90_put_var( fid, varid(9),   (/real(WCaMg)/),          start=(/i,j,n/), count=(/1,1,1/) )

                    end do

                    ! Parameters:
                    ierr = nf90_put_var( fid, varid(4),   (/n/),                    start=(/n/), count=(/1/) )
                    ierr = nf90_put_var( fid, varid(12),  (/real(Ea)/),             start=(/n/), count=(/1/) )
                    ierr = nf90_put_var( fid, varid(13),  (/real(T0)/),             start=(/n/), count=(/1/) )
                    ierr = nf90_put_var( fid, varid(19),  (/real(Kwest)/),          start=(/n/), count=(/1/) )
                    ierr = nf90_put_var( fid, varid(20),  (/real(kw)/),             start=(/n/), count=(/1/) )
                    ierr = nf90_put_var( fid, varid(21),  (/real(sigma)/),          start=(/n/), count=(/1/) )
                    ierr = nf90_put_var( fid, varid(22),  real(CaMg_rock),          start=(/1,n/), count=(/nlith,1/) )

                    ! WSI
                    ierr = nf90_put_var( fid, varid(23),  (/real(WSI)/),            start=(/n/), count=(/1/) )
  
                  end if

                end do
              end do
            end do
            print '(I3,A6)',  ( 100 * ((i1-1)*Nkw + i2)*Nsigma*NCaMgBulk*NCaMgSed )  &
                              /  (NKwest*Nkw*Nsigma*NCaMgBulk*NCaMgSed)              &
                           ,'% done'
          end do
        end do

        !*******************************************************************************************************!
        !*******************************************************************************************************!


       
        ! Close output file
        ierr = nf90_close( fid )
        call nf90_check( ierr, 'Error while closing output file' )




end program
