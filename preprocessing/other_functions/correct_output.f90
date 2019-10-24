program correct_output

        use io_module, only: nf90check
        use netcdf
        implicit none



        !================================!
        !=========  DECLARATION =========!
        !================================!

        integer:: m, n, i, j, k

        ! Geographic fields, input files and variables names

        integer, parameter:: ncorr = 5
        character(len=*), dimension(ncorr), parameter:: output_file_name = &
                (/ &
                'output/parameter_exploration_30min.nc',       &
                'output/parameter_exploration_30min_krp2.nc',  &
                'output/parameter_exploration_30min_krp3.nc',  &
                'output/parameter_exploration_30min_krp5.nc',  &
                'output/parameter_exploration_30min_krp10.nc'  &
                /)
        character(len=*), parameter:: W_varname = 'weathering'

        character(len=*), parameter:: litho_fname = '../../../calibration/present_day_30min/lith_mask.nc'
        character(len=*), parameter:: litho_varname = 'frac'
        character(len=*), parameter:: litho_dimname = 'lith'

        ! Missing-value
        double precision:: fillval

        ! Variables:
        double precision, dimension(:), allocatable:: lon, lat
        double precision, dimension(:,:), allocatable:: W, corr
        double precision, dimension(:,:,:), allocatable:: litho_frac
        integer, dimension(:), allocatable:: list_cont_i, list_cont_j
        double precision, dimension(1):: Wscal
        integer:: nlon, nlat, ncontpxl, nlith, nparam

        integer:: fid, dimid, varid, ierr









        !==================================================================================================!
        !==========================   LOAD CLIMATIC AND GEOGRAPHIC CONDITIONS   ===========================!
        !==================================================================================================!

        ierr = nf90_open(litho_fname, NF90_NOWRITE, fid)
        call nf90check(ierr, 'Error while openning file '//litho_fname)
        ierr = nf90_inq_dimid( fid, 'lon', dimid )
        call nf90check(ierr, 'Error while inquiring ID of dimension lon in file '//litho_fname)
        ierr = nf90_inquire_dimension( fid, dimid, len=nlon )
        call nf90check(ierr, 'Error while inquiring length of dimension lon in file '//litho_fname)
        ierr = nf90_inq_dimid( fid, 'lat', dimid )
        call nf90check(ierr, 'Error while inquiring ID of dimension lat in file '//litho_fname)
        ierr = nf90_inquire_dimension( fid, dimid, len=nlat )
        call nf90check(ierr, 'Error while inquiring length of dimension lat in file '//litho_fname)
        ierr = nf90_inq_dimid( fid, litho_dimname, dimid )
        call nf90check(ierr, 'Error while inquiring ID of dimension '//litho_dimname//' in file '//litho_fname)
        ierr = nf90_inquire_dimension( fid, dimid, len=nlith )
        call nf90check(ierr, 'Error while inquiring length of dimension '//litho_dimname//' in file '//litho_fname)


        ! Allocate variables
        allocate( lon(nlon) )
        allocate( lat(nlat) )
        allocate( W(nlon,nlat) )
        allocate( litho_frac(nlon,nlat,nlith) )
        allocate( corr(nlon,nlat) )
        allocate( list_cont_i(nlon*nlat) )
        allocate( list_cont_j(nlon*nlat) )

        ! load litho:
        ierr = nf90_inq_varid( fid, litho_varname, varid )
        call nf90check(ierr, 'Error while inquiring ID of variable '//litho_varname//' in file '//litho_fname)
        ierr = nf90_get_var( fid, varid, litho_frac )
        call nf90check(ierr, 'Error while getting variable '//litho_varname//' in file '//litho_fname)
        ierr = nf90_close(fid)
        call nf90check(ierr, 'Error while closing file '//litho_fname)


        do m = 1,ncorr

          print *, output_file_name(m)

          ! open output file
          ierr = nf90_open( output_file_name(m), NF90_WRITE, fid )
          call nf90check( ierr, 'Error while openning existing output file'//output_file_name(m) )

          ierr = nf90_inq_dimid( fid, 'parameterization', dimid )
          call nf90check(ierr, 'Error while inquiring ID of dimension "parameterization" in file'//output_file_name(m))
          ierr = nf90_inquire_dimension( fid, dimid, len=nparam )
          call nf90check(ierr, 'Error while getting length of dimension "parameterization" in file'//output_file_name(m))

          ierr = nf90_inq_varid( fid, W_varname, varid )
          call nf90check( ierr, 'Error while inquiring variable '//W_varname//' ID in output file'//output_file_name(m) )

          ! get W fillvalue
          ierr = nf90_get_att( fid, varid, '_FillValue', fillval )
          call nf90check( ierr, 'Error while inquiring variable '//W_varname//' attribute "_FillValue" in output file'//output_file_name(m) )

          ! get first param slice
          ierr = nf90_get_var( fid, varid, W, start=(/1,1,1/), count=(/nlon,nlat,1/) )
          call nf90check( ierr, 'Error while gettingg variable '//W_varname//' in output file'//output_file_name(m) )

          ! Test:
          ierr = nf90_get_var( fid, varid, Wscal, start=(/1,1,1/), count=(/1,1,1/) )
          call nf90check( ierr, 'Error while gettingg variable '//W_varname//' in output file'//output_file_name(m) )
          ierr = nf90_put_var( fid, varid, Wscal, start=(/1,1,1/), count=(/1,1,1/) )
          call nf90check( ierr, 'Error while putting variable '//W_varname//' in output file'//output_file_name(m) )


          ! loop to get continental points
          k = 0
          do j = 1,nlat
            do i = 1,nlon
              if (W(i,j)/=fillval) then
                if ( litho_frac(i,j,1) >= 1.d+0-1.d-13 ) then ! check if the point is not fully water/ice
                  print *, 'ERROR: fully water/ice point with a weathering value'
                  print *, 'point:', i, j
                  print *, 'weathering:', W(i,j)
                  print *, 'lithology fractions:'
                  print *, litho_frac(i,j,:)
                  stop
                end if
                k = k + 1
                list_cont_i(k) = i
                list_cont_j(k) = j
                !***************************************!
                corr(i,j) = 1.d+0 / (1-litho_frac(i,j,1)) ! 1 / continental fraction
                !***************************************!
              end if
            end do
          end do
          ncontpxl = k


          ! CHECKING: draw a map of the world in ascii file, to see if the orientation is good
          open(unit=1,file='toto',status='replace',action='write')
          do j = nlat,nlat-200,-1
            do i = 1,300
              if (W(i,j)==fillval) then
                write(unit=1,fmt='(A1)',advance='no') ' '
              else
                write(unit=1,fmt='(A1)',advance='no') '#'
              end if
            end do
            write(1,*)
          end do
          close(unit=1)


          ! loop on parameterizations and correction:
          do n = 1,nparam

            do k = 1,ncontpxl

              i = list_cont_i(k)
              j = list_cont_j(k)

              !*************************************************************************!
              ierr = nf90_get_var( fid, varid, Wscal, start=(/i,j,n/), count=(/1,1,1/) )
              Wscal(1) = Wscal(1)*corr(i,j)
              ierr = nf90_put_var( fid, varid, Wscal, start=(/i,j,n/), count=(/1,1,1/) )
              !*************************************************************************!

            end do

          end do

          ! close output file
          ierr = nf90_close(fid)
          call nf90check( ierr, 'Error while closing output file'//output_file_name(n) )



        end do




end program
