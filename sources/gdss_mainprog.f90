program gdss_mainprog

        use io_module, only: outinfo, make_input_output
        use netcdf_io_functions, only: close_file, put_var_real0D, put_var_real1D, put_var_real2D, put_var_real3D
        use dynsoil_steady_state_module, only: dynsoil_geographic_loop, CaMg_weathering, litho_average
        use climate_module, only: get_climate, find_weathering_interval, interpolate_climate, inverse_interpolation

        implicit none

        ! name of input-output interface file
        character(len=*), parameter:: IO_FNAME = '../IO_INTERFACE'

        ! dimension size
        integer:: nlon, nlat, nlith, nrun, ncont

        ! main variables
        double precision:: CO2, Fvol
        double precision, dimension(:), allocatable:: CO2_levels, forcing
        double precision, dimension(:,:), allocatable:: cell_area, land_area, curr_temp, curr_runf, slope, W_all
        double precision, dimension(:,:,:), allocatable:: temperature, runoff, lith_frac, reshp_lith_frac, h, E, xs, W
        logical, dimension(:,:), allocatable:: cont_points

        ! other variables
        integer, dimension(:), allocatable:: list_cont_i, list_cont_j
        integer:: ForwBckw
        integer, dimension(8):: computer_time

        ! parameters
        double precision, dimension(:,:,:), allocatable:: params ! nparams x nlith x nrun

        ! output netCDF file
        type(outinfo):: output_info

        ! indices
        integer:: i, j, k, n, nstep




        print *
        print *
        print *
        print *, '<><><><><><><><><><><><><><><><><>'
        print *, 'GEOCLIM-DynSoil-steady-state model'
        print *, '<><><><><><><><><><><><><><><><><>'

        !====================================================================!
        ! INPUT/OUTPUT: read interface file, load inputs, create output file !
        !====================================================================!

        call make_input_output( IO_FNAME, &
                                cell_area, land_area, temperature, runoff, slope, lith_frac, CO2_levels, &
                                curr_temp, curr_runf, h, E, xs, W, W_all, &
                                nlon, nlat, nlith, nrun, &
                                ForwBckw, forcing, params, output_info )




        print *
        print *
        print *, 'Pre-computations...'

        !========================!
        ! Set continental points !
        !========================!

        allocate(cont_points(nlon,nlat))
        allocate(list_cont_i(nlon*nlat))
        allocate(list_cont_j(nlon*nlat))

        k = 0
        do j = 1,nlat
          do i = 1,nlon
            if (land_area(i,j)>0) then
              k = k+1
              list_cont_i(k) = i
              list_cont_j(k) = j
              cont_points(i,j) = .true.
            else
              cont_points(i,j) = .false.
            end if
          end do
        end do

        ! number of continental points:
        ncont = k




        ! Reshape lithology fraction variable to increase computation efficiency (put lithology as first dimension)
        ! ---------------------------------------------------------------------------------------------------------

        allocate( reshp_lith_frac(nlith, nlon, nlat) )
        reshp_lith_frac = reshape(lith_frac, (/nlith, nlon, nlat/), order=(/2,3,1/))




        !=========================================================================================================================!
        !                                                      Run the model                                                      !
        !=========================================================================================================================!



        print *
        print *
        print *, 'Start run'
        call date_and_time(values=computer_time)
        write(*, fmt='(I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') &
            computer_time(1),'/',computer_time(2),'/',computer_time(3),' ', &
            computer_time(5),':',computer_time(6),':',computer_time(7),':',computer_time(8)
        print *


        !-------------------------------------------------------------------------------------------------!
        ! Forward run, model forced by atmospheric CO2 level, compute equilibrium volcanic degassing Fvol !
        !-------------------------------------------------------------------------------------------------!
        if (ForwBckw==1) then



          ! parameterization loop
          do n = 1,nrun

            ! Forcing
            CO2 = forcing(n)

            ! Computation
            !*********************************************************************************************************!
            call get_climate( CO2, CO2_levels, temperature, runoff, list_cont_i(1:ncont), list_cont_j(1:ncont), &
                              curr_temp, curr_runf )

            call dynsoil_geographic_loop( list_cont_i(1:ncont), list_cont_j(1:ncont),                         &
                                          land_area, curr_temp, curr_runf, slope, reshp_lith_frac, &
                                          h, E, xs, W, Fvol, params(:,:,n)                                    )
            !*********************************************************************************************************!


            print *, '#######################################'
            print *, n, '/', nrun
            print *, 'Sil. wth.:', Fvol


            ! write run-dependent output
            !---------------------------

            k = 1 ! AREA:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(land_area),  begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 2 ! LITHOLOGY FRACTION:
            if (output_info%variables(k)%unlimited) &
                call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(lith_frac),  begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))

            k = 3 ! ATMOSPHERIC CO2:
            if (output_info%variables(k)%unlimited) &
                call put_var_real0D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(CO2),        begin=(/n/), length=(/1/))

            k = 4 ! VOLCANIC DEGASSING:
            if (output_info%variables(k)%unlimited) &
                call put_var_real0D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(Fvol),       begin=(/n/), length=(/1/))

            k = 5 ! TEMPERATURE:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(curr_temp),  begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 6 ! RUNOFF:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(curr_runf),  begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 7 ! SLOPE:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(slope),      begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 8 ! EROSION:
            if (output_info%variables(k)%unlimited) &
                call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(reshape(E, (/nlon, nlat, nlith/), order=(/3,1,2/))), &
            !                            reshape: put lithology as last dim: ------^
                                    begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))

            k = 9 ! REGOLITH THICKNESS:
            if (output_info%variables(k)%unlimited) &
                call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(reshape(h, (/nlon, nlat, nlith/), order=(/3,1,2/))), &
            !                            reshape: put lithology as last dim: ------^
                                    begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))

            k = 10 ! WEATHERING IN M/YR:
            if (output_info%variables(k)%unlimited) &
                call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(reshape(W, (/nlon, nlat, nlith/), order=(/3,1,2/))), &
            !                            reshape: put lithology as last dim: ------^
                                    begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))

            k = 11 ! WEATHERING OF Ca-Mg (mol/m2/y):
            if (output_info%variables(k)%unlimited) then
              call CaMg_weathering(list_cont_i(1:ncont), list_cont_j(1:ncont), W, params(:,:,n))
              call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                  real(reshape(W, (/nlon, nlat, nlith/), order=(/3,1,2/))), &
            !                          reshape: put lithology as last dim: ------^
                                  begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))
            end if

            k = 12 ! PRIMARY PHASES PROPORTION AT SURFACES:
            if (output_info%variables(k)%unlimited) &
                call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(reshape(xs, (/nlon, nlat, nlith/), order=(/3,1,2/))), &
            !                            reshape: put lithology as last dim: -------^
                                    begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))

            k = 13 ! WEATHERING OF Ca-Mg (mol/m2/y) -- ALL LITHOLOGY:
            if (output_info%variables(k)%unlimited) then
              call litho_average(list_cont_i(1:ncont), list_cont_j(1:ncont), reshp_lith_frac, W, W_all)
              call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                  real(W_all), begin=(/1,1,n/), length=(/nlon,nlat,1/))
            end if


          end do



        !---------------------------------------------------------------------------------------------------------!
        ! Backward run, model forced by volcanic CO2 degassing, inverse for the equilibrium atmospheric CO2 level !
        !---------------------------------------------------------------------------------------------------------!
        elseif (ForwBckw==-1) then



          ! parameterization loop
          do n = 1,nrun


            print *, '#######################################'
            print *, n, '/', nrun
            write(unit=*,fmt='(A26)',advance='no') 'number of iteration steps:'


            ! Forcing
            Fvol = forcing(n)

            ! Inverse computation
            !**********************************************************************************************!
            call inverse_interpolation(Fvol, list_cont_i(1:ncont), list_cont_j(1:ncont), params(:,:,n),    &
                                       land_area, temperature, runoff, slope, reshp_lith_frac, CO2_levels, &
                                       CO2, curr_temp, curr_runf, h, E, xs, W, nstep                       )
            !**********************************************************************************************!

            print *, nstep
            print *, 'CO2 level:', CO2


            ! write run-dependent output
            !---------------------------

            k = 1 ! AREA:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(land_area),  begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 2 ! LITHOLOGY FRACTION:
            if (output_info%variables(k)%unlimited) &
                call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(lith_frac),  begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))

            k = 3 ! ATMOSPHERIC CO2:
            if (output_info%variables(k)%unlimited) &
                call put_var_real0D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(CO2),        begin=(/n/), length=(/1/))

            k = 4 ! VOLCANIC DEGASSING:
            if (output_info%variables(k)%unlimited) &
                call put_var_real0D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(Fvol),       begin=(/n/), length=(/1/))

            k = 5 ! TEMPERATURE:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(curr_temp),  begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 6 ! RUNOFF:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(curr_runf),  begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 7 ! SLOPE:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(slope),      begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 8 ! EROSION:
            if (output_info%variables(k)%unlimited) &
                call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(reshape(E, (/nlon, nlat, nlith/), order=(/3,1,2/))), &
            !                            reshape: put lithology as last dim: ------^
                                    begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))

            k = 9 ! REGOLITH THICKNESS:
            if (output_info%variables(k)%unlimited) &
                call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(reshape(h, (/nlon, nlat, nlith/), order=(/3,1,2/))), &
            !                            reshape: put lithology as last dim: ------^
                                    begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))

            k = 10 ! WEATHERING IN M/YR:
            if (output_info%variables(k)%unlimited) &
                call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(reshape(W, (/nlon, nlat, nlith/), order=(/3,1,2/))), &
            !                            reshape: put lithology as last dim: ------^
                                    begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))

            k = 11 ! WEATHERING OF Ca-Mg (mol/m2/y):
            if (output_info%variables(k)%unlimited) then
              call CaMg_weathering(list_cont_i(1:ncont), list_cont_j(1:ncont), W, params(:,:,n))
              call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                  real(reshape(W, (/nlon, nlat, nlith/), order=(/3,1,2/))), &
            !                          reshape: put lithology as last dim: ------^
                                  begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))
            end if

            k = 12 ! PRIMARY PHASES PROPORTION AT SURFACES:
            if (output_info%variables(k)%unlimited) &
                call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(reshape(xs, (/nlon, nlat, nlith/), order=(/3,1,2/))), &
            !                            reshape: put lithology as last dim: -------^
                                    begin=(/1,1,1,n/), length=(/nlon,nlat,nlith,1/))

            k = 13 ! WEATHERING OF Ca-Mg (mol/m2/y) -- ALL LITHOLOGY:
            if (output_info%variables(k)%unlimited) then
              call litho_average(list_cont_i(1:ncont), list_cont_j(1:ncont), reshp_lith_frac, W, W_all)
              call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                  real(W_all), begin=(/1,1,n/), length=(/nlon,nlat,1/))
            end if


          end do


        end if



        !------------------------------!
        ! write run-independent output !
        !------------------------------!


        k = 1 ! AREA:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                               real(land_area))

        k = 2 ! LITHOLOGY FRACTION:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                               real(lith_frac))

        k = 3 ! ATMOSPHERIC CO2:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real0D(output_info%file_id, output_info%variables(k)%varid, &
                               real(CO2))

        k = 4 ! VOLCANIC DEGASSING:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real0D(output_info%file_id, output_info%variables(k)%varid, &
                               real(Fvol))

        k = 5 ! TEMPERATURE :
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                               real(curr_temp))

        k = 6 ! RUNOFF:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                               real(curr_runf))

        k = 7 ! SLOPE:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                               real(slope))

        k = 8 ! EROSION:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                               real(reshape(E, (/nlon, nlat, nlith/), order=(/3,1,2/))))
        !                           reshape: put lithology as last dim: ------^

        k = 9 ! REGOLITH THICKNESS:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                               real(reshape(h, (/nlon, nlat, nlith/), order=(/3,1,2/))))
        !                           reshape: put lithology as last dim: ------^

        k = 10 ! WEATHERING IN M/YR:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                               real(reshape(W, (/nlon, nlat, nlith/), order=(/3,1,2/))))
        !                           reshape: put lithology as last dim: ------^

        k = 11 ! WEATHERING OF Ca-Mg (mol/m2/y):
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) then
          call CaMg_weathering(list_cont_i(1:ncont), list_cont_j(1:ncont), W, params(:,:,n))
          call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                              real(reshape(W, (/nlon, nlat, nlith/), order=(/3,1,2/))))
          !                        reshape: put lithology as last dim: ------^
        end if

        k = 12 ! PRIMARY PHASES PROPORTION AT SURFACES:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real3D(output_info%file_id, output_info%variables(k)%varid, &
                               real(reshape(xs, (/nlon, nlat, nlith/), order=(/3,1,2/))))
        !                           reshape: put lithology as last dim: -------^

        k = 13 ! WEATHERING OF Ca-Mg (mol/m2/y) -- ALL LITHOLOGY:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) then
          call litho_average(list_cont_i(1:ncont), list_cont_j(1:ncont), reshp_lith_frac, W, W_all)
          call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                              real(W_all), begin=(/1,1,n/), length=(/nlon,nlat,1/))
        end if




        print *
        print *
        print *, 'End run'
        call date_and_time(values=computer_time)
        write(*, fmt='(I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') &
            computer_time(1),'/',computer_time(2),'/',computer_time(3),' ', &
            computer_time(5),':',computer_time(6),':',computer_time(7),':',computer_time(8)
        print *




        !=============!
        ! Close files !
        !=============!

        call close_file(output_info%file_id)



end program
