program gdss_mainprog

        use io_module, only: outinfo, IPARAM, IFORC, make_input_output
        use netcdf_io_functions, only: close_file, put_var_real0D, put_var_real1D, put_var_real2D, put_var_real3D
        use dynsoil_steady_state_module, only: dynsoil_geographic_loop
        use climate_module, only: get_climate, find_weathering_interval, interpolate_climate

        implicit none

        ! name of input-output interface file
        character(len=*), parameter:: IO_FNAME = '../IO_INTERFACE'

        ! MAXIMUM ALLOWED RELATIVE IMBALANCE (|Fsilw-Fvol|/Fvol):
        ! For Backward run only
        double precision, parameter:: PREC = 1.d-6

        ! dimension size
        integer:: nlon, nlat, nlith, nrun, ncont

        ! main variables
        double precision:: CO2, Fvol, Fsilw
        double precision, dimension(:), allocatable:: CO2_levels
        double precision, dimension(:,:), allocatable:: cell_area, land_area, curr_temp, curr_runf, slope, h, E, xs, W, WCaMg
        double precision, dimension(:,:,:), allocatable:: temperature, runoff, lith_frac, reshp_lith_frac
        logical, dimension(:,:), allocatable:: cont_points

        ! other variables
        integer, dimension(:), allocatable:: list_cont_i, list_cont_j
        double precision:: lsum, Fsilw0, Fsilw1, CO2_0, CO2_1, dFsilw_dCO2, CO2_00, CO2_11, xi
        integer:: ForwBckw, k0, k1, k00, k11
        integer, dimension(8):: computer_time

        ! parameters
        include 'dynsoil_physical_parameters.inc'
        double precision, dimension(:), allocatable:: CaMg_rock

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
                                cell_area, land_area, temperature, runoff, slope, lith_frac, CaMg_rock, CO2_levels, &
                                curr_temp, curr_runf, h, E, xs, W, WCaMg, &
                                nlon, nlat, nlith, nrun, &
                                ForwBckw, output_info )




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




        !===========================================================================================================!
        ! Reshape lithology fraction variable to increase computation efficiency (put lithology as first dimension) !
        ! + redefine litho to get rid of the oceanic part (it is define as the fraction of the CELL covered by each !
        !   litho class and we want the fraction of the LAND part of the cell covered by each litho class)          !
        !===========================================================================================================!

        allocate(reshp_lith_frac(nlith,nlon,nlat))

        do j = 1,nlat
          do i =1,nlon
            lsum = sum(lith_frac(i,j,:))
            if (lsum>0) then
              !+++++++++ litho redefinition ++++++++++!
              lith_frac(i,j,:) = lith_frac(i,j,:)/lsum
              !+++++++++++++++++++++++++++++++++++++++!
            end if
            reshp_lith_frac(:,i,j) = lith_frac(i,j,:)
          end do
        end do




        !==========================================================================================================================!
        !--------------------------------------------------------------------------------------------------------------------------!
        !==========================================================================================================================!



        print *
        print *
        print *, 'Start run'
        call date_and_time(values=computer_time)
        write(*, fmt='(I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') &
            computer_time(1),'/',computer_time(2),'/',computer_time(3),' ', &
            computer_time(5),':',computer_time(6),':',computer_time(7),':',computer_time(8)
        print *




        !=========================================================================================================================!
        !                                                      Run the model                                                      !
        !=========================================================================================================================!


        !-------------------------------------------------------------------------------------------------!
        ! Forward run, model forced by atmospheric CO2 level, compute equilibrium volcanic degassing Fvol !
        !-------------------------------------------------------------------------------------------------!
        if (ForwBckw==1) then



          ! parameterization loop
          do n = 1,nrun


            ! read parameters
            read(unit=IPARAM, fmt=*)  ke, a, b, krp, Ea_rp, T0_rp, h0, kd, kw, Ea, T0, sigma, CaMg_rock

            ! read forcings
            read(unit=IFORC, fmt=*) CO2


            ! Computation
            !*********************************************************************************************************!
            call get_climate( CO2, CO2_levels, temperature, runoff, list_cont_i(1:ncont), list_cont_j(1:ncont), &
                              curr_temp, curr_runf )

            call dynsoil_geographic_loop( list_cont_i(1:ncont), list_cont_j(1:ncont),                         &
                                          CaMg_rock, land_area, curr_temp, curr_runf, slope, reshp_lith_frac, &
                                          h, E, xs, W, WCaMg, Fvol                                            )
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
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(E),          begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 9 ! REGOLITH THICKNESS:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(h),          begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 10 ! WEATHERING IN M/YR:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(W),          begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 11 ! WEATHERING:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(WCaMg),      begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 12 ! PRIMARY PHASES PROPORTION AT SURFACES:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(xs),         begin=(/1,1,n/), length=(/nlon,nlat,1/))


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


            ! read parameters
            read(unit=IPARAM, fmt=*)  ke, a, b, krp, Ea_rp, T0_rp, h0, kd, kw, Ea, T0, sigma, CaMg_rock


            ! read forcings
            read(unit=IFORC, fmt=*) Fvol

            ! find interval
            call find_weathering_interval( Fvol,                                                              &
                                           list_cont_i(1:ncont), list_cont_j(1:ncont),                        &
                                           CaMg_rock, land_area, temperature, runoff, slope, reshp_lith_frac, &
                                           h, E, xs, W, WCaMg,                                                &
                                           k0, k1, Fsilw0, Fsilw1                                             )

            CO2_00 = CO2_levels(k0)
            CO2_11 = CO2_levels(k1)
            k00 = k0
            k11 = k1

            ! initialisation
            CO2_0 = CO2_levels(k0)
            CO2_1 = CO2_levels(k1)

            Fsilw = Fsilw0
            CO2 = CO2_0

            ! Iteration to inverse equilibrium CO2
            ! Use the two levels found k0 and k1 to approximate  the derivative d(Fsilw) / d(CO2) find a new "guess" of CO2, and
            ! iterate again until the wished precision is reached

            nstep = 0

            do while ( abs((Fsilw-Fvol)/Fvol) > PREC )

              if (Fsilw>Fvol) then
                Fsilw1 = Fsilw
                CO2_1 = CO2
              else
                Fsilw0 = Fsilw
                CO2_0 = CO2
              end if

              dFsilw_dCO2 = (Fsilw1 - Fsilw0) / (CO2_1 - CO2_0)

              !++++++++++++++++++++++++++++++++++++++!
              CO2  =  CO2  +  (Fvol-Fsilw)/dFsilw_dCO2
              !++++++++++++++++++++++++++++++++++++++!

              ! Climate interpolation
              xi = (CO2-CO2_00)/(CO2_11-CO2_00)
              call interpolate_climate( list_cont_i(1:ncont), list_cont_j(1:ncont), temperature, runoff, k00, k11, xi, &
                                        curr_temp, curr_runf )

              ! Dynsoil computation at new climate => new guess of Fsilw
              call dynsoil_geographic_loop( list_cont_i(1:ncont), list_cont_j(1:ncont),                         &
                                            CaMg_rock, land_area, curr_temp, curr_runf, slope, reshp_lith_frac, &
                                            h, E, xs, W, WCaMg, Fsilw                                           )

              nstep = nstep + 1

            end do


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
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(E),          begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 9 ! REGOLITH THICKNESS:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(h),          begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 10 ! WEATHERING IN M/YR:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(W),          begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 11 ! WEATHERING:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(WCaMg),      begin=(/1,1,n/), length=(/nlon,nlat,1/))

            k = 12 ! PRIMARY PHASES PROPORTION AT SURFACES:
            if (output_info%variables(k)%unlimited) &
                call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                                    real(xs),         begin=(/1,1,n/), length=(/nlon,nlat,1/))


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
           call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                               real(E))

        k = 9 ! REGOLITH THICKNESS:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                               real(h))

        k = 10 ! WEATHERING IN M/YR:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                               real(W))

        k = 11 ! WEATHERING:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                               real(WCaMg))

        k = 12 ! PRIMARY PHASES PROPORTION AT SURFACES:
        if (output_info%variables(k)%write_var .and. (.not. output_info%variables(k)%unlimited)) &
           call put_var_real2D(output_info%file_id, output_info%variables(k)%varid, &
                               real(xs))



        ! close parameter file
        close(unit=IPARAM)




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

        if (nrun>0) then
          close(unit=IPARAM)
        end if

        close(unit=IFORC)



end program
