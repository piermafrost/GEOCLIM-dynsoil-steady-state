module climate_module
  implicit none



  !====================!
  ! Common parameters: !
  !====================!

  character(len=*), parameter:: INTERPOLATION_MODE = 'log' ! 'linear' or 'log'

  ! MAXIMUM ALLOWED RELATIVE IMBALANCE (|Fsilw-Fvol|/Fvol):
  ! For CO2 inversion (backward run)
  double precision, parameter:: PREC = 1.d-6



  contains

    !============================!
    ! Functions and subroutines: !
    !============================!


    function interp_coeff(x, x0, x1)
        double precision:: interp_coeff
        double precision, intent(in):: x, x0, x1
        select case (INTERPOLATION_MODE)
          case ('linear')
            interp_coeff = (x - x0) / (x1 - x0)
          case ('log')
            interp_coeff = log(x/x0) / log(x1/x0)
          case default
            print *
            print *, 'ERROR in "climate_module": illegal interpolation mode "'//INTERPOLATION_MODE//'"'
            print *, 'Expect "linear" or "log"'
            stop
        end select
    end function


  !--------------------------------------------------------------------------------------------------------------------------------!


    subroutine get_climate( CO2, CO2_levels, GMST, temperature, runoff, list_i, list_j, curr_temp, curr_runf, curr_GMST )
      double precision, intent(in):: CO2
      double precision, dimension(:), intent(in):: CO2_levels, GMST
      double precision, dimension(:,:,:), intent(in):: temperature, runoff
      integer, dimension(:), intent(in):: list_i, list_j
      double precision, dimension(:,:), intent(inout)::  curr_temp, curr_runf
      double precision, intent(out):: curr_GMST
      double precision:: xi
      integer:: k0, k1

      call find_co2_interval( CO2, CO2_levels, k0, k1, xi )

      ! Climate interpolation
      call interpolate_climate( list_i, list_j, temperature, runoff, GMST, k0, k1, xi, curr_temp, curr_runf, curr_GMST )

    end subroutine


    !------------------------------------------------------------------------------------------------------------------------------!


    subroutine interpolate_climate( list_i, list_j, temp_lvl, runoff_lvl, gmst_lvl, k0, k1, xi, temp, runoff, gmst )
      integer, dimension(:), intent(in):: list_i, list_j
      double precision, dimension(:,:,:), intent(in):: temp_lvl, runoff_lvl
      double precision, dimension(:),     intent(in):: gmst_lvl
      integer, intent(in):: k0, k1
      double precision, intent(in):: xi
      double precision, dimension(:,:), intent(inout):: temp, runoff
      double precision, intent(out):: gmst
      integer:: i, j, k, ncont

      ncont = size(list_i)

      do k = 1,ncont

        i = list_i(k)
        j = list_j(k)

        ! "linear" interpolation
        temp(i,j)   = (1-xi)*temp_lvl(i,j,k0) + xi*temp_lvl(i,j,k1)
        runoff(i,j) = (1-xi)*runoff_lvl(i,j,k0) + xi*runoff_lvl(i,j,k1)

      end do

      gmst = (1-xi)*gmst_lvl(k0) + xi*gmst_lvl(k1)

    end subroutine


    !------------------------------------------------------------------------------------------------------------------------------!


    subroutine find_co2_interval( CO2, CO2_levels, kinf, ksup, xi )
      ! Find Between which CO2 levels (CO2_inf and CO2_sup) is a given CO2 concentration and compute its relative position between
      ! these level: xi = ( CO2 - CO2_inf ) / ( CO2_sup - CO2_inf )

      double precision, intent(in)::  CO2
      double precision, dimension(:), intent(in):: CO2_levels
      integer, intent(out):: kinf, ksup
      double precision, intent(out):: xi
      integer:: k, nlev

      nlev = size(CO2_levels)

      if (nlev==1) then

        kinf = 1
        ksup = 1
        xi = 0

      else

        kinf = 1
        ksup = nlev
        ! find by dichotomy:
        do while ( ksup-kinf > 1 )
          k = kinf + (ksup-kinf)/2
          if (CO2 >= CO2_levels(k)) then
            kinf = k
          else
            ksup = k
          end if
        end do

        ! interpolation coefficient:
        xi = interp_coeff(CO2, CO2_levels(kinf), CO2_levels(ksup))

    end if

    end subroutine


    !------------------------------------------------------------------------------------------------------------------------------!


    subroutine find_weathering_interval( imposed_Fsilw, &
                    list_i, list_j, param, area, temp, runoff, slope, lith_frac, h, E, xs, W, &
                    kinf, ksup, Fsilw_inf, Fsilw_sup )
      ! Find for a given total CO2 consumpotion FSILW the two adjacent levels CO2_inf and CO2_sup such as:
      !     Fsilw(CO2_inf) <= FSILW < Fsilw(CO2_sup)

      use dynsoil_steady_state_module, only: dynsoil_geographic_loop

      double precision, intent(in):: imposed_Fsilw
      integer, dimension(:), intent(in):: list_i, list_j
      double precision, dimension(:,:), intent(in):: param ! [param x litho]
      double precision, dimension(:,:), intent(in):: area, slope
      double precision, dimension(:,:,:), intent(in):: temp, runoff, lith_frac
      double precision, dimension(:,:,:), intent(out):: h, E, xs, W
      double precision, intent(out):: Fsilw_inf, Fsilw_sup
      double precision:: Fsilw_curr
      integer, intent(out):: kinf, ksup
      integer:: k, nlev

      nlev = size(temp,3)
      kinf = 1
      ksup = nlev

      ! Initialisation (two levels)
      call dynsoil_geographic_loop( list_i, list_j, &
                                    area, temp(:,:,kinf), runoff(:,:,kinf), slope, lith_frac, &
                                    h, E, xs, W, Fsilw_inf, param )
      call dynsoil_geographic_loop( list_i, list_j, &
                                    area, temp(:,:,ksup), runoff(:,:,ksup), slope, lith_frac, &
                                    h, E, xs, W, Fsilw_sup, param )

      ! find by dichotomy:
      do while ( ksup-kinf > 1 )

        k = kinf + (ksup-kinf)/2

        call dynsoil_geographic_loop( list_i, list_j, &
                                      area, temp(:,:,k), runoff(:,:,k), slope, lith_frac, &
                                      h, E, xs, W, Fsilw_curr, param )

        if ( imposed_Fsilw >= Fsilw_curr ) then
          kinf = k
          Fsilw_inf = Fsilw_curr
        else
          ksup = k
          Fsilw_sup = Fsilw_curr
        end if

      end do


    end subroutine


    !------------------------------------------------------------------------------------------------------------------------------!


    subroutine inverse_interpolation(imposed_Fsilw, list_i, list_j, param, area, temp, runoff, slope, lith_frac, CO2_levels, GMST, &
                                     CO2, curr_GMST, curr_temp, curr_runf, h, E, xs, W, nstep)

      use dynsoil_steady_state_module, only: dynsoil_geographic_loop

      double precision, intent(in):: imposed_Fsilw
      integer, dimension(:), intent(in):: list_i, list_j
      double precision, dimension(:,:),   intent(in):: param ! [param x litho]
      double precision, dimension(:,:),   intent(in):: area, slope
      double precision, dimension(:,:,:), intent(in):: temp, runoff, lith_frac
      double precision, dimension(:),     intent(in):: CO2_levels, GMST
      double precision, dimension(:,:),   intent(out):: curr_temp, curr_runf
      double precision, dimension(:,:,:), intent(out):: h, E, xs, W
      double precision, intent(out):: CO2, curr_GMST
      integer, intent(out):: nstep
      double precision:: Fsilw, Fsilw0, Fsilw1, CO2_00, CO2_11, CO2_0, CO2_1, dFsilw_dCO2, xi
      integer:: k00, k11, k0, k1


      ! find interval
      call find_weathering_interval( imposed_Fsilw, list_i, list_j,                            &
                                     param, area, temp, runoff, slope, lith_frac, h, E, xs, W, &
                                     k0, k1, Fsilw0, Fsilw1                                    )

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

      do while ( abs((Fsilw - imposed_Fsilw) / imposed_Fsilw)  >  PREC )

        if (Fsilw > imposed_Fsilw) then
          Fsilw1 = Fsilw
          CO2_1 = CO2
        else
          Fsilw0 = Fsilw
          CO2_0 = CO2
        end if

        select case (INTERPOLATION_MODE)

          case ('linear')
            dFsilw_dCO2 = (Fsilw1 - Fsilw0) / (CO2_1 - CO2_0)
            !++++++++++++++++++++++++++++++++++++++++++++++++++++++!
            CO2  =  CO2_0  +  (imposed_Fsilw - Fsilw0) / dFsilw_dCO2
            !++++++++++++++++++++++++++++++++++++++++++++++++++++++!

          case ('log')
            dFsilw_dCO2 = (Fsilw1 - Fsilw0) / log(CO2_1/CO2_0)
            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
            CO2  =  CO2_0  *  exp((imposed_Fsilw - Fsilw0) / dFsilw_dCO2)
            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        end select

        ! Climate interpolation
        xi = interp_coeff(CO2, CO2_00, CO2_11)

        call interpolate_climate( list_i, list_j, temp, runoff, GMST, k0, k1, xi, curr_temp, curr_runf, curr_GMST )

        ! Dynsoil computation at new climate => new guess of Fsilw
        call dynsoil_geographic_loop( list_i, list_j,                               &
                                      area, curr_temp, curr_runf, slope, lith_frac, &
                                      h, E, xs, W, Fsilw, param                     )

        nstep = nstep + 1

      end do

    end subroutine


    !------------------------------------------------------------------------------------------------------------------------------!


end module
