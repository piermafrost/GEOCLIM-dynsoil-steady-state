module climate_module
  implicit none

  contains

    subroutine get_climate( CO2, CO2_levels, temperature, runoff, list_i, list_j, curr_temp, curr_runf )
      double precision, intent(in):: CO2
      double precision, dimension(:), intent(in):: CO2_levels
      double precision, dimension(:,:,:), intent(in):: temperature, runoff
      integer, dimension(:), intent(in):: list_i, list_j
      double precision, dimension(:,:), intent(inout)::  curr_temp, curr_runf
      double precision:: xi
      integer:: k0, k1

      call find_co2_interval( CO2, CO2_levels, k0, k1, xi )

      ! Climate interpolation
      call interpolate_climate( list_i, list_j, temperature, runoff, k0, k1, xi, curr_temp, curr_runf )

    end subroutine


    !------------------------------------------------------------------------------------------------------------------------------!


    subroutine interpolate_climate( list_i, list_j, temp_lvl, runoff_lvl, k0, k1, xi, temp, runoff )
      integer, dimension(:), intent(in):: list_i, list_j
      double precision, dimension(:,:,:), intent(in):: temp_lvl, runoff_lvl
      integer, intent(in):: k0, k1
      double precision, intent(in):: xi
      double precision, dimension(:,:), intent(inout):: temp, runoff
      integer:: i, j, k, ncont

      ncont = size(list_i)

      do k = 1,ncont

        i = list_i(k)
        j = list_j(k)

        ! Linear interpolation
        temp(i,j) = (1-xi)*temp_lvl(i,j,k0) + xi*temp_lvl(i,j,k1)
        runoff(i,j) = (1-xi)*runoff_lvl(i,j,k0) + xi*runoff_lvl(i,j,k1)

      end do

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
      xi = ( CO2 - CO2_levels(kinf) )  /  ( CO2_levels(ksup) - CO2_levels(kinf) )

    end subroutine


    !------------------------------------------------------------------------------------------------------------------------------!


    subroutine find_weathering_interval( imposed_Fsilw, &
                    list_i, list_j, CaMg_rock, area, temp, runoff, slope, lith_frac, h, E, xs, W, WCaMg, &
                    kinf, ksup, Fsilw_inf, Fsilw_sup )
      ! Find for a given total CO2 consumpotion FSILW the two adjacent levels CO2_inf and CO2_sup such as:
      !     Fsilw(CO2_inf) <= FSILW < Fsilw(CO2_sup)

      use dynsoil_steady_state, only: dynsoil_geographic_loop

      double precision, intent(in):: imposed_Fsilw
      integer, dimension(:), intent(in):: list_i, list_j
      double precision, dimension(:), intent(in):: CaMg_rock
      double precision, dimension(:,:), intent(in):: area, slope
      double precision, dimension(:,:,:), intent(in):: temp, runoff, lith_frac
      double precision, dimension(:,:), intent(inout):: h, E, xs, W, WCaMg
      double precision, intent(out):: Fsilw_inf, Fsilw_sup
      double precision:: Fsilw_curr
      integer, intent(out):: kinf, ksup
      integer:: k, nlev

      nlev = size(temp,3)
      kinf = 1
      ksup = nlev

      ! Initialisation (two levels)
      call dynsoil_geographic_loop( list_i, list_j, &
                                    CaMg_rock, area, temp(:,:,kinf), runoff(:,:,kinf), slope, lith_frac, &
                                    h, E, xs, W, WCaMg, Fsilw_inf )
      call dynsoil_geographic_loop( list_i, list_j, &
                                    CaMg_rock, area, temp(:,:,ksup), runoff(:,:,ksup), slope, lith_frac, &
                                    h, E, xs, W, WCaMg, Fsilw_sup )

      ! find by dichotomy:
      do while ( ksup-kinf > 1 )

        k = kinf + (ksup-kinf)/2

        call dynsoil_geographic_loop( list_i, list_j, &
                                      CaMg_rock, area, temp(:,:,k), runoff(:,:,k), slope, lith_frac, &
                                      h, E, xs, W, WCaMg, Fsilw_curr )

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


end module
