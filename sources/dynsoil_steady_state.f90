module dynsoil
implicit none

contains


  function steady_state_weathering(E, x_p_surf)
    double precision, intent(in):: E, x_p_surf
    double precision:: steady_state_weathering
    !++++++++++++++++++++++++++++++++++++++!
    steady_state_weathering = E*(1-x_p_surf)
    !++++++++++++++++++++++++++++++++++++++!
  end function


  !--------------------------------------------------------------------------------------------------------------------------------!


  subroutine dynsoil_steady_state( temp,runoff,slope, h, E, xs, W )

    use dynsoil_empirical_laws, only: reg_prod_opt, erosion, eq_reg_thick, dissolution_constant, eq_x_p_surface
    !
    double precision, intent(in):: temp, runoff, slope
    double precision, intent(out):: h, E, xs, W
    double precision:: RPopt, Kmain


      ! test for non-null runoff points. If runoff is null, nothing is done
      ! because all fluxes (erosive, chemical...) are null.
      !
      if (runoff>0) then

        ! optimal regolith production rate:
        RPopt = reg_prod_opt(temp,runoff)
        ! regolith erosion rate:
        E = erosion(temp,runoff,slope)
        ! regolith thickness
        h = eq_reg_thick(RPopt,E)

        ! dissolution constant:
        Kmain = dissolution_constant(temp,runoff)

        ! Primary phases proportion at surface:
        xs = eq_x_p_surface(h, E, Kmain)

        ! Weathering rate:
        W = steady_state_weathering(E, xs)

      else

        ! null runoff: set all fluxes to 0
        h = 0
        E = 0
        xs = 1
        W  = 0
    
      end if


  end subroutine


  !--------------------------------------------------------------------------------------------------------------------------------!


  subroutine dynsoil_geographic_loop( list_i, list_j, CaMg_rock, area, temp, runoff, slope, lith_frac, h, E, xs, W, WCaMg, Fsilwth )
    integer, dimension(:), intent(in):: list_i, list_j
    double precision, dimension(:), intent(in):: CaMg_rock
    double precision, dimension(:,:), intent(in):: area, temp, runoff, slope
    double precision, dimension(:,:,:), intent(in):: lith_frac
    double precision, dimension(:,:), intent(inout):: h, E, xs, W, WCaMg
    double precision, intent(out):: Fsilwth
    integer:: i, j, k, ncont

    ncont = size(list_i)

    Fsilwth = 0

    do k = 1,ncont

      i = list_i(k)
      j = list_j(k)

      ! DynSoil (volumetric weathering)
      call dynsoil_steady_state( temp(i,j), runoff(i,j), slope(i,j), h(i,j), E(i,j), xs(i,j), W(i,j) )

      ! Lithology-dependent CaMg weathering
      WCaMg(i,j) = sum( W(i,j)*lith_frac(:,i,j)*CaMg_rock )

      ! Total CO2 consumption:
      Fsilwth = Fsilwth + area(i,j)*WCaMg(i,j)

    end do

  end subroutine


end module
