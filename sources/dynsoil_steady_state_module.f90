module dynsoil_steady_state_module
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


  subroutine dynsoil_steady_state( temp, runoff, slope, h, E, xs, W, param )

    use dynsoil_empirical_laws, only: reg_prod_opt, erosion, eq_reg_thick, dissolution_constant, eq_x_p_surface
    !
    double precision, intent(in):: temp, runoff, slope
    double precision, dimension(:), intent(out):: h, E, xs, W
    double precision, dimension(:,:), intent(in):: param ! [param x litho]
    double precision:: RPopt, Kmain
    integer:: k, nlith

    nlith = size(param, 2)

    do k = 1,nlith

      ! regolith erosion rate:
      E(k) = erosion(temp, runoff, slope, param(:,k))

      ! test for non-null erosion. If erosion is null (for instance, because of null runoff,
      ! or erosion lower than numerical precision), nothing is done, all fluxes
      ! (erosive, chemical...) are set to 0, as reg thick, and x_p_surf to 1.
      !
      if (E(k)>0) then

        ! optimal regolith production rate:
        RPopt = reg_prod_opt(temp, runoff, param(:,k))
        ! regolith thickness
        h(k) = eq_reg_thick(RPopt, E(k), param(:,k))

        ! dissolution constant:
        Kmain = dissolution_constant(temp, runoff, param(:,k))

        ! Primary phases proportion at surface:
        xs(k) = eq_x_p_surface(h(k), E(k), Kmain, param(:,k))

        ! Weathering rate:
        W(k) = steady_state_weathering(E(k), xs(k))

      else

        ! null runoff: set all fluxes to 0
        h(k) = 0
        E(k) = 0
        xs(k) = 1
        W(k)  = 0
    
      end if

    end do


!    !****************************************************************!
!    !
!    ! GEOCLIM "old" weathering empirical laws (from Oliva et al., 2003)
!    ! Designed for 2 LITHOLOGIES: Basalts (#1) and Granit (#2).
!    ! These formulas are not part of DynSoil model.
!    !
!    ! Basalt weathering
!    W(1)   =   6d-3  *  dexp((-42300./8.134)*(1./(temp+273.15) - 1./288.15))  *  100*runoff  /  param(13,1)
!    !                                                                            ^^^^        ^^^^^^^^^^^^^^
!    ! Granit weathering
!    W(2)   =   8d-4  *  dexp((-48200./8.134)*(1./(temp+273.15) - 1./288.15))  *  100*runoff  /  param(13,2)
!    !                                                                            ^^^^        ^^^^^^^^^^^^^^
!    !                                                    Convert runoff in cm/yr |           |
!    !                     Divide by [CaMg] amount to cancel it (will be multiplied by later) |


  end subroutine


  !--------------------------------------------------------------------------------------------------------------------------------!


  subroutine dynsoil_geographic_loop(list_i, list_j, area, temp, runoff, slope, lith_frac, h, E, xs, W, Fsilwth, param)
    integer, dimension(:), intent(in):: list_i, list_j
    double precision, dimension(:,:), intent(in):: area, temp, runoff, slope
    double precision, dimension(:,:,:), intent(in):: lith_frac
    double precision, dimension(:,:,:), intent(out):: h, E, xs, W
    double precision, intent(out):: Fsilwth
    double precision, dimension(:,:), intent(in):: param ! [param x litho]
    integer:: i, j, k, ncont

    ncont = size(list_i)

    Fsilwth = 0

    do k = 1,ncont

      i = list_i(k)
      j = list_j(k)

      ! DynSoil (volumetric weathering)
      call dynsoil_steady_state( temp(i,j), runoff(i,j), slope(i,j), h(:,i,j), E(:,i,j), xs(:,i,j), W(:,i,j), param )

      ! Total CO2 consumption:
      Fsilwth   =   Fsilwth  +  area(i,j) * sum( W(:,i,j)*lith_frac(:,i,j)*param(13, :) )
      !   note: [CaMg] = param(13, i_litho)

    end do

  end subroutine


  !--------------------------------------------------------------------------------------------------------------------------------!


  subroutine CaMg_weathering(list_i, list_j, W, param)
    integer, dimension(:), intent(in):: list_i, list_j
    double precision, dimension(:,:,:), intent(inout):: W
    double precision, dimension(:,:), intent(in):: param
    integer:: i, j, k, ncont

    ncont = size(list_i)

    do k = 1,ncont

      i = list_i(k)
      j = list_j(k)

      W(:,i,j) = param(13,:)*W(:,i,j)
      !   note: [CaMg] = param(13, i_litho)

    end do

  end subroutine


  !--------------------------------------------------------------------------------------------------------------------------------!


  subroutine litho_average(list_i, list_j, lith_frac, varin, varout)
    integer, dimension(:), intent(in):: list_i, list_j
    double precision, dimension(:,:,:), intent(in):: lith_frac, varin
    double precision, dimension(:,:), intent(out):: varout
    integer:: i, j, k, ncont

    ncont = size(list_i)

    do k = 1,ncont

      i = list_i(k)
      j = list_j(k)

      varout(i,j) = sum(lith_frac(:,i,j)*varin(:,i,j))

    end do

  end subroutine


end module
