module dynsoil_empirical_laws
implicit none

contains


!----------------------------------------------------------------------------------------------------------------------------------!

  function erosion(temp,runoff,slope)
    include 'dynsoil_physical_parameters.inc'
    double precision, intent(in):: temp, runoff, slope ! runoff (m/y) ; temperature (°C) !
    double precision:: erosion
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !erosion = ke * (runoff**a) * (slope**b) * max( 2. , temp ) ! BQART exponents
    erosion = ke * (runoff**a) * (slope**b)
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function reg_prod_opt(temp,runoff) 
    include 'dynsoil_physical_parameters.inc'
    double precision, intent(in):: temp, runoff ! runoff (m/y) ; temperature (°C) !
    double precision:: reg_prod_opt
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    reg_prod_opt = krp * runoff * exp(-(Ea_rp/Rgas)*(1./(temp+273.15)-1./T0_rp))
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   Optimal soil production rate. Should be multiplied by the soil production function to get the actual soil production rate
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function soil_prod_func(h_soil)
    include 'dynsoil_physical_parameters.inc'
!    double precision, parameter:: knorm = (k1*d1/d2)**(-d2/(d1-d2)) - k1*((k1*d1/d2)**(-d1/(d1-d2)))  ! maximum of non-normalized SPF
    double precision, intent(in):: h_soil
    double precision:: soil_prod_func
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    soil_prod_func = ( exp(-1*h_soil/d1) - k1*exp(-1*h_soil/d2) ) / knorm ! normalized HUMPED FUNCTION
    soil_prod_func = h0/h_soil                                            ! inverse function
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   Soil production function 
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function eq_reg_thick(RPopt,E)
    include 'dynsoil_physical_parameters.inc'
    double precision, intent(in):: RPopt,E
    double precision:: eq_reg_thick
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    eq_reg_thick = d1*log(RPopt/E) ! for humped function
!    eq_reg_thick = h0*RPopt/E      ! for inverse law
    eq_reg_thick = max( h0*log(RPopt/E), 0.d+0 ) ! for exponential function
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   Estimated soil equilibrium thickness by approximating the SPF.
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function dissolution_constant(temp,runoff)
    include 'dynsoil_physical_parameters.inc'
    double precision, intent(in):: temp, runoff ! runoff (m/y) ; temperature (°C) !
    double precision:: dissolution_constant
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    dissolution_constant = kd * (1-exp(-kw*runoff)) * exp( (Ea/Rgas) * (1/T0-1/(temp+273.15)) )
!   rock dissolution constant, computed from climatic variables
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function eq_x_p_surface(h,E,K)
    include 'dynsoil_physical_parameters.inc'
    double precision, intent(in):: h, E, K
    double precision:: eq_x_p_surface
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    eq_x_p_surface = exp( -1*K * ((h/E)**(sigma+1)) / (sigma+1) )
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   Primary phases proportion at surface at steady-state
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

!  function SP_fraction(temp,runoff,h_soil)
!    include 'dynsoil_physical_parameters.inc'
!    double precision, intent(in):: temp, runoff, h_soil ! temperature (°C) ; runoff (m/y) ; h_soil (m)
!    double precision:: TauWat, TauSP ! water residence time (y) ; SP charcteristic time (y)
!    double precision:: SP_fraction
!    if (runoff > 0) then
!!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      TauWat = reg_porosity * h_soil / runoff
!      TauSP = TauSP0 !* exp((Ea_SP/Rgas)*(1./(temp+273.15)-1./T0))
!!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      SP_fraction = xSP_min  +  (xSP_max-xSP_min) / ( 1 + (TauSP/TauWat) )
!!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    else
!      SP_fraction = xSP_max
!    end if
!  end function SP_fraction
!
!!----------------------------------------------------------------------------------------------------------------------------------!
!
!  function Li_fractionation(temp)
!   include 'dynsoil_physical_parameters.inc'
!   double precision, intent(in):: temp ! temperature (°C)
!   double precision:: Li_fractionation
!   Li_fractionation = aDland/(temp+273.15)**2 + bDland
!  end function


end module
