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
    reg_prod_opt = krp * runoff * dexp(-(Ea_rp/Rgas)*(1./(temp+273.15)-1./T0_rp))
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
!    soil_prod_func = ( dexp(-1*h_soil/d1) - k1*dexp(-1*h_soil/d2) ) / knorm ! normalized HUMPED FUNCTION
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
!    eq_reg_thick = d1*dlog(RPopt/E) ! for humped function
!    eq_reg_thick = h0*RPopt/E      ! for inverse law
    eq_reg_thick = max( h0*dlog(RPopt/E), 0.d+0 ) ! for exponential function
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   Estimated soil equilibrium thickness by approximating the SPF.
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function dissolution_constant(temp,runoff)
    include 'dynsoil_physical_parameters.inc'
    double precision, intent(in):: temp, runoff ! runoff (m/y) ; temperature (°C) !
    double precision:: dissolution_constant
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    dissolution_constant = kd * (1-dexp(-kw*runoff)) * dexp( (Ea/Rgas) * (1/T0-1/(temp+273.15)) )
!   rock dissolution constant, computed from climatic variables
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function eq_x_p_surface(h,E,K)
    include 'dynsoil_physical_parameters.inc'
    double precision, intent(in):: h, E, K
    double precision:: eq_x_p_surface
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    eq_x_p_surface = dexp( -1*K * ((h/E)**(sigma+1)) / (sigma+1) )
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   Primary phases proportion at surface at steady-state
  end function

!----------------------------------------------------------------------------------------------------------------------------------!


end module
