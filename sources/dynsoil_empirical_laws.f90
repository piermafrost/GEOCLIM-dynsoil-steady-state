module dynsoil_empirical_laws
implicit none

  double precision, parameter:: Rgas = 8.314

contains


  ! PARAMETERS:
  ! nparams [x nlitho x nrun]
  ! #1:  ke
  ! #2:  a
  ! #3:  b
  ! #4:  krp
  ! #5:  EA_rp
  ! #6:  T0_rp
  ! #7:  h0
  ! #8:  kd
  ! #9:  kw
  ! #10: EA
  ! #11: T0
  ! #12: sigma
  ! #13: CaMg


!----------------------------------------------------------------------------------------------------------------------------------!

  function erosion(temp, runoff, slope, param)
    ! Physical erosion rate
    double precision:: erosion
    double precision, intent(in):: temp, runoff, slope ! runoff (m/y) ; temperature (°C) !
    double precision, dimension(13):: param
    double precision:: ke, a, b
    ke = param(1)
    a  = param(2)
    b  = param(3)
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !erosion = ke * (runoff**a) * (slope**b) * max( 2. , temp ) ! BQART-like
    erosion = ke * (runoff**a) * (slope**b)
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function reg_prod_opt(temp, runoff, param) 
!   Optimal soil production rate. Should be multiplied by the soil production function to get the actual soil production rate
    double precision:: reg_prod_opt
    double precision, intent(in):: temp, runoff ! runoff (m/y) ; temperature (°C) !
    double precision, dimension(13):: param
    double precision:: krp, Ea_rp, T0_rp
    krp   = param(4)
    Ea_rp = param(5)
    T0_rp = param(6)
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    reg_prod_opt = krp * runoff * dexp(-(Ea_rp/Rgas)*(1./(temp+273.15)-1./T0_rp))
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

!  function soil_prod_func(h_soil, param)
!    double precision:: soil_prod_func
!!    double precision, parameter:: knorm = (k1*d1/d2)**(-d2/(d1-d2)) - k1*((k1*d1/d2)**(-d1/(d1-d2)))  ! maximum of non-normalized SPF
!    double precision, intent(in):: h_soil
!    double precision, dimension(13):: param
!    double precision:: h0
!    h0 = param(7)
!!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!    soil_prod_func = ( dexp(-1*h_soil/d1) - k1*dexp(-1*h_soil/d2) ) / knorm ! normalized HUMPED FUNCTION
!    soil_prod_func = h0/h_soil                                            ! inverse function
!!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!   Soil production function 
!  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function eq_reg_thick(RPopt, E, param)
    double precision:: eq_reg_thick
    double precision, intent(in):: RPopt,E
    double precision, dimension(13):: param
    double precision:: h0
    h0 = param(7)
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    eq_reg_thick = d1*dlog(RPopt/E) ! for humped function (estimation by approximatig the SPF)
!    eq_reg_thick = h0*RPopt/E      ! for inverse law
    eq_reg_thick = max( h0*dlog(RPopt/E), 0.d+0 ) ! for exponential function
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function dissolution_constant(temp, runoff, param)
!   rock dissolution constant, computed from climatic variables
    double precision:: dissolution_constant
    double precision, intent(in):: temp, runoff ! runoff (m/y) ; temperature (°C) !
    double precision, dimension(13):: param
    double precision:: kd, kw, Ea, T0
    kd = param(8)
    kw = param(9)
    Ea = param(10)
    T0 = param(11)
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    dissolution_constant = kd * (1-dexp(-kw*runoff)) * dexp( (Ea/Rgas) * (1/T0-1/(temp+273.15)) )
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function eq_x_p_surface(h, E, K, param)
!   Primary phases proportion at surface at steady-state
    double precision:: eq_x_p_surface
    double precision, intent(in):: h, E, K
    double precision, dimension(13):: param
    double precision:: sigma
    sigma = param(12)
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    eq_x_p_surface = dexp( -1*K * ((h/E)**(sigma+1)) / (sigma+1) )
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  end function

!----------------------------------------------------------------------------------------------------------------------------------!


end module
