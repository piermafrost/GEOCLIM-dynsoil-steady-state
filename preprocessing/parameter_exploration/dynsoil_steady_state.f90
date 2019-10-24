module dynsoil
implicit none

contains


subroutine dynsoil_steady_state( temp,runoff,slope,h,E, FPdiss )
use dynsoil_empirical_laws, only: dissolution_constant

include 'dynsoil_physical_parameters.inc'
!!
!!!!!!!!!!!!!!!!!!!!!!! INPUT-OUTPUT VARIABLES: !!!!!!!!!!!!!!!!!!!!!!
double precision, intent(in):: temp, runoff, slope, h, E            !!
double precision, intent(out):: FPdiss                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERNAL VARIABLES: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision:: Kmain,xs                                                        !!
!                                                                                  !!
integer:: k                                                                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  ! test for non-null runoff points. If runoff is null, nothing is done
  ! because all fluxes (erosive, chemical...) are null.
  ! * * * * * * * * * !
  if (runoff>0) then  !
  ! * * * * * * * * * !



    !====================================================!
    !========= climate-only dependent variables =========!
    !====================================================!

    ! dissolution constant:
    Kmain = dissolution_constant(temp,runoff)




    !===============================================================================================================!
    !==============                                    WEATHERING                                    ===============!
    !===============================================================================================================!
    !                                                                                                               !
    !           Solve the equations  ' dz/dt = RP + dissrate*dz/dx ' & ' dtau/dt = 1 + dissrate*dtau/dx '           !
    !           AT STEADY-STATE (West 2012 weathering model)                                                        !
    !                                                                                                               !
    !===============================================================================================================!


    !*********************************!
    ! REGOLITH SURFACE (xs and taus): !
    !*********************************!

    xs = exp( -1*Kmain * ((h/E)**(sigma+1)) / (sigma+1) )

    FPdiss = E*(1-xs)


  ! * * !
  else  !
  ! * * !

    ! null runoff: set all fluxes to 0
    FPdiss  = 0
    

  ! * * * !
  end if  !
  ! * * * !




end subroutine


end module
