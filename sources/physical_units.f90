module physical_units
!===========================================================================!
! States the physical units of the variables used by GEOCLIM (%reference),  !
! and accepted alternative units with their conversion method:              !
!     conversion = (/factor, offset/)                                       !
!     converted_unit = factor*native_units + offset                         !
! Users can add here as many alternative units as they want by editing the  !
! subroutines that define them (def_area_units, def_temperature_units, ...) !
! Be carefull to always declare the number of accepted alternative units!   !
! If needing, the maximum number of alternative units (N_ACCEPTED_MAX) can  !
! be increased, as well as the maximum number of characters (UNITS_CHARLEN) !
!===========================================================================!

implicit none

    integer, parameter, private:: N_ACCEPTED_MAX = 100
    ! Maximum number of accepted alternative units.
    ! All variables with type "units" have the same length of alternative units array,
    ! the variable 'Naccepted' tells the actual number of alternative units

    integer, parameter, private:: UNITS_CHARLEN = 100

    double precision, parameter, private:: YRLEN = 365.2422 ! number of days in a year
    double precision, parameter, private:: RHOWAT = 1d3 ! water density, in kg/m3


    type alternative_units
        character(len=UNITS_CHARLEN):: string
        double precision, dimension(2):: conversion ! (/factor, offset/)
    end type

    !<><><><><><><><><><><><><><><><><><><><><>!
    type units
        character(len=UNITS_CHARLEN):: reference
        integer:: naccepted
        type(alternative_units), dimension(N_ACCEPTED_MAX):: accepted
    end type
    !<><><><><><><><><><><><><><><><><><><><><>!




    contains


    function not_attributed()
        type(alternative_units):: not_attributed
        not_attributed%string = '//_NOT_ATTRIBUTED_//'
        not_attributed%conversion = (/ 0d0, -1d36 /)
    end function


    !----------------------------------------------------------------------!


    function area_units()
        type(units):: area_units
        character(len=UNITS_CHARLEN):: reference_units
        integer:: naccepted
        type(alternative_units), dimension(N_ACCEPTED_MAX):: accept_units

        reference_units = 'm2'

        naccepted = 15
        !                        units name                      conversion: (/factor, offset/)
        accept_units(1)%string = 'm^2';          accept_units(1)%conversion = (/1d0, 0d0/)
        accept_units(2)%string = 'm**2';         accept_units(2)%conversion = (/1d0, 0d0/)
        accept_units(3)%string = 'km2';          accept_units(3)%conversion = (/1d6, 0d0/)
        accept_units(4)%string = 'km^2';         accept_units(4)%conversion = (/1d6, 0d0/)
        accept_units(5)%string = 'km**2';        accept_units(5)%conversion = (/1d6, 0d0/)
        accept_units(6)%string = 'Mkm2';         accept_units(6)%conversion = (/1d12, 0d0/)
        accept_units(7)%string = 'Mkm^2';        accept_units(7)%conversion = (/1d12, 0d0/)
        accept_units(8)%string = '1e6km2';       accept_units(8)%conversion = (/1d12, 0d0/)
        accept_units(9)%string = '1e6 km2';      accept_units(9)%conversion = (/1d12, 0d0/)
        accept_units(10)%string = '1e6km^2';     accept_units(10)%conversion = (/1d12, 0d0/)
        accept_units(11)%string = '1e6 km^2';    accept_units(11)%conversion = (/1d12, 0d0/)
        accept_units(12)%string = '1e12m2';      accept_units(12)%conversion = (/1d12, 0d0/)
        accept_units(13)%string = '1e12 m2';     accept_units(13)%conversion = (/1d12, 0d0/)
        accept_units(14)%string = '1e12m^2';     accept_units(14)%conversion = (/1d12, 0d0/)
        accept_units(15)%string = '1e12 m^2';    accept_units(15)%conversion = (/1d12, 0d0/)

        accept_units( naccepted+1 : N_ACCEPTED_MAX ) = not_attributed()

        area_units%reference = reference_units
        area_units%naccepted = naccepted
        area_units%accepted = accept_units

    end function


    !----------------------------------------------------------------------!


    function fraction_units()
        type(units):: fraction_units
        character(len=UNITS_CHARLEN):: reference_units
        integer:: naccepted
        type(alternative_units), dimension(N_ACCEPTED_MAX):: accept_units

        reference_units = 'dimensionless'

        naccepted = 12
        !                        units name                      conversion: (/factor, offset/)
        accept_units(1)%string = '';             accept_units(1)%conversion = (/1d0, 0d0/)
        ! Note: an empty string is what is returned if the units is not found
        accept_units(2)%string = '-';            accept_units(2)%conversion = (/1d0, 0d0/)
        accept_units(3)%string = '1';            accept_units(3)%conversion = (/1d0, 0d0/)
        accept_units(4)%string = '0-1';          accept_units(4)%conversion = (/1d0, 0d0/)
        accept_units(5)%string = '(0-1)';        accept_units(5)%conversion = (/1d0, 0d0/)
        accept_units(6)%string = '0 - 1';        accept_units(6)%conversion = (/1d0, 0d0/)
        accept_units(7)%string = '(0 - 1)';      accept_units(7)%conversion = (/1d0, 0d0/)
        accept_units(9)%string = 'unitless';     accept_units(8)%conversion = (/1d0, 0d0/)
        accept_units(9)%string = '%';            accept_units(9)%conversion = (/1d-2, 0d0/)
        accept_units(10)%string = 'percent';     accept_units(10)%conversion = (/1d-2, 0d0/)
        accept_units(11)%string = 'per cent';    accept_units(11)%conversion = (/1d-2, 0d0/)
        accept_units(12)%string = 'fraction';    accept_units(12)%conversion = (/1d0, 0d0/)

        accept_units( naccepted+1 : N_ACCEPTED_MAX ) = not_attributed()

        ! record in output:
        fraction_units%reference = reference_units
        fraction_units%naccepted = naccepted
        fraction_units%accepted = accept_units

    end function


    !----------------------------------------------------------------------!


    function slope_units()
        type(units):: slope_units
        character(len=UNITS_CHARLEN):: reference_units
        integer:: naccepted
        type(alternative_units), dimension(N_ACCEPTED_MAX):: accept_units

        reference_units = 'm/m'

        naccepted = 6
        !                        units name                      conversion: (/factor, offset/)
        accept_units(1)%string = 'dimensionless'; accept_units(1)%conversion = (/1d0, 0d0/)
        accept_units(2)%string = 'unitless';      accept_units(2)%conversion = (/1d0, 0d0/)
        accept_units(3)%string = '-';             accept_units(3)%conversion = (/1d0, 0d0/)
        accept_units(4)%string = '%';             accept_units(4)%conversion = (/1d-2, 0d0/)
        accept_units(5)%string = 'percent';       accept_units(5)%conversion = (/1d-2, 0d0/)
        accept_units(6)%string = 'per cent';      accept_units(6)%conversion = (/1d-2, 0d0/)

        accept_units( naccepted+1 : N_ACCEPTED_MAX ) = not_attributed()

        ! record in output:
        slope_units%reference = reference_units
        slope_units%naccepted = naccepted
        slope_units%accepted = accept_units

    end function


    !----------------------------------------------------------------------!


    function temperature_units()
        type(units):: temperature_units
        character(len=UNITS_CHARLEN):: reference_units
        integer:: naccepted
        type(alternative_units), dimension(N_ACCEPTED_MAX):: accept_units

        reference_units = 'degrees_celsius'

        naccepted = 22
        !                        units name                             conversion: (/factor, offset/)
        accept_units(1)%string = '°';                   accept_units(1)%conversion = (/1d0, 0d0/)
        accept_units(2)%string = '°C';                  accept_units(2)%conversion = (/1d0, 0d0/)
        accept_units(3)%string = 'C';                   accept_units(3)%conversion = (/1d0, 0d0/)
        accept_units(4)%string = 'Celsius';             accept_units(4)%conversion = (/1d0, 0d0/)
        accept_units(5)%string = 'celsius';             accept_units(5)%conversion = (/1d0, 0d0/)
        accept_units(6)%string = 'degC';                accept_units(6)%conversion = (/1d0, 0d0/)
        accept_units(7)%string = 'deg C';               accept_units(7)%conversion = (/1d0, 0d0/)
        accept_units(8)%string = 'deg_C';               accept_units(8)%conversion = (/1d0, 0d0/)
        accept_units(9)%string = 'degree C';            accept_units(9)%conversion = (/1d0, 0d0/)
        accept_units(10)%string = 'degrees C';          accept_units(10)%conversion = (/1d0, 0d0/)
        accept_units(11)%string = 'degree_C';           accept_units(11)%conversion = (/1d0, 0d0/)
        accept_units(12)%string = 'degrees_C';          accept_units(12)%conversion = (/1d0, 0d0/)
        accept_units(13)%string = 'degree_celsius';     accept_units(13)%conversion = (/1d0, 0d0/)
        accept_units(14)%string = 'degree_Celsius';     accept_units(14)%conversion = (/1d0, 0d0/)
        accept_units(15)%string = 'degrees_Celsius';    accept_units(15)%conversion = (/1d0, 0d0/)
        accept_units(16)%string = 'K';                  accept_units(16)%conversion = (/1d0, -273.15d0/)
        accept_units(17)%string = 'Kelvin';             accept_units(17)%conversion = (/1d0, -273.15d0/)
        accept_units(18)%string = 'kelvin';             accept_units(18)%conversion = (/1d0, -273.15d0/)
        accept_units(19)%string = 'F';                  accept_units(19)%conversion = (/5d0/9d0, -32d0*5d0/9d0/)
        accept_units(20)%string = '°F';                 accept_units(20)%conversion = (/5d0/9d0, -32d0*5d0/9d0/)
        accept_units(21)%string = 'Fahrenheit';         accept_units(21)%conversion = (/5d0/9d0, -32d0*5d0/9d0/)
        accept_units(22)%string = 'fahrenheit';         accept_units(22)%conversion = (/5d0/9d0, -32d0*5d0/9d0/)

        accept_units( naccepted+1 : N_ACCEPTED_MAX ) = not_attributed()

        ! record in output:
        temperature_units%reference = reference_units
        temperature_units%naccepted = naccepted
        temperature_units%accepted = accept_units

    end function


    !----------------------------------------------------------------------!


    function runoff_units()
        type(units):: runoff_units
        character(len=UNITS_CHARLEN):: reference_units
        integer:: naccepted
        type(alternative_units), dimension(N_ACCEPTED_MAX):: accept_units

        reference_units = 'm/y'

        naccepted = 32
        !                        units name                     conversion: (/factor, offset/)
        !------------------------------------------------------------------------------------------!
        accept_units(1)%string = 'mm/d';        accept_units(1)%conversion = (/1d-3*YRLEN, 0d0/)
        accept_units(2)%string = 'mm/day';      accept_units(2)%conversion = (/1d-3*YRLEN, 0d0/)
        accept_units(3)%string = 'cm/d';        accept_units(3)%conversion = (/1d-2*YRLEN, 0d0/)
        accept_units(4)%string = 'cm/day';      accept_units(4)%conversion = (/1d-2*YRLEN, 0d0/)
        accept_units(5)%string = 'm/d';         accept_units(5)%conversion = (/YRLEN, 0d0/)
        accept_units(6)%string = 'm/day';       accept_units(6)%conversion = (/YRLEN, 0d0/)
        accept_units(7)%string = 'mm/y';        accept_units(7)%conversion = (/1d-3, 0d0/)
        accept_units(8)%string = 'mm/yr';       accept_units(8)%conversion = (/1d-3, 0d0/)
        accept_units(9)%string = 'mm/a';        accept_units(9)%conversion = (/1d-3, 0d0/)
        accept_units(10)%string = 'cm/y';       accept_units(10)%conversion = (/1d-2, 0d0/)
        accept_units(11)%string = 'cm/yr';      accept_units(10)%conversion = (/1d-2, 0d0/)
        accept_units(12)%string = 'cm/a';       accept_units(11)%conversion = (/1d-2, 0d0/)
        accept_units(13)%string = 'm/yr';       accept_units(13)%conversion = (/1d0, 0d0/)
        accept_units(14)%string = 'm/a';        accept_units(14)%conversion = (/1d0, 0d0/)
        accept_units(15)%string = 'mm/h';       accept_units(15)%conversion = (/1d-3*YRLEN*24, 0d0/)
        accept_units(16)%string = 'mm/hr';      accept_units(16)%conversion = (/1d-3*YRLEN*24, 0d0/)
        accept_units(17)%string = 'mm/hour';    accept_units(17)%conversion = (/1d-3*YRLEN*24, 0d0/)
        accept_units(18)%string = 'cm/h';       accept_units(18)%conversion = (/1e-2*YRLEN*24, 0d0/)
        accept_units(19)%string = 'cm/hr';      accept_units(19)%conversion = (/1e-2*YRLEN*24, 0d0/)
        accept_units(20)%string = 'cm/hour';    accept_units(20)%conversion = (/1e-2*YRLEN*24, 0d0/)
        accept_units(21)%string = 'm/h';        accept_units(21)%conversion = (/YRLEN*24, 0d0/)
        accept_units(22)%string = 'm/hr';       accept_units(22)%conversion = (/YRLEN*24, 0d0/)
        accept_units(23)%string = 'm/hour';     accept_units(23)%conversion = (/YRLEN*24, 0d0/)
        accept_units(24)%string = 'mm/s';       accept_units(24)%conversion = (/1d-3*YRLEN*24*60*60, 0d0/)
        accept_units(25)%string = 'mm/sec';     accept_units(25)%conversion = (/1d-3*YRLEN*24*60*60, 0d0/)
        accept_units(26)%string = 'cm/s';       accept_units(26)%conversion = (/1d-2*YRLEN*24*60*60, 0d0/)
        accept_units(27)%string = 'cm/sec';     accept_units(27)%conversion = (/1d-2*YRLEN*24*60*60, 0d0/)
        accept_units(28)%string = 'm/s';        accept_units(28)%conversion = (/YRLEN*24*60*60, 0d0/)
        accept_units(29)%string = 'm/sec';      accept_units(29)%conversion = (/YRLEN*24*60*60, 0d0/)
        accept_units(30)%string = 'kg/m2/s';    accept_units(30)%conversion = (/YRLEN*24*60*60/RHOWAT, 0d0/)
        accept_units(31)%string = 'kg m-2 s-1'; accept_units(31)%conversion = (/YRLEN*24*60*60/RHOWAT, 0d0/)
        accept_units(32)%string = 'kg/(s*m2)';  accept_units(32)%conversion = (/YRLEN*24*60*60/RHOWAT, 0d0/)

        accept_units( naccepted+1 : N_ACCEPTED_MAX ) = not_attributed()

        ! record in output:
        runoff_units%reference = reference_units
        runoff_units%naccepted = naccepted
        runoff_units%accepted = accept_units

    end function


end module
