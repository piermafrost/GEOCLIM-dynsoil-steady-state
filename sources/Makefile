# This is a makefile I designed with my preferred flags for the two fortran compilers I mostly use:
# gfortran and ifort.
# It can be use in 3 modes: "standard", "optim" (for optimized execution, with no security check)
# and "debug" (with extra debugging options)
# The default mode is "standard" (stated by the variable 'MODE') and the default compiler is
# "gfortran" (stated by the variable 'FC').
# The other configuration variable is 'USE_NETCDF', = "yes", by default, as netCDF library is
# needed by the program. This variable shouldn't be changed.
# The path of the netCDF library is stated by the variable 'ncpath' lines 127-135. This is likely
# computer-dependent, change it if necessary.
# The name of the executable is 'gdss'
#
# Use this Makefile as follows:
#   `make executable_name MODE=your_chose_mode FC=your_chose_compiler USE_NETCDF=yes_or_no`
# or simply
#   `make`
# that will do the same as
#   `make gdss MODE=standard FC=gfortran USE_NETCDF=yes`
#
# `make clean` will delete all the intermediate compilation files
# `make clc` will delete them and delete the executable
# `make check` will display the compilation command without executing anything
# eg: `make check FC=ifort MODE=debug`




# Compililation mode:
MODE ?= standard


# Use netCDF library:
USE_NETCDF ?= yes


# Macros:
FC ?= gfortran
DUMMYFC = $(FC)
# If you directly do `ifeq ($(FC), gfortran)`, for some osbcure reasons, it doesn't work.
# Whereas `ifeq ($(DUMMYFC), gfortran)` works.


# compiler-dependent flags:

                            ##########################
ifeq ($(DUMMYFC), gfortran) # GNU compiler: gfortran #
                            ##########################

        # Fortran language options:
        #lang_flags = -ffree-form -std=f95 # Basic fortran standard
        lang_flags = -ffree-form -std=f2003 # some syntax of my programs are not legal in fortran 95  standard
                                            # (like the function 'get_command_argument')
        #lang_flags = -ffree-form -std=f2008 # for using 'kind' syntax

        # Error and warning options:
        ifeq ($(MODE), debug)
                warn_flags = -pedantic-errors -Wall -Wconversion-extra -Wimplicit-interface -Wunderflow -Wextra -Wunreachable-code
        else
                warn_flags = -Wpedantic -Waliasing -Wampersand -Wline-truncation -Wcharacter-truncation -Wtabs -Wconversion -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wunused-dummy-argument -Wunused -Wunreachable-code
        endif

        # Debugging options:
        ifneq ($(MODE), optim)
                debug_flags = -ffpe-trap=invalid,zero,overflow -g3
        endif

        # Code generation options:
        ifneq ($(MODE), optim)
                code_flags = -fbounds-check
        endif

        # Optimization options:
        ifeq ($(MODE), debug)
                optim_flags = -O0 -fstack-protector-all
        endif
        ifeq ($(MODE), optim)
                optim_flags = -O3
        endif
        # if none of previous: assume standard mode:
        optim_flags ?= -O1

endif

                         #########################
ifeq ($(DUMMYFC), ifort) # Intel compiler: ifort #
                         #########################

        # Fortran langauge options:
        lang_flags = -free -132

        # Error and warning options (these are the ones I think necessary, some of them are already set by default):
        ifeq ($(MODE), debug)
                warn_flags = -warn all
        else
                warn_flags = -warn general,usage,declaration,unused,truncated_source,interfaces,ignore_loc
        endif

        # Debugging options:
        ifeq ($(MODE), debug)
                debug_flags = -debug full
        endif

        # Options at run time
        ifeq ($(MODE), debug)
                code_flags = -check all -fp-stack-check -traceback
        endif

        # Optimization options:
        ifeq ($(MODE), debug)
                optim_flags = -O0 -fstack-protector-all -fstack-security-check
        endif
        ifeq ($(MODE), optim)
                optim_flags = -Ofast
        endif
        # if none of previous: assume standard mode:
        optim_flags ?= -O1

endif



########################
# LIBRARY and INCLUDE: #
########################


# NetCDF library:
#!!!!!!!!!!!!!!!!!!!!!!!!!!!#
ncpath = /usr
#!!!!!!!!!!!!!!!!!!!!!!!!!!!#
ifeq ($(USE_NETCDF), yes)
        inc_flags += -I$(ncpath)/include
        lib_flags += -L$(ncpath)/lib
        netcdf_flags = -lnetcdf -lnetcdff
endif



FFLAGS = $(lang_flags) $(warn_flags) $(debug_flags) $(code_flags) $(optim_flags) $(inc_flags) $(lib_flags) $(netcdf_flags)


#############################################################################

ext = f90


# Specific executable making rule:
gdss: ascii_io_functions.o climate_module.o dynsoil_empirical_laws.o dynsoil_steady_state_module.o gdss_mainprog.o io_module.o netcdf_io_functions.o 
	$(FC) $(FFLAGS) $^ -o $@


# General object making rule
%.o: %.$(ext)
	$(FC) -c $(FFLAGS) $< -o $@



# Additional dependencies

io_module.o: netcdf_io_functions.o ascii_io_functions.o
dynsoil_steady_state_module.o: dynsoil_empirical_laws.o  
climate_module.o: dynsoil_steady_state_module.o
gdss_mainprog.o: io_module.o netcdf_io_functions.o dynsoil_steady_state_module.o climate_module.o




#############################################################################


.PHONY: check clean clc

check:
	@echo "$(FC) $(FFLAGS)"

clean:
	rm -f *.o *.mod

clc: clean
	rm -f gdss
