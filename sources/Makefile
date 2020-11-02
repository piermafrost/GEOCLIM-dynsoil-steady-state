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


ifeq ($(findstring pgfortran, $(FC)), pgfortran)

#=========================#
# PGI compiler: pgfortran #
#=========================#

# Fortran language option
lang_flags = -Mfreeform

# Error and warning options:
ifeq ($(MODE), optim)
warn_flags = -Minform=severe
else ifeq ($(MODE), debug)
warn_flags = -Minform=inform
else # assume standard mode
warn_flags = -Minform=warn
endif

# Debugging options:
ifeq ($(MODE), debug)
debug_flags = -g
endif

# Code generation options:
ifneq ($(MODE), optim)
code_flags = -Mbounds
endif
ifeq ($(MODE), debug)
code_flags += -traceback
endif

# Optimization options:
ifeq ($(MODE), debug)
optim_flags = -O0
else ifeq ($(MODE), optim)
optim_flags = -O4 -fast
else
# if none of previous: assume standard mode:
optim_flags = -O2
endif



else ifeq ($(findstring gfortran, $(FC)), gfortran)

#========================#
# GNU compiler: gfortran #
#========================#

# Fortran language options:
lang_flags = -ffree-form -std=gnu
# '-std=f2003', '-std=f2008'

# Error and warning options:
ifeq ($(MODE), debug)
warn_flags = -Wall -Wconversion-extra -Wimplicit-interface -Wextra -Wunreachable-code
else
warn_flags = -Wpedantic -Waliasing -Wampersand -Wline-truncation -Wcharacter-truncation -Wconversion -Wimplicit-interface -Wimplicit-procedure -Wunused-dummy-argument -Wunreachable-code
endif

# Debugging options:
ifneq ($(MODE), optim)
debug_flags = -ffpe-trap=invalid,zero,overflow
endif

# Code generation options:
ifneq ($(MODE), optim)
code_flags = -fbounds-check
endif
ifeq ($(MODE), debug)
code_flags += -fbacktrace -g3
endif

# Optimization options:
ifeq ($(MODE), debug)
optim_flags = -O0 -fstack-protector-all
else ifeq ($(MODE), optim)
optim_flags = -O3
else
# if none of previous: assume standard mode:
optim_flags = -O1
endif



else ifeq ($(findstring ifort, $(FC)), ifort)

#========================
# Intel compiler: ifort #
#========================

# Fortran language options:
lang_flags = -free -132

# Error and warning options:
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
else ifeq ($(MODE), normal)
code_flags = -check bounds
endif

# Optimization options:
ifeq ($(MODE), debug)
optim_flags = -O0 -fstack-protector-all -fstack-security-check
else ifeq ($(MODE), optim)
optim_flags = -Ofast
else
# if none of previous: assume standard mode:
optim_flags = -O1
endif



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
