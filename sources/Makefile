# This is a makefile for GEOCLIM-DynSoil-steady-state. Support both intel and GNU compiler.
#
# This makefile includes one optional argument variable: MODE, defining the compilation model
# By default, MODE = standard, but two other modes can be specified while calling make: "debug" and "optim"
#     make MODE=standard    (equivalent to just "make")
#     make MODE=optim       (for fast use)
#     make MODE=debug       (extra debugging options)
#
# The path for netCDF library is defined using the shell command `nc-config --prefix` (line 36)
# If you encounter issues with netCDF library, try and redefine the path in the `make` command
#
# Use this Makefile as follows:
#   `make executable_name MODE=your_chosen_mode FC=your_chosen_compiler ncpath=netcdf_libray_path`
# or simply
#   `make`
# that will do the same as
#   `make gdss MODE=standard FC=gfortran`
#
# `make clean` will delete all the intermediate compilation files
# `make clc` will delete them and delete the executable
# `make check` will display the compilation command without executing anything
# eg: `make check FC=ifort MODE=debug`


#################################################
# Main customization variables (default values) #
#################################################

# Fortran compiler
FC = gfortran

# Compilation mode (default):
MODE = standard

# Path of netCDF library (must contains subrepertories lib/ and include/)
ncpath = $(shell nc-config --prefix)


#############################
# Compiler-dependent flags: #
#############################


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

# Error and warning options:
ifeq ($(MODE), debug)
warn_flags = -Wall -Wconversion-extra -Wimplicit-interface -Wextra -Wunreachable-code
else
warn_flags = -Wpedantic -Waliasing -Wampersand -Wline-truncation -Wcharacter-truncation -Wconversion -Wimplicit-interface -Wimplicit-procedure -Wunused-dummy-argument -Wunreachable-code
endif

# Error and warning options:
ifeq ($(MODE), debug)
warn_flags = -Wall -Wconversion-extra -Wimplicit-interface -Wextra -Wunreachable-code
else
warn_flags = -Waliasing -Wampersand -Wline-truncation -Wcharacter-truncation -Wconversion -Wimplicit-interface -Wimplicit-procedure -Wunused-dummy-argument -Wunreachable-code
ifeq ($(MODE), standard)
warn_flags += -Wpedantic
endif
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
inc_flags += -I$(ncpath)/include
lib_flags += -L$(ncpath)/lib
netcdf_flags = -lnetcdf -lnetcdff


##############################################
# Complete list of fortran compilation flags #
##############################################

main_flags = $(lang_flags) $(warn_flags) $(debug_flags) $(code_flags) $(optim_flags)
FFLAGS = $(main_flags) $(inc_flags) $(lib_flags) $(netcdf_flags)


#############################################################################

#############
# Root path #
#############

# Check if path assignment uses correct path
# If not, add 'updatepath' to the executable dependencies

# get root path (remove 'sources' from current directory path)
pathcommand = character(len=*), parameter:: gdss_root_path = "$(CURDIR:sources=)"
# get path assignment command (last line of 'path.inc' file)
currpathcommand = $(shell tail -1 path.inc)

ifneq ($(pathcommand), $(currpathcommand))
pathsource = updatepath path.inc
else
pathsource = path.inc
endif



#############################################################################


########################
# Compilation commands #
########################


# Executable file:
execut = ../executables/gdss


# Specific executable making rule:
$(execut): ascii_io_functions.o climate_module.o dynsoil_empirical_laws.o dynsoil_steady_state_module.o gdss_mainprog.o io_module.o netcdf_io_functions.o physical_units.o
	$(FC) $^ -o $@ $(FFLAGS)


# General object making rule
%.o: %.f90
	$(FC) $< -o $@ -c $(FFLAGS)


# File dependencies

dynsoil_steady_state_module.o: dynsoil_empirical_laws.o
climate_module.o: dynsoil_steady_state_module.o
io_module.o: ascii_io_functions.o climate_module.o netcdf_io_functions.o physical_units.o $(pathsource)
gdss_mainprog.o: climate_module.o dynsoil_steady_state_module.o io_module.o netcdf_io_functions.o



#############################################################################


############################
# Update root path command #
############################

# create path.inc file if it doesn't exist
path.inc:
	touch $@

.PHONY: updatepath # make sure to always execute that command when called
updatepath:
	echo '! File automatically generated by `make`' > path.inc
	echo '$(pathcommand)' >> path.inc


####################
# Special commands #
####################

.PHONY: clean clc echo echo_FC

echo:
	@echo "$(FC) $(FFLAGS)"

echo_FC:
	@echo $(FC)

clean:
	rm -f *.o *.mod *__gen_mod* *__genmod*
clc:
	rm -f *.o *.mod *__gen_mod* *__genmod* $(execut)

