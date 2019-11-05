# GEOCLIM-dynsoil-steady-state

GEOCLIM - DynSoil - Steady-State computes the geographically-distributed weathering rate (+ erosion, regolith thickness...) at steady-state, according to the given climatology (temperature and runoff) for several CO2 levels (2 at least!), the topographic slope and the lithology fraction in each grid cell.
The Total weathering flux is assumed to equals the CO2 degassing by the solid Earth.

The model has two mode:
	1) Forward mode: impose a CO2 level to the model, it get the corresponding climate field and computes weathering
	2) Backward mode: impose a volcanic degassing to the model, it finds by inversion the CO2 level for which the total CO2 consumption equals the degassing. The precision level can be set in the main program "gdss_mainprog.f90".

The continental climate (ie: temperature and runoff) is interpolated linearly between the different CO2 levels. At this time, the algorithm use for the interpolation does not allow extrapolation beyond the lowest and highest level.
The climate interpolation is done by the subroutines of the module "climate_module.f90". DynSoil integration is done by subroutines of the module "dynsoil_steady_state_module.f90". The inversion (in backward mode) is done in the main program.

In addition the model can be run either in single run mode, for which it computes weathering 1 time for 1 given set of parameter, or in multiple runs mode, for which an netCDF unlimited dimension is created, and the model is run as many time as number of given sets of parameters.

## Repository architecture
The input/output information is stated in the file "IO_INTERFACE.txt", as well as the flags for model mode (Forward, Backward, Single run, Multirun). This file is in the root directory.
It includes the names of input files to use of continent area, lithology, temperature, runoff, slope. The model check the consistency of this input dataset (size, coordinates, units) and DynSoil is run accordingly this geographic-climatic-lithologic settings.
These input files should be at netCDF format. The names of the dimensions on which are recorded the input variable in these netCDF files has also to be stated in "IO_INTERFACE.txt", section INPUT CONDITION.

The executable files is meant to be run in the directory "executable/". For this reason, the program looks for the relative path "../IO_INTERFACE.txt", as specified line 11 of the main program "gdss_mainprog.f90". For this reason, the code will crash if it is not run directory one level after the root directory. Every other paths can be changed in the file "../IO_INTERFACE.txt" without recompiling the code.

The physical parameters of DynSoil and the forcings (CO2 or degassing) are read in two separated text files whose names are stated in "IO_INTERFACE.txt". See the documentation in "IO_INTERFACE.txt" for the format of these files.

Conventionally, the parameter file is in the repertory "parameters/", the forcings file is in the repertory "forcings/" and the outputs are stored in the repertory "outputs/". Templates can be found in those repertories. Note that the reading format of the parameter file is not the same in the single run or the multirun mode.

"IO_INTERFACE.txt" is read by the subroutine 'make_input_output' in the file "io_module.f90", it creates the netCDF output file (conventionnaly in the "output/" directory) plus a fortran scratch text file with the ID of the variables to record (unit=IOUT), and potentially two others scratch files (unit=IFORC and unit=IPARAM) recording the forcings and the parameters to be used in the main program.

## Compilance checks
The program conducts a certain number of compliance checks in the subroutine "make_input_output" (file "io_module.f90"), that includes checking the forcing variables units (area, temperature, runoff, slope...), the consistency of the geographic axes of all forcing variables and the consistency of their continental mask. For this last check only, to program is not automatically aborted in case of failure, but the user is asked for continuing to run it (0 to abort, 1 to continue). This argument can be enter interactively, or the user can enter it as an argument of the executable:
./executable 1   => run the program whatever the case
./executable 0   => abort in case of continental mask inconsistency

## Compilation tips:
This program uses syntax allowed in fortran 2003 or later (subroutines with allocatable dummy arguments, and the use of the subroutine "get_command_argument"). Make sure your compiler support this syntax (with gfortran: "-std=f2003")

## Parameter exploration:
This repository also contains a version of the model where the user can vary the parameters for a unique given CO2 level. It can be found in the repertory "preprocessing/parameter_exploration/". See the README in this last repertory for more information.
