# GEOCLIM-dynsoil-steady-state

GEOCLIM - DynSoil - Steady-State computes the geographically-distributed weathering rate (+ erosion, regolith thickness...) at steady-state, according to the given climatology (temperature and runoff) for several CO2 levels, the topographic slope and the lithology fraction in each grid cell.
The Total weathering flux is assumed to equals the CO2 degassing by the solid Earth.

The model has two mode:
	1) Forward mode: impose a CO2 level to the model, it get the corresponding climate field and computes weathering
	2) Backward mode: impose a volcanic degassing to the model, it finds by inversion the CO2 level for which the total CO2 consumption equals the degassing. The precision level can be set in the main program "dynsoil_climate_mainprog.f90".

In addition the model can be run either in single run mode, for which it computes weathering 1 time for 1 given set of parameter, or in multiple runs mode, for which an netCDF unlimited dimension is created, and the model is run as many time as number of given sets of parameters.

The input/output information is stated in the file "IO_INTERFACE.txt", as well as the flags for model mode (Forward, Backward, Single run, Multirun).
It includes the names of input files to use of continent area, lithology, temperature, runoff, slope. The model check the consistency of this input dataset (size, coordinates, units) and DynSoil is run accordingly this geographic-climatic-lithologic settings.
These input files should be at netCDF format. The names of the dimensions on which are recorded the input variable in these netCDF files has also to be stated in "IO_INTERFACE.txt", section INPUT CONDITION.

The physical parameters of DynSoil and the forcings (CO2 or degassing) are read in two separated text files whose names are stated in "IO_INTERFACE.txt". See the documentation in "IO_INTERFACE.txt" for the format of these files.

Conventionally, the parameter file is in the repertory "parameters/", the forcings file is in the repertory "forcings/" and the outputs are stored in the repertory "outputs/". Templates can be found in those repertories. Note that the reading format of the parameter file is not the same in the single run or the multirun mode.

"IO_INTERFACE.txt" is read by the subroutine 'make_input_output' in the file "io_module.f90", it creates the netCDF output file plus a fortran scratch text file with the ID of the variables to record (unit=IOUT), and potentially two others scratch files (unit=IFORC and unit=IPARAM) recording the forcings and the parameters to be used in the main program.

The climate interpolation is done by the subroutines of the module "climate_module.f90". DynSoil integration is done by subroutines of the module "dynsoil_steady_state.f90". The inversion (in backward mode) is done in the main program.
