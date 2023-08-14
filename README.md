# GEOCLIM-dynsoil-steady-state

## Foreword: specificity of the current Github branch
The branch "NG" (standing for "New Guinea") corresponds to the code and the input data as it was used for the study "The rise of New Guinea and the fall of Neogene global temperatures" (Martin et al., PNAS, 2023).

The version of the Fortran source code is the same as the version 2.0 (master branch of January 3rd 2023).

#### Quick summary: how to reproduce the GEOCLIM simulations from the study
1. Download the current repository (current branch)
2. Make sure you have installed a Fortran compiler and a netCDF-Fortran library associated.
3. Compile the code (in "sources/") with `make MODE=optim` (see section "compilation tips").
4. Execute the following linking and running loop from the root directory:
```
for f in config_templates/IO_INTERFACE.NG_*
do
  ln -s -f $f IO_INTERFACE
  cd executables/
  ./gdss 0 1 2
  cd ../
done
```
This generates the 19 outputs files, using ~900 Mo of disk space.

Note: to generate lighter outputs, one can choose not to output any geographically-distributed variables (weathering rates, lithology...).
The variable "degassing" is the area-integral of the weathering field.
To do so, put ".false." for those variables in the last block of the IO\_INTERFACE files.

## Update from version 1.1 (previous release)
The version 2 is a major update from version 1, with a large restructuring of the architecture and of the configuration and input files.
This version is still under development.

## Presentation
GEOCLIM - DynSoil - Steady-State computes geographically-distributed chemical weathering rates (along with associated parameters) at steady-state, according to the given climatology (temperature and runoff) for several CO<sub>2</sub> levels (2 at least!), the topographic slope and the lithology fraction in each grid cell. The climatology is typically taken from the output of a general circulation climate model. The total weathering flux is assumed to equal the CO<sub>2</sub> degassing by the solid Earth.

The DynSoil component of the model was developed by Pierre Maffre during his PhD research with Yves Goddéris and is detailed in his thesis: Interactions entre tectonique, érosion, altération des roches silicatées et climat à l'échelle des temps géologiques: rôle des chaînes de montagnes. Océan, Atmosphère. Université Toulouse III Paul Sabatier, 2018. Français. https://tel.archives-ouvertes.fr/tel-02059359
The overall framework of integrating a silicate weathering model with varying climatology to model ancient pCO<sub>2</sub> levels is that of the GEOCLIM suite of models developed by Yves Goddéris: https://geoclimmodel.wordpress.com.

The model has two modes:
	1) Forward mode: impose a CO<sub>2</sub> level to the model, it gets the corresponding climate field and computes weathering.
	2) Backward mode: impose a volcanic degassing to the model, it finds by inversion the CO<sub>2</sub> level for which the total CO<sub>2</sub> consumption equals the degassing. The precision level can be set in the main program "gdss_mainprog.f90".

The continental climate (i.e. temperature and runoff) is interpolated linearly between the different CO<sub>2</sub> levels. At this time, the algorithm used for the interpolation does not allow extrapolation beyond the lowest and highest CO<sub>2</sub> level of any of the climate model runs.
The climate interpolation is done by the subroutines of the module "climate_module.f90". DynSoil integration is done by subroutines of the module "dynsoil_steady_state_module.f90". The inversion (in backward mode) is done in the main program.

In addition, the model can be run either in single run mode, for which it computes weathering one time for one given set of parameters, or in multiple runs mode, for which an netCDF unlimited dimension is created, and the model is run as many times as the number of a given set of parameters.

## Repository architecture
The input/output information is stated in the file "IO_INTERFACE.txt", as well as the flags for model mode (Forward, Backward, Single run, Multirun). This file is in the root directory.
It includes the names of input files to use of continent area, lithology, temperature, runoff, slope. The model checks the consistency of this input dataset (size, coordinates, units) and DynSoil is run accordingly using these geographic-climatic-lithologic settings.
These input files should be in netCDF format. The names of the dimensions for the input variables in these netCDF files also need to be stated in "IO_INTERFACE.txt" in the INPUT CONDITION section.

The executable files are meant to be run in the directory "executables/". For this reason, the program looks for the relative path "../IO_INTERFACE.txt", as specified on line 11 of the main program "gdss_mainprog.f90". Hence, the code will crash if it is not run directory one level after the root directory. All other paths can be changed in the file "IO_INTERFACE.txt" without recompiling the code.

The physical parameters of DynSoil and the forcings (CO<sub>2</sub> or degassing) are read in two separated text files whose names are stated in "IO_INTERFACE.txt". See the documentation in "IO_INTERFACE.txt" for the format of these files.

Conventionally, the parameter file is in the folder "parameters/", the forcings file is in the folder "forcings/" and the outputs are stored in the folder "outputs/". Templates can be found in those folders. Note that the reading format of the parameter file is not the same in the single run or the multirun mode.

"IO_INTERFACE.txt" is read by the subroutine 'make_input_output' in the file "io_module.f90", it creates the netCDF output file (conventionnaly in the "output/" directory) plus a fortran scratch text file with the ID of the variables to record (unit=IOUT), and potentially two others scratch files (unit=IFORC and unit=IPARAM) recording the forcings and the parameters to be used in the main program.

## Compilance checks
The program conducts a certain number of compliance checks in the subroutine "make_input_output" (file "io_module.f90"), that includes checking the forcing variables units (area, temperature, runoff, slope...), the consistency of the geographic axes of all forcing variables and the consistency of their continental mask. For this last check only, the program is not automatically aborted in case of failure, but the user is offered the following options: 0: abort the program, 1: remove the erratic cells (set their land area to 0), 2: remove only the erratics cells having missing-value (ie: ignore lithological fraction issue), 3: ignore the check and run the program as it is (not recommended). This argument can be enter interactively, or the user can enter it as an argument of the executable: "./executable 0", "./executable 1"...

## Compilation tips:
This program uses syntax allowed in fortran 2003 or later (subroutines with allocatable dummy arguments, and the use of the subroutine "get_command_argument"). Make sure your compiler supports this syntax (with gfortran: "-std=f2003").
A Makefile is also available in the "sources" folder.

## Parameter exploration:
This repository also contains a version of the model where the user can vary the parameters for a unique given CO<sub>2</sub> level. It can be found in the folder "preprocessing/parameter_exploration/". See the README in this folder for more information.

## Run example and experiments
The model is currently configure for using the input files present in the folder "input/" (see file "IO_INTERFACE.txt"). It corresponds to an example of boundary conditions consisting of present-day continent postion, topography, lithology and climate (at 3 CO2 levels: 286ppm, 572ppm and 1144ppm) excluding all the ice-sheets. The file prescribing the land area of each cell ("input/land_area_360_720_PD_match_clim_slope_lith.nc") is set to exactly match the cells having the full data (temperature, runoff, slope and lithology).
This example configuration is set in multirun, backward mode. It uses a subset of 10 parameter combinations (in the file "parameters/test_params_REDUCED.csv") of the ensemble of combinations used in the study "Emergence of the Southeast Asian islands as a driver for Neogene cooling" by Park et al. The full list of parameter combinations is available in "parameters/test_params.csv". The degassing forcing file corresponding to these parameter combinations (and used in the current configuration) is "forcings/degassing_test_params_REDUCED.txt" ("forcings/degassing_test_params.txt" for the full ensemble).
With this configuration, *all the runs should yield an equilibrium CO2 level of 286ppm*.

The rest of the input files corresponding to all the scenarios tested in the study of Park et al. can be found in the repository https://github.com/Swanson-Hysell-Group/GEOCLIM_Modern, as well as the code for generating the inputs and analyzing the outputs.
