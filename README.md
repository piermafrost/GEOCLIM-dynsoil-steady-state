# GEOCLIM-dynsoil-steady-state

## Foreword: specificity of the current Github branch
The branch "SEAI" (standing for "South-East Asian Island") corresponds to the code as it was used for the study "The Role of Southeast Asian Island Topography on Indo-Pacific Climate and Silicate Weathering " (Chiang et al., submitted to Paleoceanography & Paleoclimatology).

The version of the Fortran source code is the same as the version 2.0 (master branch of January 3rd 2023).
The only difference is the added option "fill_missing=.true." line 412 of "source/io_module.f90" to allow for inputs of global temperature field without missing-value.

#### Summary: how to reproduce the GEOCLIM simulations from the study
1. Download the current repository (current branch)
2. Make sure you have installed a Fortran compiler and a netCDF-Fortran library associated.
3. Download the input files from the Dryad dataset ???.
   It contains the CESM individual year climatology, that cannot be stored on the GitHub repository for reason of memory space.
   Put the archive file "ann_climo.tar.gz" in "input/GCM_ANN_clim/", and expand it with `tar xfz ann_climo.tar.gz`
4. The only configuration files currently stored in "config_templates/" are for the GEOCLIM simulations using the y41-70 climatology as input.
   To generate the configuration files for each individual 1-year climatology, go to "config_templates/" and launch: `./generate_yby_config_files.sh IO_INTERFACE.E1850C5_y41-70_*`.
   This will create 360 more configuration files.
5. Compile the code (in "sources/") with `make MODE=optim` (see section "compilation tips").
6. Run the successive 372 simulations with the following command (linking-and-running loop) **from the root directory**:
```
for f in config_templates/IO_INTERFACE.E1850C5_*
do
  ln -s -f $f IO_INTERFACE
  cd executables/
  ./gdss 0 1 2
  cd ../
done
```
This generates the 1 output netCDF file per simulations in "outputs/", using ~11 GB of total disk space.

#### Notes:
* About the size of the ouputs:  
The simulations "PIctrl", "noSEAI_1xCO2" and "noSEAItopo_flatSEAIslope_1xCO2", for each individual year of the y41-70 time-series, will output the 2D weathering field for each 573 parameter combinations.
If you do not want to output the weathering fields, write ".false." (instead of ".true.") for the output variable "weathering_lithmean", at the last line of "IO_INTERFACE.E1850C5_y41-70_PIctrl", "IO_INTERFACE.E1850C5_y41-70_PIctrl" and "IO_INTERFACE.E1850C5_y41-70_PIctrl", in "config_templates/" *before generating the year-by-year config files with "./generate_yby_config_files.sh"*.
The outputs would then use only ~32 MB of disk space.
Note that only the average weathering field over the 573 parameterizations were presented in Chiang et al.

* The variable "degassing" is the area-integral of the weathering field (i.e., *Global weathering rate*, as reported in Chiang et al.).

* Not all of the simulations of this repository were presented in Chiang et al.
Only "PIctrl", "noSEAItopo_flatSEAIslope_1xCO2", "noSEAI_1xCO2", "noSEAItopo_flatSEAIslope_equil" and "noSEAI_equil" (for each individual year of the y41-70 time-series) were actually used.


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
