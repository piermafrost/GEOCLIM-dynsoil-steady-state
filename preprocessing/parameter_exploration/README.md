# PARAMETER EXPLORATION

This sub-program runs GEOCLIM with a fixed climatology several times, automatically changing parameters each time. This sub-program is useful when exploring the sensitivity of the model to changes in under-constrained parameters.

## Set up

The sub-program could be modified to explore a wider range of parameters, but in its current implementation it is most directly set-up to do the following:

* `./parameter_exploration.f90`
    * Line 20
        * Set `TEST_SINGLE_PARAMETERIZATION` to `.false.` if you would like to explore several parameter combinations. Otherwise, set to `.true.` to explore only a single parameter combination (the first one).
    * Line 26
        * Set `extend_output_file` to `.false.` to create a new output file, or `.true.` to extend an existing output file.
    * Line 27
        * Set the output file name and location.
    * Lines 31-38
        * Select the files that store the boundary conditions here.
    * Lines 83-92
        * Set values for the fixed parameters here (i.e. those parameters that are not varied from run to run).
        * If you would like to vary these parameters too, the simplest way to do so would be to run the sub-program by varying the default parameters (as described below), then manually changing one or more of these fixed parameters, then rerunning the sub-program.
* `./parameter_list/`
    * This directory stores the text files that list the values of the various parameters to be explored.
    * Simply list the values of the parameters that you wish to test in these text files - all permutations of the values within these text files will be explored.
