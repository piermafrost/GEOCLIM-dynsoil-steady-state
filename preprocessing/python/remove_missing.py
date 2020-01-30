"""
This python script read the forcing files temperature, runoff, slope and lithology
(files defined in the "FORCINGS VARIABLES" section) and set area = 0 where missing
values are found (list of area files to modify defined in section "AREA VARIABLES")
It uses the packages "netCDF4" and "numpy"
"""

import netCDF4 as nc
import numpy as np


######################
##  AREA VARIABLES  ##
######################
##  Enter here all the netCDF files containing the variable "area" that you want to modify.

farea = list()
farea.append( nc.Dataset('../../input/land_area_noIA_match_clim_slope_lith.nc', mode='r+') )
farea.append( nc.Dataset('../../input/land_area_redIA10_match_clim_slope_lith.nc', mode='r+') )
farea.append( nc.Dataset('../../input/land_area_redIA5_match_clim_slope_lith.nc', mode='r+') )
farea.append( nc.Dataset('../../input/land_area_redIA_match_clim_slope_lith.nc', mode='r+') )

AREA_VARNAME = 'area'

#######################


#########################
##  FORCING VARIABLES  ##
#########################
##  Enter here all the netCDF files containing the "forcing variables" whose missing-values
##  should correspond to area = 0
##  There are currently 4 "forcing" variables: temperature, runoff, slope and lithology
##    Note: Lithology must have a shape ( nx , ny , (nlith+1) ), the first lithology class
##          standing for 'water and ice body', and if the fraction of this lithology is 1,
##          the point if considered as a missing point.

f0 = list()
#
# Temperature:
f0.append( nc.Dataset('../../input/GFDL_temp_360_720_PD.nc') )
TEMP_VARNAME = 'tmp'
#
# Runoff
f0.append( nc.Dataset('../../input/GFDL_runoff_360_720_PD.nc') )
RUNF_VARNAME = 'rnf'
#
# Slope
f0.append( nc.Dataset('../../input/slope_360_720_PD.nc') )
SLOP_VARNAME = 'slope'
#
# Lithology fraction
f0.append( nc.Dataset('../../input/lith_mask_360_720_PD.nc') )
LITH_VARNAME = 'frac'

#########################




# Define missing-value mask ('True' if missingpoint, 'False' elsewhere)

mask = np.logical_or( f0[0][TEMP_VARNAME][:,:,:].mask.any(0), f0[1][RUNF_VARNAME][:,:,:].mask.any(0) )
mask = np.logical_or( mask, f0[2][SLOP_VARNAME][:,:].mask )
mask = np.logical_or( mask, f0[3][LITH_VARNAME][0,:,:]>=(1-1e-6))


# Close forcing files

for f in f0:
    f.close()


#++++++++++++++++++++++++++++++++#
# Set area = 0 in missing-points #
#++++++++++++++++++++++++++++++++#

for f in farea:
    a = f[AREA_VARNAME][:,:]
    a[mask] = 0
    f['area'][:,:] = a
    f.close()
