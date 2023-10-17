#!/bin/bash

# Replicate config files "IO_INTERFACE' into same config file using single-year climatology, everything else constant
# i.e., replace pattern "y41-70" by single year from 41 to 70 (e.g., "y66"), and link to the corresponding files in input/GCM_ANN_clim/yby/

# CALL EXEMPLE:
#   `./generate_yby_config_files.sh IO_INTERFACE`
#   `./generate_yby_config_files.sh IO_INTERFACE.E1850C5*`


for fin in "$@" # -> loop on input arguments
do
    for y in $(seq 41 70) # -> loop on years
    do
        # duplicate config file and change climatology year
	fout=${fin/y41-70/y$y}
	cp $fin $fout
	sed -i -e "s/y41-70/y$y/" $fout
	sed -i -e "s/GCM_ANN_clim\//GCM_ANN_clim\/yby\/y-$y\//" $fout
    done
done

