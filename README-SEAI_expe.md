# How to run the GEOCLIM-DynSoil-steady-state experiments

### Compilation

(from the GDSS root directory)

`cd sources/`
`make MODE=optim`
`cd ..`

### Execution

(from the GDSS root directory)

##### Choose which experiment to run:

From the GDSS root directory:
Select the configuration file ("IO\_INTERFACE") corresponding to the experiment

* ERA5 30min x 30min CONTROL run:                           `ln -s -f config_templates/IO_INTERFACE.ERA5-ctrl IO_INTERFACE`
* GFDL 30min x 30min Pre-Industrial CONTROL run:            `ln -s -f config_templates/IO_INTERFACE.GFDL-PIctrl IO_INTERFACE`
* GFDL "Park et al" Pre-Industrial CONTROL run:             `ln -s -f config_templates/IO_INTERFACE.GFDL-PIctrl IO_INTERFACE`
* CESM 0.9x1.25 slab ocean Pre-Indusctrial CONTROL run:     `ln -s -f config_templates/IO_INTERFACE.CESM-PIctrl IO_INTERFACE`
* GFDL 30min x 30min run with reduced SEAI "5Ma scenario":  `ln -s -f config_templates/IO_INTERFACE.GFDL-redSEAI-5Ma IO_INTERFACE`
* GFDL 30min x 30min run with reduced SEAI "10Ma scenario": `ln -s -f config_templates/IO_INTERFACE.GFDL-redSEAI-10Ma IO_INTERFACE`
* GFDL 30min x 30min run with reduced SEAI "15Ma scenario": `ln -s -f config_templates/IO_INTERFACE.GFDL-redSEAI-15Ma IO_INTERFACE`
* GFDL 30min x 30min run with removed SEAI:                 `ln -s -f config_templates/IO_INTERFACE.GFDL-noSEAI IO_INTERFACE`
* GFDL "Park et al" run with reduced SEAI "5Ma scenario":   `ln -s -f config_templates/IO_INTERFACE.GFDL-redSEAI-5Ma-Parketal IO_INTERFACE`
* GFDL "Park et al" run with reduced SEAI "10Ma scenario":  `ln -s -f config_templates/IO_INTERFACE.GFDL-redSEAI-10Ma-Parketal IO_INTERFACE`
* GFDL "Park et al" run with reduced SEAI "15Ma scenario":  `ln -s -f config_templates/IO_INTERFACE.GFDL-redSEAI-15Ma-Parketal IO_INTERFACE`
* GFDL "Park et al" run with removed SEAI:                  `ln -s -f config_templates/IO_INTERFACE.GFDL-noSEAI-Parketal IO_INTERFACE`
* ERA5 CONTROL run: `ln -s -f config_templates/ IO_INTERFACE`
* ERA5 CONTROL run: `ln -s -f config_templates/ IO_INTERFACE`
* ERA5 CONTROL run: `ln -s -f config_templates/ IO_INTERFACE`

##### Run the experiment:

(from the GDSS root directory)

`cd executables`
`./gdss 0 1 2`
`cd ..`


### Notes

The difference between 'GFDL 30min x 30min \*' and 'GFDL "Park et al" \*' is that the later use the land area file from
Park et al. 2020 for the CONTROL run, and the degassing flux from this CONTROL run. 
Whereas the former use the "new" land fraction file. This results in an ~1% difference in weathering flux (ie, degassing)
