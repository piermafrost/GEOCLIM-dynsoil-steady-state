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
* CESM 0.9x1.25 slab ocean Pre-Indusctrial CONTROL run:     `ln -s -f config_templates/IO_INTERFACE.CESM-PIctrl IO_INTERFACE`
* GFDL 30min x 30min run with reduced SEAI "5Ma scenario":  `ln -s -f config_templates/IO_INTERFACE.GFDL-redSEAI-5Ma IO_INTERFACE`
* GFDL 30min x 30min run with reduced SEAI "10Ma scenario": `ln -s -f config_templates/IO_INTERFACE.GFDL-redSEAI-10Ma IO_INTERFACE`
* GFDL 30min x 30min run with reduced SEAI "15Ma scenario": `ln -s -f config_templates/IO_INTERFACE.GFDL-redSEAI-15Ma IO_INTERFACE`
* GFDL 30min x 30min run with removed SEAI:                 `ln -s -f config_templates/IO_INTERFACE.GFDL-noSEAI IO_INTERFACE`
* ERA5 CONTROL run: `ln -s -f config_templates/ IO_INTERFACE`
* ERA5 CONTROL run: `ln -s -f config_templates/ IO_INTERFACE`
* ERA5 CONTROL run: `ln -s -f config_templates/ IO_INTERFACE`

##### Run the experiment:

(from the GDSS root directory)

`cd executables`
`./gdss 0 1 2`
`cd ..`

