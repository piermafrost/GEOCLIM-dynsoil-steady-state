import numpy as np
import netCDF4 as nc

## First degassing estimates
#degrel = np.loadtxt('/home/piermafrost/Documents/GET/Ordovicien/degassing_relat_Chloe.txt', skiprows=1)
## New degassing estimates (June 2021)
degrel = np.loadtxt('/home/piermafrost/Documents/GET/Ordovicien/degassing_relat_Chloe_new.txt', skiprows=1)
age = np.int_(degrel[:,0])
fact = degrel[:,1:]
mean_fact = fact[:,0].mean()

CTRL_DATA = {'IPSL': '../output/Ordovician/gdss-output_IPSL-FOAM-SST_CTRL.nc',
             'FOAM-CPLD': '../output/Ordovician/gdss-output_FOAM_coupled_CTRL.nc',
             'FOAM-SLAB': '../output/Ordovician/gdss-output_FOAM_slab_CTRL.nc',
             'FOAM-OLD': '../output/Ordovician/gdss-output_FOAM_slab_old_CTRL.nc'}

#OUTROOT = {'IPSL': 'degassing_IPSL-FOAM-SST_',
#           'FOAM-CPLD': 'degassing_FOAM-coupled_',
#           'FOAM-SLAB': 'degassing_FOAM-slab_',
#           'FOAM-OLD': 'degassing_FOAM-slab_old_'}
OUTROOT = {'IPSL': 'degassing_new_IPSL-FOAM-SST_',
           'FOAM-CPLD': 'degassing_new_FOAM-coupled_',
           'FOAM-SLAB': 'degassing_new_FOAM-slab_',
           'FOAM-OLD': 'degassing_new_FOAM-slab_old_'}

def make(expe):

    fin = nc.Dataset(CTRL_DATA[expe])
    degass = fin['degassing'][:].data

    for k,a in enumerate(age):

        f = open(OUTROOT[expe]+'{:}Ma.txt'.format(a), mode='w')
        f.write('# CO2 degassing from solid Earth (mol/y)\n')
        #for a in [1, mean_fact] + list(fact[k,:]):
        for a in [mean_fact] + list(fact[k,:]): # do not use CTRL degassing for "new" degassing forcings
            f.write('{:}\n'.format(degass*a))

    f.close()
    fin.close()


make('IPSL')
make('FOAM-CPLD')
make('FOAM-SLAB')
make('FOAM-OLD')
