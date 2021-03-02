import numpy as np
import netCDF4 as nc

degrel = np.loadtxt('/home/piermafrost/Documents/GET/Ordovicien/degassing_relat_Chloe.txt', skiprows=1)
age = np.int_(degrel[:,0])
fact = degrel[:,1:]

CTRL_DATA = {'IPSL': '../output/Ordovician/gdss-output_IPSL-FOAM-SST_CTRL.nc',
             'FOAM-CPLD': '../output/Ordovician/',
             'FOAM=SLAB': '../output/Ordovician/'}

OUTROOT = {'IPSL': 'degassing_IPSL-FOAM-SST_',
           'FOAM-CPLD': 'degassing_FOAM-coupled_',
           'FOAM-SLAB': 'degassing_FOAM-slab_'}

def make(expe):

    fin = nc.Dataset(CTRL_DATA[expe])
    degass = fin['degassing'][:].data

    for k,a in enumerate(age):

        f = open(OUTROOT[expe]+'{:}Ma.txt'.format(a), mode='w')
        f.write('# CO2 degassing from solid Earth (mol/y)\n')
        for a in [1] + list(fact[k,:]):
            f.write('{:}\n'.format(degass*a))

    f.close()
    fin.close()

make('IPSL')