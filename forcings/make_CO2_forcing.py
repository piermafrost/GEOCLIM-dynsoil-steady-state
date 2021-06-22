'''
Generate 1 forcing file of CO2 for each experimental set-up (i.e., each GCM) and
time slice, corresponding the mean tropical SST from Chloé (Song et al., or
Goldberg et al.)
'''

import sys
sys.path.append('/home/piermafrost/Documents/GET/Ordovicien/results_GDSS')
from sst_script import inverse_SST
import numpy as np

# Mean tropical SST estimates
sst = np.loadtxt('/home/piermafrost/Documents/GET/Ordovicien/SST-tropical_Ordovicien_Song_Goldberg.txt',
                 delimiter='\t', skiprows=1)
age = np.int_(sst[:,0])
sst = sst[:,1:]

OUTROOT = {'lmdz': 'CO2-from-trop-SST_IPSL-FOAM-SST_',
           'coupled': 'CO2-from-trop-SST_FOAM-coupled_',
           'slab': 'CO2-from-trop-SST_FOAM-slab_',
           'slab_old': 'CO2-from-trop-SST_FOAM-slab_old_'}

def make(expe):

    for k,a in enumerate(age):

        f = open(OUTROOT[expe]+'{:}Ma.txt'.format(a), mode='w')
        f.write('# atmospheric CO2 (ppm) yielding the expected mean tropical SST\n')
        f.write('#  (#1 Goldbert et al. 2021, #2 Song et al. 2019)\n')
        for t in sst[k,:]:
            co2 = inverse_SST(expe, str(a), t, interp='log', lat_range=(-21,21))
            f.write('{:}\n'.format(co2))

        f.close()


make('lmdz')
make('coupled')
make('slab')
make('slab_old')
