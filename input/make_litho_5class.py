import netCDF4 as nc
import numpy as np

###########################################################################################
NLITH = 5
litho_class = ['metamorphics', 'felsics', 'intermediates', 'mafics', 'siliclastic sediments', 'carbonates']
litho_levels = [[12,15], [2,6,9], [7,10], [8,11], [0,1,3], [4,5]]
###########################################################################################
# Open files:
fin = nc.Dataset('/home/piermafrost/data/litho_Hartmann2012/LiMW_30minx30min.nc', mode='r')
f = nc.Dataset('lithology_fraction_PI_{}class_30minx30min.nc'.format(NLITH), mode='w', format='NETCDF3_CLASSIC')
###########################################################################################


# Global attributes
f.title = 'Area fraction of {} lithological classes on a regular 30min x 30min grid'.format(NLITH)
f.earth_model = 'WGS84'
f.source = 'Hartmann & Moodsorf, G3, 2013. doi:10.1029/2012GC004370'


# Copy dimensions, except lithology:
for dim in fin.dimensions.keys():
    if dim=='lithology':
        f.createDimension('lithology', NLITH)
    else:
        f.createDimension(dim, fin.dimensions[dim].size)


# Copy variables:
for var, newvar in zip(['longitude', 'latitude', 'lithology', 'polygon_fraction'], ['longitude', 'latitude', 'lithology', 'lithfrac']):
    if hasattr(fin[var], 'fill_value'):
        fill_value = fin[var].fill_value
    else:
        if var=='polygon_fraction':
            fill_value = -1e36
        else:
            fill_value = None

    v = f.createVariable(newvar, datatype=fin[var].datatype, dimensions=fin[var].dimensions, fill_value=fill_value)

    # Variables attributes:
    for attname in ['axis', 'long_name', 'standard_name', 'units']:
        if hasattr(fin[var], attname):
            setattr(v, attname, getattr(fin[var], attname))

    # Lithology variable attributes:
    if var=='lithology':
        for k in range(NLITH):
            setattr(v, 'class_{}'.format(k+1), litho_class[k])

    # Put variables
    if var=='lithology':
        v[:] = range(NLITH)

    elif var=='polygon_fraction':
        lithfrac = np.ma.zeros(v.shape, dtype=v.dtype)
        for k in range(NLITH):
            lithfrac[k,:,:] = fin[var][litho_levels[k], :, :].sum(0)

        sumlithfrac = lithfrac.sum(0) + fin[var][litho_levels[-1], :, :].sum(0)
        # Add carbonates in total sum ^^^^^^^^^^^****************^^^^^^^^^^^^^^
        mask = (sumlithfrac == 0)
        lithfrac[:,mask] = fill_value
        lithfrac[:,mask].mask = True
        lithfrac /= sumlithfrac

        v[:,:,:] = lithfrac

    else:
        if v.ndim == 0:
            v.assignValue(fin[var].getValue())
        elif v.ndim == 1:
            v[:] = fin[var][:]
        elif v.ndim == 2:
            v[:,:] = fin[var][:,:]
        elif v.ndim == 3:
            v[:,:,:] = fin[var][:,:,:]
        elif v.ndim == 4:
            v[:,:,:,:] = fin[var][:,:,:,:]
        elif v.ndim == 5:
            v[:,:,:,:,:] = fin[var][:,:,:,:,:]
        elif v.ndim == 6:
            v[:,:,:,:,:,:] = fin[var][:,:,:,:,:,:]
        elif v.ndim == 7:
            v[:,:,:,:,:,:,:] = fin[var][:,:,:,:,:,:,:]
        else:
            raise ValueError('Too many dimensions')


# Close files
fin.close()
f.close()

