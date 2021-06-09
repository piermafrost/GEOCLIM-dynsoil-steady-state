import netCDF4 as nc
import numpy as np

###########################################################################################
NLITH = 7
litho_class = ['water/ice', 'metamorphic', 'felsic', 'intermediate', 'mafic', 'carbonate', 'siliclastic sediments']
litho_levels = [[13,14], [12,15], [2,6,9], [7,10], [8,11], [4,5], [0,1,3]]
NOUTLITH = 5
OUTPUT_LITHO_LEV = [1,2,3,4,6]
###########################################################################################
# Open files:
fin = nc.Dataset('/home/piermafrost/data/litho_Hartmann2012/LiMW_30minx30min.nc', mode='r')
f = nc.Dataset('lithology_fraction_{}class_30minx30min.nc'.format(NOUTLITH), mode='w', format='NETCDF3_CLASSIC')
landfrac = fin['polygon_fraction'][:,:,:].sum(0)
FILLVALUE = -9.96921e+36
###########################################################################################


# Global attributes
f.title = 'Area fraction of {} lithological classes on a regular 30min x 30min grid'.format(NOUTLITH)
f.earth_model = 'WGS84'
f.source = 'Hartmann & Moodsorf, G3, 2013. doi:10.1029/2012GC004370'


# Copy dimensions, except lithology:
for dim in fin.dimensions.keys():
    if dim=='lithology':
        f.createDimension('lithology', NOUTLITH)
    else:
        f.createDimension(dim, fin.dimensions[dim].size)


# Copy variables:
for var in fin.variables.keys():
    if hasattr(fin[var], 'fill_value'):
        fill_value = fin[var].fill_value
    elif var=='polygon_fraction':
        fill_value = FILLVALUE
    else:
        fill_value = None

    v = f.createVariable(var, datatype=fin[var].datatype, dimensions=fin[var].dimensions, fill_value=fill_value)

    # Variables attributes:
    for attname in ['axis', 'long_name', 'standard_name', 'units']:
        if hasattr(fin[var], attname):
            setattr(v, attname, getattr(fin[var], attname))

    # Lithology variable attributes:
    if var=='lithology':
        for k,klith in enumerate(OUTPUT_LITHO_LEV):
            setattr(v, 'class_{}'.format(k+1), litho_class[klith])

    # Put variables
    if var=='lithology':
        v[:] = range(NOUTLITH)

    elif var=='polygon_fraction':
        for k,klith in enumerate(OUTPUT_LITHO_LEV):
            tmpv = np.ma.masked_invalid(fin[var][litho_levels[klith], :, :].sum(0)  /  (landfrac - fin[var][litho_levels[0], :, :].sum(0)))
            tmpv.set_fill_value(FILLVALUE)
            tmpv.data[tmpv.mask] = FILLVALUE
            v[k,:,:] = tmpv

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

