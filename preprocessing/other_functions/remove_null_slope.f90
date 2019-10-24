program remove_null_slope
use netcdf

double precision, dimension(360,180):: slope
integer:: ID, VID, ierr

ierr = nf90_open( 'slope1x1_nonzero.nc', NF90_WRITE , ID )
print *, nf90_strerror(ierr)
ierr = nf90_inq_varid( ID, 'slope', VID )
print *, nf90_strerror(ierr)
ierr = nf90_get_var( ID, VID, slope )
print *, nf90_strerror(ierr)
print *,  minval(minval( slope , 1 , slope>0 ),1)
where (slope==0) slope = minval(minval( slope , 1 , slope>0 ),1)
ierr = nf90_put_var( ID, VID, slope )
print *, nf90_strerror(ierr)
ierr = nf90_close(ID)
print *, nf90_strerror(ierr)

end program
