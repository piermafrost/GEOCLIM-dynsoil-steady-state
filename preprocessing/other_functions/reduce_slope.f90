program reduce_slope

use netcdf
implicit none


  character(len=*), parameter:: filein='/Users/yves/fortran/GEOCLIM3/calibration/365Ma/Slope365Ma.nc'
  character(len=*), parameter:: fileout='/Users/yves/fortran/GEOCLIM3/calibration/365Ma/Slope365Ma_allshifted.nc'
  integer, parameter:: nx=128, ny=128 !nx=360 , ny=180

  real:: missval, Smin, Smax
  real:: lon(nx),lat(ny)
  real, dimension(nx,ny):: slope

  integer:: fileid, varid, ierr, Odimids(2), Ovarids(4)

! +++++++++++++++++++++++++++++ !
  real, parameter:: factor=0.5
! +++++++++++++++++++++++++++++ !

!=========================================================================!

! File openning:
  ierr = nf90_open(filein, NF90_NOWRITE, fileid)
  call isitok(ierr,'Error while openning input file')

! Loading variable
  ierr = nf90_inq_varid( fileid, 'lon', varid )
  call isitok(ierr,'Error while getting "lon" variable identifiant')
  ierr = nf90_get_var( fileid, varid, lon )
  call isitok(ierr,'Error while getting variable "lon"')
  ierr = nf90_inq_varid( fileid, 'lat', varid )
  call isitok(ierr,'Error while getting "lat" variable identifiant')
  ierr = nf90_get_var( fileid, varid, lat )
  call isitok(ierr,'Error while getting variable "lat"')
  ierr = nf90_inq_varid( fileid, 'slope', varid )
  call isitok(ierr,'Error while getting "slope" variable identifiant')
  ierr = nf90_get_var( fileid, varid, slope )
  call isitok(ierr,'Error while getting variable "slope"')

! Loading attributes "missing-value"
  missval = 9.96921e+36
  ierr = nf90_get_att( fileid, varid, '_FillValue', missval)
  call isitok(ierr,'Error while getting variable attribute "Missing-value"')
  if (ierr/=0) missval = 9.96921e+36

! File closing:
  ierr = nf90_close(fileid)
  call isitok(ierr,'Error while closing input file')

!=========================================================================!

  Smin = minval(minval(slope,1,(slope/=missval.and.slope>0)),1)
  Smax = maxval(maxval(slope,1,(slope/=missval)),1)

  where (slope==0) slope = Smin

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  where (slope/=missval)
!    slope = slope * factor**( (log(slope)-log(Smin))/(log(Smax)-log(Smin)) ) ! modif high values and not low values
!     slope = slope * factor**( (log(Smax)-log(slope))/(log(Smax)-log(Smin)) ) ! modif low values and not high values
     slope = slope * factor ! modif all
  end where
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!=========================================================================!

! Output file creation:
  ierr = nf90_create(fileout, NF90_CLOBBER, fileid)
  call isitok(ierr,'Error while output file creation')

! Dimensions:
  ierr = nf90_def_dim( fileid, 'lon' , nx, Odimids(1) )
  call isitok(ierr,'Error while defining dimension "lon"')
  ierr = nf90_def_dim( fileid, 'lat' , ny, Odimids(2) )
  call isitok(ierr,'Error while defining dimension "lat"')

! Variables:
  ierr = nf90_def_var(fileid, 'lon', NF90_FLOAT, Odimids(1), Ovarids(1))
  call isitok(ierr,'Error while defining variable "lon"')
  ierr = nf90_def_var(fileid, 'lat', NF90_FLOAT, Odimids(2), Ovarids(2))
  call isitok(ierr,'Error while defining variable "lat"')
  ierr = nf90_def_var(fileid, 'slope', NF90_FLOAT, Odimids, Ovarids(3))
  call isitok(ierr,'Error while defining variable "slope"')

! Variables attributes:
  ! Names:
  ierr = nf90_put_att(fileid, Ovarids(1) , 'name', 'longitude')
  call isitok(ierr,'Error while putting variable "lon" attribute "name"')
  ierr = nf90_put_att(fileid, Ovarids(2) , 'name', 'latitude')
  call isitok(ierr,'Error while putting variable "lat" attribute "name"')
  ierr = nf90_put_att(fileid, Ovarids(3) , 'name', 'cell-mean slope')
  call isitok(ierr,'Error while putting variable "slope" attribute "name"')
  ! Units:
  ierr = nf90_put_att(fileid, Ovarids(1) , 'units', 'degrees_east')
  call isitok(ierr,'Error while putting variable "lon" attribute "units" ')
  ierr = nf90_put_att(fileid, Ovarids(2) , 'units', 'degrees_north')
  call isitok(ierr,'Error while putting variable "lat" attribute "units" ')
  ierr = nf90_put_att(fileid, Ovarids(3) , 'units', 'm/m')
  call isitok(ierr,'Error while putting variable "slope" attribute "units" ')
  ! Missing-values:
  ierr = nf90_put_att(fileid, Ovarids(3) , '_FillValue', missval)
  call isitok(ierr,'Error while putting variable "reg_thck" attribute "Missing-value"')

  ierr = nf90_enddef(fileid)
  call isitok(ierr,'Error while end of definition')

! Writting:
  ierr = nf90_put_var(fileid, Ovarids(1) , lon )
  call isitok(ierr,'Error while putting variable "lon"')
  ierr = nf90_put_var(fileid, Ovarids(2) , lat )
  call isitok(ierr,'Error wile putting variable "lat"')
  ierr = nf90_put_var(fileid, Ovarids(3) , slope )
  call isitok(ierr,'Error wile putting variable "slope"')

  ! Output file closing
  ierr = nf90_close(fileid)
  call isitok(ierr,'Error while output file closing')

end program reduce_slope


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine isitok(ierr,message)
use netcdf
implicit none
  integer, intent(in):: ierr
  character(len=*), intent(in):: message
  if (ierr/=0) then
    print *, message
    print *, nf90_strerror(ierr)
  end if
end subroutine

