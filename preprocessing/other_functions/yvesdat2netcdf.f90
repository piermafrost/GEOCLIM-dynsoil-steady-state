module check_m
implicit none
contains
subroutine check(ierr,message)
use netcdf
integer, intent(in):: ierr
character(len=*), intent(in):: message
if (ierr/=NF90_NOERR) then
  print *, message//' '//nf90_strerror(ierr)
  call abort
end if
end subroutine
end module

program yvesdat2netcdf
use netcdf
use check_m, only: check
implicit none

character(len=*), parameter:: ifname= '../calibration/present_day/runoff1x1.out' !'pangea2014/Aire_200MaVegDyn.dat' !Aire_tot.dat'
character(len=*), parameter:: ofname= '../calibration/present_day/runoff1x1.nc' !'slope_200Ma_test_FOAM-48x40.nc'
character(len=*), parameter:: gridfname = '-' !'/Users/piermafrost/Documents/data/GCM/FOAM-48x40_grid.nc' ! '-' for no file
character(len=*), parameter:: varname='runoff'
character(len=*), parameter:: varunits='m/y'
character(len=*), parameter:: varfillvalname='_FillValue'
double precision, parameter:: ifillval = -50!0.d+0 !-9999.9999d+36
double precision, parameter:: ofillval = 9.96921d+36
logical, parameter:: new_output_file = .true.

integer, parameter:: nlon=360, nlat=180

double precision, dimension(nlon):: lon
double precision, dimension(nlat):: lat
double precision, dimension(nlon,nlat):: var
double precision, dimension(nlon*nlat):: varline
double precision:: CO2lev

integer:: i,j,k, nco2, ierr, fid, dimid(2), varid(3)


! read input text file:
open(unit=1,file=ifname,status='old',action='read')

! read in lines (one line per CO2 level)
!nco2 = 0
!CO2lev = 0.d+0
!do while (CO2lev/=280.0d+0)
!  read(unit=1,fmt=*) CO2lev, varline
!end do
!! ravelling:
!k = 0
!do j = 1,nlat
!  do i = 1,nlon
!    k = k+1
!    if (varline(k) == ifillval) then
!      var(i,j) = ofillval
!    else
!      var(i,j) = varline(k)
!    end if
!  end do
!end do

! read in column
do j = 1,nlat
  do i = 1,nlon
    read(unit=1,fmt=*) var(i,j)
    if (var(i,j) == ifillval) then
      var(i,j) = ofillval
    else ! CONVERSION:
      var(i,j) = 1e-3*var(i,j) !0.005
    end if
  end do
end do

close(unit=1)

! longitude and latitude
if (gridfname=='-') then
  lon = (/ ( (dble(i)-0.5)*360/dble(nlon) - 180 ,i=1,nlon) /)
  lat = (/ ( (dble(j)-0.5)*180/dble(nlat) - 90  ,j=1,nlat) /)
else
  ierr = nf90_open( gridfname , NF90_NOWRITE , fid ) ;  call check(ierr,'error while openning file '//gridfname)
  ierr = nf90_inq_varid( fid , 'lon' , varid(1) ) ;  call check(ierr,'error while inquiring variable "lon" ID')
  ierr = nf90_inq_varid( fid , 'lat' , varid(2) ) ;  call check(ierr,'error while inquiring variable "lat" ID')
  ierr = nf90_get_var( fid , varid(1) , lon ) ;  call check(ierr,'error while getting variable "lon"')
  ierr = nf90_get_var( fid , varid(2) , lat ) ;  call check(ierr,'error while getting variable "lat"')
  ierr = nf90_close( fid )
end if

! write output netcdf file
if (new_output_file) then

  ierr = nf90_create( ofname , NF90_CLOBBER , fid )  ;  call check(ierr,'error while creating file "'//ofname//'":')

  ierr = nf90_def_dim( fid , 'lon' , nlon , dimid(1) )  ;  call check(ierr,'error while defining dimension "lon":') 
  ierr = nf90_def_dim( fid , 'lat' , nlat , dimid(2) )  ;  call check(ierr,'error while defining dimension "lat":') 

  ierr = nf90_def_var( fid , 'lon' , NF90_FLOAT , dimid(1) , varid(1) )
  call check(ierr,'error while defining variable "lon":') 
  ierr = nf90_def_var( fid , 'lat' , NF90_FLOAT , dimid(2) , varid(2) )
  call check(ierr,'error while defining variable "lat":') 

  ierr = nf90_put_att( fid , varid(1) , 'axis' , 'X' )
  ierr = nf90_put_att( fid , varid(1) , 'name' , 'longitude' )
  ierr = nf90_put_att( fid , varid(1) , 'units' , 'degrees_east' )
  ierr = nf90_put_att( fid , varid(2) , 'axis' , 'Y' )
  ierr = nf90_put_att( fid , varid(2) , 'name' , 'latitude' )
  ierr = nf90_put_att( fid , varid(2) , 'units' , 'degrees_north' )

  ierr = nf90_enddef( fid )

  ierr = nf90_put_var( fid , varid(1) , lon )  ;  call check(ierr,'error while putting variable "lon":')
  ierr = nf90_put_var( fid , varid(2) , lat )  ;  call check(ierr,'error while putting variable "lat":')

  ierr = nf90_redef( fid )

else

  ierr = nf90_open( ofname , NF90_WRITE , fid )  ;  call check(ierr,'error while openning file "'//ofname//'":')
  ierr = nf90_redef( fid )
  ierr = nf90_inq_dimid( fid , 'lon' , dimid(1) )  ;  call check(ierr,'error while inquiring dimension "lon" ID:') 
  ierr = nf90_inq_dimid( fid , 'lat' , dimid(2) )  ;  call check(ierr,'error while inquiring dimension "lat" ID:') 

end if


ierr = nf90_def_var( fid , varname , NF90_FLOAT , dimid , varid(3) )
call check(ierr,'error while defining variable "'//varname//'":') 

ierr = nf90_put_att( fid , varid(3) , 'name' , varname )
call check(ierr,'error while putting variable "'//varname//'" attribute "name" ('//varname//'):') 
ierr = nf90_put_att( fid , varid(3) , 'units' , varunits )
call check(ierr,'error while putting variable "'//varname//'" attribute "units" ('//varunits//'):') 
ierr = nf90_put_att( fid , varid(3) , varfillvalname , real(ofillval) )
call check(ierr,'error while putting variable "'//varname//'" attribute '//varfillvalname//':') 

ierr = nf90_enddef( fid )

ierr = nf90_put_var( fid , varid(3) , var )  ;  call check(ierr,'error while putting variable "'//varname//'":')

ierr = nf90_close( fid )


end program
