program create_flat_slope
use netcdf
implicit none

  character(len=*), parameter::filein='/Users/yves/fortran/GEOCLIM4_thea/preproc/dynsoil/dynsoil_48x40x10_nullreg.nc'
  character(len=*), parameter::fileout='/slope_200Ma_48x40.nc'
  integer, parameter:: nx=48, ny=40
! =======================================
  real, parameter:: NOISE_AMPLITUDE = 5.
  real, parameter:: MEANSLOPE = 0.03354 & ! 0.001 0.002 0.005 0.01234 0.02 0.03354
                               *  2.*log(NOISE_AMPLITUDE)/(NOISE_AMPLITUDE-1./NOISE_AMPLITUDE)
! =======================================

  integer:: i,j, ierr, ifid,ofid,idimid(2),odimid(2),ivarid(3),ovarid(3)
  real:: lon(nx), lat(ny)
  real, dimension(nx,ny):: Slp0, Slp
  real:: noise, missval
  
  character(len=*), parameter:: missvalname='_FillValue'


  ierr = nf90_open(   filein  , NF90_NOWRITE , ifid )  ;  call check(ierr)
  ierr = nf90_create( fileout , NF90_CLOBBER , ofid )  ;  call check(ierr)

  ierr = nf90_inq_varid( ifid , 'lon'   , ivarid(1) )  ;  call check(ierr)
  ierr = nf90_inq_varid( ifid , 'lat'   , ivarid(2) )  ;  call check(ierr)
  ierr = nf90_inq_varid( ifid , 'slope' , ivarid(3) )  ;  call check(ierr)

  ierr = nf90_def_dim( ofid , 'lon' , nx , odimid(1) )  ;  call check(ierr)
  ierr = nf90_def_dim( ofid , 'lat' , ny , odimid(2) )  ;  call check(ierr)

  ierr = nf90_def_var( ofid , 'lon' , NF90_REAL , odimid(1) , ovarid(1) )  ;  call check(ierr)
  ierr = nf90_copy_att( ifid , ivarid(1) , 'name'  , ofid , ovarid(1) )    ;  call check(ierr)
  ierr = nf90_copy_att( ifid , ivarid(1) , 'units' , ofid , ovarid(1) )    ;  call check(ierr)
  ierr = nf90_def_var( ofid , 'lat' , NF90_REAL , odimid(2) , ovarid(2) )  ;  call check(ierr)
  ierr = nf90_copy_att( ifid , ivarid(2) , 'name'  , ofid , ovarid(2) )    ;  call check(ierr)
  ierr = nf90_copy_att( ifid , ivarid(2) , 'units' , ofid , ovarid(2) )    ;  call check(ierr)
  ierr = nf90_def_var( ofid , 'slope' , NF90_REAL , odimid  , ovarid(3) )  ;  call check(ierr)
  ierr = nf90_copy_att( ifid , ivarid(3) , 'name'  , ofid , ovarid(3) )    ;  call check(ierr)
  !ierr = nf90_copy_att( ifid , ivarid(3) , 'units' , ofid , ovarid(3) )    ;  call check(ierr)
  ierr = nf90_put_att( ofid , ovarid(3) , 'units' , 'm/m' )                ;  call check(ierr)
  ierr = nf90_get_att( ifid , ivarid(3) , missvalname  , missval )         ;  call check(ierr)
  ierr = nf90_put_att( ofid , ovarid(3) , missvalname  , missval )         ;  call check(ierr)

  ierr = nf90_enddef( ofid )

  ierr = nf90_get_var( ifid , ivarid(1) , lon )   ;  call check(ierr)
  ierr = nf90_put_var( ofid , ovarid(1) , lon )   ;  call check(ierr)
  ierr = nf90_get_var( ifid , ivarid(2) , lat )   ;  call check(ierr)
  ierr = nf90_put_var( ofid , ovarid(2) , lat )   ;  call check(ierr)
  ierr = nf90_get_var( ifid , ivarid(3) , Slp0 )  ;  call check(ierr)

  !---------------------------------------------------------------------!

  call init_random_seed()

  do j = 1,ny
    do i = 1,nx

      if (Slp0(i,j)/=missval) then
         call mult_noise(noise,NOISE_AMPLITUDE)
!        ****************************
         Slp(i,j) = MEANSLOPE * noise
!        ****************************
      else
        Slp(i,j) = missval
      end if

    end do
  end do

  !---------------------------------------------------------------------!

  ierr = nf90_put_var( ofid , ovarid(3) , Slp )  ;  call check(ierr)

  ierr = nf90_close( ifid )  ;  call check(ierr)
  ierr = nf90_close( ofid )  ;  call check(ierr)


end program create_flat_slope


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine mult_noise(x,ampl)
implicit none
!
  real,intent(out):: x
  real,intent(in):: ampl
  real:: logampl
!
  logampl = log(ampl)
  call random_number(x)
  x = (2*x - 1)  *  logampl
  x = exp(x)
!
end subroutine mult_noise


subroutine init_random_seed()
implicit none
!
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed
!
  call random_seed(size = n)
  allocate(seed(n))
!
  call system_clock(count=clock)
!
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)
!
  deallocate(seed)
!
end subroutine init_random_seed


subroutine check(ierr)
  use netcdf
  implicit none
  integer, intent(in):: ierr
  if (ierr/=NF90_NOERR) then
    print *, nf90_strerror(ierr)
    call abort()
  end if
end subroutine
