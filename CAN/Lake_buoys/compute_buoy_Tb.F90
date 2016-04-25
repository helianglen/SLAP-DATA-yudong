
! compute lake Tb from water temperature 
program lismemtest 


real :: fghz, theta, tsurf 
real :: em(2), tb0cm(2), tb70cm(2), t0cm, t70cm
integer :: i, iagrc, ierr
character*100 :: infile

 if (iargc() .ne. 1) then
   write(*, *) "usage: compute_buoy_tb inputfile" 
   stop 
 end if 

 call getarg(1, infile) 

theta = 0.0    ! zenith angle in deg,  looking straight down 
fghz = 1.4

  open(15, file=trim(infile), form='formatted', status='old')
  Do i=1, 20000
     read(15, *, iostat=ierr) t0cm, t70cm
     if (ierr .eq. 0 ) then 
       call cmem_water(fghz, theta, t0cm+273.16, em, tb0cm)
       call cmem_water(fghz, theta, t70cm+273.16, em, tb70cm)
       ! tbV and tbH identical 
       write(*, '(F12.3, 3(",", F12.3) )') t0cm, tb0cm(1), t70cm, tb70cm(1) 
    else 
      exit
    end if 
  End Do 
  close(15) 

End 
