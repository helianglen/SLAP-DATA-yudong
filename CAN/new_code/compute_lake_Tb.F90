
! test desert module 
program lismemtest 


real :: fghz, theta, tsurf 
real :: em(2), tb(2) 
integer :: i 

theta = 0.0    ! zenith angle in deg,  looking straight down 
fghz = 1.4

! Header 
write(*, *) "fghz          emH            emV"
 Do i = 1, 100
  tsurf= i*0.3
  call cmem_water(fghz, theta, tsurf+273.16, em, tb) 
  write(*, '(3F10.3)') tsurf, tb(1), tb(2)
 End Do 

End 
