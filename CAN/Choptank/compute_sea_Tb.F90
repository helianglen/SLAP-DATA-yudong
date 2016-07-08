
! compute lake Tb from water temperature 
program lismemtest 

!average readings
real, parameter :: tsurf=19.0, sal=15.5   ! C and PPT

real :: fghz, theta
real :: em(2), tb(2)

theta = 0.0    ! zenith angle in deg,  looking straight down 
fghz = 1.4

  call sea_water(fghz, theta, tsurf+273.15, sal, em, tb)

  write(*, '(3F12.3)') tsurf, tb(1), em(1)  

End 
