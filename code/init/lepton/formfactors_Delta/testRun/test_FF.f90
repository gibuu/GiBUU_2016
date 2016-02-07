
program testFF
use particleProperties, only: initParticleProperties
use FF_Delta_production
implicit none

real, parameter :: dq=0.005
real, parameter :: dw=0.1
integer :: i,j
real :: Qs,W,c3v,c4v,c5v,c6v,c3a,c4a,c5a,c6a

call initParticleProperties

open(10,file='deltaFF.dat')

do j=0,10
   W=dw*j+1.1
   do i=0,400
      Qs=dq*i
      call formfactors_Delta(Qs,W,0,c3v,c4v,c5v,c6v,c3a,c4a,c5a,c6a)
      write(10,'(10F12.5)') Qs,W,c3v,c4v,c5v,c6v,c3a,c4a,c5a,c6a
   end do
end do

write(*,*) 'end of test FF delta'

end program testFF
