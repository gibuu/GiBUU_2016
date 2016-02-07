program testBosted
  use ParamEP
  use inputGeneral

  implicit none

  real :: W,Q2,eps,XS
  integer :: i

  call readinputGeneral
  call init_database


  eps = 0.99
  Q2 = 1.0

  call CalcParamEP(1.5,Q2,eps, XS) ! dummy in order to read input

  do i=110,300
     W = i*0.01
     call CalcParamEP(W,Q2,eps, XS)
     write(119,'(3f12.4,1P,e12.4)') W,Q2,eps, XS 
  end do

end program testBosted
