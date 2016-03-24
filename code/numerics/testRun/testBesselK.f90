program testBesselK
  use besselK

  implicit none

  call test2

contains
  subroutine test2

    use fgsl, only: kn => fgsl_sf_bessel_kcn

    integer :: ix
    real :: x

    do ix=1,500
       x = ix * 0.01
       write(112,*) x, besselK2(x), kn(2,x)
    end do

  end subroutine test2

end program testBesselK
