program testMaxwell
  use randomMaxwell

  implicit none

  call check1

contains

  subroutine check1
    use histf90

    integer :: ir, nr=1000000
    real :: r
    type(histogram) :: HH

    call CreateHist(HH, 'vals', 0.0, 10.0, 0.1)

    call initMaxwell(0.135, 0.500)

    do ir=1,nr

       r = rnMaxwell()
       call AddHist(HH, r, 1.0)
    end do

    call WriteHist(HH, 114, mul=1.0/nr)

  end subroutine check1


end program testMaxwell
