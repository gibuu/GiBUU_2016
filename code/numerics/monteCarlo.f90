!***************************************************************************
!****m* /monteCarlo
! NAME
! module monteCarlo
!
! PURPOSE
! This module includes routines, which are useful for Monte Carlo decisions.
!***************************************************************************
module monteCarlo

  implicit none
  private

  public :: MonteCarloChoose, MonteCarloChoose2Dim

contains

  !***************************************************************************
  !****m* monteCarlo/MonteCarloChoose
  ! NAME
  ! integer function MonteCarloChoose (a, total)
  ! PURPOSE
  ! Monte-Carlo decision according to eq. 5.43-5.44 of Oliver Buss' thesis:
  ! * Given an array of weigths (which do not have to be normalized or positive),
  !   this routine chooses one channel and calculates the weight which this event should
  !   then be assigned.
  ! INPUTS
  ! real, dimension(:), intent(in) :: a    -- Array of weights
  ! RESULT
  ! real, intent(out), optional  :: total   -- Perweight which should be assigned to the chosen channel
  ! integer            :: channel -- The chosen channel
  !***************************************************************************
  integer function MonteCarloChoose (a, total)
    use random, only: rn
    use callStack, only: traceback

    real, dimension(:), intent(in) :: a
    real, intent(out),  optional   :: total

    real :: tot,x,r
    integer :: i

    tot=sum(abs(a))

    if (present(total)) total = 0.
    MonteCarloChoose = 0 ! indicating failure
    if (tot==0.) return ! ==> failure

    x=0.
    r=rn()*tot
    do i=lbound(a,dim=1),ubound(a,dim=1)
       x=x+abs(a(i))
       if (r<=x) then
          if (present(total)) total = sign(tot,a(i))
          MonteCarloChoose = i
          return
       end if
    end do

    write(*,*) 'Error in MonteCarloChoose: No decision has been performed! ', a(lbound(a,dim=1):ubound(a,dim=1))
    call traceback()

  end function MonteCarloChoose


  !***************************************************************************
  !****m* monteCarlo/MonteCarloChoose2Dim
  ! NAME
  ! function MonteCarloChoose2Dim (a, total_out) result(channel)
  ! PURPOSE
  ! Monte-Carlo decision according to eq. 5.43-5.44 of Oliver Buss' thesis:
  ! * Given an array of weigths (which do not have to be normalized or positive),
  !   this routine chooses one channel and calculates the weight which this event should
  !   then be assigned.
  ! INPUTS
  ! real, dimension(1:,1:), intent(in) :: a    -- Array of weights
  ! RESULT
  ! real, intent(out), optional   :: total_out   -- Perweight which should be assigned to the chosen channel
  ! integer, dimension(1:2)      :: channel -- The chosen channel (attention: first channel =1 !!!)
  !***************************************************************************
  function MonteCarloChoose2Dim (a, total_out) result(channel)
    use random, only: rn
    use callStack, only: traceback

    real, dimension(1:,1:), intent(in) :: a
    real, intent(out), optional        :: total_out
    integer, dimension(1:2)            :: channel

    real :: total,x,r
    integer :: i,j

    total=sum(abs(a))

    x=0.
    r=rn()*total
    do i=lbound(a,dim=1),ubound(a,dim=1)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          x=x+abs(a(i,j))
          if (r<=x) then
             if (present(total_out)) total_out = sign(total,a(i,j))
             channel=(/i,j/)
             return
          end if
       end do
    end do

    write(*,*) 'Error in MonteCarloChoose2Dim: No decision has been performed'
    call traceback()

  end function MonteCarloChoose2Dim



end module monteCarlo
