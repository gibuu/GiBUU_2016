module deltaTest

  real, parameter :: mom=0.1
  real, parameter ::  mass=1.2
  integer, parameter ::  charge=-1

contains
  subroutine initDeltaTest(part)
    ! Initialisiert Deltas, die in der x-y Ebene mit z=0 starten und dann mit gegebenem p_z
    ! in die z-Richtung loslaufen.
    use particleDefinition
    use random
    use idTable, only : delta
    use deltaWidth
    implicit none
    integer :: j
    type(particle), dimension(:,:) :: part
    real :: imsig2,imsig3,imsigq
    integer :: i
    real :: rho


    open(100,file="width_oset.dat")
    do i=0,100
       rho=0.002*float(i)
       call  deloset(mass,rho,imsig2,imsig3,imsigq)
       write(100,'(5G18.4)') rho,imsig2,imsig3, imsigq,2*(imsig2+imsig3+imsigq)
    end do
    close(100)

    do i=lbound(part,dim=1),ubound(part,dim=1)
       call setToDefault(part(i,:))
       do j=1,50
          part(i,j)%ID=delta
          part(i,j)%charge=charge
          part(i,j)%mass=mass
          part(i,j)%momentum=0.
          part(i,j)%momentum(3)=mom
          part(i,j)%position(3)=0.
          part(i,j)%position(1)=rn()*2.
          part(i,j)%position(2)=rn()*2.
       end do
    end do
  end subroutine initDeltaTest


  subroutine countDeltas(part,time)
    ! Zählt alle Deltas, die noch nicht gestossen haben und schreibt die Zahl nach "numDeltas.dat".
    use particleDefinition
    use random
    use idTable, only : delta
    implicit none
    integer :: i,j
    type(particle), dimension(:,:) :: part
    real :: time
    integer :: numDeltas

    numDeltas=0
    do i=lbound(part,dim=1),ubound(part,dim=1)
       do j=lbound(part,dim=2),ubound(part,dim=2)
          if(part(i,j)%ID.eq.delta) then
             if(part(i,j)%event(1).eq.0) then
                numDeltas=numDeltas+1
             end if
          end if
       end do
    end do
    open(120,file="numDeltas.dat",position='append')
    write(120,*) time, numdeltas
    close(120)


  end subroutine countDeltas



end module deltaTest
