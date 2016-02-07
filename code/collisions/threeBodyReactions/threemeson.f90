!***************************************************************************
!****m* /ThreeMeson
! NAME
! module ThreeMeson
!
! PURPOSE
! This module implements the back reaction 3 -> of mesonic decay channels
!
! INPUTS
! (none)
!***************************************************************************
module ThreeMeson


  
  implicit none

  PRIVATE

  PUBLIC :: DoThreeMeson
  
contains
  
  subroutine DoThreeMeson()
    use constants, only: mPi,pi
    use decayChannels
    use IdTable, only: nMes, pion, rho, omegaMeson, phi
    use mediumDefinition
    use mesonPotentialModule, only: vecMes_massShift
    use mesonWidthMedium, only: decayWidthMesonMedium, WidthMesonMedium
    use particleDefinition
    use ParticleProperties
    use twoBodyTools, only: pCM
    use VolumeElements

    integer :: iDecay,ID,dID,iChannel !,i
    logical :: doInit
    integer, dimension(1:3) :: iEns,iInd
    type(particle), POINTER :: Part1, Part2, Part3
    real :: srts
    integer, dimension(1:4,1:20,1:2) :: ArrChannel
    integer, dimension(1:4) :: nChannel
    real :: gammaIn,gammaTot,Phi3,m0
    real, dimension(nDecays) :: decayWidth
    logical :: pauliFlag
    type(medium)        :: mediumATcollision
    real,dimension(0:3) :: momentumLRF

    ! indicate same particles in every decay channel.
    ! bool coding: 1,2,4
    ! * part 1 and part 2: '1 or 2' = 3
    ! * part 1 and part 2 and part 3: '1 or 2 or 4' = 7
    ! Please keep up to date if changed in decayChannels
    integer, parameter, dimension(1:4) :: isSameBoolArr = (/ 3, 0, 7, 0 /)


!!$    do i=1,200
!!$       srts = 2*mPi + 0.02d0*i
!!$       Phi3 = calcPhi3(srts,mPi)
!!$       write(*,*) srts,phi3,pCM(srts,mPi,mPi)/(4*pi*srts)
!!$    enddo
!!$    stop

    
    ArrChannel = 0
    nChannel = 0
    
    do ID=pion,pion+nMes-1
       do iDecay=1,nDecays
          dID=-hadron(ID)%decaysID(iDecay) ! 3Body are negative
          if (dID>0) then
             if (hadron(ID)%decays(iDecay) > 0.0) then
                nChannel(dID) = nChannel(dID)+1
                ArrChannel(dID,nChannel(dID),1) = ID
                ArrChannel(dID,nChannel(dID),2) = iDecay
             end if
          end if
       end do
    end do

    write(*,*) 'decay channels to consider:'
    do dID=1,nDecay3bodyMeson
       write(*,'(50i4)') dID,nChannel(dID),ArrChannel(dID,1:nChannel(dID),1)
       write(*,'(50i4)') dID,nChannel(dID),ArrChannel(dID,1:nChannel(dID),2)
    end do
    
    
    do dID=1,nDecay3bodyMeson
       write(*,*) 'dID = ',dID

       
       doInit = .true.
       do
          if (.not.VolumeElements_3Body(Decay3BodyMeson(dID), isSameBoolArr(dID), doInit, iEns,iInd, Part1,Part2,Part3)) exit

          ! 3 particles found

          srts = sqrtS(Part1,Part2,Part3)
          select case(dID)
          case(2,3)
             Phi3 = CalcPhi3(srts,mPi)
          case(1,4)
             Phi3 = CalcPhi3(srts,hadron(102)%mass)
          end select

          write(*,*) iEns,iInd,srts,Phi3
          
          do iChannel=1,nChannel(dID)
             ID = ArrChannel(dID,iChannel,1)
             decayWidth = decayWidthMesonMedium(ID,srts,0,pauliFlag)
             gammaIn  = decayWidth(ArrChannel(dID,iChannel,2))
             gammaTot = WidthMesonMedium(ID,srts, momentumLRF, mediumATcollision)
             m0 = hadron(ID)%mass
             if (ID==rho .or. ID==omegaMeson .or. ID==phi) then
                ! take into account a possible in-medium mass shift of the vector mesons
                m0 = m0 + vecMes_massShift(ID, mediumATcollision%density)
             end if


             
             
             write(*,*) ID,gammaTot,gammaIn
             
!             write(*,*) srts, gammaTotal
          end do

          
       end do
       
    end do ! iDecay
    
    stop
    
  end subroutine DoThreeMeson

  ! there is threeBodyPhaseSpace/Integrate_3bodyPS, which does the same,
  ! but just by simply adding 100 values.
  ! This here is closer to mesonWidthVacuum/threepi
  ! return value in GeV^2
  real function calcPhi3(srts,m3)

    use gauss_integration, only: sg20r, rg20r
    use twoBodyTools, only: pCM
    use constants, only: mPi,pi

    real, intent(in) :: srts     ! sqrt(s) = mass of decaying particle
    real, intent(in) :: m3       ! mass of third particle
    integer :: ns
    integer, parameter :: n=3
    real, dimension(1:20*n) :: absi,orde
    real, parameter :: mm = mPI**2
    real :: mm3,s,resu

    if (srts>2.*mPi+m3) then

       call sg20r(2.*mPi,srts-m3,n,absi,ns)
       absi = absi**2
       mm3 = m3**2
       s = srts**2
       
       orde = sqrt( ((absi-s-mm3)**2-4*s*mm3) * (absi-4*mm) )

       call rg20r(2.*mPi,srts-m3,n,orde,resu)
       calcPhi3 = resu/(64*pi**3*s)

    else
       calcPhi3 = 0
    end if

  end function calcPhi3

  
end module ThreeMeson
