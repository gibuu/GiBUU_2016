!****************************************************************************
!****m* /pionDelta_resonance
! NAME
! module pionDelta_resonance
! PURPOSE
! Includes the cross sections for pion-delta scattering in the resonance regime
! Implemented are the following reactions:
!   * pion delta -> X
! Public routines:
!   * pionDelta
! NOTE
! pi Delta -> K Kbar N not yet implemented
!****************************************************************************
module pionDelta_resonance

  implicit none
  Private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.


  Public :: pionDelta

contains


  !****************************************************************************
  !****m* pionDelta_resonance/pionDelta
  ! NAME
  ! subroutine pionDelta(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,plotFlag)
  !
  ! PURPOSE
  ! Evaluates pion Delta -> anything cross sections and returns also a "preevent"
  !
  ! RESULT
  ! * real, intent(out)                                        :: sigmaTot         ! total Xsection
  ! * real, intent(out)                                        :: sigmaElast       ! elastic Xsection
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross sections to decide on a final state with
  ! maximal 3 final state particles. These are returned in the vector teilchenOut. The kinematics of these teilchen is
  ! only fixed in the case of a single produced resonance. Otherwise the kinematics still need to be established. The
  ! result is:
  ! * type(preEvent),dimension(1:3), intent(out)               :: teilchenOut     ! colliding particles
  !
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : baryon Resonances
  ! * 2-particle : Kaon Lambda , Sigma Kaon
  ! NOTE
  ! pi Delta -> K Kbar N not yet implemented
  !****************************************************************************

  subroutine pionDelta(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,plotFlag)

    use idTable, only: pion, Delta, nbar, kaon, Lambda, sigmaResonance
    use particleDefinition
    use mediumDefinition
    use preEventDefinition, only : preEvent
    use twoBodyTools, only : velocity_correction, convertToAntiParticles, pcm,searchInInput
    use RMF, only : getRMF_flag

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Input
    real, intent(in)                              :: srts                  ! sqrt(s) in the process
    type(particle),dimension(1:2), intent(in)     :: teilchenIn            ! colliding particles
    type(medium), intent(in)                      :: mediumATcollision     ! Medium informations at the position of the collision

    logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections
    real, intent(in) ,dimension(0:3)              :: momentumLRF           ! Total Momentum in LRF

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Output
    type(preEvent),dimension(1:3), intent(out) :: teilchenOut      ! colliding particles
    real, intent(out)                          :: sigmaTot         ! total Xsection
    real, intent(out)                          :: sigmaElast       ! elastic Xsection


    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! CrossSections named according to final states
    real ::  lambdaKaon
    real, dimension  (0:1)  :: sigmaKaon             ! index = charge of final state kaon
    real, dimension (-2:1)  :: kaonKaonBarN          ! index= final state
!    real ,  dimension(lBound(baryon,dim=1):uBound(baryon,dim=1)) :: sigmaRes      ! pi N -> R crosssection

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Field to store the resonance masses
!    real ,  dimension(lBound(baryon,dim=1):uBound(baryon,dim=1)) :: massRes       !  Resonance masses


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real :: fluxCorrector        ! Correction of the fluxfactor due to different velocities
                                 ! in the medium compared to the vacuum
    type(particle) :: pion_particle, delta_particle
    logical :: antiParticleInput, failFlag


    ! partial cross sections for pi Delta -> R
    real , dimension(Delta:nbar) :: sigmaRes
    ! Field to store the resonance masses
    real , dimension(Delta:nbar) :: massRes      ! Resonance masses

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    teilchenOut(:)%ID=0                    ! ID of produced particles
    teilchenOut(:)%charge=0                ! Charge of produced particles
    teilchenOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    teilchenOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(teilchenIn,pion,delta,pion_particle,delta_particle,failFlag)
    If (failFlag) then
       Write(*,*) 'Wrong input in KaonNuc', teilchenIn%ID
    end if

    If(pion_particle%antiParticle) then
          ! This case is not considered yet
          write(*,*) 'pion is antiparticle in "pionnuc"!!!',teilchenIN%ID,teilchenIN%antiparticle
          stop
    end if

    If(delta_particle%antiParticle) then
   ! Invert all particles in antiparticles
       delta_particle%antiParticle=.false.
       delta_particle%charge=-delta_particle%charge
       pion_particle%charge=-pion_particle%charge
       antiParticleInput=.true.
    else
          antiParticleInput=.false.
    end if

    ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
    if( .not.getRMF_flag() ) then
      fluxCorrector=velocity_correction(teilchenIn)
    else
      fluxCorrector=1.
    end if

    ! (2) Evaluate the cross sections
    call evaluateXsections

    ! Cutoff to kick the case out, that the cross section is zero
    if(sigmaTot.lt.1E-12) then
       sigmatot=0.
       sigmaElast=0.
       return
    end if

    ! (3) Plot them if wished
    If(Present(PlotFlag).or.debugFlag) then
       If (plotFlag.or.debugFlag)  call makeOutput
    end if

    ! (4) Define final state
    call MakeDecision

    ! (5) Check Output
    If (Sum(teilchenOut(:)%Charge).ne.pion_particle%Charge+delta_particle%Charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error' ,pion_particle%Charge,delta_particle%Charge,&
            &  teilchenOut(:)%Charge,teilchenOut(:)%ID
       stop
    end if

    ! (6) Invert particles in antiParticles if input included antiparticles
    If(antiParticleInput) then
       IF(debugFlagAnti) write(*,*) teilchenOut
       call convertToAntiParticles(teilchenOut)
       IF(debugFlagAnti) write(*,*) teilchenOut
    end if

  contains

    subroutine evaluateXsections
      !****s* pionDelta/evaluateXsections
      ! NAME
      ! subroutine evaluateXsections
      ! PURPOSE
      ! Evaluates all Xsections for pi N scattering
      !***
      use particleDefinition
!       use mediumDefinition, only: vacuum
      use particleProperties, only: hadron
      use parametrizationsBarMes, only: huangd, huanglamd
      use resonanceCrossSections, only: barMes_R_barMes, barMes2resonance
      use clebschGordan, only: clebschSquared
      use constants, only: mK

      real :: sigmaHuangDelta
      !logical :: background, propagated
      !real,dimension(0:3)                 :: momentum_vacuum        ! Total Momentum in vacuum
      real, dimension(1:3) :: position
      logical :: perturbative

      integer :: kaoncharge
      real :: isoZ_sigma

      position=0.5*(teilchenIN(1)%position+teilchenIN(2)%position)
      if(teilchenIN(1)%perturbative.or.teilchenIN(2)%perturbative) then
         perturbative=.true.
      else
         perturbative=.false.
      end if

      !momentum_vacuum(1:3)=teilchenIn(1)%momentum(1:3)+teilchenIn(2)%momentum(1:3)
      !momentum_vacuum(0)=FreeEnergy(teilchenIn(1))+FreeEnergy(teilchenIn(2))


      !#################################################################
      ! Evaluate partial cross sections
      !######################################################################


      !*****************************************************************************************************
      ! -> Kaon KaonBar N
      !*****************************************************************************************************
      kaonKaonBarN(:)=0.

      ! Not yet implemented


      !*****************************************************************************************************
      ! -> Lambda Kaon
      !*****************************************************************************************************
      lambdaKaon=0.
      If (srts > hadron(Lambda)%mass + mK) then
         If((Delta_particle%charge+pion_particle%charge.eq.0).or.(Delta_particle%charge+pion_particle%charge.eq.1)) then
            ! hunaglamd gives : delta^{++} pi- -> Lambda Kaon^{+}
            lambdaKaon = 2.*huangLamd(srts)*clebschSquared(1.5,1.0,0.5,real(delta_particle%Charge)-0.5,real(pion_particle%Charge))
            ! divide out "delta^{++} pi- -> Lambda Kaon^{+}" clebsch and multiply be real clebsch
         end if
      end if

      !*****************************************************************************************
      ! -> Sigma Kaon
      !*****************************************************************************************
      sigmaKaon=0.
      If((pion_particle%charge+delta_particle%charge.le.2).and.(pion_particle%charge+delta_particle%charge.ge.-1)) then

         ! divide out Clebsch Gordan coefficient for (pi- D++ -> S0 k) and multiply by initial clebsch Gordan
         sigmaHuangDelta = 6.*huangd(srts)*clebschSquared(1.0,1.5,0.5,real(pion_particle%charge),real(delta_particle%charge)-0.5)

         do kaonCharge=0,1
            isoZ_Sigma=pion_particle%charge+delta_particle%charge-kaonCharge
            if(abs(isoZ_sigma).lt.1.01) then
               sigmaKaon(kaonCharge)=sigmaHuangDelta*clebschSquared(1.0,0.5,0.5,isoZ_Sigma,real(kaonCharge)-0.5)
            end if
         end do
      end if

      !*******************************************************************************************
      ! -> R
      !*****************************************************************************************

      ! Full resonance contribution in the medium
      sigmaRes = barMes2resonance (pion,delta,pion_particle%charge,delta_particle%charge,.true.,mediumAtCollision, &
                                   momentumLRF,massRes,pion_particle%Mass,delta_particle%Mass,position,perturbative,srts)

      !###################################################################################################
      ! evaluate elastic Xsection
      !###################################################################################################

      sigmaElast=barMes_R_barMes(pion,delta,pion,delta,&
           & pion_particle%Charge,delta_particle%Charge,pion_particle%Charge,delta_particle%Charge, &
           & .false.,.false.,MediumAtCollision,momentumLRF,pion_particle%Mass,delta_particle%Mass, &
           & position,perturbative,srts)

      !###################################################################################################
      ! Do the flux correction for each channel
      !###################################################################################################

      ! We do this for each channel since they might show up seperately in the output if makeoutput is called
      If(fluxCorrector_flag) then
         ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
         sigmaElast=sigmaElast*fluxcorrector
         lambdaKaon=lambdaKaon *fluxcorrector
         sigmaKaon=sigmaKaon*fluxcorrector
         kaonKaonBarN=kaonKaonBarN *fluxcorrector
         sigmaRes=sigmaRes *fluxcorrector
      end if


      !###################################################################################################
      ! Sum up everything for the total cross section
      !###################################################################################################
      ! Be careful since sigma elast is already included in the partial cross sections, therefore it is not
      ! included in the total cross section

      sigmaTot=lambdaKaon +Sum ( sigmaKaon ) + sum ( kaonKaonBarN ) + sum (sigmaRes )


    end subroutine evaluateXsections


    !********************************************************************************************************************************
    !****s* pionDelta/makeDecision
    ! NAME
    ! subroutine MakeDecision
    ! PURPOSE
    ! Decides on the final state which is returned via teilchenOut by Monte-Carlo.
    !  * Assigns charges and ID's.
    !  * Only for resonance-production also the mass is assigned, since the mass of the resonance needed to be calculated earlier.
    ! The Monte-Carlo routine is adding up channels until the sum is exceeding x*sigma(total). x has a flat distribution in [0,1].
    ! The last added channel is then the one which is chosen for the event. After choosing a channel, the subroutine is returning to
    ! the calling routine.
    !********************************************************************************************************************************
    subroutine MakeDecision
      use random, only: rn

      real :: summe,cut,cut2
      integer :: totalCharge,charge,resID

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=pion_particle%Charge+delta_particle%Charge
      !############################################################
      ! (1) Resonance production
      !############################################################
      If (sum(sigmaRes)>=cut) then
         summe=0.
         cut2=rn()*sum(sigmaRes)
         Do resId=Delta,nbar
            summe=summe+sigmaRes(resID)
            If (summe>=cut2) exit
         End do
         teilchenOut(1)%Id=resID
         teilchenOut(1)%Charge=totalCharge
         teilchenOut(1)%Mass=massRes(resID)
         return
      end if
      cut=cut-sum(sigmaRes)

      !############################################################
      ! (2) Two -body final states
      !############################################################

      If ( lambdaKaon .ge.cut) then
         teilchenOut(1)%Id=kaon
         teilchenOut(2)%Id=Lambda
         teilchenOut(1)%Charge=totalCharge
         teilchenOut(2)%Charge=0
         return
      end if
      cut=cut-lambdaKaon

      If (Sum ( sigmaKaon ) .ge.cut) then
         teilchenOut(1)%Id=kaon
         teilchenOut(2)%Id=sigmaResonance
         cut2=rn()*sum(sigmaKaon)
         Do charge=0,1
            If(sum(sigmaKaon(0:charge)).ge.cut2) then
               teilchenOut(1)%Charge=charge
               teilchenOut(2)%Charge=totalCharge-charge
               exit
            end if
         End do
         return
      end if
      cut=cut-sum(sigmaKaon)


      !############################################################
      ! (3) Three body final state
      !############################################################

      !      Not yet implemented

      !############################################################
      ! Error message if no channel is chosen
      !############################################################
      write (*,*) 'Error in makeDecision : No decision made', &
           & cut,    sum(kaonKaonBarN) , teilchenOut(:)%ID,teilchenOut(:)%Charge,sigmaTot

      Stop

    end subroutine MakeDecision


    !*****************************************************************************************
    !****s* pionDelta/makeOutput
    ! NAME
    ! subroutine makeOutput
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV].
    ! Filenames:
    ! * 'piDelta_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'piDelta_strangeProd.dat'        : strangeness production
    ! * 'piDelta_resProd.dat'            : Baryon resonance production
    !*****************************************************************************************
    subroutine makeOutput
      logical, save :: initFlag=.true.
      real :: plab
      character(len=30), parameter :: outputFile(1:3) = (/ 'piDelta_sigTotElast.dat', &
                                                           'piDelta_strangeProd.dat', &
                                                           'piDelta_resProd.dat    ' /)

      plab=SQRT(((srts**2-pion_particle%Mass**2-delta_particle%Mass**2)/2./delta_particle%Mass)**2-pion_particle%Mass**2)

      If (initFlag) then
         Open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         Open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         Open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         write (101,*) '# srts, plab, sigmaTot, sigmaElast '
         write (102,*) '# srts, plab, lambdaKaon, sigmaKaon(0:1), kaonKaonBarN(-2:1)'
         write (103,*) '# srts, plab, sigmaRes(2:40)'
         initFlag=.false.
      else
         Open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
      end If
      write (101,'(4F9.3)') srts, plab,sigmaTot, sigmaElast
      write (102,'(9F9.3)') srts, plab,lambdaKaon, sigmaKaon, kaonKaonBarN
      write (103,'(41F9.3)') srts, plab, sigmaRes(2:40)
      Close(101)
      Close(102)
      Close(103)
    end subroutine makeOutput

  end subroutine pionDelta


end module pionDelta_resonance
