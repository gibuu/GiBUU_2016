!****************************************************************************
!****m* /sigmaNucleon
! NAME
! module sigmaNucleon
! PURPOSE
! Includes the cross sections for sigma-nucleon scattering in the resonance regime.
! Public routines:
! * sigmaNuc
!****************************************************************************
module sigmaNucleon

  implicit none
  Private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  Public :: sigmaNuc

contains


  !**************************************************************************************************
  !****s* sigmaNucleon/sigmaNuc
  ! NAME
  ! subroutine sigmaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,useHiEnergy,HiEnergySchwelle,plotFlag)
  !
  ! PURPOSE
  ! Evaluates sigma Nucleon -> anything cross sections and returns also a "preevent"
  !
  ! INPUTS
  ! * real, intent(in)                              :: srts                  ! sqrt(s) in the process
  ! * type(particle),dimension(1:2), intent(in)     :: teilchenIn            ! colliding particles
  ! * type(medium), intent(in)                      :: mediumATcollision     ! Medium informations at the position of the collision
  ! * real, intent(in) ,dimension(0:3)              :: momentumLRF           ! Total Momentum in LRF
  !
  ! High energy matching:
  ! * logical,intent(in)                            :: useHiEnergy
  ! * .true. if High-Energy cross sections are given by paramBarMesHE
  ! * real,intent(in)                               :: HiEnergySchwelle
  ! * threshold sqrt(s) for paramBarMesHE, i.e. at which energy the cross sections of paramBarMesHE are used
  !
  ! Debugging:
  ! * logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections
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
  ! The cross sections are based upon a parametrization by Golubeva. See routine golub in
  ! parametrizationBarMes.
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : baryon Resonances
  ! * 2-particle : pi N, sigma N
  !**************************************************************************************************
  subroutine sigmaNuc (srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,&
                       useHiEnergy,HiEnergySchwelle,plotFlag)

    use idTable
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
    logical,intent(in)                            :: useHiEnergy            ! .true. if High-Energy cross sections are given by paramBarMesHE
    real,intent(in)                               :: HiEnergySchwelle      ! threshold sqrt(s) for paramBarMesHE
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Output
    type(preEvent),dimension(1:3), intent(out) :: teilchenOut      ! colliding particles
    real, intent(out)                          :: sigmaTot         ! total Xsection
    real, intent(out)                          :: sigmaElast       ! elastic Xsection

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Cross sections
    real, dimension(Delta:nbar) :: sigmaRes    ! sigma N -> R crosssection
    real,dimension(-1:1) ::     piN            ! -> pi N, index denotes pion charge
    real ::     sigmaN                         ! -> sigma N
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Field to store the resonance masses
    real, dimension(Delta:nbar) :: massRes     !  Resonance masses
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real :: fluxCorrector        ! Correction of the fluxfactor due to different velocities
                                 ! in the medium compared to the vacuum
    type(particle) :: sigma_particle, nucleon_particle
    logical :: antiParticleInput, failFlag

    ! Field to store the resonance masses
          ! Resonance masses

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    teilchenOut(:)%ID=0                    ! ID of produced particles
    teilchenOut(:)%charge=0                ! Charge of produced particles
    teilchenOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    teilchenOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(teilchenIn,sigmaMeson,nucleon,sigma_particle,nucleon_particle,failFlag)
    If (failFlag) then
       Write(*,*) 'Wrong input in SigmaNuc', teilchenIn%ID
    end if

    If(sigma_particle%antiParticle) then
       ! This case is not considered yet
       write(*,*) 'sigma is antiparticle in "sigmaNuc"!!!',teilchenIN%ID,teilchenIN%antiparticle
       stop
    end if

    If(nucleon_particle%antiParticle) then
       ! Invert all particles in antiparticles
       nucleon_particle%Charge        =  -nucleon_particle%Charge
       nucleon_particle%antiparticle  = .false.
       sigma_particle%Charge          =  -sigma_particle%Charge
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

    ! Cutoff to kick those case out, that the cross section is zero
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
    If (Sum(teilchenOut(:)%Charge).ne.nucleon_particle%charge+sigma_particle%charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error' ,sigma_particle%Charge, &
            & nucleon_particle%Charge, teilchenOut(:)%Charge,teilchenOut(:)%ID
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
      use resonanceCrossSections, only: barMes_R_barMes, barMes2resonance
      use mediumDefinition, only : vacuum
      !use parametrizationsBarMes, only : golub
      use idTable, only : nucleon, pion, sigmaMeson
      use parBarMes_HighEnergy, only : paramBarMesHE
      use clebschGordan,only : clebschSquared
      use particleProperties, only: hadron
      use constants, only: mN, mPi

      real, dimension(1:3) ::  position
      logical :: perturbative
      !real, dimension(1:4) :: sigmaGolub
      !real, dimension(-1:1) :: piN_Golub
      !real :: pFinal,pInitial,detailedBalanceFactor
      integer :: pionCharge, nucCharge
      real :: elastic_Vacuum
      real, dimension(0:3) :: momentum_vacuum

      real :: sigmaTotal_HE,sigmaElast_HE ! High energy matchin

      position=0.5*(teilchenIN(1)%position+teilchenIN(2)%position)
      if(teilchenIN(1)%perturbative.or.teilchenIN(2)%perturbative) then
         perturbative=.true.
      else
         perturbative=.false.
      end if

      momentum_vacuum(1:3)=teilchenIn(1)%momentum(1:3)+teilchenIn(2)%momentum(1:3)
      momentum_vacuum(0)=FreeEnergy(teilchenIn(1))+FreeEnergy(teilchenIn(2))

      !#################################################################
      ! Evaluate partial cross sections
      !######################################################################


      !call golub(srts,1,sigmaGolub) ! Results by Golubeva


      !*******************************************************************************************
      ! sigma N -> pi N
      !*****************************************************************************************

      ! Full resonance contribution in the vacuum
      sigmaRes = barMes2resonance (sigmaMeson,nucleon,sigma_particle%charge,nucleon_particle%charge,.true.,vacuum, &
                                   momentum_vacuum,massRes,sigma_particle%Mass,nucleon_particle%Mass,position,perturbative,srts)

      ! High energy matching
      piN=0.
      If (srts > mN + mPi) then
         If(useHiEnergy) then
            call paramBarMesHE(HiEnergySchwelle,sigmaMeson,nucleon,sigma_particle%charge,nucleon_particle%charge,&
                 & mediumAtCollision,sigmaTotal_HE,sigmaElast_HE)
            do pionCharge=-1,1
               nucCharge=sigma_particle%charge+nucleon_particle%charge-pionCharge
               if((nucCharge.eq.0).or.(nucCharge.eq.1)) then
                  piN(pionCharge)=Max(sigmaTotal_HE-sum(sigmaRes),0.)*clebschSquared(1.,0.5,0.5, &
                       & real(pionCharge),real(nucCharge)-0.5)
               end if
            end do
         end if
      end if


      !*******************************************************************************************
      ! sigma N -> sigma N
      !*****************************************************************************************


      elastic_vacuum=barMes_R_barMes(sigmaMeson,nucleon,sigmaMEson,nucleon,&
           & sigma_particle%Charge,nucleon_particle%Charge,sigma_particle%Charge,nucleon_particle%Charge,  &
           & .false.,.false.,Vacuum,momentum_vacuum,sigma_particle%Mass,nucleon_particle%Mass, &
           & position,perturbative,srts)

      if (srts > mN + hadron(omegaMeson)%minmass) then
        sigmaN=Max(0., sigmaElast_HE - elastic_vacuum)
      else
        sigmaN=0.
      end if


      !*******************************************************************************************
      ! sigma N -> R
      !*****************************************************************************************

      ! Full resonance contribution in the medium
      sigmaRes = barMes2resonance(sigmaMeson,nucleon,sigma_particle%charge,nucleon_particle%charge,.true., &
                                  mediumAtCollision,momentumLRF,massRes, &
                                  sigma_particle%Mass,nucleon_particle%Mass,position,perturbative,srts)

      !###################################################################################################
      ! evaluate elastic Xsection
      !###################################################################################################

      sigmaElast=barMes_R_barMes(sigmaMeson,nucleon,sigmaMeson,nucleon,&
           & sigma_particle%Charge,nucleon_particle%Charge,sigma_particle%Charge,nucleon_particle%Charge, &
           & .false.,.false.,MediumAtCollision,momentumLRF,&
           & sigma_particle%Mass,nucleon_particle%Mass,position,perturbative,srts)+sigmaN


      !###################################################################################################
      ! Do the flux correction for each channel
      !###################################################################################################

      If(fluxCorrector_flag) then
         ! We do this for each channel since they might show up seperately in the output if makeoutput is called
         sigmaElast=sigmaElast*fluxcorrector
         sigmaN=sigmaN*fluxcorrector
         piN=piN*fluxcorrector
         sigmaRes=sigmaRes *fluxcorrector
      end if

      !###################################################################################################
      ! Sum up everything for the total cross section
      !###################################################################################################
      ! Be careful since sigma elast is already included in the partial cross sections, therefore it is not
      ! included in the total cross section

      sigmaTot=sigmaN + sum( piN ) + sum (sigmaRes )

    end subroutine evaluateXsections



    subroutine makeDecision
      use random, only : rn

      real :: summe, cut, cut2
      integer :: resID, totalCharge, pionCharge

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=sigma_particle%Charge+nucleon_particle%Charge
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

      ! sigma N production
      If(sigmaN.ge.cut) then
         teilchenOut(1)%Id=sigmaMeson
         teilchenOut(2)%Id=nucleon

         teilchenOut(1)%Charge=0
         teilchenOut(2)%Charge=totalCharge
         return
      end if
      cut=cut-sigmaN

      ! piN production
      Do pionCharge=-1,1
         If(piN(pionCharge).ge.cut) then
            teilchenOut(1)%Id=pion
            teilchenOut(2)%Id=nucleon

            teilchenOut(1)%Charge=pionCharge
            teilchenOut(2)%Charge=totalCharge-pionCharge
            return
         end if
         cut=cut-piN(pionCharge)
      end do

      ! Not event was generated:
      write(*,*) 'Error in makedecision of sigmaNuc', piN, cut
      stop

    end subroutine makeDecision

    !****************************************************************************
    !****s* sigmaNuc/makeOutput
    ! NAME
    ! subroutine makeOutput
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV].
    ! Filenames:
    ! * 'sigmaN_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'sigmaN_resProd.dat'            : Baryon resonance production
    ! * 'sigmaN_nonStrange_nuk.dat'     : non-strange meson with nucleon in final state
    !****************************************************************************
    subroutine makeOutPut
      logical, save :: initFlag=.true.
      real :: plab
      character(len=30), parameter :: outputFile(1:3) = (/ 'sigmaN_sigTotElast.dat   ', &
                                                           'sigmaN_resProd.dat       ', &
                                                           'sigmaN_nonStrange_nuk.dat' /)

      plab=SQRT(((srts**2-sigma_particle%mass**2-nucleon_particle%mass**2)/2./nucleon_particle%mass)**2-sigma_particle%mass**2)

      If (initFlag) then
         Open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         Open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         Open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         write (101,*) '# srts, plab, sigmaTot, sigmaElast '
         write (102,*) '# srts, plab, sigmaRes(2:40) '
         write (103,*) '# srts, plab, piN(-1:1), sigmaN '
         initFlag=.false.
      else
         Open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
      end If
      write (101,'(4F9.3)') srts, plab,sigmaTot, sigmaElast
      write (102,'(41F9.3)') srts, plab, sigmaRes(2:40)
      write (103,'(6F9.3)') srts, plab,piN, sigmaN
      Close(101)
      Close(102)
      Close(103)

    end subroutine makeOutPut

  end subroutine sigmaNuc


end module sigmaNucleon
