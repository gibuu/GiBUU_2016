!****************************************************************************
!****m* /phiNucleon
! NAME
! module phiNucleon
! PURPOSE
! Includes the cross sections for phi-nucleon scattering in the resonance regime.
! Public routines:
! * phiNuc
!****************************************************************************
module phiNucleon
  implicit none
  Private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  Public :: phiNuc

contains

  !****************************************************************************
  !****s* phiNucleon/phiNuc
  ! NAME
  ! subroutine phiNuc (srts,teilchenIN,mediumATcollision,teilchenOUT,sigmaTot,sigmaElast,useHiEnergy,HiEnergySchwelle,plotFlag)
  !
  ! PURPOSE
  ! Evaluates phi Nucleon -> anything cross sections and returns also a "preevent"
  !
  ! INPUTS
  ! * real, intent(in)                              :: srts                  ! sqrt(s) in the process
  ! * type(particle),dimension(1:2), intent(in)     :: teilchenIn            ! colliding particles
  ! * type(medium), intent(in)                      :: mediumATcollision     ! Medium informations at the position of the collision
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
  ! OUTPUT
  ! * real, intent(out)                                        :: sigmaTot         ! total Xsection
  ! * real, intent(out)                                        :: sigmaElast       ! elastic Xsection
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross sections to decide on a final state with
  ! maximal 3 final state particles. These are returned in the vector teilchenOut. The kinematics of these teilchen is
  ! only fixed in the case of a single produced resonance. Otherwise the kinematics still need to be established. The
  ! result is:
  ! * type(preEvent),dimension(1:3), intent(out)               :: teilchenOut     ! colliding particles
  !
  ! The cross sections are based upon a parametrization by Golubeva. See routine golub_phi in
  ! parametrizationBarMes.
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : baryon Resonances
  ! * 2-particle : pi N, phi N, pi pi N
  !****************************************************************************
  subroutine phiNuc (srts, teilchenIN, mediumATcollision, teilchenOUT, sigmaTot, sigmaElast, &
                     useHiEnergy, HiEnergySchwelle, plotFlag)

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
    logical,intent(in)                            :: useHiEnergy            ! .true. if High-Energy cross sections are given by paramBarMesHE
    real,intent(in)                               :: HiEnergySchwelle      ! threshold sqrt(s) for paramBarMesHE

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Output
    type(preEvent),dimension(1:3), intent(out) :: teilchenOut      ! colliding particles
    real, intent(out)                          :: sigmaTot         ! total Xsection
    real, intent(out)                          :: sigmaElast       ! elastic Xsection

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Cross sections
    real,dimension(-1:1) ::     piN            ! -> pi N, index denotes pion charge
    real ::     phiN         ! -> phi N
    real ::     pipiN          ! -> pi pi N

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real :: fluxCorrector        ! Correction of the fluxfactor due to different velocities
                                 ! in the medium compared to the vacuum
    type(particle) :: phi_particle, nucleon_particle
    logical :: antiParticleInput, failFlag



    ! partial cross sections for phi N -> R
    !real , dimension(:),Allocatable :: sigmaRes

    ! Field to store the resonance masses
    !real , dimension(:),Allocatable ::massRes      ! Resonance masses
    !Allocate(sigmaRes(lBound(baryon,dim=1):uBound(baryon,dim=1)))
    !Allocate(massRes(lBound(baryon,dim=1):uBound(baryon,dim=1)))


    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    teilchenOut(:)%ID=0                    ! ID of produced particles
    teilchenOut(:)%charge=0                ! Charge of produced particles
    teilchenOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    teilchenOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(teilchenIn,phi,nucleon,phi_particle,nucleon_particle,failFlag)
    If (failFlag) then
       Write(*,*) 'Wrong input in PhiNuc', teilchenIn%ID
    end if

    If(phi_particle%antiParticle) then
       ! This case is not considered yet
       write(*,*) 'phi is antiparticle in "phiNuc"!!!',teilchenIN%ID,teilchenIN%antiparticle
       stop
    end if

    If(nucleon_particle%antiParticle) then
       ! Invert all particles in antiparticles
       nucleon_particle%Charge        =  -nucleon_particle%Charge
       nucleon_particle%antiparticle  = .false.
       phi_particle%Charge          =  -phi_particle%Charge
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
    If (Sum(teilchenOut(:)%Charge).ne.nucleon_particle%charge+phi_particle%charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error' ,phi_particle%Charge, &
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

    !**************************************************************
    !****s* phiNuc/evaluateXsections
    ! NAME
    ! subroutine evaluateXsections
    ! PURPOSE
    ! Evaluates phi Nucleon -> anything cross sections
    ! NOTES
    ! There are no resonance contributions to phi N scattering.
    ! The contributions are given by Golubeva (see golub_phi).
    !**************************************************************
    subroutine evaluateXsections

      use parametrizationsBarMes, only : golub_phi
      use parBarMes_HighEnergy, only : paramBarMesHE
      use idTable, only : nucleon, pion, phi
      use output, only : writeparticle
      use constants, only: mN, mPi
      use particleProperties, only: hadron

      real, dimension(1:4) :: sigmaGolub
      real :: pFinal, pInitial,detailedBalanceFactor

      real :: sigmaTotal_HE,sigmaElast_HE ! High energy matchin


      !#################################################################
      ! Evaluate partial cross sections
      !######################################################################

      sigmaGolub = golub_phi (srts, phi_particle%mass)   ! results by Golubeva for phi meson scattering


      !*******************************************************************************************
      ! phi N -> pi N
      !*****************************************************************************************
      ! piN = cross section by Golubeva (pi^- p-> phi N) by detailed balance - vacuum cross section by resonances

      pFinal=pcm(srts,mPi,mN)
      pInitial=pcm(srts,phi_particle%mass,nucleon_particle%mass)
      If(pinitial.lt.1E-12) then
         write(*,*) 'WARNING: pInitial is zero in phiNuc', pinitial
         write(*,*) 'phi meson:'
         call writeparticle(6,0,0,phi_Particle)
         write(*,*) 'nucleon:'
         call writeparticle(6,0,0,nucleon_Particle)
         detailedBalanceFactor= 0.
      else
         detailedBalanceFactor= 1./3.*(pFinal/pInitial)**2
         ! given by detailed balance: factor 1/3 due to (2j+1)-Terms in cross section and different
         ! spins in initial and final state
      end if

      sigmaGolub(1)=3./2.* sigmaGolub(1)  ! 3./2. since Golub returns cross section pi^- p -> phi N, with 3./.2 we divide by the isospin-clebsch

      If(nucleon_particle%charge.eq.0) then
         piN(-1)=  2./3.*sigmaGolub(1) * detailedBalanceFactor
         piN(0)=   1./3.*sigmaGolub(1) * detailedBalanceFactor
         piN(1)=   0.
      else If(nucleon_particle%charge.eq.1) then
         piN(-1)= 0.
         piN(0)=  1./3.*sigmaGolub(1)  * detailedBalanceFactor
         piN(1)=  2./3.*sigmaGolub(1)  * detailedBalanceFactor
      else
         write(*,*) 'Error in phiNuc', nucleon_particle
      end if

      !*******************************************************************************************
      ! phi N -> phi N
      !*****************************************************************************************

      if (srts>mN+hadron(phi)%minmass) then
        phiN=Max(0., sigmaGolub(3))
      else
        phiN=0.
      end if

      !*******************************************************************************************
      ! -> pi pi N
      !*******************************************************************************************

      ! Evaluate Phi N -> N pi pi by taking the total inelastic crossection by Golub and subtracting then the
      ! inelastic cross sections for "phi N -> N pi"

      If(srts.gt.(2*mPi+mN)) then
         pipiN=Max(0.,sigmaGolub(4)-Sum(piN))
         ! Matching to High energy region
         if(useHiEnergy) then
            call paramBarMesHE(HiEnergySchwelle,phi,nucleon,phi_particle%charge,nucleon_particle%charge,&
                 & mediumAtCollision,sigmaTotal_HE,sigmaElast_HE)
            pipiN=max(pipiN,sigmaTotal_HE-Sum(piN)-phiN)
         end if
      else
         piPiN=0.
      end if


      !###################################################################################################
      ! evaluate elastic Xsection
      !###################################################################################################

      sigmaElast=phiN

      !###################################################################################################
      ! Do the flux correction for each channel
      !###################################################################################################

      If(fluxCorrector_flag) then
         ! We do this for each channel since they might show up seperately in the output if makeoutput is called
         sigmaElast=sigmaElast*fluxcorrector
         phiN=phiN*fluxcorrector
         piN=piN*fluxcorrector
         pipiN=pipiN*fluxcorrector
      end if

      !###################################################################################################
      ! Sum up everything for the total cross section
      !###################################################################################################

      sigmaTot=phiN + sum( piN ) + pipiN

    end subroutine evaluateXsections


    !333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
    !333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
    !333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333


    subroutine makeDecision
      use random, only: rn, ranCharge

      real :: cut!,cut2
      integer :: totalCharge!,resID

      integer,dimension (1:3)  :: izmin,izmax,izout     ! needed for ranCharge
      logical :: ranChargeFlag

      integer :: pionCharge

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=phi_particle%Charge+nucleon_particle%Charge

      !############################################################
      ! (1) Two -body final states
      !############################################################

      ! phi N production
      If(phiN.ge.cut) then
         teilchenOut(1)%Id=phi
         teilchenOut(2)%Id=nucleon

         teilchenOut(1)%Charge=0
         teilchenOut(2)%Charge=totalCharge
         return
      end if
      cut=cut-phiN

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

      !############################################################
      ! (2) Three body final state
      !############################################################
      ! pi pi N production

      If(pipiN.ge.cut) then
         teilchenOut(1)%Id=pion
         teilchenOut(2)%Id=pion
         teilchenOut(3)%Id=nucleon
         izmin=(/-1,-1,0/)
         izmax=(/1,1,1/)
         call rancharge(izmin,izmax,totalCharge,izout,ranChargeFlag)
         If (.not.ranChargeFlag) write(*,*) 'Error in rancharge :',izmin,izmax,totalCharge,izout,ranChargeFlag
         teilchenOut(1)%Charge=izout(1)
         teilchenOut(2)%Charge=izout(2)
         teilchenOut(3)%Charge=izout(3)
         return
      end if

      ! No event was generated:
      write(*,*) 'Error in makedecision of phiNuc'
      write(*,*) srts
      write(*,*) sigmatot
      write(*,*) phiN
      write(*,*) piN
      write(*,*) pipiN
      write(*,*) cut
      stop

    end subroutine makeDecision


    !*************************************************************************************************
    !****s* phiNuc/makeOutput
    ! NAME
    ! subroutine makeOutput
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV].
    ! Filenames:
    ! * 'phiN_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'phiN_nonStrange_nuk.dat'     : non-strange meson with nucleon in final state
    !*************************************************************************************************
    subroutine makeOutPut

      logical, save :: initFlag=.true.

      ! The output files
      character(30), dimension(1:2), parameter :: outputFile = (/ 'phiN_sigTotElast.dat   ', 'phiN_nonStrange_nuk.dat' /)
      real :: plab

      plab=SQRT(((srts**2-phi_particle%mass**2-nucleon_particle%mass**2)/2./nucleon_particle%mass)**2-phi_particle%mass**2)

      If (initFlag) then
         Open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         Open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         write (101,*) '# srts, plab, sigmaTot, sigmaElast '
         write (102,*) '# srts, plab, piN(-1:1), phiN, pipiN '
         initFlag=.false.
      else
         Open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
      end If
      write (101,'(4F9.3)') srts, plab,sigmaTot, sigmaElast
      write (102,'(7F9.3)') srts, plab,piN, phiN, pipiN
      Close(101)
      Close(102)
    end subroutine makeOutPut

  end subroutine phiNuc


end module phiNucleon
