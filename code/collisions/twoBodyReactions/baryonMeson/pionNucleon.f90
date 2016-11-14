!****************************************************************************
!****m* /pionNucleon
! NAME
! module pionNucleon
! PURPOSE
! Includes the cross sections for pion-nucleon scattering in the resonance regime.
! Implemented are the following reactions:
! * pion nucleon -> X
! Public routines:
! * pionNuc
! Public variables:
! * matrixDeltaEta
!****************************************************************************
module pionNucleon
  implicit none
  Private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  !****************************************************************************
  !****g* pionNucleon/matrixDeltaEta
  ! SOURCE
  !
  real, parameter :: matrixDeltaEta=7.
  ! PURPOSE
  ! Matrix Element for pi N -> eta Delta ; see old "commcoll"->matdeta
  !****************************************************************************

  Public :: pionNuc, matrixDeltaEta

contains


  !****************************************************************************
  !****s* pionNucleon/pionNuc
  ! NAME
  ! subroutine pionNuc (srts, teilchenIn, mediumATcollision, momentumLRF, teilchenOut, sigmaTot, sigmaElast, plotFlag)
  !
  ! PURPOSE
  ! Evaluates pion Nucleon -> anything cross sections and returns also a "preevent"
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross
  ! sections to decide on a final state with maximal 3 (or 4) final state particles.
  ! These are returned in the vector teilchenOut.
  ! The kinematics of these teilchen is only fixed in the case of a single
  ! produced resonance.
  ! Otherwise the kinematics still need to be established.
  !
  ! INPUTS
  ! * real                          :: srts        -- sqrt(s) in the process (free one!!)
  ! * type(particle),dimension(1:2) :: teilchenIn  -- colliding particles
  ! * type(medium) :: mediumATcollision  --  Medium information at the position of the collision
  ! * real, dimension(0:3) :: momentumLRF -- total momentum in LRF
  ! * logical, OPTIONAL :: plotFlag -- Switch on plotting of the  Xsections
  !
  ! RESULT
  ! * real    :: sigmaTot         -- total Xsection
  ! * real    :: sigmaElast       -- elastic Xsection
  ! * type(preEvent),dimension(1:4) :: teilchenOut -- colliding particles
  !
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : baryon Resonances
  ! * 2-particle : pi N, omega N, phi N, eta Delta,  Kaon Lambda , Sigma Kaon,
  !   pion Delta, rho N, rhoDelta, sigma N, pion P11_1440, rho Delta
  ! * 3-particle : omega pi N, phi pi N, K KBar N, Lambda Kaon Pion, Sigma Kaon Pion
  ! * Not yet implemented : pi N -> pi pi N which is used for the matching to
  !   the high energy region
  !****************************************************************************
  subroutine pionNuc (srts, teilchenIn, mediumATcollision, momentumLRF, teilchenOut, sigmaTot, sigmaElast, plotFlag)
    use monteCarlo, only: MonteCarloChoose2Dim
    use idTable
    use particleDefinition
    use mediumDefinition, only: medium, vacuum
    use preEventDefinition, only : preEvent
    use twoBodyTools, only : velocity_correction, convertToAntiParticles, searchInInput
    use RMF, only : getRMF_flag

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real, intent(in)                            :: srts
    type(particle), dimension(1:2), intent(in)  :: teilchenIn
    type(medium), intent(in)                    :: mediumATcollision
    real, intent(in), dimension(0:3)            :: momentumLRF
    type(preEvent), dimension(1:4), intent(out) :: teilchenOut
    real, intent(out)                           :: sigmaTot, sigmaElast
    logical, intent(in), optional               :: plotFlag

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! cross sections named according to final states
    real :: omegaN, phiN, etaDelta                                         ! exactly one possible final state
    real :: omegaPiN, phiPiN                                               ! same Xsection for all isospin channels
    real :: lambdaKaon, sigmaN, etaN, piPiN, piPiPiN
    real, dimension (-1:1)  :: piN, pionDelta, piP11_1440, rhoN, rhoDelta  ! index = charge of outgoing meson
    real, dimension  (0:1)  :: sigmaKaon                                   ! index = charge of final state kaon
    real, dimension (-2:1)  :: kaonKaonBarN                                ! index= final state
    real, dimension (0:1,-1:1) :: lambdaKaonPion, SigmaKaonPion            ! 1st index: kaon charge, 2nd index: pion charge
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    logical      :: antiParticleINPUT     ! set to .true. if antiparticle in the input
    real         :: fluxCorrector         ! Correction due to incoming flux factor
    real         :: pionMass, nukMass
    integer      :: pionCharge, nukCharge
    logical      :: failFlag
    type(particle) :: pion_particle, nucleon_particle
    real, dimension(Delta:nbar) :: sigmaRes   ! partial cross sections for pion N -> R
    real, dimension(Delta:nbar) :: massRes    ! Resonance masses

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    teilchenOut(:)%ID=0                    ! ID of produced particles
    teilchenOut(:)%charge=0                ! Charge of produced particles
    teilchenOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    teilchenOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(teilchenIn,pion,nucleon,pion_particle,nucleon_particle,failFlag)
    If (failFlag) then
       Write(*,*) 'Wrong input in pionnuc', teilchenIn%ID
       stop
    end if

    pionCharge = pion_particle%charge
    pionMass   = pion_particle%mass
    nukCharge  = nucleon_particle%charge
    nukMass    = nucleon_particle%mass

    If (pion_particle%antiParticle) then
       ! This case should be impossible
       write(*,*) 'pion is antiparticle in "pionnuc"!!!',teilchenIN%ID,teilchenIN%antiparticle
       stop
    end if

    If (nucleon_particle%antiParticle) then ! Invert all particles in antiparticles
       nukCharge = -nukcharge
       pionCharge= -pionCharge
       antiParticleInput=.true.
    else
       antiParticleInput=.false.
    end if


    ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
    if (.not. getRMF_flag()) then
       fluxCorrector=velocity_correction(teilchenIn)
    else
       fluxCorrector=1.
    end if

    ! (2) Evaluate the cross sections
    call evaluateXsections

    ! Cutoff to kick the case out, that the cross section is zero
    if (sigmaTot<1E-12) then
       sigmatot=0.
       sigmaElast=0.
       return
    end if


    ! (3) Plot them if wished
    If (Present(PlotFlag).or.debugFlag) then
       If (plotFlag.or.debugFlag)  call makeOutput
    end if

    ! (4) Define final state
    call MakeDecision

    ! (5) Check Output
    If (Sum(teilchenOut(:)%Charge).ne.pionCharge+nukCharge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error' ,pionCharge, &
            & nukCharge, teilchenOut(:)%Charge,teilchenOut(:)%ID
       stop
    end if

    ! (6) Re-Invert particles if antiparticles in input
    If (antiParticleInput) then
       if (debugFlagAnti) write(*,*) teilchenOut
       call convertToAntiParticles(teilchenOut)
       if (debugFlagAnti) write(*,*) teilchenOut
    end if

  contains

    !***********************************************************************
    !****s* pionNuc/evaluateXsections
    ! NAME
    ! subroutine evaluateXsections
    ! PURPOSE
    ! Evaluates all Xsections for pi N scattering
    !***********************************************************************
    subroutine evaluateXsections
      use particleProperties, only : hadron
      use parametrizationsBarMes, only: piN_elastic, piN_chargeExchange, huangLam, sibirtpi, golub_omega, golub_phi, &
                                        huang, pin_to_strangebaryon_kaon_pion_matrix
      use parBarMes_HighEnergy, only : paramBarMesHE_pion
      use resonanceCrossSections, only: barMes_R_barMes, barMes2resonance
      use clebschGordan, only: clebschSquared
      use constants, only: mN, mPi, mK
      use twoBodyTools, only: pCM
      use twoBodyPhaseSpace, only: Integrate_2bodyPS_resonance

      real, dimension(1:4) :: sigmaGolub, sigmaHuang
      real :: elastic, sigmaHuangLam, sibirtpiCross, sigma_pitot, sigma_HiEl, isoZ_nuk
      logical :: background, propagated, perturbative, successElastic
      integer :: deltaCharge,rhoCharge,new_pionCharge,new_nukCharge,new_deltaCharge,new_p11Charge
      real, dimension(0:3) :: momentum_vacuum        ! Total Momentum in vacuum
      real, dimension(1:3) :: position
      real, dimension(-1:1) :: sigma_total_param, sigma_elast_param
      real, dimension(1:5) :: ps

      position=0.5*(teilchenIN(1)%position+teilchenIN(2)%position)

      perturbative = (teilchenIN(1)%perturbative.or.teilchenIN(2)%perturbative)

      momentum_vacuum(1:3)=teilchenIn(1)%momentum(1:3)+teilchenIn(2)%momentum(1:3)
      momentum_vacuum(0)=FreeEnergy(teilchenIn(1))+FreeEnergy(teilchenIn(2))

      isoZ_nuk=nukCharge-0.5

      !##############################################################
      ! Evaluate partial cross sections
      !##############################################################


      !**************************************************************
      ! pi N -> pi N
      ! Index : (-1:1)  Charge of outgoing pion
      !**************************************************************

      piN=0.

      ! (1) contribution due to non-propagated resonances
      background=.true.
      propagated=.false.

      Do new_pionCharge=-1,1
         new_nukCharge=pionCharge+nukCharge-new_pionCharge
         if (DebugFlag) write(*,*) 'charges=' , new_nukCharge,new_PionCharge
         if ((new_NukCharge==1).or.(new_NukCharge==0)) then
            piN(new_PionCharge) = barMes_R_barMes(pion,nucleon,pion,nucleon,pionCharge,nukCharge, &
                                                  new_pionCharge,new_nukCharge,background,propagated,Vacuum,momentum_vacuum, &
                                                  pionMass,nukMass,position,perturbative,srts)
            if (DebugFlag) write(*,*) piN(new_PionCharge)
         end if
      End do

      ! Elastic : Subtract resonance contribution from data
      if (pionCharge /= 0) then
         If (nukcharge == 0) then
            elastic = piN_elastic(-pionCharge,srts,successElastic)
            if (successElastic) &
               piN(pionCharge) = max(0.,elastic-barMes_R_barMes(pion,nucleon,pion,nucleon, &
                                                                pionCharge,nukcharge,pionCharge,nukcharge,.false.,.true., &
                                                                Vacuum,momentum_vacuum,pionMass,nukMass,position,perturbative,srts))
         else if (nukcharge == 1) then
            elastic = piN_elastic(pionCharge,srts,successElastic)
            if (successElastic) &
               piN(pionCharge) = max(0.,elastic-barMes_R_barMes(pion,nucleon,pion,nucleon, &
                                                                pionCharge,nukcharge,pionCharge,nukcharge,.false.,.true., &
                                                                Vacuum,momentum_vacuum,pionMass,nukMass,position,perturbative,srts))
         else
            stop 'error in pionNUC'
         end if
      else
        ! pi0 = (pi+ + pi-)/2
        elastic = 0.5 * (piN_elastic(+1,srts,successElastic) + piN_elastic(-1,srts,successElastic))
        if (successElastic) &
          piN(pionCharge) = max(0.,elastic-barMes_R_barMes(pion,nucleon,pion,nucleon, &
                                                           pionCharge,nukcharge,pionCharge,nukcharge,.false.,.true., &
                                                           Vacuum,momentum_vacuum, pionMass,nukMass,position,perturbative,srts))
      end if

      sigmaElast = barMes_R_barMes(pion,nucleon,pion,nucleon,&
                                   pionCharge,nukCharge,pionCharge,nukCharge, .false.,.false.,&
                                   MediumAtCollision,momentumLRF,pionMass,nukMass,position,perturbative,srts) &
                   + piN(pionCharge)


      ! Charge Exchange : Subtract resonance contribution from data
      If (nukcharge==0 .and. pionCharge==0) then       ! pi n -> pi p
         piN(pionCharge-1)= max(0.,piN_chargeExchange(srts)-barMes_R_barMes(pion,nucleon,pion,nucleon, &
                            pionCharge,nukcharge,pionCharge-1,nukcharge+1,.false.,.true.,Vacuum,momentum_vacuum,&
                            pionMass,nukMass,position,perturbative,srts) )
      else if (nukcharge==1 .and. pionCharge<=0) then ! pi  p -> pi n
         piN(pionCharge+1)= max(0.,piN_chargeExchange(srts)-barMes_R_barMes(pion,nucleon,pion,nucleon, &
                            pionCharge,nukcharge,pionCharge+1,nukcharge-1,.false.,.true.,Vacuum,momentum_vacuum,&
                            pionMass,nukMass,position,perturbative,srts))
      end if

      !        write(*,*) 'piN mit background' , piN

      !**************************************************************
      ! -> omega N
      !**************************************************************
      sigmaGolub(1:2) = golub_omega (srts)     ! Golub returns the cross section for pi- p -> omega n
      ! There is always a I=1/2 resonance as intermediate state :
      omegaN=clebschSquared(1.,0.5,0.5,real(pionCharge),isoZ_nuk) &
             * 3./2.* (sigmaGolub(1)-barMes_R_barMes(pion,nucleon,omegaMeson,nucleon,-1,1,0,0,.false.,.true.,Vacuum, &
                                                     momentum_vacuum,pionMass,nukMass,position,perturbative,srts))
      ! factor 3./2. to divide by the isospin clebsch for pi- p -> omega n
      If (omegaN<0.) then
         omegaN=0.
         If(debugFlag) Write(*,*) 'Problem in pionNuc : pion N -> Omega N &
              &    resonance Contribution greater than elementary crossSection'
      end if

      ! -> omega pi N
      ! Assume same cross section for all channels
      omegaPiN=sigmaGolub(2)

      !**************************************************************
      ! -> phi N
      !**************************************************************
      sigmaGolub = golub_phi (srts, 0.)   ! Golub returns the cross section for pi- p -> phi n
      phiN = clebschSquared(1.,0.5,0.5,real(pionCharge),isoZ_nuk) * 3./2. * sigmaGolub(1)

      !**************************************************************
      ! -> phi  pi N
      !**************************************************************
      ! Assume same cross section for all channels
      phiPiN=sigmaGolub(2)

      !**************************************************************
      ! -> Kaon KaonBar N
      !**************************************************************
      kaonKaonBarN(:)=0.
      sibirtpiCross = sibirtpi(srts)
      ! pion proton ->Kaon KaonBar N
      ! Index :
      ! 1  outgoing K^+   kbar^0                 ! =-1
      ! 0  outgoing K^0   kbar^0                 ! =-2
      ! -1 outgoing K^0  Kbar^-                  ! =1
      ! -2 outgoing K^+  kbar^-                  ! =0
      If (nukCharge==1) then
         Select Case(pionCharge)
         Case(-1)
            kaonKaonBarN(1)=0.
            kaonKaonBarN(0)=sibirtpiCross
            kaonKaonBarN(-1)=sibirtpiCross/2.
            kaonKaonBarN(-2)=sibirtpiCross
         Case(0)
            kaonKaonBarN(1)=sibirtpiCross
            kaonKaonBarN(0)=sibirtpiCross/4.
            kaonKaonBarN(-1)=0.
            kaonKaonBarN(-2)=sibirtpiCross/4.
         Case(1)
            kaonKaonBarN(1)=sibirtpiCross/2.
            kaonKaonBarN(0)=0.
            kaonKaonBarN(-1)=0.
            kaonKaonBarN(-2)=0.
         End select
      else
         ! Change piMinus <-> piPlus and outgoing channels acording to isospin rotation scheme above
         ! 1  outgoing K^+   kbar^0    =>  -1 after iso-rotation
         ! 0  outgoing K^0   kbar^0    => -2
         ! -1 outgoing K^0  Kbar^-     =>  1
         ! -2 outgoing K^+  kbar^-      => 0
         Select Case(pionCharge)
         Case(1)
            kaonKaonBarN(-1)=0.
            kaonKaonBarN(-2)=sibirtpiCross
            kaonKaonBarN(1)=sibirtpiCross/2.
            kaonKaonBarN(0)=sibirtpiCross
         Case(0)
            kaonKaonBarN(-1)=sibirtpiCross
            kaonKaonBarN(-2)=sibirtpiCross/4.
            kaonKaonBarN(1)=0.
            kaonKaonBarN(0)=sibirtpiCross/4.
         Case(-1)
            kaonKaonBarN(-1)=sibirtpiCross/2.
            kaonKaonBarN(-2)=0.
            kaonKaonBarN(1)=0.
            kaonKaonBarN(0)=0.
         End select
      end if


      !**************************************************************
      ! -> Lambda Kaon
      !**************************************************************
      lambdaKaon=0.
      If (srts>(hadron(Lambda)%mass+mK)) then
         ! hunaglam gives : pi^{-} p -> Lambda Kaon^{0}
         sigmaHuangLam = huangLam (srts)
         Select Case(pionCharge)
         Case(-1)
            If (nukCharge==1) &
               lambdaKaon=sigmaHuangLam-barMes_R_barMes(pion,nucleon,kaon,Lambda, &
                                                        -1,1,0,0,.false.,.true.,Vacuum,momentum_vacuum, &
                                                        pionMass,nukMass,position,perturbative,srts)
         Case(0)
            lambdaKaon=0.5* (sigmaHuangLam-barMes_R_barMes(pion,nucleon,kaon,Lambda, &
                                                           -1,1,0,0,.false.,.true.,Vacuum,momentum_vacuum, &
                                                           pionMass,nukMass,position,perturbative,srts) )
         Case(1)
            If (nukCharge==0) &
               lambdaKaon=sigmaHuangLam-barMes_R_barMes(pion,nucleon,kaon,Lambda, &
                                                        -1,1,0,0,.false.,.true.,Vacuum,momentum_vacuum, &
                                                        pionMass,nukMass,position,perturbative,srts)
         End Select

         If (lambdaKaon<0.) then
            lambdaKaon=0.
            Write(*,*) 'Problem in pionNuc : pion N -> Lambda Kaon resonance', &
                       '  Contribution greater than elementary crossSection'
         end if
      end if

      !**************************************************************
      ! -> Sigma Kaon
      !**************************************************************
      ! sigmaKaon(0:1) : Index is charge of final state kaon
      ! sigmaHuang(1) = pi^{+}  p  ->   K^{+}  Sigma+      !
      ! sigmaHuang(2) = pi^{0}  p  ->   K^{+}  Sigma0       !
      ! sigmaHuang(3) = pi^{-}  p  ->   K^{0}  Sigma0      !
      ! sigmaHuang(4) = pi^{-}  p  ->   K^{+}   Sigma-     !
      sigmaKaon(:)=0.
      If (srts>(hadron(SigmaResonance)%mass+mK)) then
         sigmaHuang = huang(srts)
         If (nukCharge==1) then
            Select Case(pionCharge)
            Case(-1)
               sigmaKaon(1)=sigmaHuang(4)
               sigmaKaon(0)=sigmaHuang(3)
            Case(0)
               sigmaKaon(1)=sigmaHuang(2)
               sigmaKaon(0)=sigmaHuang(3)   ! by isospin consideration
            Case(1)
               sigmaKaon(1)=sigmaHuang(1)
               sigmaKaon(0)=0.
            End select
         else ! neutron Xsections by charge conjugation
            Select Case(pionCharge)
            Case(-1)
               sigmaKaon(1)=0.
               sigmaKaon(0)=sigmaHuang(1)
            Case(0)
               sigmaKaon(1)=sigmaHuang(3)
               sigmaKaon(0)=sigmaHuang(2)
            Case(1)
               sigmaKaon(1)=sigmaHuang(3)
               sigmaKaon(0)=sigmaHuang(4)
            End select
         end if
      end if

      !**************************************************************
      ! -> eta Delta
      !**************************************************************
      etaDelta=0.
      If (srts>(hadron(Delta)%minmass+hadron(eta)%mass)) then
        ps = Integrate_2bodyPS_resonance (Delta, srts, hadron(eta)%mass, 0.)
        etaDelta = matrixDeltaEta * ps(1) / (pCM(srts,mN,mPi)*srts**2) &
                      * clebschSquared(1.,0.5,1.5,real(pionCharge),isoZ_nuk)
      end if

      !**************************************************************
      ! -> pion Delta
      !**************************************************************
      pionDelta=0.
      If (srts>(hadron(Delta)%minmass+mPi)) then
         Do new_PionCharge=-1,1
            deltaCharge=pionCharge+nukCharge-new_PionCharge
            If ((deltaCharge>=-1).and.(deltaCharge<=2)) then
               pionDelta(new_PionCharge)=barMes_R_barMes(pion,nucleon,pion,Delta,pionCharge,nukCharge,&
                                                         new_PionCharge,deltaCharge,background,propagated,Vacuum,momentum_vacuum,&
                                                         pionMass,nukMass,position,perturbative,srts)
            end If
         End Do
      End If

      !**************************************************************
      ! -> rho N
      !**************************************************************
      rhoN=0.
      If (srts>(hadron(rho)%minmass+mN)) then
         Do new_Nukcharge=0,1
            rhoCharge=pionCharge+nukCharge-new_NukCharge
            If ((rhoCharge>=-1).and.(rhoCharge<=1)) then
               rhoN(rhoCharge)=barMes_R_barMes(pion,nucleon,rho,nucleon,pionCharge,nukCharge,rhoCharge, &
                                               new_nukCharge,background,propagated,Vacuum,momentum_vacuum, &
                                               pionMass,nukMass,position,perturbative,srts)
            end If
         end Do
      end If

      !**************************************************************
      !  -> rhoDelta
      !**************************************************************
      ! Evaluate backGrounds by subtracting Resonance cross section
      ! of propagated resonances off the cross section for all resonances:
      rhoDelta=0.
      If (srts>(hadron(delta)%minmass+hadron(rho)%minmass)) then
         Do new_DeltaCharge=-1,2
            rhoCharge=pionCharge+nukCharge-new_DeltaCharge
            If ((rhoCharge>=-1).and.(rhoCharge<=1)) then
               rhoDelta(rhoCharge)=barMes_R_barMes(pion,nucleon,rho,Delta,pionCharge,nukCharge,&
                                                   rhoCharge,new_DeltaCharge,background,propagated,Vacuum,momentum_vacuum,&
                                                   pionMass,nukMass,position,perturbative,srts)
            end if
         end do
      end if

      !**************************************************************
      ! -> pion P11_1440
      !**************************************************************
      piP11_1440=0.
      If (srts>(hadron(P11_1440)%minmass+mPi)) then
         Do new_p11Charge=0,1
            new_pionCharge=pionCharge+nukCharge-new_p11Charge
            If ((new_pionCharge>=-1).and.(new_pionCharge<=1)) then
               piP11_1440(new_pionCharge)=barMes_R_barMes(pion,nucleon,pion,P11_1440,pionCharge,nukCharge, &
                                          new_pionCharge,new_p11Charge,background,propagated,Vacuum,momentum_vacuum, &
                                          pionMass,nukMass,position,perturbative,srts)
            end If
         end Do
      end If

      !**************************************************************
      ! -> sigma N
      !**************************************************************
      sigmaN=0.
      If ((nukCharge+pionCharge==0).or.(nukCharge+pionCharge==1)) then
         If (srts>(hadron(sigmaMeson)%minmass+mN)) then
            sigmaN=barMes_R_barMes(pion,nucleon,sigmaMeson,nucleon,pionCharge,nukCharge, &
                                   0,nukCharge+pionCharge,background,propagated,Vacuum,momentum_vacuum, &
                                   pionMass,nukMass,position,perturbative,srts)
         end If
      end If

      !**************************************************************
      ! -> eta N
      !**************************************************************
      etaN=0.
      If ((nukCharge+pionCharge==0).or.(nukCharge+pionCharge==1)) then
         If (srts>(hadron(eta)%mass+mN)) then
            etaN=barMes_R_barMes(pion,nucleon,eta,nucleon,pionCharge,nukCharge, &
                                 0,nukCharge+pionCharge,background,propagated,Vacuum,momentum_vacuum, &
                                 pionMass,nukMass,position,perturbative,srts)
         end If
      end If

      !**************************************************************
      ! pi N -> R
      !**************************************************************
      sigmaRes = barMes2resonance (pion,nucleon,pionCharge,nukCharge,.true.,mediumAtCollision, &
                                   momentumLRF,massRes,pionMass,nukMass,position,perturbative,srts)

      !**************************************************************
      ! pi N -> Sigma Kaon Pion
      ! pi N -> Lambda Kaon Pion
      !**************************************************************
      call piN_to_strangeBaryon_kaon_pion_matrix(srtS,pionCharge,nukCharge, LambdaKaonPion, SigmaKaonPion)

      !**************************************************************
      ! -> pi pi N, pi pi pi N
      !**************************************************************
      ! We use these channels to put missing strength of the Xsection
      ! at higher energies into this channel.
      ! But first we put missing strength into the elastic channel.

      piPiN  =0.
      piPiPiN=0.

      sigma_HiEl=0. ! This is the hi energy elastic contribution

      If (srts>1.75) then

         call paramBarMesHE_pion(srts,sigma_total_param,sigma_elast_param)

         if (nukCharge==1) then
            sigma_piTot=sigma_total_param( pionCharge)
            sigma_HiEl =sigma_elast_param( pionCharge) ! in vacuum !!!
         else
            sigma_piTot=sigma_total_param(-pionCharge)
            sigma_HiEl =sigma_elast_param(-pionCharge) ! in vacuum !!!
         end if

         sigma_HiEl = max(0.0, sigma_HiEl - barMes_R_barMes(pion,nucleon,pion,nucleon,pionCharge,nukCharge, &
                                                            pionCharge,nukCharge,.false.,.false.,Vacuum,momentum_vacuum, &
                                                            pionMass,nukMass,position,perturbative,srts) & ! in vacuum !!!
                               - piN(pionCharge))

         piN(pionCharge)= piN(pionCharge)+ sigma_HiEl
         sigmaElast     = sigmaElast     + sigma_HiEl


         ! pipiN=sigTot-all partial cross sections evaluated before
         pipiN = max(0.,sigma_piTot - (omegaN +  phiN + etaDelta  +   Sum(piN) +lambdaKaon +Sum ( sigmaKaon ) + &
                 sum ( kaonKaonBarN ) + &
                 sum (pionDelta) +sum ( piP11_1440 ) + sum (  rhoN ) + sum (rhoDelta ) + &
                 sigmaN + etaN + sum (sigmaRes ) + omegaPiN + phiPiN  + &
                 sum( sigmaKaonPion  ) + sum( LambdaKaonPion )))

         if (srts<2.1) then
            if (srts<1.95) then
               piPiN = 0.0
            else
               piPiN=piPiN*(srts-1.95)/0.15 ! Smooth transition
            endif
         end if

         !-----------------------------------------
         ! now we select between 2 pion and 3 pion:
         !-----------------------------------------
         !
         pipipiN = 0.66 * pipiN
         pipiN = pipiN - pipipiN
         !
         !-----------------------------------------

      end if


      !##############################################################
      ! Do the flux correction for each channel
      !##############################################################
      ! We do this for each channel since they might show up
      ! separately in the output, if makeoutput is called.
      If (fluxCorrector_flag) then
         sigmaElast=sigmaElast*fluxcorrector
         omegaN=omegaN*fluxcorrector
         phiN=phiN *fluxcorrector
         etaDelta=etaDelta*fluxcorrector
         piN=piN*fluxcorrector
         lambdaKaon=lambdaKaon *fluxcorrector
         sigmaKaon=sigmaKaon*fluxcorrector
         kaonKaonBarN=kaonKaonBarN *fluxcorrector
         pionDelta=pionDelta*fluxcorrector
         piP11_1440=piP11_1440 *fluxcorrector
         rhoN=rhoN*fluxcorrector
         rhoDelta=rhoDelta*fluxcorrector
         sigmaN=sigmaN *fluxcorrector
         etaN=etaN *fluxcorrector
         sigmaRes=sigmaRes *fluxcorrector
         omegaPiN=omegaPiN*fluxcorrector
         phiPiN=phiPiN  *fluxcorrector
         LambdaKaonPion=LambdaKaonPion * fluxcorrector
         SigmaKaonPion =SigmaKaonPion  * fluxcorrector
         pipiN=pipiN*fluxcorrector
         pipipiN=pipipiN*fluxcorrector
      end if

      !##############################################################
      ! Sum up everything for the total cross section
      !##############################################################
      ! Be careful since sigma elast is already included in the
      ! partial cross sections, therefore it is not included in the
      ! total cross section

      sigmaTot = omegaN + phiN + etaDelta + Sum(piN) + lambdaKaon &
               + sum(sigmaKaon ) + sum (kaonKaonBarN) &
               + sum(pionDelta)  + sum(piP11_1440) &
               + sum(rhoN) + sum(rhoDelta) &
               + sigmaN + etaN + sum (sigmaRes ) + omegaPiN + phiPiN  &
               + sum(LambdaKaonPion) + sum(SigmaKaonPion) + pipiN + pipipiN


    end subroutine evaluateXsections


    !****************************************************************************
    !****s* pionNuc/makeDecision
    ! NAME
    ! subroutine MakeDecision
    ! PURPOSE
    ! Decides on the final state which is returned via teilchenOut by Monte-Carlo.
    !  * Assigns charges and ID's.
    !  * Only for resonance-production also the mass is assigned, since the mass
    !    of the resonance needed to be calculated earlier.
    ! The Monte-Carlo routine is adding up channels until the sum is exceeding
    ! x*sigma(total). x has a flat distribution in [0,1].
    ! The last added channel is then the one which is chosen for the event.
    !****************************************************************************
    subroutine MakeDecision
      use random, only: rn, ranCharge

      real :: summe,cut,cut2
      integer :: totalCharge,charge,resID
      integer,dimension (1:3)  :: izmin,izmax,izout     ! needed for ranCharge
      integer,dimension (1:4)  :: izmin4,izmax4,izout4     ! needed for ranCharge
      logical :: ranChargeFlag
      integer, dimension (1:2) :: channel_index

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=pionCharge+nukCharge

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

      If (omegaN>=cut) then
         teilchenOut(1:2)%Id = (/ omegaMeson, nucleon /)
         teilchenOut(1:2)%Charge = (/ 0, totalCharge /)
         return
      end if
      cut=cut-omegaN

      if (phiN>=cut) then
         teilchenOut(1:2)%Id = (/ phi, nucleon /)
         teilchenOut(1:2)%Charge = (/ 0, totalCharge /)
         return
      end if
      cut=cut-phiN

      If (etaDelta>=cut) then
         teilchenOut(1:2)%Id = (/ eta, Delta /)
         teilchenOut(1:2)%Charge = (/ 0, totalCharge /)
         return
      end if
      cut=cut-etaDelta

      If (Sum(piN)>=cut) then
         teilchenOut(1:2)%Id = (/ pion, nucleon /)
         cut2=rn()*sum(piN)
         Do charge=-1,1
            If(sum(piN(-1:charge)).gt.cut2) then
               teilchenOut(1:2)%Charge = (/ charge, totalCharge-charge /)
               exit
            end if
         End do
         return
      end if
      cut=cut-Sum(piN)

      If (lambdaKaon>=cut) then
         teilchenOut(1:2)%Id = (/ kaon, Lambda /)
         teilchenOut(1:2)%Charge = (/ totalCharge, 0 /)
         return
      end if
      cut=cut-lambdaKaon

      If (sum(sigmaKaon)>=cut) then
         teilchenOut(1:2)%Id = (/ kaon, sigmaResonance /)
         cut2=rn()*sum(sigmaKaon)
         Do charge=0,1
            If (sum(sigmaKaon(0:charge))>=cut2) then
               teilchenOut(1:2)%Charge = (/ charge, totalCharge-charge /)
               exit
            end if
         End do
         return
      end if
      cut=cut-sum(sigmaKaon)

      If (sum(pionDelta)>=cut) then
         teilchenOut(1:2)%Id = (/ pion, Delta /)
         cut2=rn()*sum(pionDelta)
         Do charge=-1,1
            If (sum(pionDelta(-1:charge))>=cut2) then
               teilchenOut(1:2)%Charge = (/ charge, totalCharge-charge /)
               exit
            end if
         End do
         return
      end if
      cut=cut-sum(pionDelta)

      If (sum(piP11_1440)>=cut) then
         teilchenOut(1:2)%Id = (/ pion, P11_1440 /)
         cut2=rn()*sum(piP11_1440)
         Do charge=-1,1
            If (sum(piP11_1440(-1:charge))>=cut2) then
               teilchenOut(1:2)%Charge = (/ charge, totalCharge-charge /)
               exit
            end if
         End do
         return
      end if
      cut=cut-sum(piP11_1440)

      If (sum(rhoN)>=cut) then
         teilchenOut(1:2)%Id = (/ rho, nucleon /)
         cut2=rn()*sum(rhoN)
         Do charge=-1,1
            If (sum(rhoN(-1:charge))>=cut2) then
               teilchenOut(1:2)%Charge = (/ charge, totalCharge-charge /)
               exit
            end if
         End do
         return
      end if
      cut=cut-sum(rhoN)

      If (sum(rhoDelta)>=cut) then
         teilchenOut(1:2)%Id = (/ rho, Delta /)
         cut2=rn()*sum(rhoDelta)
         Do charge=-1,1
            If (sum(rhoDelta(-1:charge))>=cut2) then
               teilchenOut(1:2)%Charge = (/ charge, totalCharge-charge /)
               exit
            end if
         End do
         return
      end if
      cut=cut-sum(rhoDelta)

      If (sigmaN>=cut) then
         teilchenOut(1:2)%Id = (/ sigmaMeson, nucleon /)
         teilchenOut(1:2)%Charge = (/ 0, totalCharge /)
         return
      end if
      cut=cut-sigmaN

      If (etaN>=cut) then
         teilchenOut(1:2)%Id     = (/ eta, nucleon /)
         teilchenOut(1:2)%Charge = (/ 0, totalCharge /)
         return
      end if
      cut=cut-etaN


      !############################################################
      ! (3) Three body final state
      !############################################################

      If (omegaPiN>=cut) then
         teilchenOut(1:3)%Id = (/ omegaMeson, pion, nucleon /)
         izmin=(/0,-1,0/)
         izmax=(/0,1,1/)
         call rancharge(izmin,izmax,totalCharge,izout,ranChargeFlag)
         If (.not.ranChargeFlag) write(*,*) 'Error in rancharge :',izmin,izmax,totalCharge,izout,ranChargeFlag
         teilchenOut(1:3)%Charge = izout(1:3)
         return
      end if
      cut=cut-omegaPiN

      if (phiPiN>=cut) then
         teilchenOut(1:3)%Id = (/ phi, pion, nucleon /)
         izmin=(/0,-1,0/)
         izmax=(/0,1,1/)
         call rancharge(izmin,izmax,totalCharge,izout,ranChargeFlag)
         If (.not.ranChargeFlag) write(*,*) 'Error in rancharge :',izmin,izmax,totalCharge,izout,ranChargeFlag
         teilchenOut(1:3)%Charge=izout(1:3)
         return
      end if
      cut=cut-phiPiN

      if (piPiN>=cut) then
         teilchenOut(1:3)%Id = (/ pion, pion, nucleon /)
         izmin=(/-1,-1,0/)
         izmax=(/1,1,1/)
         call rancharge(izmin,izmax,totalCharge,izout,ranChargeFlag)
         If (.not.ranChargeFlag) write(*,*) 'Error in rancharge :',izmin,izmax,totalCharge,izout,ranChargeFlag
         teilchenOut(1:3)%Charge=izout(1:3)
         return
      end if
      cut=cut-piPiN

      If (sum(kaonKaonBarN)>=cut) then
         teilchenOut(1:3)%Id = (/ kaon, kaonBar, nucleon /)
         cut2=rn()*sum(kaonKaonBarN)
         Do charge=-2,1
            If (sum(kaonKaonBarN(-2:charge))>=cut2) then
               Select case(charge)
               case(-2)  ! outgoing K^+  kbar^-
                  teilchenOut(1:3)%Charge = (/ 1, -1, totalCharge /)
               case(-1)  ! outgoing K^0  Kbar^-
                  teilchenOut(1:3)%Charge = (/ 0, -1, totalCharge+1 /)
               case(0)  ! outgoing K^0   kbar^0
                  teilchenOut(1:3)%Charge = (/ 0, 0, totalCharge /)
               case(1)  ! outgoing K^+   kbar^0
                  teilchenOut(1:3)%Charge = (/ 1, 0, totalCharge-1 /)
               End Select
               exit
            end if
         End do
         return
      end if
      cut=cut-sum(kaonKaonBarN)

      If (sum(SigmaKaonPion)>=cut) then
         ! pi N-> Sigma Kaon Pion
         ! Choose channel
         teilchenOut(1:3)%ID = (/ kaon, pion, sigmaResonance /)
         channel_index=MonteCarloChoose2Dim(SigmaKaonPion)
         teilchenOut(1)%Charge = 0+channel_index(1)-1
         teilchenOut(2)%Charge = -1+channel_index(2)-1
         teilchenOut(3)%Charge = totalCharge-sum(teilchenOut(1:2)%Charge)
         return
      end if
      cut=cut-sum(SigmaKaonPion)

      If (sum(LambdaKaonPion)>=cut) then
         ! pi N-> Lambda Kaon Pion
         ! Choose charge channel
         teilchenOut(1:3)%ID = (/ kaon, pion, Lambda /)
         channel_index=MonteCarloChoose2Dim(LambdaKaonPion)
         teilchenOut(1)%Charge = 0+channel_index(1)-1
         teilchenOut(2)%Charge = -1+channel_index(2)-1
         teilchenOut(3)%Charge = 0
         return
      end if
      cut=cut-sum(LambdaKaonPioN)


      !############################################################
      ! (4) Four body final state
      !############################################################

      if (pipipiN>=cut) then
         teilchenOut(1:4)%Id = (/ pion, pion, pion, nucleon /)
         izmin4=(/-1,-1,-1,0/)
         izmax4=(/ 1, 1, 1,1/)
         call rancharge(izmin4,izmax4,totalCharge,izout4,ranChargeFlag)
         If (.not.ranChargeFlag) write(*,*) 'Error in rancharge :',izmin4,izmax4,totalCharge,izout4,ranChargeFlag
         teilchenOut(1:4)%Charge=izout4(1:4)
         return
      end if
      cut=cut-pipipiN


      !############################################################
      ! Error message if no channel is chosen
      !############################################################
      write (*,*) 'Error in makeDecision : No decision made', &
                   cut, sum(kaonKaonBarN), teilchenOut(:)%ID, teilchenOut(:)%Charge, sigmaTot
      Stop

    end subroutine MakeDecision


    !****************************************************************************
    !****s* pionNuc/makeOutput
    ! NAME
    ! subroutine makeOutput
    ! PURPOSE
    ! Writes several cross sections to file as function of srts and plab [GeV].
    ! Filenames:
    ! * 'piN_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'piN_strangeProd.dat'        : strangeness production
    ! * 'piN_resProd.dat'            : Baryon resonance production
    ! * 'piN_nonStrange_nuk.dat'     : non-strange meson with nucleon in final state
    ! * 'piN_nonStrange_delta.dat '  : non-strange meson with delta in final state
    !****************************************************************************
    subroutine makeOutput
      logical, save :: initFlag=.true.
      real :: plab
      character(len=24), parameter :: outputFile(1:5) = (/ 'piN_sigTotElast.dat     ', &
                                                           'piN_strangeProd.dat     ', 'piN_resProd.dat         ', &
                                                           'piN_nonStrange_nuk.dat  ', 'piN_nonStrange_delta.dat' /)

      plab=SQRT(((srts**2-pionMass**2-nukMass**2)/2./nukMass)**2-pionMass**2)

      If (initFlag) then
         Open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         Open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         Open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         Open(file=outputFile(4),UNIT=104,Status='Replace',Action='Write')
         Open(file=outputFile(5),UNIT=105,Status='Replace',Action='Write')
         write (101,*) '# srts, plab, sigmaTot, sigmaElast '
         write (102,'(2A)') '# srts, plab, lambdaKaon, sigmaKaon(0:1), kaonKaonBarN(-2:1), sigmaKaonPion(0,-1:1)', &
                            'sigmaKaonPion(1,-1:1), lambdaKaonPion(0,-1:1),lambdaKaonPion(1,-1:1)'
         write (103,*) '# srts, plab, sum(sigmaRes), sigmaRes(2:40)'
         write (104,*) '# srts, plab, piN(-1:1),  phiN, omegaN, piP11_1440(-1:1) '
         write (104,*) '# , omegaPiN, phiPiN ,rhoN(-1:1) ,sigmaN,etaN, pipiN,pipipiN'
         write (105,*) '# srts, plab,  pionDelta(-1:1) , etaDelta, rhoDelta(-1:1)   '
         initFlag=.false.
      else
         Open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(4),UNIT=104,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(5),UNIT=105,Status='old',Position='Append',Action='Write')
      end If

      write (101, '(4F9.4)') srts, plab, sigmaTot, sigmaElast
      write (102,'(21F9.4)') srts, plab, lambdaKaon, sigmaKaon, kaonKaonBarN, sigmaKaonPion(0,-1:1), &
                                         sigmaKaonPion(1,-1:1), lambdaKaonPion(0,-1:1), lambdaKaonPion(1,-1:1)
      write (103,'(42F9.4)') srts, plab, sum(sigmaRes), sigmaRes(2:40)
      write (104,'(19F9.4)') srts, plab, piN(-1:1), phiN, omegaN, piP11_1440(-1:1), omegaPiN, phiPiN, &
                                         rhoN(-1:1), sigmaN, etaN, pipiN, pipipiN
      write (105, '(9F9.4)') srts, plab, pionDelta, etaDelta, rhoDelta

      Close(101)
      Close(102)
      Close(103)
      Close(104)
      Close(105)

    end subroutine makeOutput

  end subroutine pionNuc


end module pionNucleon
