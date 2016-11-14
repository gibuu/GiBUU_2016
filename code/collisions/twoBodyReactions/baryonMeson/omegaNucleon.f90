!****************************************************************************
!****m* /omegaNucleon
! NAME
! module omegaNucleon
! PURPOSE
! Includes the cross sections for omega-nucleon scattering in the resonance regime.
! Public routines:
! * omegaNuc
!****************************************************************************
module omegaNucleon
  implicit none
  Private

  ! Debug-flags
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! To decide wether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  Public :: omegaNuc

contains

  !**************************************************************************
  !****s* omegaNucleon/omegaNuc
  ! NAME
  ! subroutine omegaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,
  ! teilchenOUT,sigmaTot,sigmaElast,useHiEnergy,HiEnergySchwelle,plotFlag)
  !
  ! PURPOSE
  ! Evaluates omega Nucleon -> anything cross sections and returns also a "preevent"
  !
  ! INPUTS
  ! * real                           :: srts               --- sqrt(s) in the process
  ! * type(particle),dimension(1:2)  :: teilchenIn         --- colliding particles
  ! * type(medium)                   :: mediumATcollision  --- Medium informations at the position of the collision
  ! * real ,dimension(0:3)           :: momentumLRF        --- Total Momentum in LRF
  !
  ! High energy matching:
  ! * logical                        :: useHiEnergy ---
  !   .true. if High-Energy cross sections are given by paramBarMesHE
  ! * real                           :: HiEnergySchwelle ---
  !   threshold sqrt(s) for paramBarMesHE, i.e. at which energy the cross
  !   sections of paramBarMesHE are used
  !
  ! Debugging:
  ! * logical,optional               :: plotFlag --- Switch on plotting of the  Xsections
  !
  ! RESULT
  ! * real       :: sigmaTot         --- total Xsection
  ! * real       :: sigmaElast       --- elastic Xsection
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross
  ! sections to decide on a final state with maximal 3 final state particles.
  ! These are returned in the vector teilchenOut. The kinematics of these
  ! teilchen is only fixed in the case of a single produced resonance.
  ! Otherwise the kinematics still need to be established. The
  ! result is:
  ! * type(preEvent),dimension(1:3) :: teilchenOut --- particles
  !
  ! The cross sections are based upon a parametrization by Golubeva.
  ! See routine golub_omega in parametrizationBarMes.
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : baryon Resonances
  ! * 2-particle : pi N, omega N, K Lambda, K Sigma
  ! * 3-particle : pi pi N
  !**************************************************************************
  subroutine omegaNuc (srts, teilchenIn, mediumATcollision, momentumLRF, teilchenOUT, sigmaTot, &
                       sigmaElast, useHiEnergy, HiEnergySchwelle, plotFlag, K_Factor)
    use idTable
    use particleDefinition
    use particleProperties, only : hadron
    use mediumDefinition
    use preEventDefinition, only : preEvent
    use twoBodyTools, only : velocity_correction, convertToAntiParticles, pcm,searchInInput
    use RMF, only : getRMF_flag
    use constants, only: mN, mPi, mK
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Input
    real, intent(in)                            :: srts
    type(particle), dimension(1:2), intent(in)  :: teilchenIn
    type(medium), intent(in)                    :: mediumATcollision
    real, intent(in) ,dimension(0:3)            :: momentumLRF
    logical, intent(in)                         :: useHiEnergy
    real, intent(in)                            :: HiEnergySchwelle
    logical, intent(in), optional               :: plotFlag
    real, intent(in)                            :: K_Factor
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Output
    type(preEvent), dimension(1:3), intent(out) :: teilchenOut
    real, intent(out)                           :: sigmaTot, sigmaElast
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Cross sections
    real,dimension(-1:1) :: piN, piN_R            ! -> pi N, index denotes pion charge
    real :: omegaN, pipiN, lambdaKaon
    real, dimension(0:1)  :: sigmaKaon            ! index = charge of final state kaon
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Local variables
    real :: fluxCorrector    ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
    real :: s
    type(particle) :: omega_particle, nucleon_particle
    logical :: antiParticleInput, failFlag
    real, dimension(Delta:nbar) :: sigmaRes, massRes    ! partial cross sections for omega N -> R, and resonance masses

    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    teilchenOut(:)%ID=0                    ! ID of produced particles
    teilchenOut(:)%charge=0                ! Charge of produced particles
    teilchenOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
    teilchenOut(:)%mass=0                  ! Mass of produced particles

    ! (1) Check  Input
    call searchInInput(teilchenIn,omegaMeson,nucleon,omega_particle,nucleon_particle,failFlag)
    If (failFlag) then
       Write(*,*) 'Wrong input in OmegaNuc', teilchenIn%ID
    end if

    If (omega_particle%antiParticle) then
       ! This case is not considered yet
       write(*,*) 'omega is antiparticle in "omegaNuc"!!!',teilchenIN%ID,teilchenIN%antiparticle
       stop
    end if

    If (nucleon_particle%antiParticle) then
       ! Invert all particles in antiparticles
       nucleon_particle%Charge        =  -nucleon_particle%Charge
       nucleon_particle%antiparticle  = .false.
       omega_particle%Charge          =  -omega_particle%Charge
       antiParticleInput              = .true.
    else
       antiParticleInput=.false.
    end if


    ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
    if (.not.getRMF_flag()) then
      fluxCorrector=velocity_correction(teilchenIn)
    else
      fluxCorrector=1.
    end if

    s=srts**2

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
    If (Sum(teilchenOut(:)%Charge).ne.nucleon_particle%charge+omega_particle%charge) then
       write(*,*) 'No charge conservation in pionNuc!!! Critical error', omega_particle%Charge, nucleon_particle%Charge, &
                                                                         teilchenOut(:)%Charge, teilchenOut(:)%ID
       stop
    end if

    ! (6) Invert particles in antiParticles if input included antiparticles
    If (antiParticleInput) then
       IF(debugFlagAnti) write(*,*) teilchenOut
       call convertToAntiParticles(teilchenOut)
       IF(debugFlagAnti) write(*,*) teilchenOut
    end if

  contains

    subroutine evaluateXsections
      use resonanceCrossSections, only: barMes_R_barMes, barMes2resonance
      use mediumDefinition, only : vacuum
      use parametrizationsBarMes, only : golub_omega, omegaN_lykasov, huang, huanglam
      use parBarMes_HighEnergy, only : paramBarMesHE
      use output, only : writeparticle

      real, dimension(1:3) ::  position
      logical :: perturbative
      real, dimension(1:2) :: sigmaGolub
      real, dimension(1:4) :: sigmaHuang
      real, dimension(-1:1) :: piN_Golub
      real :: p_piN, p_omegaN, ratio
      real :: pFinal,pInitial,detailedBalanceFactor,elast_R
      integer :: pionCharge, nucCharge
      real, dimension(0:3) :: momentum_vacuum
      real :: sigmaTotal_HE,sigmaElast_HE ! High energy matchin
      logical :: lDummy

      position=0.5*(teilchenIN(1)%position+teilchenIN(2)%position)
      if (teilchenIN(1)%perturbative.or.teilchenIN(2)%perturbative) then
         perturbative=.true.
      else
         perturbative=.false.
      end if

      momentum_vacuum(1:3)=teilchenIn(1)%momentum(1:3)+teilchenIn(2)%momentum(1:3)
      momentum_vacuum(0)=FreeEnergy(teilchenIn(1))+FreeEnergy(teilchenIn(2))

      !#################################################################
      ! Evaluate partial cross sections
      !######################################################################

      !*******************************************************************************************
      ! omega N -> pi N
      !*****************************************************************************************
      ! piN = [cross section by Golubeva (pi^- p-> omega N) by detailed balance] - [vacuum cross section via resonances]

      If (srts > mPi + mN) then
        sigmaGolub = golub_omega (srts)            ! Results by Golubeva

        pFinal = pCM(srts, mPi, mN)
        pInitial = pCM(srts, omega_Particle%mass, nucleon_Particle%mass)

        If (pinitial<1E-12) then
          write(*,*) 'WARNING: pInitial is zero in omegaNuc', pinitial
          write(*,*) 'omega meson:'
          call writeparticle(6,0,0,omega_Particle)
          write(*,*) 'nucleon:'
          call writeparticle(6,0,0,nucleon_Particle)
          detailedBalanceFactor= 0.
        else
          detailedBalanceFactor= 1./3.*(pFinal/pInitial)**2
          ! given by detailed balance: factor 1/3 due to (2j+1)-Terms in cross section and different
          ! spins in initial and final state
        end if

        sigmaGolub(1) = 3./2. * sigmaGolub(1)  ! 3./2. since Golub returns cross section pi^- p -> omega N, with 3./.2 we divide by the isospin-clebsch

        If (nucleon_particle%charge==0) then
          piN_golub(-1) = 2./3. * sigmaGolub(1) * detailedBalanceFactor
          piN_golub (0) = 1./3. * sigmaGolub(1) * detailedBalanceFactor
          piN_golub (1) = 0.
        else If (nucleon_particle%charge==1) then
          piN_golub(-1) = 0.
          piN_golub (0) = 1./3. * sigmaGolub(1) * detailedBalanceFactor
          piN_golub (1) = 2./3. * sigmaGolub(1) * detailedBalanceFactor
        else
          write(*,*) 'Error in omegaNuc', nucleon_particle
        end if

        Do pionCharge=-1,1
          nucCharge=omega_particle%charge+nucleon_particle%charge-pionCharge
          If ((nucCharge==0).or.(nucCharge==1)) then
              piN_R(pionCharge) = barMes_R_barMes(omegaMeson, nucleon, pion, nucleon, &
                                      omega_particle%Charge, nucleon_particle%Charge, pionCharge, nucCharge, &
                                      .false., .false., Vacuum, momentum_vacuum, omega_particle%Mass, nucleon_particle%Mass, &
                                      position, perturbative, srts)
              piN(pionCharge) = max(0., piN_golub(pionCharge) - piN_R(pionCharge))
          else
              piN(pionCharge) = 0.
          end if
        end do
      else
        piN = 0.
      end if

      !*******************************************************************************************
      ! omega N -> omega N
      !*****************************************************************************************

      elast_R = barMes_R_barMes(omegaMeson,nucleon,omegaMeson,nucleon,omega_particle%Charge,nucleon_particle%Charge, &
                                omega_particle%Charge,nucleon_particle%Charge,.false.,.false.,MediumAtCollision,momentumLRF, &
                                omega_particle%Mass,nucleon_particle%Mass,position,perturbative,srts)

      ! subtract resonance contribution from (direct) omega N channel
      if (srts>mN+hadron(omegaMeson)%minmass) then
        omegaN=Max(0., omegaN_lykasov(srts,omega_Particle%mass,1)  - elast_R)
      else
        omegaN=0.
      end if

      !*******************************************************************************************
      ! omega N -> R
      !*****************************************************************************************

      ! Full resonance contribution in the medium
      sigmaRes = barMes2resonance (omegaMeson,nucleon,omega_particle%charge,nucleon_particle%charge,.true.,mediumAtCollision, &
                                   momentumLRF,massRes,omega_particle%Mass,nucleon_particle%Mass,position,perturbative,srts)

      !*******************************************************************************************
      ! -> pi pi N
      !*******************************************************************************************

      ! Evaluate Omega N -> N pi pi by taking the total inelastic cross section by Lykasov (times a K-factor)
      ! and subtracting then all included inelastic cross sections.
      ! Inelastic resonance contribution = full contribution - elastic contribution

      If (srts > 2*mPi + mN) then
         pipiN=max(0.,omegaN_lykasov(srts,omega_Particle%mass,2)*K_Factor-(Sum(sigmaRes)-elast_R)-Sum(piN))
         ! Matching to High energy region
         if(useHiEnergy) then
            call paramBarMesHE (HiEnergySchwelle,omegaMeson,nucleon,omega_particle%charge,nucleon_particle%charge, &
                                mediumAtCollision,sigmaTotal_HE,sigmaElast_HE)
            sigmaTotal_HE = (sigmaTotal_HE-sigmaElast_HE)*K_factor + sigmaElast_HE
            pipiN=max(pipiN,sigmaTotal_HE-Sum(sigmaRes)-Sum(piN)-omegaN)
         end if
      else
         piPiN=0.
      end if

      !###################################################################################################
      ! evaluate elastic Xsection
      !###################################################################################################

      sigmaElast = omegaN + elast_R

      ! Correction factor to the pion-nucleon cross sections by detailed balance
      ! (see J. Cugnon et al, PRC 40, 1822 (1989))
      ! c.m. momenta of pion-nucleon and omega-nucleon
      p_piN    = pCM(s,mPi,nucleon_particle%Mass, lDummy)
      p_omegaN = pCM(s,omega_particle%Mass,nucleon_particle%Mass, lDummy)
      if(p_omegaN.gt.1.e-06) then
         ratio=p_piN/p_omegaN
      else
         ratio=1.
      end if

      !**************************************************************
      ! -> Lambda Kaon
      !**************************************************************
      lambdaKaon=0.
      If(srts.gt.(hadron(Lambda)%mass+mK)) then
         ! huanglam gives : pi^{-} p -> Lambda Kaon^{0}
         lambdaKaon = 0.5 * huangLam(srts) * ratio  ! assume that sigma(omega p -> Lambda K^+) = sigma(pi^0 p -> Lambda K^+) * p_piN/p_omegaN
                                             ! No resonance contribution
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
      If(srts.gt.(hadron(SigmaResonance)%mass+mK)) then
         sigmaHuang = huang(srts) * ratio       ! correction due to detailed balance
         If(nucleon_particle%Charge.eq.1) then
            sigmaKaon(1)=sigmaHuang(2)     ! assume that sigma(omega p -> K^+ Sigma^0) = sigma(pi^0 p -> K^+ Sigma^0)
            sigmaKaon(0)=2.*sigmaKaon(1)   ! by isospin consideration
         else                              ! neutron Xsections by charge conjugation
            sigmaKaon(0)=sigmaHuang(2)
            sigmaKaon(1)=2.*sigmaKaon(0)
         end if
      end if

      !###################################################################################################
      ! Do the flux correction for each channel
      !###################################################################################################

      If(fluxCorrector_flag) then
         ! We do this for each channel since they might show up seperately in the output if makeoutput is called
         sigmaElast=sigmaElast*fluxcorrector
         omegaN=omegaN*fluxcorrector
         piN=piN*fluxcorrector
         pipiN=pipiN*fluxcorrector
         sigmaRes=sigmaRes *fluxcorrector
         lambdaKaon=lambdaKaon*fluxcorrector
         sigmaKaon=sigmaKaon*fluxcorrector
      end if

      !###################################################################################################
      ! Sum up everything for the total cross section
      !###################################################################################################
      ! Be careful since sigma elast is already included in the partial cross sections (omegaN and sigmaRes),
      ! therefore it is not included in the total cross section.

      sigmaTot = omegaN + sum(piN) + pipiN + sum (sigmaRes) + lambdaKaon + sum(sigmaKaon)

    end subroutine evaluateXsections


    subroutine makeDecision
      use random, only: rn, ranCharge

      real :: summe, cut, cut2
      integer :: resID, totalCharge, pionCharge, charge
      integer,dimension (1:3)  :: izmin,izmax,izout     ! needed for ranCharge
      logical :: ranChargeFlag

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=omega_particle%Charge+nucleon_particle%Charge
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

      ! omega N production
      If(omegaN.ge.cut) then
         teilchenOut(1)%Id=omegaMeson
         teilchenOut(2)%Id=nucleon

         teilchenOut(1)%Charge=0
         teilchenOut(2)%Charge=totalCharge
         return
      end if
      cut=cut-omegaN

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

      ! Kaon Lambda production
      If (lambdaKaon .ge.cut) then
         teilchenOut(1)%Id=kaon
         teilchenOut(2)%Id=Lambda
         teilchenOut(1)%Charge=totalCharge
         teilchenOut(2)%Charge=0
         return
      end if
      cut=cut-lambdaKaon

      ! Kaon Sigma production
      If (sum(sigmaKaon) .ge.cut) then
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
      write(*,*) 'Error in makedecision of omegaNuc', srts, sigmaTot, omegaN, sum(piN), pipiN, sum(sigmaRes), cut
      stop

    end subroutine makeDecision


    !******************************************************************************************
    !****s* omegaNucleon/makeOutput
    ! NAME
    ! subroutine makeOutput
    !
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV]
    ! .
    ! Filenames:
    ! * 'omegaN_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'omegaN_resProd.dat'            : Baryon production (Resonances with ID's 1:40)
    ! * 'omegaN_nonStrange_nuk.dat'     : non-strange meson with nucleon in final state
    ! * 'omegaN_strangeProd.dat'        : Kaon and hyperon in final state
    !******************************************************************************************

    subroutine makeOutPut

      logical, save :: initFlag=.true.
      real :: plab
      character(30), parameter :: outputFile(1:4) = (/ 'omegaN_sigTotElast.dat   ', 'omegaN_resProd.dat       ', &
                                                       'omegaN_nonStrange_nuk.dat', 'omegaN_strangeProd.dat   ' /)

      plab=SQRT(((s-omega_particle%mass**2-nucleon_particle%mass**2)/2./nucleon_particle%mass)**2-omega_particle%mass**2)

      If (initFlag) then
         Open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         Open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         Open(file=outputFile(3),UNIT=103,Status='Replace',Action='Write')
         Open(file=outputFile(4),UNIT=104,Status='Replace',Action='Write')
         write (101,*) '# srts, plab, sigmaTot, sigmaElast '
         write (102,*) '# srts, plab, sum(sigmaRes), sigmaRes(2:40) '
         write (103,*) '# srts, plab, piN(-1:1), piN_R(-1:1), omegaN, pipiN'
         write (104,*) '# srts, plab, lambdaKaon, sigmaKaon(0:1)'
         initFlag=.false.
      else
         Open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(3),UNIT=103,Status='old',Position='Append',Action='Write')
         Open(file=outputFile(4),UNIT=104,Status='old',Position='Append',Action='Write')
      end If
      write (101, '(4F10.4)') srts, plab, sigmaTot, sigmaElast
      write (102,'(42F10.4)') srts, plab, sum(sigmaRes), sigmaRes(2:40)
      write (103,'(10F10.4)') srts, plab, piN, piN_R, omegaN, pipiN
      write (104, '(5F10.4)') srts, plab, lambdaKaon, sigmaKaon
      Close(101)
      Close(102)
      Close(103)
      Close(104)

    end subroutine makeOutPut
  end subroutine omegaNuc


end module omegaNucleon
