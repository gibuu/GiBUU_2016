!*******************************************************************************
!****p* /GiBUU
! NAME
! program GiBUU
!
! PURPOSE
! This is the main file of the Giessen Boltzmann-Uehling-Uhlenbeck (GiBUU) code.
! It steers the whole simulation, including initialization, time evolution
! (propagation & collisions) and analysis.
!
! COPYRIGHT
! (C) 2005-2014 The GiBUU Team (see full list of authors below)
!
! The GiBUU code is licensed under the GNU General Public License (GPL v2).
! See accompanying LICENSE file or http://www.gnu.org/licenses/gpl-2.0.html.
!
! AUTHOR
! Academic Supervisor:
! * Ulrich Mosel
!
! Alphabetical list of authors:
! * Oliver Buss
! * Thomas Falter
! * Theodoros Gaitanos
! * Kai Gallmeister
! * David Kalok
! * Murat Kaskulov
! * Alexei Larionov
! * Olga Lalakulich
! * Ivan Lappo-Danilewski
! * Tina Leitner
! * Ulrich Mosel
! * Birger Steinmueller
! * Janus Weil
!
! Email:
! * gibuu@projects.hepforge.org
!
! Postal address:
! * Institut fuer Theoretische Physik I
! * Justus-Liebig-Universitaet Giessen
! * Heinrich-Buff-Ring 16
! * D-35392 Giessen
!
! Telephone:
! * +49 (0)641 99 33344   (U. Mosel)
!
! SEE ALSO
! The full physics included in the code is described in this review paper:
! * O. Buss et al., Phys. Rept. 512 (2012) 1,
!   http://inspirehep.net/record/912923
! Additional information and documentation can be found on the website:
! * https://gibuu.hepforge.org
!
! BUGS
! Please report bugs and suggestions for improvements to
! gibuu@projects.hepforge.org.
!*******************************************************************************
program GiBUU

  use output, only: header, chapter, subchapter, DoPR, timeMeasurement, printTime, writeParticleVector
  use version, only: printVersion
  use inputGeneral, only: eventType, delta_T, num_Energies, num_runs_SameEnergy, current_run_number, readinputGeneral
  use statistics, only : splashInfo
  use particleDefinition, only: particle, setToDefault, setNumbersToDefault
  use nucleusDefinition, only: tnucleus
  use collisionNumbering, only : writeCountedEvents
  use random, only: setRandom
  use particleProperties, only: initParticleProperties
  use checks, only: ChecksSetDefaulSwitches

  implicit none


  !*************************************************************************
  !****ig* GiBUU/realParticles
  ! SOURCE
  !
  type(particle),Allocatable :: realParticles(:,:)
  !
  ! PURPOSE
  ! The particle vector: real particles
  !
  ! NOTES
  ! * First index : ensemble
  ! * Second index : position within ensemble
  !*************************************************************************

  !*************************************************************************
  !****ig* GiBUU/pertParticles
  ! SOURCE
  !
  type(particle),Allocatable :: pertParticles(:,:)
  !
  ! PURPOSE
  ! The particle vector: perturbative particles
  !
  ! NOTES
  ! * First index : ensemble
  ! * Second index : position within ensemble
  !*************************************************************************

  !*************************************************************************
  !****ig* GiBUU/targetNuc
  ! SOURCE
  !
  type(tNucleus),pointer :: targetNuc => NULL()
  !
  ! PURPOSE
  ! Target nucleus
  !*************************************************************************

  !*************************************************************************
  !****ig* GiBUU/projectileNuc
  ! SOURCE
  !
  type(tNucleus),pointer :: projectileNuc => NULL()   ! Projectile nucleus
  !
  ! PURPOSE
  ! Projectile nucleus
  !*************************************************************************


  integer :: j,k
  logical :: raiseEnergy
  real, save :: delta_T_max   !field to store initial delta_T from inputGeneral


  ! Formats for output:
  character(200), parameter :: format1 = '(79("="),/,5("=")," ",A," ",i5," / ",i5,/,79("="),/)'
  character(200), parameter :: format2 = '(79("#"),/,5("#")," ",A," ",i8," / ",i8,/,79("#"),/)'
  character(200), parameter :: format3 = '(5("=")," ",A," ",i5,"/",i5,"  ",5("#")," ",A," ",i9,"/",i9)'

  call PrintTime("START")
  write(*,header) 'BUU simulation: start'

  !============= Output the version of the code
  call PrintVersion

  !============= Initialize the database entries, i.e. the fields "baryon" and "meson"
  write(*,chapter) 'Init database:...'
  call readInputGeneral
  call initParticleProperties
  write(*,chapter) 'Init database: finished'
  call ChecksSetDefaulSwitches(EventType)

  raiseEnergy=.false.
  delta_T_max=delta_T

  energyLoop: Do j=1,num_Energies         ! loop over different energies
     write(*,format1) 'Energy loop :',j,num_Energies

     subsequentRunsLoop: Do k=1,num_runs_SameEnergy   ! loop over subsequent runs
        write(*,format2) 'Run loop :',k,num_runs_SameEnergy

        !============= Output the present status:
        open(123,file='main.run',status='unknown')
        write(123,format3) 'Energy loop :',j,num_Energies, 'Run loop :',k,num_runs_SameEnergy
        close(123)

        !============= Take care of random numbers:
        call setRandom()
        current_run_number = current_run_number + 1

        !============= Initialize configuration
        write(*,subchapter) 'Init starting configuration'
        call initConfig(raiseEnergy)

        !============= Do some analysis after the init and before the run.
        call analysis(.false.,.true.)

        !============= Transport the hadronic matter
        write(*,subchapter) 'Do RUN'
        call run

        !============= Final analysis
        write(*,subchapter) 'Analysing output of the run'
        call analysis(k.eq.num_runs_sameEnergy,.false.)
        raiseEnergy=.false.

        !============= Print out statistical informations
        write(*,*) '(FINAL) number of counted Events:'
        call writeCountedEvents(0)

        !============= Sets info module back to normal
        call splashInfo

     End do subsequentRunsLoop
     raiseEnergy=.true.
  End do energyLoop

  call finalCleanup()

  write(*,header) 'BUU simulation: finished'
  call PrintTime("STOP")


contains


  !*************************************************************************
  !****s* GiBUU/initConfig
  ! NAME
  ! subroutine initConfig(energyRaiseFlag)
  !
  ! PURPOSE
  ! * Initializes starting configuration for particles in each run.
  ! * Here we have some "select case(eventType)":
  !   for each eventType different subroutines take care of the proper
  !   initialization.
  ! * The vectors realParticles and pertParticles are allocated.
  !
  ! INPUTS
  ! * logical :: energyRaiseFlag -- if (true) then we raise the bombarding
  !   energy (used to perform an energy scan)
  !
  !*************************************************************************
  subroutine initConfig(energyRaiseFlag)
    use inputGeneral,               only: length_perturbative, length_real, numEnsembles, printParticleVectors
    use initBox,                    only: initializeBox, BoostToEps
    use initPion,                   only: InitPionInduced
    use initPionBox,                only: initializePionBox
    use initHiPion,                 only: InitHiPionInduced
    use initLowPhoton,              only: initialize_lowPhoton, lowPhotonInit_getRealRun
    use lowElectron,                only: init_lowElectron
    use initHiLepton,               only: InitHiLeptonInduced
    use initHeavyIon,               only: initHeavyIonCollision
    use initElementary,             only: initElementaryCollision
    use initHadron,                 only: initHadronInduced
    use initNeutrino,               only: init_Neutrino, neutrinoInit_getRealRun
    use initTransportGivenParticle, only: init_transportGivenParticle
    use initExternal,               only: initializeExternal, ExternalIsPerturbative
    use nucleus,                    only: getTarget, getProjectile
    use densityModule,              only: updateDensity, updateRMF, storeFields
    use RMF, only : getRMF_flag
    use yukawa, only: updateYukawa
    use coulomb, only: updateCoulomb
    use energyCalc, only: updateEnergies
    use eventtypes
    use checks, only: ChecksCallEnergy
    use propagation, only: propagate_cascade, updateVelocity
    use insertion, only : GarbageCollection
    use PILCollected, only: PILCollected_ZERO
    use initNucleus_in_PS,only : initNucPhaseSpace, shiftBack
    use groundStateAnalysis,only : countMass
    use baryonPotentialModule, only: HandPotentialToDensityStatic
    use deuterium_PL, only : deuteriumPL_assign
    use DeltaTest, only: initDeltaTest

    logical, intent(in) :: energyRaiseFlag
    integer :: lengthReal = 0  ! maximum number of real particles per ensemble
    integer :: lengthPert = 0  ! maximum number of perturbative particles per ensemble
    logical, save :: first = .true.

    if (first) then

       if (eventType.ne.elementary) then
          targetNuc => getTarget()                 !set up target resting at 0
          lengthReal = targetNuc%mass ! Real particles fixed by nucleons in target

          if (targetNuc%ReAdjustForConstBinding) then
             write(*,*) 'we now have to readjust the density'
             call HandPotentialToDensityStatic(targetNuc)
          end if

       else if (.not. associated(targetNuc)) then
          allocate(targetNuc) ! dummy target, only %velocity is used
       end if

      if (eventType==HeavyIon) projectileNuc => getProjectile()    ! set up projectile resting at 0

      !...Defining maximal lenghts of particle vectors: (Default)

      select case (eventType)
      Case(elementary) ! 0=elementary collision
         lengthPert =   1
         lengthReal = 100
      Case(HeavyIon,ExternalSource)   ! 1=Heavy Ion collision, 22=Transport for external source
         lengthPert =   10
         lengthReal = 3000
      Case(LoPion)       ! 2=Pion Nucleus collision
         lengthPert =  300
      Case(RealPhoton)      ! 3=photon nucleus collision
         if (lowPhotonInit_getRealRun()) then
            lengthReal = targetNuc%mass*11
            lengthPert = 1
         else
            lengthPert = max(100,15*targetNuc%mass)
         end if
      Case(LoLepton)   ! 4=virtual photon collision
         lengthPert = 10000
      Case(Neutrino)   ! 5=neutrino nucleus collision
         if (neutrinoInit_getRealRun()) then
            lengthReal = targetNuc%mass*11
            lengthPert = 1
         else
            lengthPert =  max(100,10*targetNuc%mass)
         end if
      Case(HiPion)     ! 12= High energy pion induced collisions
         lengthPert = 10000
      Case(HiLepton)   ! 14= high energy lepton nucleus collision
         lengthPert = 5000
      Case(InABox)
         lengthPert =  1
         lengthReal = 5000
      Case(InABox_pion,inABox_delta)
         lengthPert =  150
         lengthReal = 5000
      Case(PionBox)
         lengthPert =  1
         lengthReal = 5000
      Case(groundState)
         lengthPert =   1
         lengthReal = 500
      Case(transportGivenParticle)
         lengthPert = 100
      Case(hadron)     !300= hadron-nucleus collision
         lengthPert =  1
         lengthReal = targetNuc%mass + 100
      end select

      !...Defining maximal lenghts of particle vectors: (From Input)

      if (length_perturbative >= 0) lengthPert = length_perturbative
      if (length_real >= 0)         lengthReal = length_real

      !...Allocate the vectors

      Allocate(realparticles(1:numEnsembles,1:lengthReal))
      Allocate(pertparticles(1:numEnsembles,1:lengthPert))

      if (DoPr(1)) then
         write(*,'(79("#"))')
         write(*,'(79("#"))')
         write(*,'("##              ","#Ensemble","*(","len(Pert)","+","len(Real)",")",'// &
              &  ' "   * ",i3," Bytes                ","##")') sizeof(realParticles(1,1))
         write(*,'("## MEMORY USAGE=",i9,"*(",i9,"+",i9,") = ",i12," <-> ",f6.1," MB ##")') &
              numEnsembles,lengthPert,lengthReal, numEnsembles*(lengthReal+lengthPert), &
              real(numEnsembles*(lengthReal+lengthPert)*sizeof(realParticles(1,1)))/(1024**2)
         write(*,'(79("#"))')
         write(*,'(79("#"))')
      end if

      first=.false.

    else

      call setToDefault(realParticles)
      call setToDefault(pertParticles)

    end if

    call setNumbersToDefault

    call GarbageCollection(pertParticles,.true.)
    call GarbageCollection(realParticles)

    call PILCollected_ZERO ! reset all particle info lists


    !...Do the Init

    Select Case(eventType)

    Case(elementary) ! 0=elementary collision

       call initElementaryCollision(realparticles,energyRaiseFlag)


    Case(HeavyIon) ! 1=Heavy Ion collision

       call initHeavyIonCollision(targetNuc,projectileNuc) !fix kinematics of the nuclei

       ! Set up test-particles which represent the nuclei:
       call initNucPhaseSpace(realparticles,projectileNuc)
       call initNucPhaseSpace(realparticles,targetNuc)
       !shift back coordinates along z-axis (see notes in initNucPhaseSpace):
       call shiftBack(realparticles)

       if( .not.getRMF_flag() ) then

          call updateDensity(realParticles)
          call updateCoulomb
          call updateYukawa(.true.)

       else

          call updateRMF(realParticles) ! This routine also updates single-particle
                                        ! energies (p(0)^*'s), velocities.
                                        ! updateDensity is now also called from updateRMF.
          if (DoPr(2)) write(*,*) 'Do delta_T propagation to get baryon 4-current'
          call propagate_cascade(realParticles,delta_T)
          call storeFields
          call updateRMF(realParticles)
          call updateCoulomb

       end if


    Case(LoPion) ! 2=Pion Nucleus collision

       ! Set up test-particles which represent the nucleus:
       call SetUpTarget(.false.)
       ! Set up test-particles which represent the pions:
       call InitPionInduced(pertParticles,energyRaiseFlag,targetNuc)

    Case(RealPhoton) ! 3= photon nucleus collision

       ! Set up test-particles which represent the nucleus:
       call SetUpTarget()
       call deuteriumPL_assign(realParticles)

       call initialize_lowPhoton(realParticles, pertParticles,energyRaiseFlag,targetNuc)

    Case(LoLepton) ! 4= virtual photon nucleus collision

       ! Set up test-particles which represent the nucleus:
       call SetUpTarget()

       call init_lowElectron(realParticles, pertParticles,energyRaiseFlag,targetNuc)

    Case(Neutrino) ! 5= neutrino nucleus collision

       ! Set up test-particles which represent the nucleus:
       call SetUpTarget()

       call init_Neutrino(realParticles,pertParticles,energyRaiseFlag,num_runs_sameEnergy,targetNuc)

    Case(HiPion) ! 12= High energy pion induced collisions

       ! Set up test-particles which represent the nucleus:
       call SetUpTarget(.false.)

       ! Set up test-particles which represent the pions:
       ! (replacing all perturbative pions by the reaction outputs)
       call InitHiPionInduced(pertParticles,realParticles,targetNuc)

    Case(HiLepton) ! 14= high energy lepton nucleus collision

       ! Set up test-particles which represent the nucleus:
       call SetUpTarget(.false.)

       call ChecksCallEnergy(-99.9,realParticles)

       ! Set up test-particles which represent the result
       ! of the high energy leptoproduction
       call InitHiLeptonInduced(realParticles,pertParticles,targetNuc)

    Case(InABox,InABox_pion,InABox_delta)  ! 31,32,33 = BOX

       ! Set up test-particles which represent the nucleus:
       call initializeBox(realParticles)  ! set up the box with nucleons

       call updateDensity(realParticles)
       call updateCoulomb
       call updateYukawa(.true.)
       call updateVelocity(realParticles)

       ! Set up test-particles which shall be propagated in the box

!       if (DoPr(1)) write(*,*)  'Vor InitBoxParticles'
       select case(eventType)
       case (inABox)
          call BoostToEps(realParticles)
       case (inABox_pion)
          call InitPionInduced(pertParticles,energyRaiseFlag)
       case (inABox_delta)
          call initDeltaTest(pertParticles)
       case DEFAULT
          write(*,*) 'error inAbox in initConfig'
          stop
       end select
!       if (DoPr(1)) write(*,*)  'Nach InitBoxParticles'

    Case(PionBox) ! 41 = Pions in a Box

       call initializePionBox(realParticles)  ! set up the box

    Case(ExternalSource) ! 22 = Transport for external hadronic source

       call initializeExternal(realParticles,pertParticles)

       if (ExternalIsPerturbative()) then
          ! Set up test-particles which represent the nucleus:
          call SetUpTarget()
       else

          call updateDensity(realParticles)
          call updateCoulomb
          if( .not.getRMF_flag() ) then
             call updateYukawa(.true.)
          else
             call updateRMF(realParticles) ! This routine also updates single-particle
             ! energies (p(0)^*'s), velocities.
             ! updateDensity is now also called from updateRMF.
             call storeFields
          end if
       endif

    Case(groundState)
       ! Set up test-particles which represent the nucleus:
       call initNucPhaseSpace(realparticles,targetNuc)

       call countMass(realparticles)

       if( .not.getRMF_flag() ) then

          call updateDensity(realParticles)
          call updateCoulomb
          call updateYukawa(.true.)

       else

          call updateRMF(realParticles) ! This routine also updates single-particle
                                        ! energies (p(0)^*'s), velocities.
                                        ! updateDensity is now also called from updateRMF.
          call storeFields
          call updateCoulomb

       end if


    Case(hadron)   !300= hadron-nucleus collision (with nonperturbative hadron)
       ! Set up test-particles which represent the nucleus:
       call initNucPhaseSpace(realparticles,targetNuc)

       if(getRMF_flag()) then
          call updateRMF(realParticles) ! This routine also updates single-particle
                                        ! energies (p(0)^*'s), velocities.
                                        ! updateDensity is now also called from updateRMF.
          call updateCoulomb
       end if
       ! Initialise a hadron:
       call initHadronInduced(realParticles,pertParticles)

       if( .not.getRMF_flag() ) then

          call updateDensity(realParticles)
          call updateCoulomb
          call updateYukawa(.true.)

       else

          call updateRMF(realParticles) ! This routine also updates single-particle
                                        ! energies (p(0)^*'s), velocities.
                                        ! updateDensity is now also called from updateRMF.
          call storeFields
          call updateCoulomb

       end if


    Case(transportGivenParticle)
       ! Set up test-particles which represent the nucleus:
       call SetUpTarget(.false.)

       call init_transportGivenParticle(realParticles,pertParticles)

    Case default
       Write(*,*) 'No valid eventtype:',eventtype,'STOP'
       Stop
    End select

    ! clean up after init routines
    call GarbageCollection(pertParticles,.true.)
    call GarbageCollection(realParticles)

    call deuteriumPL_assign (realParticles)

    if (.not. getRMF_flag()) then
       !... update momentum(0) after updating the mean fields
       call updateEnergies(realParticles)
       call updateEnergies(pertParticles)
       ! Update velocities
       call updateVelocity(pertParticles)
       call updateVelocity(realParticles)
    end if

    If (PrintParticleVectors) then
       !output initializations:
       write(*,'(A)') "The initialized particle vectors can be found in the files *_Init_*.dat"
       call writeParticleVector('RealParticles_Init',realParticles)
       call writeParticleVector('PertParticles_Init',pertParticles)
    end if

    ! For safety a final garbage collection:
    call GarbageCollection(pertParticles,.true.)
    call GarbageCollection(realParticles)

    !=== Elementary checks ===

    call ChecksCallEnergy(0.0,realParticles)

  end subroutine initConfig

  !***********************************************************************
  !****s* GiBUU/SetUpTarget
  ! NAME
  ! subroutine SetUpTarget(DoUpdateEnergies)
  !
  ! PURPOSE
  ! Initialize the test particles representing the target nucleus
  !
  ! INPUTS
  ! * logical, OPTIONAL :: DoUpdateEnergies -- Flag to indicate, whether
  !   'updateEnergies' should be called
  ! OUTPUT
  ! The particle vector realparticles is modified, also global arrays
  ! connected with this.
  !***********************************************************************
  subroutine SetUpTarget(DoUpdateEnergies)
    use initNucleus_in_PS,only : initNucPhaseSpace
    use densityModule, only: updateDensity
    use coulomb, only: updateCoulomb
    use yukawa, only: updateYukawa
    use propagation, only: updateVelocity
    use energyCalc, only: updateEnergies

    logical, intent(in), OPTIONAL :: DoUpdateEnergies

    logical :: DoUE

    DoUE = .true.
    if (present(DoUpdateEnergies)) DoUE = DoUpdateEnergies

    call initNucPhaseSpace(realparticles,targetNuc)
    call updateDensity(realParticles)
    call updateCoulomb
    call updateYukawa(.true.)
    call updateVelocity(realParticles)
    if (DoUE) call updateEnergies(realParticles)

  end subroutine SetUpTarget


  !*****************************************************************************
  !****s* GiBUU/run
  ! NAME
  ! subroutine run
  !
  ! PURPOSE
  ! Propagate a given realParticle and pertParticle vector
  ! according to the BUU equations.
  !*****************************************************************************
  subroutine run
    use povray, only: povray_output
    use RMF, only : getRMF_flag
    use yukawa, only: updateYukawa
    use densitymodule, only: updateDensity, gridSpacing, gridSize
    use coulomb, only: updateCoulomb
    use propagation, only: Propagate
    use propagation_RMF, only: Propagate_RMF
    use checks, only: ChecksCallAll, ChecksCallEnergy, evaluateTimeStep, CheckGridSize
                     !,ChecksCallOccupied
    use collisionTerm, only: collideMain, forceDecays
    use insertion, only: GarbageCollection
    use thermoDynamics, only: updateTemperature
    use energyCalc, only: updateEnergies
    use hadronFormation, only: formation
    use collisionNumbering, only: writeCountedEvents, nullCountedEvents, nulln_participants
    use pionGamma, only: count_Pions
    use DeltaTest, only: countDeltas
    use eventtypes
    use HeavyIonAnalysis,               only: DoHeavyIonAnalysisTime, HeavyIon_evol
    use sourceAnalysis,                 only: DoSourceAnalysis,getSMM_Flag,resetSMMvalues,stopGiBUU
    use HiLeptonAnalysis,               only: HiLeptonAnalysisPerTime
    use HiPionAnalysis,                 only: HiPionAnalysisPerTime
    use hadronAnalysis,                 only: DoHadronAnalysisTime
    use PionBoxAnalysis,                only: DoPionBoxAnalysisTime
    use transportGivenParticleAnalysis, only: transportGivenParticle_analyze
    use Dilepton_Analysis,              only: Dilep_Decays
    use radiativeDeltaDecay, only: DoRadiativeDeltaDecay
    use inputGeneral, only: numTimeSteps, printParticleVectorTime, variableTimeStep, time_max, povray_switch, freezeRealParticles, &
                            checkGridSize_Flag, DoFragmentNucleons
    use deuterium_PL, only: deuteriumPL_assign
    use FragmentNucleons, only: DoAddFragmentNucleons
    use BoxAnalysis, only: DoBoxAnalysisTime
    use FreezeoutAnalysis, only: DoFreezeoutAnalysisPerTime
    use EventOutputAnalysis, only: DoEventOutput

    integer :: timeStep
    real :: time,delta_T_new

    time=0.

    call nullCountedEvents(0)
    call nulln_participants

    !loop over time steps
    If (numTimeSteps==0) then
       if (eventtype==RealPhoton .or. eventtype==HiPion .or. eventtype==HiLepton) then
          call Dilep_Decays(0.,pertParticles,1)
       else if (eventType==HeavyIon .or. eventType==hadron) then
          call Dilep_Decays(0.,realParticles,1)
       end if

       if (eventtype==neutrino.or.eventtype==transportGivenParticle) call DoRadiativeDeltaDecay(0,pertParticles)

       call deuteriumPL_assign (realParticles)
    end if

    if (eventtype==Hadron .or. eventType==groundState .or. eventType==ExternalSource) then
       if ( getSMM_Flag() ) then !determine properties of fragmenting source(s)
          call DoSourceAnalysis(realParticles,time,delta_T,.false.,targetNuc)
       endif
    end if

    if (printParticleVectorTime) then
       Select Case(eventType)
       case (HeavyIon,Hadron,groundState)
          call DoHeavyIonAnalysisTime (realParticles, time)
       end select
    endif

    ! event output at initialization (before time evolution)
    call DoEventOutput(realparticles, pertParticles, 0)

    select case (eventType)
    case (inABox)
      call DoBoxAnalysisTime(realParticles, 0)
    case (inABox_pion)
      call count_Pions(pertParticles,time)
    case (inABox_delta)
      call countDeltas(pertParticles,time)
    end select

    TimeStep=0

    PhaseSpaceEvolution : Do while(time<time_max-delta_T/2.)

       if (DoPR(1)) call timeMeasurement(.true.) ! Reset stopWatch

       if(variableTimeStep) then
          Select Case(eventtype)
          case(elementary,Heavyion,hadron)
             ! real-real collisions
             call evaluateTimeStep(1,1.2,delta_T_max,time,delta_T_new)
          case default
             ! pert-real collisions
             call evaluateTimeStep(2,1.,delta_T_max,time,delta_T_new)
          end select
          delta_T=delta_T_new
       end if

       TimeStep=TimeStep+1

       ! Reduce time step in some time-window: *********************
       !if(5. <= time .and. time <= 15.) then
       !  delta_T=delta_T_max/10.
       !else
       !  delta_T=delta_T_max
       !end if
       !************************************************************

       time=time+delta_T

       write(*,'(79("*"),/,5("*")," TimeStep ",i4,": ",f7.3," fm"/,79("*"))') timeStep,time

       !=== Do some analysis/statistics on "per timestep" basis:

       if (povray_switch) call povray_output (timestep, realParticles, pertParticles)

       !=== Some "Per-Time-Step - Analysis"

       select case (eventtype)

       case (RealPhoton)
          call Dilep_Decays(time, pertParticles, timeStep/numTimeSteps)  ! Dilepton Analysis

       case (neutrino)
          call DoRadiativeDeltaDecay(TimeStep,pertParticles)

       case (HiLepton)
          call HiLeptonAnalysisPerTime(time,pertParticles)
          if ( getSMM_Flag() ) then !determine properties of fragmenting source(s)
             call DoSourceAnalysis(realParticles,time,delta_T,.false.,targetNuc)
             if (stopGiBUU) EXIT PhaseSpaceEvolution !stop BUU-run after on-set of equilibration
          endif
          call Dilep_Decays(time, pertParticles, timeStep/numTimeSteps)  ! Dilepton Analysis

       case(HiPion)
          call HiPionAnalysisPerTime(timestep,time,pertParticles)

       case(HeavyIon)
          if ( getSMM_Flag() ) then !determine properties of fragmenting source(s)
             call DoSourceAnalysis(realParticles,time,delta_T,.false.,targetNuc,projectileNuc)
             if (stopGiBUU) EXIT PhaseSpaceEvolution !stop BUU-run after on-set of equilibration
          endif
          call Dilep_Decays (time, realParticles, timeStep/numTimeSteps)  ! Dilepton Analysis

       case(Hadron,groundState,ExternalSource)
          if ( getSMM_Flag() ) then !determine properties of fragmenting source(s)
             call DoSourceAnalysis(realParticles,time,delta_T,.false.,targetNuc)
             if (stopGiBUU) EXIT PhaseSpaceEvolution !stop BUU-run after on-set of equilibration
          endif
          call Dilep_Decays (time, realParticles, timeStep/numTimeSteps)  ! Dilepton Analysis

       case(transportGivenParticle)
          call transportGivenParticle_analyze(pertParticles,timestep)
          call DoRadiativeDeltaDecay(TimeStep,pertParticles)

       case(InABox)
          call DoBoxAnalysisTime(realParticles, timestep)

       case(PionBox)
          call DoPionBoxAnalysisTime(realParticles, timestep)

       end select

       call DoFreezeoutAnalysisPerTime(timestep,time,pertParticles,realParticles)

       ! event output at particular timestep
       call DoEventOutput(realparticles, pertParticles, timestep)

       if (DoPR(1)) call timeMeasurement() ! Print stopWatch

       !=== Elementary checks ===

       call ChecksCallAll(timestep,time,realParticles,pertParticles)

       if (DoPR(1)) call timeMeasurement() ! Print stopWatch

       !=== Garbage Collection ===
       !         Write(*,*) '**Do Garbage Collection'
       call GarbageCollection(pertParticles, .true.)
       call GarbageCollection(realParticles)

       !=== propagation ===
       Write(*,*) '**Do Propagation'
       if (.not. getRMF_flag()) then
         ! Use Skyrme-like potential + MDI
         call Propagate (realParticles, pertParticles, delta_T, TimeStep)
       else
         ! Use Relativistic Mean Fields (RMF):
         call propagate_RMF (realParticles, pertParticles, delta_T, TimeStep)
       end if

       !=== updating mean fields ===
       if(.not.freezeRealParticles) then
           if( .not.getRMF_flag() ) then
              call updateDensity(realParticles)
              call updateCoulomb
              call updateYukawa(.false.)
              !=== update momentum(0) after updating the mean fields ===
              call updateEnergies(realParticles)
           end if
           call updateTemperature(realParticles)
       end if

       if( .not.getRMF_flag() ) call updateEnergies(pertParticles,.true.)

       !=== collisions ===

       Write(*,*) '**Do Collisions'
       call formation(pertParticles,realParticles,time,.false.)
       call nullCountedEvents(1)
       call nulln_participants
       call collideMain(pertParticles,realParticles,time)

       call deuteriumPL_assign (realParticles)

       if (DoPR(1)) call timeMeasurement() ! Print stopWatch

       !=== checks ===

       call ChecksCallEnergy(time,realParticles)

       If(checkGridSize_Flag) then !checks only if particles start to escape out of grid
          call CheckGridSize(realParticles,time,time_max,eventType,gridSpacing,gridSize)
       end if

       if (DoPr(2)) call writeCountedEvents(1,time)
       if (DoPr(1)) call writeCountedEvents(2)

       select case (eventtype)
       case (inABox_pion)
          call count_Pions(pertParticles,time) ! = Count all pions
       case (inABox_delta)
          call countDeltas(pertParticles,time) ! = Count all pions
       case (HeavyIon,groundState,ExternalSource)
          call HeavyIon_evol (realparticles, time)
       case (hadron)
          call HeavyIon_evol (realparticles, time)
          call DoHadronAnalysisTime (realparticles, time, .false.)
       end select


       !=== Print the particle vector frequently (timeSequence) for time > timeForOutput ===
       !    Since each eventclass uses its own format, it was convenient to introduce also here
       !    the select case statement.

       if (printParticleVectorTime) then

          Select Case(eventType)
          case (HeavyIon,Hadron,groundState,ExternalSource)
             call DoHeavyIonAnalysisTime (realParticles, time)
          end select

       endif

       if (DoPR(1)) call timeMeasurement() ! Print stopWatch

    End do PhaseSpaceEvolution

    ! === Forcing decays of all particles at the end ===
    call formation(pertParticles,realParticles,time,.true.)
    call ForceDecays(pertParticles,realParticles,time)

    if( .not.getRMF_flag() ) then
       ! === update momentum(0) after the forced decays ===
       if(.not.freezeRealParticles) call updateEnergies(realParticles)
       call updateEnergies(pertParticles,.true.)
    end if

    !=== Elementary checks ===

    call ChecksCallEnergy(time,realParticles)
    call ChecksCallAll(timestep,time,realParticles,pertParticles)

    ! === Add nucleons from Fragmentation

    if (DoFragmentNucleons) then
       call DoAddFragmentNucleons(realParticles,targetNuc)
    end if

    ! === Again "Per-Time-Step"-Analysis after the forced decays

    select case (eventtype)
    case (HeavyIon)
       if (getSMM_Flag()) then
          ! Reset variables needed for sourceAnalysis for the next subsequent run
          call DoSourceAnalysis(realParticles,time,delta_T,.true.,targetNuc,projectileNuc)
          call resetSMMvalues
       endif

    case (Hadron)
       if (getSMM_Flag()) then
          ! Reset variables needed for sourceAnalysis for the next subsequent run
          call DoSourceAnalysis(realParticles,time,delta_T,.true.,targetNuc)
          call resetSMMvalues
       endif
       call DoHadronAnalysisTime (realparticles, time, .true.)

    case (ExternalSource)
       if (getSMM_Flag()) then
          ! Reset variables needed for sourceAnalysis for the next subsequent run
          call DoSourceAnalysis(realParticles,time,delta_T,.true.,targetNuc)
          call resetSMMvalues
       endif

    case (neutrino,transportGivenParticle)
       call DoRadiativeDeltaDecay(numTimeSteps+1,pertParticles)

    case (HiLepton)
       if (getSMM_Flag()) then
          ! Reset variables needed for sourceAnalysis for the next subsequent run
          call DoSourceAnalysis(realParticles,time,delta_T,.true.,targetNuc)
          call resetSMMvalues
       endif

    end select

    call DoFreezeoutAnalysisPerTime(timestep,time,pertParticles,realParticles)

!     call ChecksCallOccupied(realParticles,pertParticles,'At the end of the run: ')
    if (DoPR(1)) call timeMeasurement() ! Print stopWatch


  end subroutine run



  !*****************************************************************************
  !****s* GiBUU/analysis
  ! NAME
  ! subroutine analysis(finalizeFlag,beforeRUN)
  ! PURPOSE
  ! Analyze output of run.
  ! INPUTS
  ! * logical :: finalizeFlag -- This flag signalises that the last run of
  !   one energy was taking place.
  ! * logical :: beforeRUN    -- Flag to indicate, whether this routine
  !   is called before or after "run" is called. Makes it possible to produce
  !   some analysis-output of the particle vector directly after its init.
  !*****************************************************************************
  subroutine analysis(finalizeFlag,beforeRUN)
    use eventtypes
    use inputGeneral, only: FinalCoulombCorrection, PrintParticleVectors, numTimeSteps
    use CoulombKorrektur, only: CoulombPropagation
    use pionXsection, only: pionXsectionAnalysis
    use HiLeptonAnalysis, only: DoHiLeptonAnalysis
    use HiPionAnalysis, only: DoHiPionAnalysis
    use pionGamma, only: evaluate_pionGamma
    use lowphotonanalysis, only: analyze_Photon
    use HeavyIonAnalysis, only: DoHeavyIonAnalysis
    use ElementaryAnalysis, only: DoElementaryAnalysis
    use lowElectronAnalysis, only: lowElectron_Analyze
    use yScalingAnalysis, only: yScaling_Analyze
    use neutrinoAnalysis, only: neutrino_Analyze
    use initLowPhoton, only: lowPhotonInit_getRealRun
    use initNeutrino, only: neutrinoInit_getRealRun
    use history, only: history_print
    use Dilepton_Analysis, only: Dilep_write_CS
    use radiativeDeltaDecay, only: radiativeDeltaDecay_write_CS
    use EventOutputAnalysis, only: DoEventOutput

    logical, intent(in) :: finalizeFlag, beforeRUN
    integer :: i,j
    logical :: first_historyPrint

    if (beforeRUN) then
       Select Case(eventType)
       Case(HiPion)
          call DoHiPionAnalysis(pertParticles,finalizeFlag,beforeRUN)
       Case(HiLepton)
          call DoHiLeptonAnalysis(realparticles,pertParticles,targetNuc%mass,finalizeFlag,beforeRUN)
       end Select

       return ! leave this routine !!!
    endif

    !********** Final Coulomb propagation

    Select Case(eventType)
    Case(LoPion,RealPhoton,neutrino,HiPion,HiLepton)
       If (FinalCoulombCorrection) call CoulombPropagation (pertParticles)
    end Select

    !********** Analysis

    Select Case(eventType)
    Case(elementary)
       write(*,*) '   Main : Calling elementary collision analysis routine'
       call DoElementaryAnalysis (realparticles, finalizeFlag)
       write(*,*) '   Main : Finished elementary collision analysis routine'
       call DoHeavyIonAnalysis(realparticles,pertParticles,finalizeFlag)

    Case(HeavyIon,hadron,ExternalSource)
       write(*,*) '   Main : Calling heavy ion induced analysis routine'
       call DoHeavyIonAnalysis(realparticles,pertParticles,finalizeFlag)
       call Dilep_write_CS()  ! write dilepton cross sections to file

    Case(LoPion)
       write(*,*) '   Main : Calling pion induced analysis routine'
       call  pionXsectionAnalysis(pertParticles,finalizeFlag)

    case(RealPhoton)
       If(PrintParticleVectors) then
          write(*,'(A)') "The particle vectors before analysis can be found in the files *_preAnalysis_*.dat"
          call writeParticleVector('RealParticles_preAnalysis',realParticles)
          call writeParticleVector('PertParticles_preAnalysis',pertParticles)
       end if
       call Dilep_write_CS()  ! write dilepton cross sections to file
       write(*,*) '   Main : Calling photon induced analysis routine'
       if(lowPhotonInit_getRealRun()) then
          call analyze_Photon(realParticles,finalizeFlag)
       else
          call analyze_Photon(pertParticles,finalizeFlag)
       end if

    Case(LoLepton)
       call lowElectron_Analyze(pertParticles,finalizeFlag)
       call yScaling_Analyze(pertParticles,finalizeFlag,eventType)

    Case(neutrino)
       write(*,*) '   Main : Calling neutrino induced analysis routine'
       if(neutrinoInit_getRealRun()) then
          call neutrino_Analyze(realParticles,finalizeFlag,num_runs_sameEnergy)
       else
          call yScaling_Analyze(pertParticles,finalizeFlag,eventType)
          call neutrino_Analyze(pertParticles,finalizeFlag,num_runs_sameEnergy)
          call radiativeDeltaDecay_write_CS(pertParticles)
       end if

    Case(HiPion)
       call Dilep_write_CS()  ! write dilepton cross sections to file
       call DoHiPionAnalysis(pertParticles,finalizeFlag)

    Case(HiLepton)
       call DoHiLeptonAnalysis(realparticles,pertParticles,targetNuc%mass,finalizeFlag)
       call Dilep_write_CS()  ! write dilepton cross sections to file

    Case(inABox_pion)
       If(finalizeFlag) call evaluate_pionGamma()

    Case(groundState)
       call DoHeavyIonAnalysis(realparticles,pertParticles,finalizeFlag)

    Case(transportGivenParticle)
       call radiativeDeltaDecay_write_CS(pertParticles)

    end Select


    If(PrintParticleVectors) then
       call writeParticleVector('RealParticles_Final',realParticles)
       call writeParticleVector('PertParticles_Final',pertParticles)

       first_historyPrint=.true.
       open(100,file="PertParticles_Final_collHistory.dat")
       Do i=1,Size(PertParticles,dim=1)
          Do j=1,Size(PertParticles,dim=2)
             If (Pertparticles(i,j)%ID > 0) then
                call history_print(i,PertParticles(i,j),100,first_historyPrint)
                first_historyPrint=.false.
             end if
          End do
       end do
       close(100)
       first_historyPrint=.true.
       open(100,file="RealParticles_Final_collHistory.dat")
       Do i=1,Size(realParticles,dim=1)
          Do j=1,Size(realParticles,dim=2)
             If (realparticles(i,j)%ID > 0) then
                call history_print(i,realParticles(i,j),100)
                first_historyPrint=.false.
             end if
          End do
       end do
       close(100)

    end if

    ! final event output after time evolution
    call DoEventOutput(realparticles, pertParticles, numTimeSteps+1)

  End subroutine analysis


  !*************************************************************************
  !****s* GiBUU/finalCleanup
  ! NAME
  ! subroutine finalCleanup
  ! PURPOSE
  ! Clean up all "global" memory. Deallocates particle vectors and projectile/target nuclei
  ! and calls cleanup routines in several modules.
  !*************************************************************************
  subroutine finalCleanup
    use eventtypes
    use baryonWidth, only: cleanupBaryon => cleanUp
    use mesonWidth, only: cleanupMeson => cleanUp
    use densityModule, only: cleanupDensity => cleanup
    use yukawa, only: cleanupYukawa => cleanup
    use thermoDynamics, only: cleanupThermo => cleanup
    use parametrizationsBarMes, only: cleanupBarMes => cleanup
    use initLowPhoton, only: cleanupLowPhoton => cleanup
    use PILCollected, only: PILCollected_DeAllocate
    use selfenergy_baryons, only: cleanupRealParts => cleanup
    use neutrinoAnalysis, only: cleanupNeutrinoAna => cleanup
    use initNeutrino, only: cleanupNeutrinoInit => cleanup
    use baryonWidthMedium_tables, only: cleanupBaryonMedium => cleanup
    use mesonWidthMedium_tables, only: cleanupMesonMedium => cleanup
    use coulomb, only: cleanupCoulomb => cleanup
    use volumeElements, only: cleanupVE => cleanup
    use hadronAnalysis, only: cleanupHadronAna => cleanup

    write (*,*)
    write (*,*) "Cleaning Up All Allocated Memory!"
    write (*,*)
    ! Particle Vectors (real & perturbative)
    DeAllocate(realParticles, pertParticles)
    ! Target & Projectile Nuclei
    if (associated(projectileNuc)) DeAllocate(projectileNuc)
    if (associated(targetNuc)) DeAllocate(targetNuc)
    ! General stuff (Particle Properties, Widths, Cross Sections, Potentials, Densities, ...
    call cleanupBaryon
    call cleanupMeson
    call cleanupDensity
    call cleanupYukawa
    call cleanupThermo
    call cleanupBarMes
    call PILCollected_DeAllocate
    call cleanupRealParts
    call cleanupBaryonMedium
    call cleanupMesonMedium
    call cleanupCoulomb
    call cleanupVE
    ! Inits, Analysis, etc (specific to different eventTypes)
    select case(eventType)
    case(RealPhoton)
      call cleanupLowPhoton
    case(neutrino)
      call cleanupNeutrinoInit
      call cleanupNeutrinoAna
    case(hadron)
      call cleanupHadronAna
    end select
  end subroutine finalCleanup


end program GiBUU
