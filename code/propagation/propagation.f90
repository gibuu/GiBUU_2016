!***************************************************************************
!****m* /propagation
! NAME
! module propagation
! PURPOSE
! Module which includes the propagation of the test-particles.
!***************************************************************************
module propagation
  implicit none

  Private

  !*************************************************************************
  !****g* propagation/predictorCorrector
  ! PURPOSE
  ! Switch for predictor-corrector method in the propagation.
  ! If .false. then simple
  ! Euler method is used.
  ! SOURCE
  !
  logical, save :: predictorCorrector=.true.
  !*************************************************************************

  !*************************************************************************
  !****g* propagation/delta_E
  ! PURPOSE
  ! Delta energy in derivatives
  ! SOURCE
  !
  real, save :: delta_E=0.01
  !*************************************************************************

  !*************************************************************************
  !****g* propagation/delta_P
  ! PURPOSE
  ! Delta Momentum in derivatives
  ! SOURCE
  !
  real, save :: delta_P=0.01
  !*************************************************************************

  !*************************************************************************
  !****g* propagation/UseCoulomb
  ! PURPOSE
  ! Whether to use coulomb force directly in propagation or not. (If switched
  ! off while coulomb is switched on in module coulomb, the effect of the
  ! coulomb potential comes in via the gradient of the potentials. With this
  ! flag you can not switch on/off coulomb, you just select, how it is
  ! treated.)
  ! SOURCE
  !
  logical,  save :: UseCoulomb=.true.
  !*************************************************************************

  !*************************************************************************
  !****g* propagation/UseHadronic
  ! PURPOSE
  ! Whether to use hadronic potentials in propagation
  ! SOURCE
  !
  logical,  save :: UseHadronic=.true.
  !*************************************************************************

  !*************************************************************************
  !****g* propagation/dh_dp0_switch
  ! PURPOSE
  ! Switch which decides whether we use dh_dp0.
  ! SOURCE
  !
  logical,  save :: dh_dp0_switch=.true.
  !*************************************************************************

  !*************************************************************************
  !****g* propagation/DerivativeType
  ! PURPOSE
  ! 1=first order Runge-Kutta, 2=second order Runge-Kutta in derivatives
  ! SOURCE
  !
  integer,  save :: DerivativeType=1
  !*************************************************************************

  !*************************************************************************
  !****g* propagation/offShellInfo
  ! PURPOSE
  ! print out offShellInfo: set to .true. automatically if offShellTransport is used
  ! SOURCE
  !
  logical,save :: offShellInfo=.false.
  !*************************************************************************

  !*************************************************************************
  !****g* propagation/offShellInfoDetail
  ! PURPOSE
  ! print out detailed offShellInfo
  ! SOURCE
  !
  logical,save :: offShellInfoDetail=.false.
  !*************************************************************************

  !*************************************************************************
  !****g* propagation/tachyonDebug
  ! PURPOSE
  ! ...
  ! SOURCE
  !
  logical,save :: tachyonDebug=.false.
  !*************************************************************************

  logical,save :: startFlag = .true.
  logical,parameter :: warnNegativeMass2 = .false.

  !***************************************************************************
  !****s* propagation/updateVelocity
  ! NAME
  ! subroutine updateVelocity(teilchen)
  ! PURPOSE
  ! Updates velocities of a single particle or a field of particles
  ! and checks for v<c.
  ! INPUTS
  ! * type(particle),intent(inOut), dimension(:,:) :: teilchen
  ! or:
  ! * type(particle),intent(inOut), dimension(:)   :: teilchen
  ! or:
  ! * type(particle),intent(inOut) :: teilchen
  !***************************************************************************
  Interface updateVelocity
     Module Procedure updateVelocity_1, updateVelocity_field,updateVelocity_matrix
  End Interface


  Public :: propagate, propagate_euler, propagate_cascade
  Public :: updateVelocity, checkVelo, gradients

contains


  !***************************************************************************
  !****s* propagation/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads out namelist "propagation"
  ! INPUTS
  ! * (none)
  ! OUTPUT
  ! * Initializes global module variables
  !***************************************************************************
  subroutine init
    use output
    use offshellPotential

    integer :: ios

    !*************************************************************************
    !****n* propagation/Propagation
    ! NAME
    ! NAMELIST propagation
    ! PURPOSE
    ! Namelist which includes the input switches:
    ! * delta_P
    ! * delta_E
    ! * UseCoulomb
    ! * UseHadronic
    ! * DerivativeType
    ! * predictorCorrector
    ! * dh_dp0_switch
    ! * offShellInfoDetail
    ! * tachyonDebug
    !*************************************************************************
    NAMELIST /propagation/ delta_P, delta_E, UseCoulomb, UseHadronic, &
                           DerivativeType, predictorCorrector, &
                           dh_dp0_switch, offShellInfoDetail, tachyonDebug

    call Write_ReadingInput('propagation',0)
    rewind(5)
    read(5,nml=propagation,iostat=ios)
    if (ios>0) then
       write(*,*)
       write(*,*) 'Maybe you use an old jobcard.'
       write(*,*) 'Variabes have been renamed:'
       write(*,*) '   coulomb --> UseCoulomb'
       write(*,*) '   hadronic --> Usehadronic'
    end if
    call Write_ReadingInput('propagation',0,ios)

    write(*,*) 'delta_P = ',delta_P
    write(*,*) 'delta_E = ',delta_E
    write(*,*) 'use coulomb force directly = ',UseCoulomb
    write(*,*) 'hadronic potential flag    = ',UseHadronic
    write(*,*) 'DerivativeType             = ',DerivativeType
    write(*,*) 'Predictor-Corrector mode   = ',predictorCorrector
    write(*,*) 'dh_dp0_switch: ',dh_dp0_switch

    if(tachyonDebug) then
       write(*,*) 'Do tachyon debug!!'
    end if

    If(get_useOffShellPotentialBaryons().or.get_useOffShellPotentialMesons()) then
       offShellInfo=.true.
       if(tachyonDebug) offShellInfoDetail=.true.
       write(*,*) 'Detailed off-shell info?', offshellInfoDetail
    end if

    call Write_ReadingInput('propagation',1)

    startFlag = .false.
  end subroutine init




  !***************************************************************************
  !****s* propagation/gradients
  ! NAME
  ! subroutine gradients(teilchen,Grad_P,Grad_R)
  ! PURPOSE
  ! determine gradients of the Hamiltonfunction
  ! INPUTS
  ! * type(particle),intent(in)  :: teilchenIn  --
  !   particle whose gradients we want to calculate
  ! OUTPUT
  ! * real, dimension(1:3) :: Grad_P -- momentum gradient
  ! * real, dimension(0:3) :: Grad_R -- space gradients =(d/dt , d/dr(1:3) )
  !   [Note that there is no minus sign!!]
  !***************************************************************************
  subroutine gradients(teilchenIn,Grad_P,Grad_R)

    use particleDefinition
    use coulomb, only : emfoca
    use densityModule
    use energyCalc
    use offShellPotential
    use minkowski
    use histf90
    use hist2Df90
    use idTable, only : photon
    use derivatives, only : finiteDifference
    use CallStack

    type(particle),intent(in)  :: teilchenIN
    real, intent(out), dimension(1:3) :: Grad_P
    real, intent(out), dimension(0:3),optional :: Grad_R

    real :: delta_R(1:3)
    real,dimension(-2:2,0:3) :: energy
    type(particle) :: teilchen
    integer :: i,j
    real :: cpot                     ! Coulomb potential
    real, dimension(1:3) :: emForce  ! electromagnetic Force
    real :: dH_dp0           ! Derivative with respect to p_0: dH/dp_0
    type(histogram),save :: dh_dp0_hist
    type(histogram2D),save :: dh_dp0_hist2D
    integer, save  :: numCalls=0
    logical,save :: DoInit=.true.
    logical :: outOfBounds!,success
    integer :: deriv_scheme

    numCalls=numCalls+1

    if (startFlag) call init

    if(get_offshell_debug()) then
       if(DoInit) then
          call createHist(dH_dp0_hist,'dH/dp_0',-2.,10.,0.1)
          call createHist2D(dH_dp0_hist2D,'Num events as function of dH/dp_0 and mass',(/0.6,-2./),(/1.3,10./),(/0.04,0.1/))
          DoInit=.false.
       end if
    end if

    ! Check Input
    if(teilchenIn%ID <= 0) then !no particle
       write(*,*) 'ID =',teilchenIn%ID
       call TRACEBACK('Error. Particle-ID<=0  in Gradient.')
    else if (teilchenIn%ID.eq.photon) then
       ! Photons do not undergo any potentials
       if (present(grad_R)) grad_R=0.
       grad_P(1:3)  =  teilchenIn%momentum(1:3)/teilchenIn%momentum(0)
       if(Dot_Product(grad_P,grad_P).gt.1.0000000001) then
          write(*,*) 'Problem with photon in gradients'
          write(*,*) '|Velo|^2=',Dot_Product(grad_P,grad_P)
          write(*,*) 'Velo=',grad_P
          write(*,*) 'Momentum=',teilchenIn%momentum
          return
       end if
    end if

    grad_P=0.

    if (present(grad_R)) then
       grad_R=0.
       If (UseCoulomb) then ! special treatment of electromagnetic part
          cpot = emfoca(teilchenIn%position(1:3),teilchenIn%momentum(1:3),teilchenIn%charge,teilchenIn%ID,emForce)
          grad_R(1:3)=grad_R(1:3)-emForce  !Grad(Potential)=-Force
       end if
    end if

    If (UseHadronic) then !with hadronic potentials

       if (present(grad_R)) then
          ! Determine Gridsize in position space
          delta_R(1:3)=getGridSpacing()

          ! Derivative with respect to r: grad_r H
          teilchen=teilchenIN
          energy=0.
          do i=1,3 !loop over x,y,z
             deriv_scheme=0
             teilchen%position=teilchenIn%position
             do j=-DerivativeType,DerivativeType
                teilchen%position(i)=teilchenIN%position(i)+float(j)*delta_R(i)
                ! determine energy due to hadronic interaction
                if(treatParticleOffShell(teilchen%Id,teilchen%offshellparameter)) then
                   energy(j,i)=hamiltonFunc_offshell(teilchen,outOFBounds)
                   if(outOfBounds) then
                      call OutOfBoundsMessage('r')
                      deriv_scheme=-j
                   end if
                else
                   call energyDetermination(teilchen,ForbidCoulomb=UseCoulomb)
                   energy(j,i)=teilchen%momentum(0)
                end if
             end do
             grad_R(i)=grad_R(i)+finiteDifference(energy(:,i),delta_R(i),derivativeType,deriv_scheme)
          end do
       end if



       ! Derivative with respect to p: grad_p H
       teilchen=teilchenIN
       energy=0.

       do i=1,3 !loop over x,y,z
          deriv_scheme=0
          teilchen%momentum=teilchenIn%momentum
          do j=-DerivativeType,DerivativeType
             teilchen%momentum(i)=teilchenIN%momentum(i)+float(j)*delta_P
             ! determine energy due to hadronic interaction
             if(treatParticleOffShell(teilchen%Id,teilchen%offshellparameter)) then
                energy(j,i)=hamiltonFunc_offshell(teilchen,outOfBounds)
                if(outOfBounds) then
                   call OutOfBoundsMessage('p')
                   deriv_scheme=-j
                end if
             else
                call energyDetermination(teilchen,ForbidCoulomb=UseCoulomb)
                energy(j,i)=teilchen%momentum(0)
             end if
          end do
          grad_P(i)=grad_P(i)+finiteDifference(energy(:,i),delta_P,derivativeType,deriv_scheme)
       End do



       !grad_R(0)=0.

       if(treatParticleOffShell(teilchenIn%Id,teilchenIn%offshellparameter).and.dh_dp0_switch) then
          !****************************************
          ! Derivative with respect to p_0: dH/dp_0
          teilchen=teilchenIN
          energy=0.
          teilchen%momentum=teilchenIn%momentum
          deriv_scheme=0
          do j=-derivativeType,derivativeType
             teilchen%momentum(0)=teilchenIN%momentum(0)+float(j)*delta_E
             ! determine energy due to hadronic interaction
             energy(j,0)=hamiltonFunc_offshell(teilchen,outofBounds)
             if(outOfBounds) then
                call OutOfBoundsMessage('E')
                deriv_scheme=-j
             end if
          end do
          dH_dp0=finiteDifference(energy(:,0),delta_E,derivativeType,deriv_scheme)

          if(offShellInfoDetail) then
             write(*,*) 'dH_dp0=', dh_dp0
             if (.not.(dH_dp0>=0. .or. dH_dp0<=0.)) then
                write(*,*) "NaN in dH_dp0   ID=",teilchenIn%Id, "     offShellParam=",teilchenIn%offshellparameter
                write(*,*) "energy=", energy(:,0)
                stop
             end if
          end if


          ! See Oliver's notes for the case of an energy dependent Hamiltonian:
          if (present(grad_R)) grad_R(1:3)=grad_R(1:3)/(1-dH_dp0)
          grad_P=grad_P/(1-dH_dp0)


          if(get_offshell_debug()) then

             call addHist(dH_dp0_hist,dH_dp0,1.)
             call AddHist2D(dH_dp0_hist2D, (/teilchen%mass,dH_dp0/),1.)

             if(mod(numCalls,1000).eq.0) then
                open(33,file='dH_dp0.dat')
                open(44,file='dH_dp0_2D.dat')
                call writeHist(dh_dp0_hist,33,mul=dH_dp0_hist%xbin)
                call writeHist2D_Gnuplot(dh_dp0_hist2D,44,mul=dH_dp0_hist2D%xbin(1)*dH_dp0_hist2D%xbin(2))
                close(33)
                close(44)
             end if
          end if
       end if
    else !no hadronic potentials
       !grad_R=0.
       grad_P=grad_P+teilchenIN%momentum(1:3)/FreeEnergy(teilchenIN)
    end if

  contains

    subroutine OutOfBoundsMessage(C)

      character*(*), intent(in) :: C
      if(offShellInfoDetail) &
           write(*,'(A,A,i4)') C,': Hamilton function off shell is out of bounds!!',j
      if(j.eq.0) then ! Out of bounds at central value
         write(*,'(A,A,I4,A,I4)') C,': Problem with outofBounds: j=',j
      else if(deriv_scheme*j.gt.0) then ! Out of bounds on both sides
         write(*,'(A,A,I4,A,I4)') C,': Problem with outofBounds: deriv_scheme='&
              & ,deriv_scheme," new deriv_scheme=",-j
      end if

    end subroutine OutOfBoundsMessage


  end subroutine gradients

  !****************************************************************************
  !****s* propagation/propagate
  ! NAME
  ! subroutine propagate (realTeilchen, pertTeilchen, delta_T, timeStep)
  ! PURPOSE
  ! This routine propagates the particle vectors.
  ! INPUTS
  ! * type(particle),intent(inOUT),dimension(:,:)  :: realTeilchen  -- real particle vector which should be propagated
  ! * type(particle),intent(inOUT),dimension(:,:)  :: pertTeilchen  -- perturbative particle vector which should be propagated
  ! * real, intent(in)   :: delta_T    -- time step size (fm/c)
  ! * integer,intent(in) :: timeStep   -- time step number
  ! OUTPUT
  ! * type(particle),intent(inOUT),dimension(:,:)  :: realTeilchen  -- real particle vector which should be propagated
  ! * type(particle),intent(inOUT),dimension(:,:)  :: pertTeilchen  -- perturbative particle vector which should be propagated
  ! NOTES
  ! * Uses predictor corrector scheme or simple Euler time stepping. See also files in "Documentation_Extra/propagation/".
  !****************************************************************************
  subroutine propagate (realTeilchen, pertTeilchen, delta_T, timeStep)
    use particleDefinition
    use inputGeneral, only: freezeRealParticles
    use densityModule, only: updateDensity
    use yukawa, only: updateYukawa
    use coulomb, only: updateCoulomb
    use offshellPotential, only: get_useOffShellPotentialMesons, get_useOffShellPotentialBaryons

    type(particle),intent(inOut),dimension(:,:) :: realTeilchen
    type(particle),intent(inOut),dimension(:,:) :: pertTeilchen
    real, intent(in) :: delta_T
    integer, intent(in) :: timeStep

    logical :: doReal
    ! Gradients of predictor step
    real, dimension(:,:,:), Allocatable :: gradR_real, gradP_real
    real, dimension(:,:,:), Allocatable :: gradR_pert, gradP_pert

    if (startFlag) call init

    doReal = (.not.freezeRealParticles) .or. (timeStep<=1)

    If (PredictorCorrector) then
       ! Allocate fields to store the gradients
       if (doReal) then
          Allocate(gradR_real(0:3, &
               lbound(realteilchen,dim=1):ubound(realteilchen,dim=1), &
               lbound(realteilchen,dim=2):ubound(realteilchen,dim=2)))
          Allocate(gradP_real(1:3, &
               lbound(realteilchen,dim=1):ubound(realteilchen,dim=1), &
               lbound(realteilchen,dim=2):ubound(realteilchen,dim=2)))
       end if
       Allocate(gradR_pert(0:3, &
            lbound(pertteilchen,dim=1):ubound(pertteilchen,dim=1), &
            lbound(pertteilchen,dim=2):ubound(pertteilchen,dim=2)))
       Allocate(gradP_pert(1:3, &
            lbound(pertteilchen,dim=1):ubound(pertteilchen,dim=1), &
            lbound(pertteilchen,dim=2):ubound(pertteilchen,dim=2)))

       !Predictor-Corrector method:
       ! (1) Predictor Step: Define the predicted value of the particle at Delta_T

       if (doReal) call predictorStep(realteilchen,GradP_real,GradR_real,delta_T)
       call predictorStep(pertteilchen,GradP_pert,GradR_pert,delta_T)

       ! (2) Update potentials
       if (doReal) then
          call updateDensity(realteilchen)
          call updateCoulomb
          call updateYukawa(.false.)
       end if

       ! (3) Do time stepping considering also the predicted gradients
       if (doReal) then
         call correctorStep(realteilchen,GradP_real,GradR_real,delta_T)
         DeAllocate(gradR_real,gradP_real)
       end if

       call correctorStep(pertteilchen,GradP_pert,GradR_pert,delta_T)
       DeAllocate(gradR_pert,gradP_pert)

    else
       if (doReal) call propagate_euler(realteilchen,delta_T)
       call propagate_euler(pertteilchen,delta_T)
    end if

    !  ** Set masses properly:
    if (get_useOffShellPotentialBaryons() .or. get_useOffShellPotentialMesons()) then
      if (doReal) call setMass(realteilchen)
      call setMass(pertteilchen)
    end if

    if (doReal) call checkVelos(realteilchen)
    call checkVelos(pertteilchen)

  end subroutine propagate


  !****************************************************************************
  !****s* propagation/setMass
  ! NAME
  ! subroutine setmass(part)
  ! PURPOSE
  ! * Resets %mass according to the offshell parameter and the actual momentum and position of the particle.
  ! * if particles is NOT meant to be treated offshell, nothing is changed.
  ! INPUTS
  ! type(particle),dimension(:,:),intent(inout) :: part
  !
  ! OUTPUT
  ! type(particle),dimension(:,:),intent(inout) :: part
  !****************************************************************************
  subroutine setMass(part)
    use particleDefinition
    use offshellPotential, only:treatParticleOffShell!, getOffShellMass
    use minkowski, only: SP
    use potentialModule, only: massDetermination, scapot

    type(particle),dimension(:,:),intent(inout) :: part
    integer :: i,j
    logical :: success


    do j=lbound(part,dim=2),ubound(part,dim=2)
       do i=lbound(part,dim=1),ubound(part,dim=1)
          if(part(i,j)%ID.le.0) cycle
          if(.not.treatParticleOffShell(part(i,j)%Id,part(i,j)%offshellparameter)) cycle

          call massDetermination(part(i,j),success=success)

          if(.not.success) then
             write(*,*) "Can't set baremass in setMass!!!"
             write(*,*) "Particle is deleted!!!"
             call protocolFailure(part(i,j),'setMass1')
             part(i,j)%ID=0
             cycle
          end if
          if (SP(part(i,j)%momentum,part(i,j)%momentum)<=0.01**2) then ! 10 MeV
             write(*,*)
             write(*,*) 'WARNING in setMass: p^mu p_mu too small! ',SP(part(i,j)%momentum,part(i,j)%momentum)
             write(*,*) 'particle is deleted!!'
             call protocolFailure(part(i,j),'setMass2')
             part(i,j)%ID = 0 ! delete particle
             cycle
          end if
          if (part(i,j)%mass<0.01) then ! 10MeV
             write(*,*)
             write(*,*) 'WARNING in setMass: baremass too small! ',part(i,j)%mass
             write(*,*) 'particle is deleted!!'
             call protocolFailure(part(i,j),'setMass3')
             part(i,j)%ID = 0 ! delete particle
             cycle
          else if (part(i,j)%mass+scapot(part(i,j))>part(i,j)%momentum(0)) then ! inv. mass > energy !
             write(*,*)
             write(*,*) 'WARNING in setMass: mass too large! ',part(i,j)%mass
             write(*,*) 'particle is deleted!!'
             call protocolFailure(part(i,j),'setMass4')
             part(i,j)%ID = 0 ! delete particle
             cycle
          end if


      end do
    end do
  end subroutine setMass


  !****************************************************************************
  !****s* propagation/predictorStep
  ! NAME
  ! subroutine predictorStep(teilchen,gradP,gradR,delta_T)
  ! PURPOSE
  ! * Routine propagates a particle vector using an Euler method and returns the gradients.
  ! INPUTS
  ! * real, intent(in) :: delta_T
  !
  ! OUTPUT
  ! type(particle),intent(inOut),dimension(:,:) :: teilchen -- vector at predicted position and momentum
  ! real, dimension(:,:,:),intent(out) :: gradR, gradP ! Gradients in R and P-Space
  !
  ! NOTES
  ! * Used as  predictor step in a predictor/corrector scheme
  !****************************************************************************
  subroutine predictorStep(teilchen,gradP,gradR,delta_T)
    use offshellPotential,only: get_useOffShellPotentialBaryons,get_useOffShellPotentialMesons
    use particleDefinition
    use minkowski, only :SP
    use output
    use IdTable, only: isBaryon

    real, dimension(1:,:,:),intent(out) :: gradP
    real, dimension(0:,:,:),intent(out) :: gradR
    type(particle),intent(inOut),dimension(:,:) :: teilchen
    real, intent(in) :: delta_T
    integer :: n,m

    Do n=lbound(teilchen,dim=1),ubound(teilchen,dim=1)
       Do m=lbound(teilchen,dim=2),ubound(teilchen,dim=2)
          if (teilchen(n,m)%ID < 0) exit
          if (teilchen(n,m)%ID == 0) cycle

          ! Evaluate dp/dt=-dH/dr and dr/dt=dH/dp
          call gradients(teilchen(n,m),GradP(:,n,m),GradR(:,n,m))

          ! Set the perturbative particle vector to the predictor values:
          Teilchen(n,m)%momentum(1:3)=teilchen(n,m)%momentum(1:3)-delta_T*gradR(1:3,n,m)
          Teilchen(n,m)%momentum( 0 )=teilchen(n,m)%momentum( 0 )+delta_T*gradR(0,n,m)
          Teilchen(n,m)%position(1:3)=teilchen(n,m)%position(1:3)+delta_T*gradP(1:3,n,m)
          Teilchen(n,m)%velocity(1:3)=gradP(1:3,n,m)

          if(SP(teilchen(n,m)%momentum,teilchen(n,m)%momentum) < -epsilon(0.0D0)) then
             if (warnNegativeMass2) then
                write(*,*)'problems in predictorstep, SP(p,p).lt.0:', SP(teilchen(n,m)%momentum,teilchen(n,m)%momentum)
                write(*,*)'gradR = ', gradR(1:3,n,m)
                call writeParticle_debug(teilchen(n,m))
             end if
             if(isBaryon(teilchen(n,m)%ID)) then
                write(*,*) 'Deleting Particle'
                call protocolFailure(teilchen(n,m),'predictor')
                teilchen(n,m)%ID=0
             end if
          end if

       end do
    end do

    if(get_useOffShellPotentialBaryons().or.get_useOffShellPotentialMesons()) call setMass(teilchen)

  end subroutine predictorStep



  !****************************************************************************
  !****s* propagation/correctorStep
  ! NAME
  ! subroutine correctorStep(teilchen,gradP,gradR,delta_T)
  !
  ! PURPOSE
  ! * Routine performs a corrector step on  a particle vector previously undergoing a simple Euler method.
  !
  ! INPUTS
  ! * real, intent(in) :: delta_T
  ! * real, dimension(:,:,:),intent(out) :: gradR, gradP ! Gradients in R and P-Space of predictor step
  ! * type(particle),intent(inOut),dimension(:,:) :: teilchen -- particle vector at predicted values
  !
  ! OUTPUT
  ! type(particle),intent(inOut),dimension(:,:) :: teilchen -- particle vector at t+Delta(t)
  !
  ! NOTES
  ! * Used as  corrector step in a predictor/corrector scheme
  !****************************************************************************
  subroutine correctorStep(teilchen,gradP,gradR,delta_T)
    use particleDefinition
    use densityModule, only : gridsize
    use inputGeneral, only : continousBoundaries

    real, dimension(1:,:,:),intent(in) :: gradP
    real, dimension(0:,:,:),intent(in) :: gradR
    type(particle),intent(inOut),dimension(:,:) :: teilchen
    real, intent(in) :: delta_T
    real, dimension(1:3) :: Grad_P_Predictor
    real, dimension(0:3) :: Grad_R_Predictor
    integer :: i
    integer :: n,m
    logical :: checkVelo_flag

    ! (3) Corrector Step: Evaluate Gradients at predicted value

    Do n=lbound(teilchen,dim=1),ubound(teilchen,dim=1)
       Do m=lbound(teilchen,dim=2),ubound(teilchen,dim=2)
          if (teilchen(n,m)%ID < 0) exit
          if (teilchen(n,m)%ID == 0) cycle

          call gradients(teilchen(n,m),Grad_P_Predictor,Grad_R_Predictor) ! Evaluate dH/dr and dH/dp

          ! (3) Do time stepping considering also the predicted gradients
          !
          !       momentum(1:3,t+delta t)=momentum(1:3,t)-delta_T*(grad_R+grad_R_Predictor)/2.
          !       momentum_predicted(1:3,t+delta t)=momentum(1:3,t)-delta_T*grad_R
          !
          ! => momentum(1:3,t+delta t)=momentum_Predicted(1:3)-delta_T*(-grad_R+grad_R_Predictor)/2.
          ! And be reminded that "teilchen" is the predictor (therefore the - sign in front of gradR and gradP).
          !
          ! The same argument also holds for the position.
          !
          teilchen(n,m)%momentum(1:3)=teilchen(n,m)%momentum(1:3)-delta_T*(-gradR(1:3,n,m)+grad_R_Predictor(1:3))/2.
          teilchen(n,m)%momentum( 0 )=teilchen(n,m)%momentum( 0 )+delta_T*(-gradR( 0, n,m)+grad_R_Predictor( 0 ))/2.
          teilchen(n,m)%position(1:3)=teilchen(n,m)%position(1:3)+delta_T*(-gradP(1:3,n,m)+grad_P_Predictor(1:3))/2.

          ! (4) Save velocity
          teilchen(n,m)%velocity(1:3)=(gradP(1:3,n,m)+grad_P_Predictor(1:3))/2.

          If(continousBoundaries) then
             ! Implement continous boundary conditions
             do i=1,3
                If(teilchen(n,m)%position(i).gt.gridsize(i)) then
                   teilchen(n,m)%position(i)= teilchen(n,m)%position(i)-2*gridsize(i)
                else if(teilchen(n,m)%position(i).lt.-gridsize(i)) then
                   teilchen(n,m)%position(i)= teilchen(n,m)%position(i)+2*gridsize(i)
                end if
             end do
          end if
          if(tachyonDebug) checkVelo_flag=checkVelo(teilchen(n,m))
       end do
    End do

  end subroutine correctorStep



  !****************************************************************************
  !****s* propagation/propagate_euler
  ! NAME
  ! subroutine propagate_euler(teilchen,delta_T)
  ! PURPOSE
  ! * Routine propagates a particle vector
  ! INPUTS
  ! * type(particle),intent(inOUT),dimension(:,:)  :: teilchen  -- particle vector which should be propagated
  ! * real, intent(in) :: delta_T ! time step (fm/c)
  ! OUTPUT
  ! * type(particle),intent(inOUT),dimension(:,:)  :: teilchen  -- particle vector which should be propagated
  ! NOTES
  ! * Uses simple Euler time stepping. See also files in "Documentation_Extra/propagation/".
  !****************************************************************************
  subroutine propagate_euler(teilchen,delta_T)

    use particleDefinition
    use nucleusDefinition
    use densityModule, only : gridsize
    use inputGeneral, only : continousBoundaries,eventtype
    use eventtypes, only: HeavyIon, hadron
    use nucleus, only : getTarget, getProjectile

    type(particle),intent(inOut),dimension(:,:) :: teilchen
    real, intent(in) :: delta_T
    type(tNucleus), pointer :: proj, targ
    integer :: n,m

    real,dimension(0:3) :: grad_R
    real,dimension(1:3) :: grad_P
    integer :: i
    logical :: checkVelo_flag

    if (startFlag) call init

    if (.not.UseHadronic .and. (eventtype==HeavyIon .or. eventtype==hadron)) then
      proj => getProjectile()
      targ => getTarget()
    end if

    Do n=lbound(teilchen,dim=1),ubound(teilchen,dim=1)
       Do m=lbound(teilchen,dim=2),ubound(teilchen,dim=2)
          if (teilchen(n,m)%ID < 0) exit
          if (teilchen(n,m)%ID == 0) cycle

          ! Simple Euler Time stepping:
          if(UseHadronic) then

             call gradients(teilchen(n,m),Grad_P,Grad_R) ! Evaluate dH/dr and dH/dp
             teilchen(n,m)%momentum(1:3)=teilchen(n,m)%momentum(1:3)-delta_T*grad_R(1:3)
             teilchen(n,m)%momentum( 0 )=teilchen(n,m)%momentum( 0 )+delta_T*grad_R( 0 )
             teilchen(n,m)%position(1:3)=teilchen(n,m)%position(1:3)+delta_T*grad_P(1:3)
             teilchen(n,m)%velocity(1:3)=grad_P(1:3)
             if(tachyonDebug) checkVelo_Flag=checkVelo(teilchen(n,m))

          else if (teilchen(n,m)%perturbative .or. teilchen(n,m)%event(1)>=1000000 &
                   .or. (eventtype/=HeavyIon .and. eventtype/=hadron) ) then ! Cascade

             call gradients(teilchen(n,m),Grad_P,Grad_R) ! Evaluate dH/dr and dH/dp
             teilchen(n,m)%position(1:3)=teilchen(n,m)%position(1:3)+delta_T*grad_P(1:3)
             teilchen(n,m)%velocity(1:3)=grad_P(1:3)


          else  ! "Frozen" cascade (presently for eventtypes HeavyIon and hadron only !)

             if(mod(teilchen(n,m)%event(1),2)==1) then

                  if(teilchen(n,m)%ID.ne.1) then
                     write(*,*)' (1) wrong particle propagated:'
                     write(*,*)teilchen(n,m)%ID,teilchen(n,m)%event
                     stop
                  endif

                  teilchen(n,m)%position(1:3) = teilchen(n,m)%position(1:3) + delta_T * targ%velocity(1:3)
                  teilchen(n,m)%velocity(1:3) = targ%velocity(1:3)

             else if( mod(teilchen(n,m)%event(1),2)==0) then

                  !if(teilchen(n,m)%ID.ne.1) then
                  !   write(*,*)' (2) wrong particle propagated:'
                  !   write(*,*)teilchen(n,m)%ID,teilchen(n,m)%event
                  !   stop
                  !endif

                  teilchen(n,m)%position(1:3) = teilchen(n,m)%position(1:3) + delta_T * proj%velocity(1:3)
                  teilchen(n,m)%velocity(1:3) = proj%velocity(1:3)

             else

                  call gradients(teilchen(n,m),Grad_P,Grad_R) ! Evaluate dH/dr and dH/dp
                  teilchen(n,m)%position(1:3)=teilchen(n,m)%position(1:3)+delta_T*grad_P(1:3)
                  teilchen(n,m)%velocity(1:3)=grad_P(1:3)

             end if

          end if

          If(continousBoundaries) then
             ! Implement continous boundary conditions
             do i=1,3
                If(teilchen(n,m)%position(i).gt.gridsize(i)) then
                   teilchen(n,m)%position(i)= teilchen(n,m)%position(i)-2*gridsize(i)
                else if(teilchen(n,m)%position(i).lt.-gridsize(i)) then
                   teilchen(n,m)%position(i)= teilchen(n,m)%position(i)+2*gridsize(i)
                end if
             end do
          end if

       end do
    end do

  end subroutine propagate_euler

  !****************************************************************************
  !****s* propagation/propagate_cascade
  ! NAME
  ! subroutine propagate_cascade(teilchen,delta_T)
  ! PURPOSE
  ! Routine propagates particles in case of cascade mode. Useful also for doing initial step in the RMF mode.
  ! INPUTS/Results
  ! type(particle), dimension(:,:) :: teilchen  ! particles which should be propagated
  ! real, intent(in) :: delta_T                 ! time step (fm/c)
  ! NOTES
  ! Straight line trajectories, no momentum change, since mean fields are switched off.
  !****************************************************************************
  subroutine propagate_cascade(teilchen,delta_T)
    use nucleusDefinition
    use particleDefinition
    use inputGeneral, only : eventtype
    use eventtypes, only: HeavyIon
    use nucleus, only : getTarget, getProjectile

    type(particle), dimension(:,:) :: teilchen  ! particles which should be propagated
    real, intent(in) :: delta_T                 ! time step (fm/c)
    type(tNucleus), pointer :: proj, targ
    integer :: i, j

    if (startFlag) call init

    if (eventtype==HeavyIon) then
      proj => getProjectile()
      targ => getTarget()
    end if

    Loop_over_ensembles : Do i=1,Size(Teilchen,dim=1)
      Loop_over_particles : Do j=1,Size(Teilchen,dim=2)

         If ( teilchen(i,j)%id == 0 ) then
            cycle Loop_over_particles
         else If( teilchen(i,j)%id < 0 ) then
            exit Loop_over_particles
         end If

         if (eventtype==HeavyIon .and. teilchen(i,j)%event(1)<1000000) then

           ! In case of HIC the particle is assumed to be "frozen"
           ! until it collides with other particle, i.e.
           ! it propagates with the speed of either target or projectile.

           if( mod(teilchen(i,j)%event(1),2) == 1 ) then

             if(teilchen(i,j)%ID.ne.1) then
               write(*,*)' In propagate_cascade (1) wrong particle propagated:'
               write(*,*) teilchen(i,j)%ID, teilchen(i,j)%event
               stop
             endif

             teilchen(i,j)%velocity(1:3) = targ%velocity(1:3)

           else if( mod(teilchen(i,j)%event(1),2) == 0 ) then

             if(teilchen(i,j)%ID.ne.1) then
               write(*,*)' In propagate_cascade (2) wrong particle propagated:'
               write(*,*)teilchen(i,j)%ID,teilchen(i,j)%event
               stop
             endif

             teilchen(i,j)%velocity(1:3) = proj%velocity(1:3)

           end if

         else

           teilchen(i,j)%velocity(1:3) = teilchen(i,j)%momentum(1:3) / FreeEnergy( teilchen(i,j) )

         end if

         teilchen(i,j)%position(1:3) = teilchen(i,j)%position(1:3) + delta_T * teilchen(i,j)%velocity(1:3)

      end do Loop_over_particles
    end do Loop_over_ensembles

  end subroutine propagate_cascade


  !****************************************************************************
  !****s* propagation/checkVelo
  ! NAME
  ! subroutine checkVelo(part)
  ! PURPOSE
  ! * Checks the velocity
  ! INPUTS
  ! type(particle),intent(inout) :: part
  !****************************************************************************
  logical function checkVelo(part)
    use IdTable, only: rho,omegaMeson,phi,photon
    use particleDefinition
    use minkowski, only: abs4
    use offshellPotential, only: treatParticleOffShell, get_offshell_debug, get_useOffShellPotentialMesons
    use hist2Df90
    use output
    use callStack,only : traceback
    use mediumDefinition
    use mediumModule, only: mediumAt

    type(particle),intent(inout) :: part
    !integer :: n,m

    logical, save :: first=.true.
    integer, save :: failures=0
    integer, save :: numtry=0
    type(histogram2D),save :: hist2D_rho,hist2D_omega,hist2D_phi
    type(medium) :: med

    checkVelo=.true.
    if(part%ID.le.0 .or. part%ID==photon) return
    numtry=numtry+1

    if (startFlag) call init

    if (first .and. get_offshell_debug() .and. get_useOffShellPotentialMesons()) then
      call createHist2D(hist2D_rho,'number of tachyons as function of mass and momentum',(/0.,0./),(/2.,2./),(/0.02,0.02/))
      call createHist2D(hist2D_omega,'number of tachyons as function of mass and momentum',(/0.,0./),(/2.,2./),(/0.02,0.02/))
      call createHist2D(hist2D_phi,'number of tachyons as function of mass and momentum',(/0.,0./),(/2.,2./),(/0.02,0.02/))
      first = .false.
    end if

    if (get_offshell_debug() .and. get_useOffShellPotentialMesons()) then
      select case (part%ID)
      case (rho)
        call AddHist2D(hist2D_rho,(/abs4(part%momentum),absMom(part)/),0.,1.)
      case (omegaMeson)
        call AddHist2D(hist2D_omega,(/abs4(part%momentum),absMom(part)/),0.,1.)
      case (phi)
        call AddHist2D(hist2D_phi,(/abs4(part%momentum),absMom(part)/),0.,1.)
      end select
    end if

    if(1. - Dot_Product(part%velocity(1:3),part%velocity(1:3)) .le. 0.) then

       failures=failures+1
       checkVelo=.false.

       med = mediumAt(part%position)

       write(*,*)
       write(*,'(A,G12.5)')'Problems in CheckVelo : v = ',sqrt( Dot_Product(part%velocity(1:3),part%velocity(1:3)))
       write(*,'(A,G12.5,A,G12.5,A,G12.5)') 'm = ',abs4(part%momentum),'; pabs = ',absMom(part), &
                                            '; rho = ',med%densityproton+med%densityneutron
       call WriteParticle(6,0,0,part)
       write(*,*)'Number of failures :',failures,' (',float(failures)/float(numtry)*100., '%)'
       call protocolFailure(part,'checkVelo')

       if (get_offshell_debug() .and. get_useOffShellPotentialMesons()) then
         select case (part%ID)
         case (rho)
           call AddHist2D(hist2D_rho,(/abs4(part%momentum),absMom(part)/),1.)
         case (omegaMeson)
           call AddHist2D(hist2D_omega,(/abs4(part%momentum),absMom(part)/),1.)
         case (phi)
           call AddHist2D(hist2D_phi,(/abs4(part%momentum),absMom(part)/),1.)
         end select

         if(mod(failures,100)==0) then
           call writeHist2D_Gnuplot(hist2D_rho,44,file='Tachyons2D_rho.dat')
           call writeHist2D_Gnuplot(hist2D_omega,45,file='Tachyons2D_omega.dat')
           call writeHist2D_Gnuplot(hist2D_phi,46,file='Tachyons2D_phi.dat')
         end if
       end if

       if(treatParticleOffShell(part%Id,part%offshellparameter)) then
          !ACCEPTABLE ONLY FOR OFFSHELL TRANSPORT, OTHERWISE SERIOUS ERROR!!!!!!!!!!
          write(*,*) 'this particle is now deleted!'
          part%ID=0
       else
          write(*,*) 'checkvelo: stop'
          call traceback()
          stop
       end if
    end if

  end function checkvelo


  !****************************************************************************
  !****s* propagation/checkVelos
  ! NAME
  ! subroutine checkVelos(part)
  ! PURPOSE
  ! * Checks the velocities
  ! INPUTS
  ! type(particle),dimension(:,:),intent(inout) :: part
  !****************************************************************************
  subroutine checkVelos(part)
    use particleDefinition

    type(particle),dimension(:,:),intent(inout) :: part
    integer :: n,m
    logical :: c

    do m=lbound(part,dim=2),ubound(part,dim=2)
       do n=lbound(part,dim=1),ubound(part,dim=1)
          if(part(n,m)%ID<=0) cycle
          c=checkVelo(part(n,m))
       end do
    end do

  end subroutine checkVelos


  !***************************************************************************
  ! cf. Interface updateVelocity
  !***************************************************************************
  subroutine updateVelocity_matrix(teilchen)
    use particleDefinition
    use output, only: DoPR

    real,dimension(1:3)  :: grad_P
    type(particle),intent(inOut), dimension(:,:) :: teilchen
    integer :: i,j
    logical :: success

    if (startFlag) call init

    if (DoPr(2)) Write(*,*) 'Updating particle velocities'

    ensLoop: Do i=lbound(teilchen, dim=1),ubound(teilchen, dim=1)
       Do j=lbound(teilchen, dim=2),ubound(teilchen, dim=2)
          if(teilchen(i,j)%ID <= 0) cycle ensLoop
          call gradients(teilchen(i,j),Grad_P) ! Evaluate dH/dp
          teilchen(i,j)%velocity(1:3)=grad_P
          success=checkVelo(teilchen(i,j))
       end do
    end do ensLoop
  end subroutine updateVelocity_matrix



  !***************************************************************************
  ! cf. Interface updateVelocity
  !***************************************************************************
  subroutine updateVelocity_field(teilchen)
    use particleDefinition

    real,dimension(1:3)  :: grad_P
    type(particle),intent(inOut), dimension(:) :: teilchen
    integer :: i
    logical :: success

    if(startFlag) call init

    Do i=lbound(teilchen, dim=1),ubound(teilchen, dim=1)
       if(teilchen(i)%ID <= 0) cycle
       call gradients(teilchen(i),Grad_P) ! Evaluate dH/dp
       teilchen(i)%velocity(1:3)=grad_P
       success=checkVelo(teilchen(i))
    end do
  end subroutine updateVelocity_field

  !***************************************************************************
  ! cf. Interface updateVelocity
  !***************************************************************************
  subroutine updateVelocity_1(teilchen,success)
    use particleDefinition
    use output, only : writeParticle_debug

    real,dimension(1:3)  :: grad_P
    type(particle),intent(inOut) :: teilchen
    logical,intent(out),optional :: success
    logical :: c

    if (startFlag) call init

    if(teilchen%ID <= 0) return
    call gradients(teilchen,Grad_P) ! Evaluate dH/dp
    teilchen%velocity(1:3)=grad_P

    c=checkVelo(teilchen)

    if (present(success)) then
       success = c
       if(tachyonDebug.and..not.success) then
          call writeParticle_debug(teilchen)
          stop
       end if
    end if
  end subroutine updateVelocity_1



  subroutine protocolFailure(part,kindOfFailure)
    use particleDefinition

    character(*)   ,intent(in) :: kindOfFailure
    type(particle) ,intent(in) :: part

    logical, save :: firstTime=.true.
    character(20) :: form
    integer, save :: numFailures=0
    form='(A40,4G13.5)'

    numFailures=numFailures+1

    if(firstTime) then
       open(22,file='propa_failures.txt')
       write(22,*)
       firstTime=.false.
    else
       open(22,file='propa_failures.txt',position='append')
    end if
    write(22,'(A)')'Problems in ', kindOfFailure
    write(22,form)'|v| = '                      , sqrt(Dot_Product(part%velocity(1:3),part%velocity(1:3)))
    write(22,form)'Particle ID: '               , part%ID
    write(22,form)'Particle mass: '             , part%mass
    write(22,form)'Particle offshellparameter: ', part%offshellparameter
    write(22,form)'Particle number: '           , part%number
    write(22,form)'Particle firstevent: '       , part%firstevent
    write(22,form)'Particle velocity: '         , part%velocity
    write(22,form)'Particle momentum: '         , part%momentum
    write(22,form)'Particle perturbative?: '    , part%perturbative
    write(22,form)'Particle charge: '           , part%Charge
    write(22,form)'Particle position: '         , part%position
    write(22,*)   'Number of failures: '        , numFailures
    write(22,*)
    close(22)
  end subroutine protocolFailure


end module propagation
