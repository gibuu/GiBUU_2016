!*************************************************************************
!****m* /initNucleus_in_PS
! NAME
! module initNucleus_in_PS
!
! PURPOSE
! Module which establishes the test-particle representation of a nucleus.
!*************************************************************************
module initNucleus_in_PS

  implicit none
  PRIVATE

  !*************************************************************************
  !****g* initNucleus_in_PS/improved_monteCarlo
  ! SOURCE
  !
  logical,save :: improved_monteCarlo=.false.
  !
  ! PURPOSE
  ! * If this flag is set to .true. then we use the information of the already
  !   initialized nucleons to decide on the position of a nucleon which has to
  !   be initialized.
  ! * This prescription does only work properly if the smearing with is really
  !   small. Therefore it is switched off by default.
  !*************************************************************************


  !*************************************************************************
  !****g* initNucleus_in_PS/improved_monteCarlo_speedup
  ! SOURCE
  !
  integer,save :: improved_monteCarlo_speedup=500
  !
  ! PURPOSE
  ! * If improved_monteCarlo  is set to .true. then this variable defines
  !   the speedup of the algorithm.
  ! * The number defines how often the density field is updated.
  ! * A large value of  this parameter yields a less accurate test-particle
  !   distribution and a faster initialization.
  !*************************************************************************


  !*************************************************************************
  !****g* initNucleus_in_PS/HiTail
  ! SOURCE
  !
  logical,save :: HiTail = .false.
  !
  ! PURPOSE
  ! If HiTail is set to .true., then a simple parametrization of n(p)
  ! is used to initialize the nucleon momenta
  ! (cf. function chooseAbsMomentum for details).
  !*************************************************************************


  !*************************************************************************
  !****g* initNucleus_in_PS/determine_Fermi_momentum_by_binding_energy
  ! SOURCE
  !
  logical, save :: determine_Fermi_momentum_by_binding_energy=.false.
  !
  ! PURPOSE
  ! If set to .true., the Fermi momentum will determined by
  ! E_B=p_f**2/(2m)+U(rho,p_F),
  ! where E_B is the binding energy per nucleon.
  !*************************************************************************


  !*************************************************************************
  !****g* initNucleus_in_PS/determine_Fermi_new_NucDLDA
  ! SOURCE
  !
  logical, save :: determine_Fermi_new_NucDLDA=.false.
  !
  ! PURPOSE
  ! If set to .true., the Fermi momentum will be set to a value
  ! such that there are no unbound nucleons at the initialisation.
  !*************************************************************************


  !*************************************************************************
  !****g* initNucleus_in_PS/useEnergySF
  ! SOURCE
  !
  logical, save :: useEnergySF=.false.
  !
  ! PURPOSE
  ! If set to .true., then a spectral function is used to choose the energy.
  !*************************************************************************


  !*************************************************************************
  !****g* initNucleus_in_PS/compressedFlag
  ! SOURCE
  !
  logical, save :: compressedFlag=.false.
  !
  ! PURPOSE
  ! If set to .true., then a spherically deformed nucleus is initialized
  ! (isotropic compression/expansion; protons & neutrons in phase).
  ! This type of deformation corresponds to a giant-monopol resonance mode.
  !*************************************************************************


  !*************************************************************************
  !****g* initNucleus_in_PS/ScaleFactor
  ! SOURCE
  !
  Real, save :: ScaleFactor=1.
  !
  ! PURPOSE
  ! If compressedFlag=.true., then rescale coordinates by ScaleFactor.
  !*************************************************************************


  !*************************************************************************
  !****g* initNucleus_in_PS/useCdA
  ! SOURCE
  !
  logical, save :: useCdA = .false.
  !
  ! PURPOSE
  ! Instead of the usual momentum distribution according a fermi gas,
  ! use the momentum parametrizations as given in:
  ! * C. Ciofi degli Ati, S. Simula, PRC 53, 1689 (1996)
  ! These exist only for 2H,3He,4He,12C,160 40Ca,56Fe,208Pb
  !*************************************************************************


  !*************************************************************************
  !****g* initNucleus_in_PS/z_shift
  ! SOURCE
  !
  real, parameter :: z_shift = 100.
  !
  ! PURPOSE
  ! Artificial shift along z-axis between projectile & target nuclei to avoid
  ! overlapping between them when the density is updated.
  ! NOTES
  ! Relevant only for heavy-ion collisions.
  !*************************************************************************


  public :: initNucPhaseSpace
  public :: shiftBack


  logical, save :: initFlag = .true.
  integer, save :: MomChooseType


contains


  !*************************************************************************
  !****s* initNucleus_in_PS/init
  ! NAME
  ! subroutine init
  !
  ! PURPOSE
  ! initalize the module
  !*************************************************************************
  subroutine init
    use output, only: Write_ReadingInput
    use baryonPotentialModule, only: getsymmetryPotFlag_baryon, getPotentialEQSType
    use RMF, only: getRMF_flag, g_rho

    !***************************************************************************
    !****n* initNucleus_in_PS/InitNucleus_in_PS
    ! NAME
    ! Namelist /InitNucleus_in_PS/
    ! PURPOSE
    ! Includes the switches:
    ! * improved_monteCarlo
    ! * improved_monteCarlo_speedup
    ! * HiTail
    ! * determine_Fermi_momentum_by_binding_energy
    ! * determine_Fermi_new_NucDLDA
    ! * useEnergySF
    ! * compressedFlag
    ! * ScaleFactor
    ! * useCdA
    !***************************************************************************
    NAMELIST /initNucleus_in_PS/ improved_monteCarlo, improved_monteCarlo_speedup, HiTail, &
                                 determine_Fermi_momentum_by_binding_energy, determine_Fermi_new_NucDLDA, &
                                 useEnergySF, compressedFlag, ScaleFactor, useCdA

    integer :: ios

    call Write_ReadingInput('initNucleus_in_PS',0)
    rewind(5)
    read(5,nml=initNucleus_in_PS,iostat=ios)
    call Write_ReadingInput('initNucleus_in_PS',0,ios)

    ! Decide how the momentum for A>2 has to be choosen:
    MomChooseType = 0
    if (.not. getRMF_flag()) then
        if (determine_Fermi_momentum_by_binding_energy) then
           MomChooseType = 1 ! = by binding energy
        else
           if (getsymmetryPotFlag_baryon()) then
              MomChooseType = 2 ! = with isospin asymmetry
           else
              if (getPotentialEQSType().eq.6.and.determine_Fermi_new_NucDLDA) then
                 MomChooseType = 4 ! = without isospin asymmetry, new_NucDLDA
              else
                 MomChooseType = 3 ! = without isospin asymmetry
              endif
           end if
        end if
        if (useCdA) MomChooseType = 5
    else
        if (g_rho/=0.) then
           MomChooseType = 2 ! = with isospin asymmetry
        else
           MomChooseType = 3 ! = without isospin asymmetry
        endif
    end if

    write(*,*) 'Type of momentum selection (for A>2):'
    select case (MomChooseType)
    case (1)
       write(*,*) '(1) fermi gas distribution, pF by binding energy'
    case (2)
       write(*,*) '(2) fermi gas distribution, p and n separate'
    case (3)
       write(*,*) '(3) fermi gas distribution, isospin symmetric'
    case (4)
       write(*,*) '(4) fermi gas distribution, isospin symmetric, new_NucDLDA'
    case (5)
       write(*,*) '(5) parametrization Ciofi degli Atti & Simula'
    end select

    if(HiTail) write(*,*) 'Nucleon momenta are initialized with a high momentum tail.'
    write(*,*)

    write(*,'(A,L8)') ' Improved Monte Carlo method?',improved_monteCarlo
    if(improved_monteCarlo) write(*,*)'# Speedup=',improved_monteCarlo_speedup

    if(compressedFlag) then
       write(*,*) 'Initialized nucleus is comressed!'
       write(*,*) 'Compression factor: ScaleFactor = ',ScaleFactor
    end if

    call Write_ReadingInput('initNucleus_in_PS',1)

    initFlag = .false.

  end subroutine init



  !***************************************************************************
  !****s* initNucleus_in_PS/initNucPhaseSpace
  ! NAME
  ! subroutine initNucPhaseSpace(teilchen,nuc)
  !
  ! PURPOSE
  ! Represents nucleus 'nuc' in phase space by testparticles which
  ! are stored in vector 'teilchen'.
  ! The ordering in the vector teilchen is choosen to be random.
  !
  ! bahaviour according mass:
  ! * A=1: Elementary event, no fermi motion.
  ! * A=2,Z=1: We call the deuterium routine, which does a special momentum
  !   distribitution.
  ! * A>2: Uses fermi gas distribution of nucleons for all nuclei.
  !
  ! NOTES
  ! For A>2 it is checked that the radius of the nucleus is properly initialized.
  ! The particles are initialized with %perturbative=.false. since all nuclei must be
  ! in the real particle vector.
  !
  ! All nucleons in a nucleus get the same "%event" number.
  ! It is given by the variable "integer, save ::eventnummer" in the code.
  ! This variable is raised by one after each nucleus-initialization. Therefore the
  ! first nucleus gets "%event=1", the second "%event=2" and so on...
  ! INPUTS
  ! * type(nucleus) :: nuc
  ! * type(particle),dimension(:,:),intent(inout) :: teilchen
  ! OUTPUT
  ! type(particle),dimension(:,:),intent(inout) :: teilchen
  !***************************************************************************
  subroutine initNucPhaseSpace (teilchen, nuc)
    use particleDefinition
    use IdTable, only: nucleon
    use nucleusDefinition
    use densityModule, only: densityAt, get_densitySwitch, set_densitySwitch, FermiMomAt, updateDensity
    use dichtedefinition
    use random, only: rn, rnOmega
    use deuterium, only : initDeuterium
    use potentialModule, only : massDetermination
    use output, only: Write_InitStatus
    use inputGeneral, only : fullensemble,numTimeSteps,eventtype
    use eventtypes, only: hadron, groundState, HeavyIon
    use energyCalc, only: energyDetermination
    use CallStack, only: TRACEBACK
    use insertion, only : GarbageCollection
    use constants, only: mN

    type(tNucleus),pointer :: nuc
    type(particle),dimension(:,:),intent(inout) :: teilchen

    integer :: producedProtons ! Counts number of produced protons
    integer :: i,k, index, offset,chargeSave,numberSave
    real :: maxDist,maxP,maxN
    integer :: densitySwitch_Save

    integer, save :: eventNummer=1
    ! All nucleons in a  nucleus get the same "%event" number. It is given
    ! by eventnummer. This variable is raised by one after each nucleus-initialization.
    ! Therefore the first nucleus gets "%event=1", the second "%event=2" and so on...

    if (initFlag) call init

    call Write_InitStatus("Nucleons in Nucleus",0)

    Write (*,*) 'Nucleus: charge=',nuc%charge,' mass=',nuc%mass, 'FermiMotion=', nuc%fermiMotion
    if (sum(Abs(nuc%velocity))>0.) write (*,'(A,4g12.5)') ' velocity=', nuc%velocity
    if (sum(Abs(nuc%position))>0.) write (*,'(A,4g12.5)') ' position=', nuc%position

    If (nuc%Mass==1) then

       write(*,*) '...Initialising elementary nucleus!!!'
       Do i=1,size(teilchen,dim=1)  !Loop over all ensembles
          index=1
          Do !Search for empty space in particle vector
             if (index>size(teilchen,dim=2)) call TRACEBACK('Real particle vector too small.')
             if (teilchen(i,index)%ID <= 0) exit
             index=index+1
          end do
          call setToDefault(teilchen(i,index))   ! set teilchen to its default values
          if (i==1) then
             call setNumber(teilchen(i,index)) ! give each particle a unique number
          else
             if (fullEnsemble) then
                teilchen(i,index)%number=teilchen(1,index)%number
                ! give each particle a unique number, has to be the same
                ! in all ensembles in full ensemble mode.
             else
                call setNumber(teilchen(i,index)) ! give each particle a unique number
             end if
          end if
          Teilchen(i,index)%event=eventNummer
          teilchen(i,index)%charge=nuc%charge
          teilchen(i,index)%position=0.
          Teilchen(i,index)%ID=Nucleon
          Teilchen(i,index)%antiparticle=.false.
          Teilchen(i,index)%perturbative=.false.
          Teilchen(i,index)%productionTime=0.
          Teilchen(i,index)%mass=mN
          Teilchen(i,index)%momentum=(/mN,0.,0.,0./)
          call boostIt (nuc, Teilchen(i,index))
       End do

    else if (nuc%Mass==2 .and. nuc%charge==1) then

       write(*,*)'...Number Ensembles=',size(teilchen(:,1))
       write(*,*)'...Number particles perEnsembles=',size(teilchen(1,:))

       call initDeuterium (teilchen, nuc, eventNummer, fullEnsemble, z_shift)

    else if (.not. nuc%DoInit) then

       if (numTimeSteps/=0 .and. ((get_densitySwitch()==2 .and. nuc%densityswitch_static==0) .or. get_densitySwitch()==0)) then
          write(*,*) 'zero static density for nuclei makes no sense for runs with FSI (numTimeSteps>0)'
          write(*,*) '-> STOP',numTimeSteps,get_densitySwitch(),nuc%densityswitch_static
          call TRACEBACK()
       end if

       if (improved_MonteCarlo) then
          call updateDensity(teilchen)
          densitySwitch_Save=get_densitySwitch()
          call set_densitySwitch(1)
          call updateDensity(teilchen)
       end if

       ! Search for MC parameters...

       maxDist = nuc%MaxDist
       maxP = nuc%MaxDens(1)
       maxN = nuc%MaxDens(2)

       write(*,*)'...Number Ensembles=',size(teilchen(:,1))
       write(*,*)'...Number particles perEnsembles=',size(teilchen(1,:))

       if (eventType==groundState) eventNummer=1

       if(fullensemble) then
          !************************************
          ! In full ensemble we take care that
          ! every particle has N (=number of
          ! ensemble) test particles, which all
          ! have the same "number", such that
          ! so called self interactions do not
          ! occur in scatterings.
          !************************************
          producedProtons=0
          Do k=1,nuc%mass   !Loop over all particles in nucleus
             Do i=1,size(teilchen,dim=1)  !Loop over all ensembles
                offset=0
                Do !Search for empty space in particle vector
                   index=k+offset
                   if (index>size(teilchen,dim=2)) call TRACEBACK('Real particle vector too small.')
                   if (teilchen(i,index)%ID <= 0) exit
                   offset=offset+1
                end do
                call setToDefault(teilchen(i,index)) !set teilchen to its default values
                Teilchen(i,index)%event=eventNummer
                if (eventType==groundState) eventNummer=eventNummer+1
                Teilchen(i,index)%ID=Nucleon
                Teilchen(i,index)%antiparticle=.false.
                Teilchen(i,index)%perturbative=.false.
                Teilchen(i,index)%productionTime=0.
                Teilchen(i,index)%mass=mN
                If (i==1) then
                   Teilchen(i,index)%charge = chooseCharge(k)
                   call setNumber(teilchen(i,index)) ! give each particle a unique number
                   chargeSave=teilchen(i,index)%charge
                   numberSave=teilchen(i,index)%number
                else
                   Teilchen(i,index)%charge= chargeSave
                   Teilchen(i,index)%number= numberSave
                end if
                Teilchen(i,index)%position = choosePosition()
                if (improved_MonteCarlo) then
                   if (mod(k*i,improved_monteCarlo_speedup)==0)  call updateDensity(teilchen)
                end if
             End do
             if (improved_MonteCarlo) then
                if (mod(k,2)==0)  write(*,*) '  nucleons', k,' of ' ,nuc%mass
             end if
          End do
          If (producedProtons/=nuc%charge) then
             Write(*,*) 'Problem in initPhaseSpace', producedProtons, nuc%charge
          end if
       else
          !************************************
          ! Parallel ensembles
          !************************************
          Do i=1,size(teilchen,dim=1)  !Loop over all ensembles
             producedProtons=0
             offset=0
             Do k=1,nuc%mass   !Loop over all particles in nucleus
                Do !Search for empty space in particle vector
                   index=k+offset
                   if (index.gt.size(teilchen,dim=2)) &
                        & call TRACEBACK('Real particle vector too small.')
                   if (teilchen(i,index)%ID <= 0) exit
                   offset=offset+1
                end do
                call setToDefault(teilchen(i,index)) !set teilchen to its default values
                call setNumber(teilchen(i,index)) ! give each particle a unique number
                Teilchen(i,index)%event=eventNummer
                if (eventType==groundState) eventNummer=eventNummer+1
                Teilchen(i,index)%ID=Nucleon
                Teilchen(i,index)%antiparticle=.false.
                Teilchen(i,index)%perturbative=.false.
                Teilchen(i,index)%productionTime=0.
                Teilchen(i,index)%mass=mN
                Teilchen(i,index)%charge = chooseCharge(k)
                Teilchen(i,index)%position = choosePosition()
                if (improved_MonteCarlo) then
                   if (mod(k*i,improved_monteCarlo_speedup)==0)  call updateDensity(teilchen)
                end if
             End do
             if (improved_MonteCarlo) then
                if (mod(i,100)==0) write(*,*) '  Ensemble',i,' of ',size(teilchen,dim=1)
             end if
             If (producedProtons/=nuc%charge) then
                Write(*,*) 'Problem in initPhaseSpace', producedProtons, nuc%charge
             end if
          End do
       end if


       !    For dynamical density the density field must be updated according
       !    to test particle positions:
       call updateDensity(teilchen)

       Do i=1,size(teilchen,dim=1)
          Do index=1,size(teilchen,dim=2)
             if (Teilchen(i,index)%ID < 0) exit
             if (Teilchen(i,index)%ID==Nucleon) then
                if (Teilchen(i,index)%event(1)==eventNummer .or. eventType==groundState) then
                   call chooseMomentum (Teilchen(i,index))
                   call boostIt (nuc, Teilchen(i,index))
                end if
             end if
          End do
       End do
       If (improved_MonteCarlo) call set_densitySwitch(densitySwitch_Save)


       if(useEnergySF) call chooseEnergySF


    else

       write(*,*) 'Problem in initNucPhaseSpace!!! Radius of nucleus is not well defined:',&
            & nuc%radius, nuc%mass, nuc%charge
       call Traceback()
    end If

    ! raise Event number such that the next nucleus which is initialized gets a higher %event :
    eventNummer=eventNummer+1
    if (eventtype==hadron) eventNummer=eventNummer+1 ! here += 2 !

    call Write_InitStatus("Nucleons in Nucleus",1)

    call GarbageCollection(teilchen)

  contains

    !*************************************************************************

    integer function chooseCharge (k)
      integer, intent(in) :: k
      real :: probabilityProton
      !Choose randomly the charge of the nucleon
      probabilityProton=float(nuc%charge-producedProtons)/float(nuc%mass-k+1)
      if (rn()<=probabilityProton) then
         chooseCharge=1
         producedProtons=producedProtons+1
      else
         chooseCharge=0
      end if
    end function chooseCharge

    !*************************************************************************

    function choosePosition() result(r)
      use dichtedefinition
      use densityStatic, only: staticDensity

      real, dimension(1:3) :: r

      real :: maxDist2
      type(dichte) :: density,densityMC

      densityMC%proton=0.
      densityMC%neutron=0.
      maxDist2 = maxDist**2

      Do
         Do !Monte Carlo distribution of position in sphere with radius maxDist
            r(1)=(1.-2.*rn())*maxDist
            r(2)=(1.-2.*rn())*maxDist
            r(3)=(1.-2.*rn())*maxDist
            If (r(1)**2+r(2)**2+r(3)**2 .le. maxDist2) exit
         End do
         density=staticdensity(r,nuc) !Evaluate density at choosen position
         if(improved_MonteCarlo) densityMC=densityAt(r)
         !Monte Carlo decision whether position is to be assigned:
         If (teilchen(i,index)%charge.eq.1) then
            if ((rn()*maxP).le.Max(density%proton(0)-densityMC%proton(0),0.)) exit
         else
            if ((rn()*maxN).le.Max(density%neutron(0)-densityMC%neutron(0),0.)) exit
         End if
      End do
      ! Prepare compressed nucleus:
      if (compressedFlag) r=(1.+ScaleFactor)*r
    end function choosePosition


    !*************************************************************************

    subroutine chooseMomentum (part)
      use constants, only: mN
      use dichteDefinition
      use determine_fmbybE, only: Teilchen_Fermi, determine_fermimomentum, determine_fermiNucDLDA
      use RMF, only : walecka, getRMF_flag

      type(particle) :: part

      real :: pFermi, pAbs
      real, dimension(1:3) :: place, p
      real :: dens, shift
      real, parameter :: cutoff=0.8 !momentum cutoff
      type(dichte) :: density

      if (.not.nuc%fermiMotion) then
         if (getRMF_flag()) call TRACEBACK('RMF and no Fermi not implemented.')
         part%momentum=(/mN,0.,0.,0./)
         part%velocity(1:3)= 0.0
         return
      end if

      place=part%position
      density = densityAt(place)

      ! Please note, only MomChooseType=2,3 respect the rescaling
      ! of the density, if ReAdjustForConstBinding

      select case(MomChooseType)
      case (1) ! = by binding energy
         dens=density%baryon(0)
         call determine_fermimomentum(dens,pFermi)

      case (2) ! = with isospin asymmetry
         pFermi = FermiMomAt(place,part%charge)

      case (3) ! = without isospin asymmetry
         pFermi = FermiMomAt(place)

      case (4) ! = without isospin asymmetry, new_NucDLDA
         dens=density%baryon(0)
         call determine_fermiNucDLDA(place,dens,pFermi)

      case (5) ! = parametrization Ciofi degli Atti & Simula
         call ChooseFermiMomCdA(nuc%mass,pAbs)
         pFermi = pAbs
      end select

      if (MomChooseType.ne.5) then
         ! Monte Carlo distribution of momentum in sphere of radius pFermi:
         do
            pAbs = pFermi * chooseAbsMomentum()
            if (pAbs .lt. cutoff) exit
         end do
      end if
      p(1:3)=pAbs*rnOmega()

      part%momentum(1:3)=p

      If ( (pAbs > pFermi) .and. (Hitail) ) then
         Teilchen_Fermi=part
         Teilchen_Fermi%momentum(1:3)=(/pFermi,0.,0./)
         call energyDetermination(Teilchen_Fermi)
         part%momentum(0)=Teilchen_Fermi%momentum(0)
         call massDetermination(part)
      end if


      if( getRMF_flag() ) then ! RMF used
         call walecka(density%baryon(0),shift)
         part%momentum(0) = Sqrt(dot_product(p,p)+(mN-shift)**2)
      else
         part%momentum(0) = Sqrt(dot_product(p,p)+mN**2)
      end if

      part%velocity(1:3)=part%momentum(1:3)/part%momentum(0)

    end subroutine chooseMomentum


    !*************************************************************************
    !****f* initNucPhaseSpace/chooseAbsMomentum
    ! NAME
    ! real function chooseAbsMomentum
    ! PURPOSE
    ! return the absolute value of the momentum in units of pFermi
    !
    ! NOTES
    ! we have to initalize according p^2*n(p) !!!
    !
    ! normally the value of p lies between 0 and 1 and n(p)=1.
    !
    ! Via the flag 'HiTail' you can switch into a mode,
    ! where the distribution n(p) is given by
    ! * 0.83             for p = 0 .. 1
    ! * h*exp(-A*p)/p^2  for p > 1
    ! .
    !
    ! The functional form of the tail is choosen to allow for simple
    ! random number generation.
    !
    ! The slope parameter of the tail is fitted to calculations
    ! of Kalok
    !
    !*************************************************************************
    function chooseAbsMomentum()
      use constants, only :rhoNull
      use dichteDefinition

      real, parameter :: v1 = 0.83         ! for lehr distribution 0.85  ! n(p) for p=0..1
      real, parameter :: v2 = 0.24 ! for lehr distribution 0.15  ! n(p) for p==1.000001
      !      real, parameter :: A = 3.1784 ! slope for exp(-A p)
      real, parameter :: A = 3.67 ! for lehr distribtution 2.3    ! slope for exp(-A p)/p^2
      real :: chooseAbsMomentum
      real :: h1, h2
      real :: h,u

      real, dimension(1:3) :: place
      type(dichte) :: density
      real :: dens
      if (.not.HiTail) then
         chooseAbsMomentum = rn()**(1./3.)
      else
         place=Teilchen(i,index)%position

         density=densityAt(place)

         dens=density%baryon(0)
         if (dens .lt. 0.00001)  then
            chooseabsMomentum=0.
            return
         end if


         h1 = v1/3*(dens/rhoNull)**(-0.03)   ! lehr(-0.04)
         h2 = h1 + v2/A*(dens/rhoNull)**(0.27)  !lehr(0.37)

         h = v2*exp(A)*(dens/rhoNull)**(0.27) ! scaling of hiTail
         u = rn()

         chooseAbsMomentum = (h2/h1*u)**(1./3.)

         if (chooseAbsMomentum.gt.1.0) then
            chooseAbsMomentum = -log(h2*A/h*(1-u))/A
         endif
      endif
    end function chooseAbsMomentum


    !*************************************************************************

    subroutine boostIt (nucl, part)
      use lorentzTrafo, only: lorentz

      type(tnucleus), intent(in) :: nucl
      type(particle) :: part

      !Do nothing if nucleus rests in calculation frame :
      If (sum(Abs(nucl%velocity))<0.000001) then
         ! Write(*,*) 'No boost necessary! Velocity of nucleus:',nucl%velocity
         ! Shift to actual position of center of mass of nucleus
         part%position(1:3)=part%position(1:3)+nucl%position
      else
         ! boost members of nucleus from restframe to calculation frame,
         ! which is moving with -beta seen from the rest frame:
         call lorentz (-nucl%velocity, part%momentum(0:3))
         If ((Abs(nucl%velocity(1))>0.1).or.(Abs(nucl%velocity(2))>0.1)) then
            Write(*,*) 'Error in initNucPhaseSpace'
            Write(*,*) 'Velocity of nucleus in transverse direction is not negligible'
            Write(*,*) 'Velocity in x,y:' , nucl%velocity(1:2)
            Write(*,*) 'Nucleus: A=', nucl%mass,'Z=', nucl%charge
            Write(*,*) 'critical error. Stop program!'
            stop
         end if
         !Length contraction of nucleus in z-direction(linear scaling)
         !Neglect length contraction in transverse direction
         part%position(3)=part%position(3)*sqrt(1-nucl%velocity(3)**2)
         !Shift to actual position of center of mass of nucleus
         part%position(1:3)=part%position(1:3)+nucl%position
      end if

      ! shift particles along z-axis to avoid overlapping between projectile & target:
      if (eventType==HeavyIon) then
         if (eventNummer==1) then
            part%position(3)=part%position(3)+z_shift
         else
            part%position(3)=part%position(3)-z_shift
         endif
      endif

    end subroutine boostIt



    !*************************************************************************
    !****s* initNucPhaseSpace/chooseEnergySF
    ! NAME
    ! subroutine chooseEnergySF
    ! PURPOSE
    ! set the energies and masses of the realparticles according to a recipe of
    ! Ankowski et al (arXiv:0711.2031v2). The energy is choosen according to a
    ! spectral function for the initial state nucleons which accounts for
    ! the shell structure.
    !
    ! NOTES
    ! * Only implemented for oxygen!
    ! * This does not affect the momentum distribution!
    !*************************************************************************
    subroutine chooseEnergySF
      use inputGeneral, only : LRF_equals_CALC_frame
      use potentialModule, only : potential_LRF
      use hist2Df90
      use histf90

      real :: Erem,pvecsq,V,effmass
      logical :: success
      type(histogram2D), save :: p_E_distri_SF,removal
      type(histogram), save :: massspectra

      call CreateHist2D(p_E_distri_SF,'p_E_distri_SF',(/0.,0./),(/.5,3./),(/0.005,0.005/))
      call CreateHist2D(removal,'removal',(/0.,0./),(/.5,.1/),(/0.001,0.001/))
      call CreateHist(massspectra,'massspectra',0.8,1.3,0.001)

      Loop_over_ensembles : Do i=1,size(teilchen,dim=1)
         Loop_over_particles : Do index=1,size(teilchen,dim=2)
            if(Teilchen(i,index)%ID==Nucleon) then
               if(Teilchen(i,index)%event(1)==eventNummer) then

                  if(.not.LRF_equals_CALC_frame) then
                     write(*,*) 'STOP: initNucPhaseSpace: chooseEnergySF works only in the "LRF=Calculation frame" assumption'
                     stop
                  end if

                  if(nuc%mass.eq.16) then
                     do
                        Erem=get_RemovalEnergy_Wr_Oxygen(Teilchen(i,index)%charge)
                        if(Erem.gt.0) exit
                     end do
                  else
                     write(*,*) 'initNucPhaseSpace: chooseEnergySF only implemented for 16-O so far -> STOP',nuc%mass
                     stop
                  end if

                  V=potential_LRF(Teilchen(i,index))
                  pvecsq=Dot_Product(Teilchen(i,index)%momentum(1:3),Teilchen(i,index)%momentum(1:3))

                  effmass=sqrt((sqrt(mN**2+pvecsq)+V)**2-pvecsq)
                  Teilchen(i,index)%momentum(0)=-Erem+effmass-V

                  !non-relativistic prescription (equivalent!)
                  !effmass=mN-Erem-pvecsq/2./mN
                  !Teilchen(i,index)%momentum(0)=sqrt(effmass**2+pvecsq)

                  call massDetermination(Teilchen(i,index),success=success)

                  if(.not.success) then
                     write(*,*) 'initNucPhaseSpace: chooseEnergySF, massDetermination success=.false.'
                  end if

                  call AddHist(massspectra,Teilchen(i,index)%mass, 1.)
                  call AddHist2D(p_E_distri_SF,(/sqrt(pvecsq),Teilchen(i,index)%momentum(0)/), 1.)
                  call AddHist2D(removal,(/sqrt(pvecsq),Erem/),1.)

               end if
            else if(Teilchen(i,index)%ID < 0) then
               exit Loop_over_particles
            end if
         End do Loop_over_particles
      End do Loop_over_ensembles

      Open(10,file='massspectra.dat')
      call WriteHist(massspectra,10)
      close(10)

      Open(10,file='p_E_distri_SF.dat')
      call WriteHist2D_Gnuplot(p_E_distri_SF,10)
      close(10)

      Open(10,file='removal.dat')
      call WriteHist2D_Gnuplot(removal,10)
      close(10)

    end subroutine chooseEnergySF

    real function get_RemovalEnergy_Wr_Oxygen(charge)
      !  Wroclaw-Fit
      use random, only: rnGauss
      integer, intent(in) :: charge
      real :: Erem,rd,stdDev,mean
      integer :: i
      real, dimension(1:3) :: E,D

      if(charge.eq.1) then
         !in MeV
         E=(/45.0,18.44,12.11/)
         D=(/70.,4.,4./)
      else if(charge.eq.0) then
         !in MeV
         E=(/47.0,21.80,15.65/)
         D=(/70.,4.,4./)
      end if

      !throw dice to choose peak
      !Wroclaw recipe: each shell comes with the same probability
      !but: possible number of nucleons in each shell depends on quantum numbers - not taken into account!
      rd=rn()
      if(rd.le.1./3.) then
         i=1
      else if(rd.le.2./3.) then
         i=2
      else
         i=3
      end if

      mean=E(i)
      stdDev=D(i)/sqrt(16.)

      Erem=rnGauss(stdDev,mean)

      !convert removal energy from MeV in GeV
      Erem=Erem/1000.

      get_RemovalEnergy_Wr_Oxygen=Erem

    end function get_RemovalEnergy_Wr_Oxygen

  end subroutine initNucPhaseSpace

  !***************************************************************************
  subroutine ChooseFermiMomCdA(nucA,k)
    use constants, only: hbarc
    use random, only :rn

    integer, intent(in) :: nucA
    real, intent(out) :: k

!     integer :: i
    real :: x,y
    real, parameter :: ymax = 0.06 ! a good dummy value
    real, parameter :: xmax = 4.0 ! in fm^-1
!     integer :: n

    do
       x = rn()*xmax * hbarc ! now in GeV
       y = x**2*FermiMomCdA(nucA,x)
       if (ymax*rn().lt.y) exit
    end do
    k = x

!!$    n = 208
!!$    ymax = 0.0
!!$    do i=1,400
!!$       x = i*0.01*0.197
!!$       y = FermiMomCdA(n,x)
!!$       write(2000+n,*) x,y,y*x**2
!!$       if (y*x**2.gt.ymax) ymax = y*x**2
!!$    end do
!!$    write(*,*) ymax
!!$    stop

  end subroutine ChooseFermiMomCdA

  !***************************************************************************
  !****f* initNucleus_in_PS/FermiMomCdA
  ! NAME
  ! real function FermiMomCdA(nucA,k)
  !
  ! PURPOSE
  ! calculate the parametrizations according:
  ! * C. Ciofi degli Ati, S. Simula, PRC 53, 1689 (1996)
  ! These exist only for 2H,3He,4He,12C,160 40Ca,56Fe,208Pb
  !
  ! INPUTS
  ! * integer :: nucA -- the size of the nucleus
  ! * real :: k -- the momentum in GeV
  !
  ! OUTPUT
  ! The function value n(k) in fm^3
  !***************************************************************************
  real function FermiMomCdA(nucA,k)
    use constants, only: hbarc

    integer, intent(in) :: nucA
    real, intent(in) :: k

    real,parameter :: A(6,7) = RESHAPE ( (/ &
         & 31.7, 1.32, 5.98, 0.00266, 0.365, 0.0, &! 3He
         & 4.33, 1.54, 0.419, 5.49, 4.90, 0.0,    &! 4He
         !--
         & 2.61, 2.66, 3.54,   0.0,  0.0,  0.0, &!  12C
         & 2.74, 3.33, 6.66,   0.0,  0.0,  0.0, &!  16O
         & 3.24, 3.72,  0.0,  11.1,  0.0,  0.0, &!  40Ca
         & 3.57, 4.97,  0.0,  19.8, 15.0,  0.0, &!  56Fe
         & 1.80, 4.77,  0.0,  25.5,  0.0, 40.3  &! 208Pb
         &/) , (/6,7/) )

    real,parameter :: B(4,7) = RESHAPE ( (/ &
         & 0.0,   0.0,  0.0,    0.0,  &!  3He
         & 0.665, 2.15, 0.0244, 0.22, &!  4He
         & 0.426, 1.60, 0.0237, 0.22, &!  12C
         & 0.326, 1.40, 0.0263, 0.22, &!  16O
         & 0.419, 1.77, 0.0282, 0.22, &!  40Ca
         & 0.230, 1.20, 0.0286, 0.22, &!  56Fe
         & 0.275, 1.01, 0.0304, 0.22  &! 208Pb
         &/) , (/4,7/) )

    integer, parameter :: iiA(8) = (/3,4,12,16,40,56,208,999/)

    integer :: iA
    real :: k2,n0,n1

    FermiMomCdA = 0.0
    k2 = (k/hbarc)**2 ! now k^2 in fm^-2
    do iA = 1,8
       if (iiA(iA).eq.nucA) exit
    end do
    if (iA.eq.8) then
       write(*,*) 'Error: CdA parametrization not available for A=',nucA
       stop
    end if

    if (iA.le.2) then
       n0 = A(1,iA)*exp(-A(2,iA)*k2)/(1.+A(3,iA)*k2)**2 &
            & + A(4,iA)*exp(-A(5,iA)*k2)/(1.+A(6,iA)*k2)**2
    else
       n0 = A(1,iA)*exp(-A(2,iA)*k2) &
            & *(1.+A(3,iA)*k2+A(4,iA)*k2**2+A(5,iA)*k2**3+A(6,iA)*k2**4)
    endif
    if (iA.le.1) then
       n1 = 7.40*exp(-1.23*k2)/(1.+3.21*k2)**2+0.0139*exp(-0.234*k2)
    else
       n1 = B(1,iA)*exp(-B(2,iA)*k2)+B(3,iA)*exp(-B(4,iA)*k2)
    endif

    FermiMomCdA = n0+n1

  end function FermiMomCdA


  subroutine shiftBack(teilchen)
    use particleDefinition

    type(particle),dimension(:,:),intent(inout) :: teilchen

    integer :: i,j

    Do i=1,size(teilchen,dim=1)
       Do j=1,size(teilchen,dim=2)
          if(Teilchen(i,j)%ID < 0) exit
          if (Teilchen(i,j)%position(3) > 0.0) then
             Teilchen(i,j)%position(3) = Teilchen(i,j)%position(3) - z_shift
          else
             Teilchen(i,j)%position(3) = Teilchen(i,j)%position(3) + z_shift
          endif
       End do
    End do

  end subroutine shiftBack


end module initNucleus_in_PS
