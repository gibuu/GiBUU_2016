!*******************************************************************************
!****m* /HeavyIonAnalysis
! NAME
! module HeavyIonAnalysis
!
! PURPOSE
! Contains output routines for heavy ion collisions.
!*******************************************************************************
module HeavyIonAnalysis

  implicit none
  PRIVATE


  !*****************************************************************************
  !****g* HeavyIonAnalysis/flag_outputReal
  ! PURPOSE
  ! If .true., then the output of the real particle vector
  ! will be written to the file 'DoHIA.dat'.
  ! SOURCE
  !
  logical, save :: flag_outputReal = .false.
  !*****************************************************************************


  !*****************************************************************************
  !****g* HeavyIonAnalysis/flag_outputPert
  ! PURPOSE
  ! If .false., then the output of the perturbative particle vector
  ! will be written to the file 'DoHIA_pert.dat'.
  ! SOURCE
  !
  logical, save :: flag_outputPert = .false.
  !*****************************************************************************


  !*****************************************************************************
  !****g* HeavyIonAnalysis/flag_outputDetailed
  ! PURPOSE
  ! Print out more detailed information at each time step
  ! from subroutine HeavyIon_evol:
  ! * rhorad_*.dat
  ! * rhoz_*.dat
  ! * rhozx_*.dat
  ! * Fields_*.dat
  ! * pauli_*.dat
  ! * dens_max.dat
  ! SOURCE
  !
  logical, save :: flag_outputDetailed = .false.
  !*****************************************************************************


  !*****************************************************************************
  !****g* HeavyIonAnalysis/pionAnalysis
  ! PURPOSE
  ! This flag generates various pion spectra (p_T, m_T, y, etc).
  ! The analysis operates under the assumption of a fixed target,
  ! and expects the collision to be performed in the CMS system
  ! (cf. cmsFlag in namelist /heavyIon/).
  ! The analysis matches the one applied to the HADES data in
  ! Agakishiev et al., Eur.Phys.J. A40 (2009) 45-49.
  ! SOURCE
  !
  logical, save :: pionAnalysis = .false.
  !*****************************************************************************


  !*****************************************************************************
  !****g* HeavyIonAnalysis/rapBinning
  ! PURPOSE
  ! Rapidity binning for the pion analysis (only used if pionAnalysis = .true.).
  ! The numbers represent the binning borders in y0. For each of the seven
  ! y0 bins, a separate mT spectrum will be generated.
  ! SOURCE
  !
  real, dimension(0:7), save :: rapBinning = (/ -0.75, -0.45, -0.15, 0.15, 0.45, 0.75, 1.05, 1.35 /)
  !*****************************************************************************


  !*****************************************************************************
  !****g* HeavyIonAnalysis/KaonAnalysis
  ! PURPOSE
  ! This flag generates various Kaon spectra and Kaon-related analyses.
  ! SOURCE
  !
  logical, save :: KaonAnalysis = .false.
  !*****************************************************************************


  logical, save :: initFlag=.true.


  PUBLIC :: DoHeavyIonAnalysis, DoHeavyIonAnalysisTime, HeavyIon_evol


contains


  !*****************************************************************************
  !****s* HeavyIonAnalysis/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads in the namelist "HICanalysis_Input"
  ! INPUTS
  ! * (none)
  ! OUTPUT
  ! * Initializes global module variables
  !*****************************************************************************
  subroutine init
    use output, only: Write_ReadingInput

    integer :: ios

    !***************************************************************************
    !****n* HeavyIonAnalysis/HICanalysis_Input
    ! NAME
    ! NAMELIST /HICanalysis_Input/
    ! PURPOSE
    ! Includes the switches:
    ! * flag_outputReal
    ! * flag_outputPert
    ! * flag_outputDetailed
    ! * pionAnalysis
    ! * rapBinning
    ! * KaonAnalysis
    !***************************************************************************
    NAMELIST /HICanalysis_Input/ flag_outputReal, flag_outputPert, flag_outputDetailed, &
                                 pionAnalysis, rapBinning, KaonAnalysis

    call Write_ReadingInput('HICanalysis_Input',0)
    rewind(5)
    read(5,nml=HICanalysis_Input,iostat=ios)
    call Write_ReadingInput('HICanalysis_Input',0,ios)

    write(*,*) 'flag_outputReal    : ', flag_outputReal
    write(*,*) 'flag_outputPert    : ', flag_outputPert
    write(*,*) 'flag_outputDetailed: ', flag_outputDetailed
    write(*,*) 'pionAnalysis       : ', pionAnalysis
    if (pionAnalysis) write(*,'(A,8F7.3)') ' rapBinning         : ', rapBinning
    write(*,*) 'KaonAnalysis       : ', KaonAnalysis


    call Write_ReadingInput('HICanalysis_Input',1)

    initFlag=.false.

  end subroutine init


  !*************************************************************************
  !****s* HeavyIonAnalysis/DoHeavyIonAnalysis
  ! NAME
  ! subroutine DoHeavyIonAnalysis(realparticles,pertParticles,finalFlag)
  !
  ! PURPOSE
  ! Makes the output of the test particle id's, charges, masses, positions and momenta
  ! on disk. Also does some statistical analysis (see subroutine histo1).
  !
  ! INPUTS
  ! * type(particle), dimension(:,:), intent(in)  :: realparticles,pertParticles  --
  !   real and perturbative particle vectors
  ! * logical, intent (in)     :: finalFlag   ! .true. if it is the last call for one specific
  !   energy, therefore final output must be made.
  !
  ! RESULT
  ! The test particle infos are printed into file 'DoHIA.dat'.
  ! The statistical results are printed to the files:
  ! * 'DoHIA1.dat'
  ! * 'DoHIA2.dat'
  !*************************************************************************
  subroutine DoHeavyIonAnalysis (realparticles, pertParticles, finalFlag)

    use IdTable, only: nucleon, pion, EOV, NOP, isMeson
    use particleDefinition
    use initHeavyIon, only : b_HI => b
    use initHadron,   only : b_had => b, particleId, antiparticle, perturbative
    use initElementary, only : b_ele => impactParameter
    use twoBodyStatistics, only : sqrts_distribution, rate
    use particleProperties, only: isStrange
    use inputGeneral, only : eventtype
    use eventtypes, only: HeavyIon, hadron, elementary
    use history, only: history_getParents

    type(particle), dimension(:,:), intent(in), target  :: realparticles
    type(particle), dimension(:,:), intent(in)          :: pertParticles
    logical,                        intent(in)          :: finalFlag

    integer :: i,j,k,numensembles,indFree,parents(1:3)
    real :: factor,stossParameter
    integer, save :: isu=0     ! counter of subsequent runs
    type(particle), dimension(1:2) :: dummy

    if (initFlag) call init

    isu = isu + 1

    if (eventtype==HeavyIon) then
       stossParameter = b_HI
    else if(eventtype==hadron) then
       stossParameter = b_had
    else if(eventtype==elementary) then
       stossParameter = b_ele    ! If negative -- actual impact parameter is randomly chosen
       ! (not shows up here)
    else
       write(*,*)' Problem in DoHeavyIonAnalysis: impact parameter is not defined'
       stossParameter=0.
    end if

    numensembles=size(realParticles,dim=1)


    !******** Real particle output: ***************************************
    if (flag_outputReal) then

       !*************************************************************************
       !****o* HeavyIonAnalysis/DoHIA.dat
       ! NAME
       ! file DoHIA.dat
       ! PURPOSE
       ! Contains the full dump of the real particle vector at the end of the simulation,
       ! including particle ID, charge, position, momentum, etc.
       !*************************************************************************
       open(30,file='DoHIA.dat',position='Append')

       do i = 1,numensembles

          do j = 1,size(realParticles,dim=2)

             if(realParticles(i,j)%ID <= 0) cycle

             if(realParticles(i,j)%antiparticle) then
                factor = -1.
             else
                factor = 1.
             endif

             if(isFree(realParticles,i,j,realParticles(i,j))) then
                indFree=1
             else
                indFree=0
             end if

             parents = history_getParents(realParticles(i,j)%history)
             do k=1,2
                if(parents(k)>200) parents(k)=200-parents(k)
             end do

             write(30,5) int(factor)*realParticles(i,j)%ID,realParticles(i,j)%charge,&
                  realParticles(i,j)%event(1),parents(1:2),realParticles(i,j)%mass,&
                  indFree,realParticles(i,j)%position(1:3),realParticles(i,j)%momentum(1:3),&
                  i,isu,stossParameter
5            format(i4,1x,i2,1x,i8,2(1x,i4),2x,f6.3,1x,i1,3(1x,f8.3),3(1x,f8.3),1x,i5,1x,i4,1x,f5.2)

          end do

       end do

       close(30)
    end if


    !******** Perturbative particle output: ***************************************
    if (flag_outputPert) then

       !*************************************************************************
       !****o* HeavyIonAnalysis/DoHIA_pert.dat
       ! NAME
       ! file DoHIA_pert.dat
       ! PURPOSE
       ! Contains the full dump of the perturbative particle vector at the end of the simulation,
       ! including particle ID, charge, position, momentum, etc.
       !*************************************************************************
       open(30,file='DoHIA_pert.dat',position='Append')

       do i = 1,numensembles

          do j = 1,size(pertParticles,dim=2)

             if (pertParticles(i,j)%ID == EOV) exit
             if (pertParticles(i,j)%ID == NOP .or. pertParticles(i,j)%ID == nucleon) cycle

             if(pertParticles(i,j)%antiparticle) then
                factor = -1.
             else
                factor = 1.
             endif

             if(isFree(realParticles,i,j,pertParticles(i,j))) then
                indFree=1
             else
                indFree=0
             end if

             parents = history_getParents(pertParticles(i,j)%history)
             do k=1,2
                if(parents(k)>200) parents(k)=200-parents(k)
             end do

             write(30,6) int(factor)*pertParticles(i,j)%ID,pertParticles(i,j)%charge,&
                  pertParticles(i,j)%event(1),parents(1:2),pertParticles(i,j)%mass,&
                  indFree, (pertParticles(i,j)%momentum(k), k=0,3),&
                  pertParticles(i,j)%perWeight,i,isu,stossParameter
6            format(i4,1x,i2,1x,i10,2(1x,i4),2x,f8.5,1x,i1,4(1x,f10.5),1x,e13.6,1x,i5,1x,i4,1x,f5.2)

          end do

       end do

       close(30)

    end if

    if (finalFlag) then
       call sqrts_distribution(dummy,0,.true.)
       call rate(dummy,dummy,0.,.true.)
    end if

    if (eventtype==hadron .and. particleId==nucleon .and. antiparticle .and. .not.perturbative ) &
         call histo1

    if (pionAnalysis) call analyze_Pions
    if (KaonAnalysis) call analyze_Kaons

  contains

    subroutine histo1

      integer, parameter :: N_max=200            ! Maximum number of pions produced in an event
      real, save :: P_Npion(0:N_max)             ! Pion multiplicity distribution
      integer, parameter :: Nmom=100             ! Number of momentum bins
      real, parameter :: dmom=0.02               ! Momentum bin (GeV/c)
      real, save, dimension(1:Nmom,0:10) :: dNpiondMom  ! Pion momentum distribution
      real, save :: pion_events

      integer :: numPions,ibin
      logical :: flag_not_only_pions
      real :: fnorm,pinumAv,momentumAbs
      type(particle), POINTER :: pPart


      if(isu.eq.1) then
         pion_events=0.
         P_Npion=0.
         dNpiondMom=0.
         open(36,file='FewPionEvents.dat')
      end if

      Ensemble_loop : do i = 1,numensembles

         numPions=0
         flag_not_only_pions=.false.

         Particle_loop1 : do j = 1,size(realParticles,dim=2)

            if(realParticles(i,j)%ID <= 0) cycle Particle_loop1

            if(realParticles(i,j)%ID.eq.pion) then
               numPions=numPions+1
            else if(isMeson(realParticles(i,j)%ID)) then
               flag_not_only_pions=.true.
            end if

         End do Particle_loop1

         if(.not.flag_not_only_pions) pion_events=pion_events+1.

         if(.not.flag_not_only_pions) then

            if(numPions.le.N_max)  P_Npion(numPions)=P_Npion(numPions)+1.

            Particle_loop2 : do j = 1,size(realParticles,dim=2)

               if(realParticles(i,j)%ID <= 0) cycle Particle_loop2

               if(numPions.le.1) then
                  if(realParticles(i,j)%antiparticle) then
                     factor = -1.
                  else
                     factor = 1.
                  endif
                  write(36,5) int(factor)*realParticles(i,j)%ID,realParticles(i,j)%charge,&
                       realParticles(i,j)%mass,&
                       (realParticles(i,j)%position(k), k=1,3),&
                       (realParticles(i,j)%momentum(k), k=1,3),&
                       i,isu
5                 format(i4,1x,i2,1x,f6.3,3(1x,f8.3),3(1x,f8.3),1x,i5,1x,i4)
               end if

               if(realParticles(i,j)%ID .ne. pion) cycle Particle_loop2

               pPart => realParticles(i,j)

               momentumAbs=sqrt(dot_product(pPart%momentum(1:3),pPart%momentum(1:3)))
               ibin=nint((momentumAbs-dmom/2.)/dmom)+1
               if(ibin.ge.1 .and. ibin.le.Nmom .and. pPart%charge.ne.0) then
                  dNpiondMom(ibin,0)=dNpiondMom(ibin,0)+1.
                  if(numPions.ge.1 .and. numPions.le.10) &
                       & dNpiondMom(ibin,numPions)=dNpiondMom(ibin,numPions)+1.
               end if

            End do Particle_loop2

         end if


      End do Ensemble_loop

      if(finalFlag) then

         ! Normas:
         !fnorm=1./float(numEnsembles)/float(isu)
         if(pion_events.gt.0.) then
            fnorm=1./pion_events
         else
            fnorm=0.
         end if

         dNpiondMom(:,:)=dNpiondMom(:,:)*fnorm/dmom
         P_Npion(:)=P_Npion(:)*fnorm


         open(36,file='DoHIA1.dat')
         write(36,*)'# Distribution of events in pion number'
         write(36,*)'# Number of events: ', numEnsembles*isu
         write(36,*)'# Npion:    P_Npion:'
         pinumAv=0.
         do j=0,N_max
            pinumAv=pinumAv+float(j)*P_Npion(j)
            write(36,*) j, P_Npion(j)
         end do
         write(36,*)'# Norma: ', sum(P_Npion(:))
         write(36,*)'# Average pion number:', pinumAv
         close(36)

         open(36,file='DoHIA2.dat')
         write(36,*)'# Charged pion momentum distribution'
         write(36,*)'# Number of events: ', numEnsembles*isu
         write(36,*)'# momentum, GeV/c:    dNpiondMom, c/GeV:'
         do j=1,Nmom
            write(36,'(f6.3,8(1x,e10.3))')  (float(j)-0.5)*dmom, dNpiondMom(j,0), &
                 & dNpiondMom(j,1:6), sum(dNpiondMom(j,7:10))
         end do
         write(36,*)'# Norma: ', sum(dNpiondMom(:,0))*dmom
         close(36)

      end if

    end subroutine histo1


    subroutine analyze_Pions
      use histMC
      use IDTable, only: EOV
      use minkowski, only: abs4, abs3
      use inputGeneral, only: numEnsembles, num_Runs_SameEnergy, num_energies
      use initHeavyIon, only: ekin_lab_Projectile
      use constants, only: mN
      use nucleusDefinition
      use nucleus, only: getProjectile
      use histMC_avg

      type(histogramMC), save :: hist_pT, hist_y, hist_y0, hist_cost
      type(histogramMC), dimension(0:7), save :: hist_mT       ! 0=total; 1:7=different rap. bins
      type(histogramMC_avg), save :: hist_v1_y, hist_v2_y, hist_v1_u, hist_v2_u
      logical, save :: first = .true.
      real, save :: w, y_cms, u_proj
      integer :: i, j, ch, rb
      real :: mom(0:3), pT, mT, m, y, y0, cost, E, p, pabs, v1, v2, ut0, beta_proj, beta(1:3), gamma
      character(40) :: str
      type(tnucleus), pointer :: proj
      integer, save :: n = 0
      real :: fac

      if (first) then
         call CreateHistMC (hist_pT,   'Pion transverse momentum spectra: dN/dpT',     0.,1.,0.01,3)
         call CreateHistMC (hist_y,    'Pion rapidity spectra: dN/dy',                -3.,3.,0.1 ,3)
         call CreateHistMC (hist_y0,   'Pion rapidity spectra: dN/dy0',               -3.,3.,0.1 ,3)
         call CreateHistMC (hist_mT,   'Pion transverse mass spectra: mT^(-2)*dN/dmT', 0.,1.,0.01,3)
         call CreateHistMC (hist_cost, 'Pion polar angle spectra: dN/dcos(theta) with 0.2<p<0.8',    -1.,1.,0.05,3)
         call CreateHistMC_avg (hist_v1_y, 'Directed flow of pions: v1(y0)',              -2.,2.,0.2 ,3)
         call CreateHistMC_avg (hist_v2_y, 'Elliptic flow of pions: v2(y0)',              -2.,2.,0.2 ,3)
         call CreateHistMC_avg (hist_v1_u, 'Directed flow of pions: v1(ut0)',              0.,4.,0.2 ,3)
         call CreateHistMC_avg (hist_v2_u, 'Elliptic flow of pions: v2(ur0)',              0.,4.,0.2 ,3)
         do i=1,7
            write(str,'(A,i1,A,f5.2,A,f5.2,A)') " (rap. bin #",i,": y0 = ",rapBinning(i-1)," ... ",rapBinning(i)," GeV)"
            hist_mT(i)%name = trim(hist_mT(i)%name) // str
         end do
         hist_pT%xDesc = 'p_T [GeV]'
         hist_pT%yDesc = (/ "pi-", "pi0", "pi+" /)
         call CopyDesc (hist_y, hist_pT)
         call CopyDesc (hist_y0, hist_pT)
         call CopyDesc (hist_cost, hist_pT)
         call CopyDesc (hist_mT, hist_pT)
         hist_y%xDesc    = 'y'
         hist_y0%xDesc   = 'y0'
         hist_cost%xDesc = 'cos(theta_cm)'
         hist_mT%xDesc   = 'm_T - m [GeV]'
         hist_v1_y%xDesc   = 'y0'
         hist_v1_y%yDesc = (/ "pi-", "pi0", "pi+" /)
         call CopyDesc_avg (hist_v2_y, hist_v1_y)
         call CopyDesc_avg (hist_v1_u, hist_v1_y)
         call CopyDesc_avg (hist_v2_u, hist_v1_y)
         hist_v1_u%xDesc   = 'ut0'
         hist_v2_u%xDesc   = 'ut0'
         w = 1. / (numEnsembles*num_Runs_SameEnergy*num_energies)   ! weight factor
         E = 2*mN + ekin_lab_Projectile
         p = sqrt(ekin_lab_Projectile**2 + 2.*mN*ekin_lab_Projectile)  ! assumption: fixed target
         y_cms = 0.5*log((E+p)/(E-p))
         print *, "analyze_Pions: y_cms = ", y_cms
         proj => getProjectile()
         beta_proj = sqrt(sum(proj%velocity**2))
         u_proj = beta_proj / sqrt(1. - beta_proj**2)
         print *, "analyze_Pions: u_proj = ", u_proj
         first = .false.
      end if

      do i=LBOUND(realParticles,1),UBOUND(realParticles,1)    ! loop over ensembles
         do j=LBOUND(realParticles,2),UBOUND(realParticles,2) ! loop over all particles
            if (realParticles(i,j)%ID == EOV) exit            ! end of current ensemble reached, jump to next one
            if (realParticles(i,j)%ID /= pion) cycle
            ! observables
            ch = realParticles(i,j)%charge+2
            mom = realParticles(i,j)%momentum
            beta = mom(1:3) / mom(0)
            gamma = 1./ sqrt( 1. - dot_product(beta,beta) )
            m = abs4(mom)                                 ! mass
            pabs = abs3(mom)                              ! absolute momentum
            y = 0.5*log((mom(0)+mom(3))/(mom(0)-mom(3)))  ! rapidity
            y0 = y/y_cms                                  ! 'normalized' rapidity
            pt = sqrt(mom(1)**2+mom(2)**2)                ! transverse momentum
            mt = sqrt(m**2+pt**2)                         ! transverse mass
            cost = mom(3)/pabs
            v1 = mom(1)/pt                                      ! v1 = < px/pT >
            v2 = (mom(1)**2 - mom(2)**2)/pt**2                  ! v2 = < (px**2 - py**2) / pT**2 >
            ut0 = gamma * sqrt(beta(1)**2+beta(2)**2) / u_proj  ! transverse comp. of scaled four-velocity
            ! histograms
            call AddHistMC (hist_pT,    pt,   ch, w)
            call AddHistMC (hist_y,     y,    ch, w)
            call AddHistMC (hist_y0,    y0,   ch, w)
            call AddHistMC (hist_mT(0), mt-m, ch, w/mt**2)
            if (pabs>0.2 .and. pabs<0.8) call AddHistMC (hist_cost,  cost, ch, w)
            if (ut0>1.0 .and. ut0<4.2) then
              call AddHistMC_avg (hist_v1_y, y0, ch, v1)
              call AddHistMC_avg (hist_v2_y, y0, ch, v2)
            end if
            if (y0>-1.8 .and. y0<0.) then
              call AddHistMC_avg (hist_v1_u, ut0, ch, v1)
              call AddHistMC_avg (hist_v2_u, ut0, ch, v2)
            end if
            ! determine rap. bin
            rb = 0
            do while (rb < 8)
               if (y0 < rapBinning(rb)) exit
               rb = rb + 1
            end do
            if (rb>0 .and. rb<8) call AddHistMC (hist_mT(rb), mt-m, ch, w/mt**2)
         end do
      end do

      n = n + 1                                        ! count runs
      fac = num_Runs_SameEnergy*num_energies/float(n)  ! multiplication factor for writing histograms

      call WriteHistMC (hist_pT,  'PionPt.dat',   mul=fac)
      call WriteHistMC (hist_y,   'PionY.dat',    mul=fac)
      call WriteHistMC (hist_y0,  'PionY0.dat',   mul=fac)
      call WriteHistMC (hist_cost,'PionCost.dat', mul=fac)
      do i=0,7
         str = ""
         if (i>0) write(str,'(A,i1)') "_rapBin",i
         call WriteHistMC (hist_mT(i), 'PionMt'//trim(str)//'.dat', mul=fac)
      end do
      call WriteHistMC_avg (hist_v1_y,'PionV1_y.dat')
      call WriteHistMC_avg (hist_v2_y,'PionV2_y.dat')
      call WriteHistMC_avg (hist_v1_u,'PionV1_u.dat')
      call WriteHistMC_avg (hist_v2_u,'PionV2_u.dat')
    end subroutine analyze_Pions


    subroutine analyze_Kaons
      use histMC
      use initHeavyIon, only: ekin_lab_Projectile
      use idTable, only: Kaon, Kaonbar, isBaryon
      use constants, only: mN
      use minkowski, only: abs4, abs3
      use inputGeneral, only: numEnsembles, num_Runs_SameEnergy, num_energies

      type(histogramMC), save :: hist_p, hist_pT, hist_y, hist_y0, hist_cost, hist_mt, hist_parents
      logical, save :: first = .true.
      real, save :: w, y_cms
      integer :: i, j, ch, parents(1:3), channel
      real :: mom(0:3), pT, mT, m, y, y0, cost, E, p, pabs

      if (first) then
        call CreateHistMC (hist_p,    'Kaon momentum spectra: dN/dp',                 0.,1.5,0.05,4)
        call CreateHistMC (hist_pT,   'Kaon transverse momentum spectra: dN/dpT',     0.,1.,0.01,4)
        call CreateHistMC (hist_y,    'Kaon rapidity spectra: dN/dy',                -3.,3.,0.1 ,4)
        call CreateHistMC (hist_y0,   'Kaon rapidity spectra: dN/dy0',               -3.,3.,0.1 ,4)
        call CreateHistMC (hist_mT,   'Kaon transverse mass spectra: mT^(-2)*dN/dmT', 0.,1.,0.01,4)
        call CreateHistMC (hist_cost, 'Kaon polar angle spectra: dN/dcos(theta)',    -1.,1.,0.05,4)
        call CreateHistMC (hist_parents, 'Kaon parent sources',                       0.5,5.5,1.,4)
        hist_p%xDesc = 'p [GeV]'
        hist_p%yDesc = (/ "K0   ", "K+   ", "K-   ", "K0bar" /)
        call CopyDesc (hist_pT,      hist_p)
        call CopyDesc (hist_y,       hist_p)
        call CopyDesc (hist_y0,      hist_p)
        call CopyDesc (hist_cost,    hist_p)
        call CopyDesc (hist_mT,      hist_p)
        call CopyDesc (hist_parents, hist_p)
        hist_pT%xDesc   = 'p_T [GeV]'
        hist_y%xDesc    = 'y'
        hist_y0%xDesc   = 'y0'
        hist_cost%xDesc = 'cos(theta_cm)'
        hist_mT%xDesc   = 'm_T - m [GeV]'
        hist_parents%xDesc = "source channel"
        w = 1. / (numEnsembles*num_Runs_SameEnergy*num_energies)   ! weight factor
        E = 2*mN + ekin_lab_Projectile
        p = sqrt(ekin_lab_Projectile**2 + 2.*mN*ekin_lab_Projectile)  ! assumption: fixed target
        y_cms = 0.5*log((E+p)/(E-p))
        print *,"analyze_Kaons: y_cms = ",y_cms
        first = .false.
      end if

      do i=LBOUND(realParticles,1),UBOUND(realParticles,1)    ! loop over ensembles
         do j=LBOUND(realParticles,2),UBOUND(realParticles,2) ! loop over all particles
            if (realParticles(i,j)%ID == EOV) exit            ! end of current ensemble reached, jump to next one
            if (realParticles(i,j)%ID == Kaon) then
              ch = realParticles(i,j)%charge+1
            else if (realParticles(i,j)%ID == Kaonbar) then
              ch = realParticles(i,j)%charge+4
            else
              cycle
            end if
            mom = realParticles(i,j)%momentum
            m = abs4(mom)                                 ! mass
            pabs = abs3(mom)                              ! absolute momentum
            y = 0.5*log((mom(0)+mom(3))/(mom(0)-mom(3)))  ! rapidity
            y0 = y/y_cms                                  ! 'normalized' rapidity
            pt = sqrt(mom(1)**2+mom(2)**2)                ! transverse momentum
            mt = sqrt(m**2+pt**2)                         ! transverse mass
            cost = mom(3)/pabs
            call AddHistMC (hist_p,    pabs, ch, w)
            call AddHistMC (hist_pT,     pt, ch, w)
            call AddHistMC (hist_y,       y, ch, w)
            call AddHistMC (hist_y0,     y0, ch, w)
            call AddHistMC (hist_mT,   mt-m, ch, w/mt**2)
            call AddHistMC (hist_cost, cost, ch, w)
            ! determine where it came from
            parents = history_getParents (realParticles(i,j)%history)
            if (parents(2) > 0) then
              ! 2-body collisions
              if (isBaryon(parents(1)) .and. isBaryon(parents(2))) then
                channel = 1  ! channel 1: BB
              else if (isMeson(parents(1)) .and. isMeson(parents(2))) then
                channel = 2  ! channel 2: mm
              else
                channel = 3  ! channel 3: mB
              end if
!               print *,"analyze_Kaons (2-body):", ch, parents(1:3)
            else if (isMeson(parents(1))) then
              channel = 4   ! channel 4: meson decay
!               print *,"analyze_Kaons (meson dec.):", ch, parents(1:3)
            else if (isBaryon(parents(1))) then
              channel = 5   ! channel 5: baryon decay
!               print *,"analyze_Kaons (baryon dec.):", ch, parents(1:3)
            else
              print *,"analyze_Kaons: unknown parents!", ch, parents(1:3), realParticles(i,j)%history
              stop
            end if
            call AddHistMC (hist_parents, float(channel), ch, w)
         end do
      end do
      call WriteHistMC (hist_p,       'KaonPlab.dat')
      call WriteHistMC (hist_pT,      'KaonPt.dat')
      call WriteHistMC (hist_y,       'KaonY.dat')
      call WriteHistMC (hist_y0,      'KaonY0.dat')
      call WriteHistMC (hist_mT,      'KaonMt.dat')
      call WriteHistMC (hist_cost,    'KaonCost.dat')
      call WriteHistMC (hist_parents, 'KaonParents.dat')
    end subroutine

  end subroutine DoHeavyIonAnalysis


  !*************************************************************************
  subroutine DoHeavyIonAnalysisTime (realPV, time)
    use particleDefinition
    use inputGeneral, only : time_max, delta_T, timeForOutput, timeSequence
    use output, only : intTochar
    use IdTable, only: EOV, NOP

    type(particle), dimension(:,:), intent(in) :: realPV    ! real particle vector
    real, intent(in) :: time

    integer, save    :: isut=0      ! number of subsequent runs
    integer, save    :: itime=0

    if (time==0.) itime=0

    if (time>= timeForOutput) then

      if (mod(itime*delta_T,timeSequence)<1E-4) then
        call writeParticleVector
        call RapiditySpectra
      endif

      itime = itime + 1

    end if

    if (abs(time-time_max) < 1E-4) isut = isut + 1

  contains

    subroutine writeParticleVector
      use minkowski, only: abs4
      use inputGeneral, only: eventtype
      use eventtypes, only: HeavyIon, hadron
      use initHeavyIon, only: b_HI => b
      use initHadron,   only: b_had => b

      integer :: i, j, fact, indFree
      real    :: stossParameter

      !*************************************************************************
      !****o* HeavyIonAnalysis/DoHIATime___.dat
      ! NAME
      ! file DoHIATime___.dat
      ! PURPOSE
      ! Contains the full dump of the real particle vector at a certain time step.
      ! The time is given in the file name and in the first line of the file.
      ! The columns have the following meaning:
      ! * 1    = particle ID
      ! * 2    = charge
      ! * 3    = vac. mass in GeV
      ! * 4    = eff. mass in GeV
      ! * 5    = isFree (particle is bound or free?)
      ! * 6-8  = position (x,y,z) in fm
      ! * 9-11 = momentum (px,py,pz) in GeV
      ! * 12   = ensemble no.
      ! * 13   = run no.
      ! * 14   = impact parameter in fm
      !*************************************************************************
      open(103,file='DoHIATime'//intToChar(nint(time/delta_T))//'.dat',position='Append')
      if (isut==0) then
         write(103,*) '# time = ',time,' fm/c'
         write(103,*)
      endif

      if (eventtype==HeavyIon) then
         stossParameter = b_HI
      else if(eventtype==hadron) then
         stossParameter = b_had
      else
         write(*,*)' Problem in DoHeavyIonAnalysisTime: impact parameter is not defined'
         stossParameter=0.
      end if

      Loop_over_ensembles: do i = 1,size(realPV,dim=1)
         Loop_over_particles : do j = 1,size(realPV,dim=2)
            If (realPV(i,j)%id == NOP) cycle Loop_over_particles
            If (realPV(i,j)%id == EOV) exit Loop_over_particles
            if (realPV(i,j)%ID > 0) then
               if (realPV(i,j)%antiparticle) then
                  fact = -1
               else
                  fact = 1
               endif
               if (isFree(realPV,i,j,realPV(i,j))) then
                  indFree=1
               else
                  indFree=0
               end if
               write(103,50) fact*realPV(i,j)%ID, realPV(i,j)%charge, &
                             realPV(i,j)%mass, abs4(realPV(i,j)%momentum)**2, indFree, &
                             realPV(i,j)%position(1:3), realPV(i,j)%momentum(1:3), &
                             i, isut+1, stossParameter
            endif
         enddo Loop_over_particles
      enddo Loop_over_ensembles
      close(103)

50    format(2i4,2f9.4,i2,6f9.3,2i6,f6.2)

    end subroutine writeParticleVector



    subroutine RapiditySpectra
      use idTable, only: nucleon
      use inputGeneral, only: numEnsembles
      use histf90

      integer, parameter :: nrap = 60
      real, parameter :: ystart=-6.0, dy=0.1

      type(histogram), allocatable, save :: dNdy(:)
      logical, save :: RapidityInitFLAG=.true.

      integer :: i, j, k, N
      real :: yb, p0, pz
      character(len=100) :: title

      if (RapidityInitFLAG) then
        N = int((time_max - timeForOutput)/timeSequence)
        ! print *,"RapiditySpectra: allocating histograms: ", N, time_max, timeForOutput, timeSequence
        allocate(dNdy(0:N))
        do k=0,N
          write(title,'(A,f7.2,A)') "nucleon rapidity distribution dN/dy at time = ", timeForOutput + k*timeSequence, " fm/c"
          call createHist (dNdy(k),  title, ystart, ystart+2*nrap*dy, dy)
        end do
        RapidityInitFLAG = .false.
      endif

      k = nint((time - timeForOutput) / timeSequence)

      Loop_over_ensembles: do i = 1,size(realPV,dim=1)
         Loop_over_particles : do j = 1,size(realPV,dim=2)
            If (realPV(i,j)%id == NOP) cycle Loop_over_particles
            If (realPV(i,j)%id == EOV) exit Loop_over_particles

            if (realPV(i,j)%ID == nucleon) then

               p0 = realPV(i,j)%momentum(0)
               pz = realPV(i,j)%momentum(3)
               yb = 0.5*log( (p0+pz)/(p0-pz) )

               call addHist (dNdy(k), yb, 1./float(numEnsembles))
            endif
         enddo Loop_over_particles
      enddo Loop_over_ensembles

      call writeHist (dNdy(k),  file='RapidityDistributions'//intTochar(nint(time/delta_T))//'.dat')

    end subroutine RapiditySpectra


  end subroutine DoHeavyIonAnalysisTime



  !*** time evolution of some global observables
  subroutine HeavyIon_evol (realparticles, time)

    use IdTable
    use particleDefinition
    use densitymodule
    use RMF, only: getRMF_flag, fourMomDen_flag, g_omega, g_rho, g_sigma, ModificationFactor
    use particleProperties, only: hadron
    use output, only : realTochar,intTochar
    use inputGeneral, only: delta_T, eventtype
    use constants, only : pi, hbarc, rhoNull, mN
    use coulomb, only : emfoca
    use PauliBlockingModule, only: pauliBlocking
    use initHadron, only : b,z,p_lab,E_bind
    use collisionNumbering, only : GetCountedEvents
    use eventtypes, only: ET_hadron => hadron
    use thermoDynamics, only: temperatureAt, muAt

    type(particle), dimension(:,:), intent(in)  :: realparticles   ! real particle vector

    real, intent(in) :: time

    real, dimension(-121:121,-2:2) :: parnum, parnum_free     ! particle numbers
    integer :: numensembles,i,j,k,id,charge,Index1,Index2,Index3,npart
    integer, save :: icall=0     ! counter of calls
    integer, save :: isut=0      ! number of subsequent run
    real :: rhobar_points,rhobar_gauss,rhorad,rhorad_n,rhorad_p,r1,r2,r,velrad1,velrad2,&
         &rhoz,rhoz_bar,rhoz_antibar,rho_bar_max,rho_antibar_max,rholrf,&
         &endens,rhobar_local,m_inv
    real :: rhoz_Bar_points,rhoz_AntiBar_points,rhoz_BoundBar_points,rhoz_TargetBar_points,&
         &energy,fnorm
    real :: mstar, pf, sigma_rad, factor !,Ef
    real, dimension(1:3) :: p2_aver
    ! Check conservation laws: *************************************
    real, dimension(0:3) :: p ! total 4-momentum of all particles
    real :: baryon_number,charge_number,strangeness
    ! *************************************************************
    ! Properties of the bound system: ****************************************
    real, parameter :: rho_min_bound=0.01*rhoNull ! minimal density
    real :: B_bound, N_bound, Z_bound ! baryon, neutron and charge numbers
    real :: rho_bound                ! mass averaged density
    real :: E_kinColl_bound          ! collective kinetic energy
    real, dimension(0:3) :: P_bound  ! total 4-momentum (w/o Coulomb contribution in P_bound(0))
    real :: E_Coul_bound             ! Coulomb energy
    !*************************************************************************
    real :: cpot
    real, dimension(1:3) :: place,impuls,velColl
    real, dimension(0:3) :: momentum

    logical :: flag_output, blockFlag

    real :: rhoz_n, rhoz_p, Ef_p, Ef_n, pf_p, pf_n, temp, mub
    real :: Ecoul   ! EcoulNorm
    real :: pion_ratio

    real, dimension(1:3)     :: number_of_baryons, Qzz_aver, r2_aver, r_aver_sq
    real, dimension(1:3)     :: number_of_neutrons,number_of_protons,np_ratio
    real, dimension(1:3,1:3) :: r_aver
    real, dimension(1:6) :: probability

    if (initFlag) call init

    if (icall==0) then
       open(32,file='evol.dat')
       write(32,*)'# Column No., quantity:'
       write(32,*)'# 1  time'
       write(32,*)'# 2  nucleon multiplicity'
       write(32,*)'# 3  delta'
       write(32,*)'# 4  Higher Nonstrange Baryonic Resonances'
       write(32,*)'# 5  pi'
       write(32,*)'# 6  K'
       write(32,*)'# 7  Kbar'
       write(32,*)'# 8  Lambda'
       write(32,*)'# 9  Lambda free'
       write(32,*)'# 10 Sigma'
       write(32,*)'# 11 Sigma free'
       write(32,*)'# 12 Xi'
       write(32,*)'# 13 Xi free'
       write(32,*)'# 14 XiStar'
       write(32,*)'# 15 XiStar free'
       write(32,*)'# 16 K^*'
       write(32,*)'# 17 Kbar^*'
       write(32,*)'# 18 Y^* (S=-1 only)'
       write(32,*)'# 19 LambdaBar'
       write(32,*)'# 20 SigmaBar'
       write(32,*)'# 21 Ybar^* (S=1 only)'
       write(32,*)'# 22 XiBar'
       write(32,*)'# 23 XiBar^*'
       write(32,*)'# 24 J/Psi'
       write(32,*)'# 25 D'
       write(32,*)'# 26 Dbar'
       write(32,*)'# 27 D^*'
       write(32,*)'# 28 Dbar^*'

       open(33,file='dens.dat',status='unknown')
       write(33,'(A)') '# time rhocen_baryons rhocen_antibaryons ptot baryon_number '// &
                       'charge strangeness temperature(central) mu_B(central)'

       open(40,file='conserv.dat')
       write(40,*)'# time: energy:  px:  py:  pz:  baryon number:  charge: strangeness:'

       if (getRMF_flag() .and. fourMomDen_flag ) then
          open(43,file='BoundSystem.dat')
          write(43,*)'# Properties of the bound system:'
          write(43,*)'# time:   B:   N:   Z:  P(0:3), AGeV:   E_kinColl, AGeV:  E_Coul, AGeV:  <rho>, fm ^-3:'
       end if

       open(44,file='collisions.dat')
       write(44,*)'# time N_2body(integrated) N_2body(timestep)'

       open(50,file='isospinRatios.dat')

    endif

    icall = icall + 1
    if( abs(time-delta_T) < 1.e-06 ) isut=isut+1
    flag_output=.false.

    parnum= 0.
    parnum_free= 0.
    rhobar_points= 0.
    endens= 0.
    p2_aver= 0.
    r2_aver=0.
    r_aver=0.0
    Qzz_aver=0.
    p(:)=0.
    baryon_number=0.
    charge_number=0.
    strangeness=0.
    number_of_baryons=0.
    number_of_neutrons=0.
    number_of_protons=0.
    E_Coul_bound=0.

    numensembles = size(realParticles,dim=1)
    ensemble_loop : do i = 1,numensembles
       particle_loop : do j = 1,size(realParticles,dim=2)

          if (realParticles(i,j)%ID == EOV) exit particle_loop
          if (realParticles(i,j)%ID == NOP) cycle particle_loop

          !     Particle numbers:
          if(realParticles(i,j)%antiparticle) then
             factor = -1.
          else
             factor = 1.
          endif
          id = realParticles(i,j)%ID*int(factor)
          charge = realParticles(i,j)%charge
          if(abs(id).le.121) then
             parnum(id,charge) = parnum(id,charge) + 1.
             if(id.eq.Lambda .or. id.eq.SigmaResonance .or. id.eq.Xi .or. id.eq.XiStar) then
                if(IsFree(realParticles,i,j,realParticles(i,j))) parnum_free(id,charge) = parnum_free(id,charge) + 1.
             end if
          end if
          p(:)=p(:)+realParticles(i,j)%momentum(:)
          if(isBaryon(realParticles(i,j)%ID)) then
             baryon_number=baryon_number+factor
             strangeness=strangeness+real(hadron(realParticles(i,j)%ID)%strangeness)*factor
          else if(isMeson(realParticles(i,j)%ID)) then
             strangeness=strangeness+real(hadron(realParticles(i,j)%ID)%strangeness)
          end if
          charge_number=charge_number+real(realParticles(i,j)%charge)

          if(abs(realParticles(i,j)%position(1)) < 0.5 .and. &
               abs(realParticles(i,j)%position(2)) < 0.5 .and. &
               abs(realParticles(i,j)%position(3)) < 0.5)  then
             !       Central baryon density:
             if(isBaryon(realParticles(i,j)%ID)) then
                rhobar_points=  rhobar_points + factor
             end if
             !       Central energy density:
             endens= endens + realParticles(i,j)%momentum(0)
          endif

          if(isBaryon(realParticles(i,j)%ID)) then

             !       Average squares of momentum components :
             p2_aver(1:3)= p2_aver(1:3) + realParticles(i,j)%momentum(1:3)**2

             if(.not.realParticles(i,j)%antiparticle) then

                !         Indexes in large grid:
                Index1=NINT(realParticles(i,j)%position(1)/gridSpacing(1))
                Index2=NINT(realParticles(i,j)%position(2)/gridSpacing(2))
                Index3=NINT(realParticles(i,j)%position(3)/gridSpacing(3))

                if(        abs(Index1).le.gridPoints(1) &
                     & .and. abs(Index2).le.gridPoints(2) &
                     & .and. abs(Index3).le.gridPoints(3)  ) then

                   if(getRMF_flag()) then
                      factor=ModificationFactor(nucleon,.true.)    ! Modification factor of the antinucleon coupling constants
                      rhobar_local=(   densityField(Index1,Index2,Index3)%baryon(0)  &
                           & + factor*totalDensity(Index1,Index2,Index3) )/(1.+factor)  ! density of the baryons
                   else
                      ! w/o RMF "totalDensity" is actually (baryon-antibaryon) density and
                      ! densityField(0,0,0)%baryon(0) is (baryon+antibaryon) density (see density.f90: addTodensityfield)
                      rhobar_local=(   densityField(Index1,Index2,Index3)%baryon(0)  &
                           & + totalDensity(Index1,Index2,Index3) )/2.                 ! density of the baryons
                   end if

                   !cut on density (0.1*rho_sat)
                   if(rhobar_local.gt.0.1*rhoNull) then
                      number_of_baryons(1)=number_of_baryons(1)+1.
                      r2_aver(1)= r2_aver(1) &
                           & + dot_product(realParticles(i,j)%position,realParticles(i,j)%position)
                      r_aver(1,:)= r_aver(1,:) + realParticles(i,j)%position(:)
                      Qzz_aver(1)=Qzz_aver(1) + 2.*realParticles(i,j)%position(3)**2 &
                           &- realParticles(i,j)%position(1)**2 - realParticles(i,j)%position(2)**2
                      if (realParticles(i,j)%ID==1 .and. realParticles(i,j)%charge==0) &
                           & number_of_neutrons(1)=number_of_neutrons(1)+1.
                      if (realParticles(i,j)%ID==1 .and. realParticles(i,j)%charge==1) &
                           & number_of_protons(1)=number_of_protons(1)+1.
                   end if

                   !cut on density (0.01*rho_sat)
                   if(rhobar_local.gt.0.01*rhoNull) then
                      number_of_baryons(2)=number_of_baryons(2)+1.
                      r2_aver(2) = r2_aver(2) &
                           & + dot_product(realParticles(i,j)%position,realParticles(i,j)%position)
                      r_aver(2,:)= r_aver(2,:) + realParticles(i,j)%position(:)
                      Qzz_aver(2)=Qzz_aver(2) + 2.*realParticles(i,j)%position(3)**2 &
                           &- realParticles(i,j)%position(1)**2 - realParticles(i,j)%position(2)**2
                      if (realParticles(i,j)%ID==1 .and. realParticles(i,j)%charge==0) &
                           & number_of_neutrons(2)=number_of_neutrons(2)+1.
                      if (realParticles(i,j)%ID==1 .and. realParticles(i,j)%charge==1) &
                           & number_of_protons(2)=number_of_protons(2)+1.
                   end if

                   !cut on binding energy (E<0)
                   if(getRMF_flag()) then
                      call Particle4Momentum_RMF(realParticles(i,j),momentum)
                      energy=momentum(0)
                   else
                      energy=realParticles(i,j)%momentum(0)
                   end if
                   place(1:3)= realParticles(i,j)%position(1:3)
                   impuls(1:3)= realParticles(i,j)%momentum(1:3)
                   energy = energy + emfoca(place,impuls,realParticles(i,j)%charge,realParticles(i,j)%ID)
                   if ( energy - realParticles(i,j)%mass < 0.) then
                      number_of_baryons(3)=number_of_baryons(3)+1.
                      r2_aver(3) = r2_aver(3) &
                           & + dot_product(realParticles(i,j)%position,realParticles(i,j)%position)
                      r_aver(3,:)= r_aver(3,:) + realParticles(i,j)%position(:)
                      if (realParticles(i,j)%ID==1 .and. realParticles(i,j)%charge==0) &
                           & number_of_neutrons(3)=number_of_neutrons(3)+1.
                      if (realParticles(i,j)%ID==1 .and. realParticles(i,j)%charge==1) &
                           & number_of_protons(3)=number_of_protons(3)+1.
                   endif

                   if(rhobar_local.gt.rho_min_bound) then
                      place(1:3)= realParticles(i,j)%position(1:3)
                      impuls(1:3)= realParticles(i,j)%momentum(1:3)
                      E_Coul_bound = E_Coul_bound + 0.5*emfoca(place,impuls,realParticles(i,j)%charge,realParticles(i,j)%ID) !see notes in evaluateTotal4Momentum_RMF
                   end if

                end if

             end if

          end if

       end do particle_loop
    enddo ensemble_loop

    parnum= parnum/float(numensembles)
    parnum_free= parnum_free/float(numensembles)
    rhobar_points= rhobar_points/float(numensembles)/1.**3
    endens= endens/float(numensembles)/1.**3
    p2_aver(:)= p2_aver(:)/float(numensembles)
    r2_aver(:)= r2_aver(:)/number_of_baryons(:)
    do i=1,3
       r_aver(i,:)  = r_aver(i,:)/number_of_baryons(i)
       r_aver_sq(i) = sqrt(dot_product(r_aver(i,:),r_aver(i,:)))
    end do
    Qzz_aver(:)=Qzz_aver(:)/number_of_baryons(:)
    p(:)=p(:)/float(numensembles)
    baryon_number=baryon_number/float(numensembles)
    charge_number=charge_number/float(numensembles)
    strangeness=strangeness/float(numensembles)
    number_of_baryons(:)=number_of_baryons(:)/float(numensembles)
    number_of_protons(:)=number_of_protons(:)/float(numensembles)
    number_of_neutrons(:)=number_of_neutrons(:)/float(numensembles)
    E_Coul_bound=E_Coul_bound/float(numensembles)

    !n/p-ratio:
    do i=1,3
       if ( number_of_protons(i) .ne. 0.0 ) then
          np_ratio(i) = number_of_neutrons(i) / number_of_protons(i)
       else
          np_ratio(i) = 0.0
       end if
    end do

    !pion-ratio pi^{-}/pi^{+}:
    if ( parnum(pion,+1).ne.0.0 ) then
       pion_ratio = parnum(pion,-1)/parnum(pion,+1)
    else
       pion_ratio = 0.0
    end if

    write(32,5) time,sum(parnum(nucleon,:)),sum(parnum(delta,:)),&
         sum(parnum(P11_1440:F37_1950,:)),&
         sum(parnum(pion,:)),sum(parnum(kaon,:)),sum(parnum(kaonBar,:)),&
         sum(parnum(Lambda,:)),sum(parnum_free(Lambda,:)),&
         sum(parnum(SigmaResonance,:)),sum(parnum_free(SigmaResonance,:)),&
         sum(parnum(Xi,:)),sum(parnum_free(Xi,:)),&
         sum(parnum(XiStar,:)),sum(parnum_free(XiStar,:)),&
         sum(parnum(kaonStar,:)),sum(parnum(kaonStarBar,:)),&
         sum(parnum(Sigma_1385:Sigma_1915,:)),&
         sum(parnum(-Lambda,:)),sum(parnum(-SigmaResonance,:)),&
         sum(parnum(-Sigma_1915:-Sigma_1385,:)),&
         sum(parnum(-Xi,:)),sum(parnum(-XiStar,:)),&
         sum(parnum(JPsi,:)),sum(parnum(dMeson,:)),sum(parnum(dBar,:)),&
         sum(parnum(dStar,:)),sum(parnum(dStarBar,:))

    write(50,5) time,(np_ratio(i),i=1,3),pion_ratio

5   format(50(1x,e13.6))

    rhobar_gauss= 0.
    do k = -2,2
       do j = -2,2
          do i = -2,2
             rhobar_gauss= rhobar_gauss + densityField(i,j,k)%baryon(0)
          enddo
       enddo
    enddo
    rhobar_gauss= rhobar_gauss/5.**3

    if(getRMF_flag()) then
       factor=ModificationFactor(nucleon,.true.)    ! Modification factor of the antinucleon coupling constants
       rhoz_bar=(densityField(0,0,0)%baryon(0)+factor*totalDensity(0,0,0))/(1.+factor)  ! density of the baryons
       rhoz_antibar=(totalDensity(0,0,0)-densityField(0,0,0)%baryon(0))/(1.+factor)  ! density of the antibaryons
    else
       ! w/o RMF "totalDensity" is actually (baryon-antibaryon) density and
       ! densityField(0,0,0)%baryon(0) is (baryon+antibaryon) density (see density.f90: addTodensityfield)
       rhoz_bar=(densityField(0,0,0)%baryon(0)+totalDensity(0,0,0))/2.
       rhoz_antibar=(densityField(0,0,0)%baryon(0)-totalDensity(0,0,0))/2.
    end if

    temp = temperatureAt ((/0.,0.,0./))
    mub = muAt (rhoz_bar, temp)
    write(33,5) time, rhoz_bar, rhoz_antibar, p(:), baryon_number, charge_number, strangeness, temp, mub

    write(40,5) time,(r_aver_sq(i), sqrt(r2_aver(i)), number_of_baryons(i), i=1,3)

    if(getRMF_flag() .and. fourMomDen_flag ) then
       P_bound(:)=0.
       E_kinColl_bound=0.
       B_bound=0.
       N_bound=0.
       Z_bound=0.
       rho_bound=0.
       do Index1=-gridpoints(1),gridpoints(1)
          do Index2=-gridPoints(2),gridPoints(2)
             do Index3=-gridPoints(3),gridPoints(3)
                rhobar_local=densityField(Index1,Index2,Index3)%baryon(0)
                if(rhobar_local.gt.rho_min_bound) then
                   P_bound(0:3)=P_bound(0:3)+fourMomentumDensity(Index1,Index2,Index3,0:3)
                   m_inv=fourMomentumDensity(Index1,Index2,Index3,0)**2 &
                        &-dot_product(fourMomentumDensity(Index1,Index2,Index3,1:3),&
                        &fourMomentumDensity(Index1,Index2,Index3,1:3))
                   m_inv=sqrt(max(0.,m_inv))
                   E_kinColl_bound=E_kinColl_bound+fourMomentumDensity(Index1,Index2,Index3,0)&
                        &-m_inv
                   B_bound=B_bound+rhobar_local
                   N_bound=N_bound+densityField(Index1,Index2,Index3)%neutron(0)
                   Z_bound=Z_bound+densityField(Index1,Index2,Index3)%proton(0)
                   rho_bound=rho_bound+rhobar_local**2
                end if
             end do
          end do
       end do
       P_bound(0:3)=P_bound(0:3)*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       E_kinColl_bound=E_kinColl_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       B_bound=B_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       N_bound=N_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       Z_bound=Z_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       rho_bound=rho_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)/B_bound

       write(43,5) time, B_bound, N_bound, Z_bound, P_bound(:)/B_bound, E_kinColl_bound/B_bound,&
            &E_Coul_bound/B_bound, rho_bound
    end if

    write(44,5) time, float(getCountedEvents(0,2,1))/float(numensembles), float(getCountedEvents(1,2,1))/float(numensembles)

    ! flush all files, so that output is actually written to disk
    flush(32)
    flush(33)
    flush(40)
    flush(44)
    flush(50)

    !-----------------------------------------------------------------------------------------------
    ! more detailed output, if flag_outputDetailed=.true.
    !-----------------------------------------------------------------------------------------------
    if ( .not. flag_outputDetailed ) return
    !-----------------------------------------------------------------------------------------------

    if( abs(time-delta_T) < 1.e-06 .or. abs(time-nint(time)) < 0.5*delta_T ) then
       flag_output=.true.
       open(34,file='rhorad_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(34,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
       write(34,'(A,f10.4)')'# time, fm/c: ', time
       write(34,'(A,i3)')'# run number: ', isut
       open(35,file='rhoz_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(35,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
       write(35,'(A,f10.4)')'# time, fm/c: ', time
       write(35,'(A,i3)')'# run number: ', isut
       open(38,file='rhozx_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(38,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
       write(38,'(A,f10.4)')'# time, fm/c: ', time
       write(38,'(A,i3)')'# run number: ', isut
       if(getRMF_flag()) then
          open(36,file='Fields_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
          write(36,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
          write(36,'(A,f10.4)')'# time, fm/c: ', time
          write(36,'(A,i3)')'# run number: ', isut
       end if
       open(37,file='pauli_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(37,'(A,f10.4)')'# time, fm/c: ', time
       write(37,'(A,i3)')'# run number: ', isut
       write(37,*)'# Momentum, GeV/c:    PauliBlock:'
    end if

    if (icall == 1) then
       open(39,file='dens_max.dat',status='unknown')
       if (eventtype == ET_hadron) &
            write(39,'(A,4(e15.8,1x))')" # hadron's b, z, p_lab, E_bind: ", b, z, p_lab, E_bind
       write(39,*)'# time  rho_bar_max  rho_antibar_max'
    end if


    Final_Output : if(flag_output) then

       write(34,*) '# r, fm:   rhorad:  rhorad_n:   rhorad_p:  velrad1/c:   velrad2/c:  PauliBlock(1-6):'
       do k= 0,min(gridpoints(1),gridpoints(2),gridpoints(3))

          r1=max(0.,(float(k)-0.5)*gridSpacing(1))
          r2=(float(k)+0.5)*gridSpacing(1)

          npart=0
          rhorad=0.
          rhorad_n=0.
          rhorad_p=0.
          velrad1=0.
          velrad2=0.
          sigma_rad=0.

          do i = 1,numensembles
             particle_loop1 : do j = 1,size(realParticles,dim=2)

                if (realParticles(i,j)%ID == EOV) exit particle_loop1
                if (realParticles(i,j)%ID == NOP) cycle particle_loop1

                place(1:3)= realParticles(i,j)%position(1:3)
                r=sqrt(dot_product(place(1:3),place(1:3)))

                if( r.gt.r1 .and. r.le.r2 ) then

                   npart=npart+1

                   !             Indexes in large grid:
                   Index1=NINT(place(1)/gridSpacing(1))
                   Index2=NINT(place(2)/gridSpacing(2))
                   Index3=NINT(place(3)/gridSpacing(3))

                   if(        abs(Index1).le.gridPoints(1) &
                        & .and. abs(Index2).le.gridPoints(2) &
                        & .and. abs(Index3).le.gridPoints(3)  ) then

                      rhobar_local=densityField(Index1,Index2,Index3)%baryon(0)
                      rhorad=rhorad+rhobar_local

                      rhorad_n=rhorad_n+densityField(Index1,Index2,Index3)%neutron(0)
                      rhorad_p=rhorad_p+densityField(Index1,Index2,Index3)%proton(0)

                      if(allocated(sigmaField)) &
                           & sigma_rad=sigma_rad+sigmaField(Index1,Index2,Index3)

                      if(rhobar_local.gt.1.e-06) then
                         velColl(1:3)=densityField(Index1,Index2,Index3)%baryon(1:3) &
                              &/rhobar_local
                         if(r.gt.1.e-06) velrad1=velrad1+dot_product(velColl(1:3),place(1:3))/r
                      end if

                      if(allocated(fourMomentumDensity)) then
                         if(fourMomentumDensity(Index1,Index2,Index3,0).gt.1.e-06) then
                            velColl(1:3)=fourMomentumDensity(Index1,Index2,Index3,1:3) &
                                 &/fourMomentumDensity(Index1,Index2,Index3,0)
                            if(r.gt.1.e-06) velrad2=velrad2+dot_product(velColl(1:3),place(1:3))/r
                         end if
                      end if

                   end if

                end if


             end do particle_loop1
          end do

          if(npart.gt.0) then
             rhorad=rhorad/float(npart)
             rhorad_n=rhorad_n/float(npart)
             rhorad_p=rhorad_p/float(npart)
             sigma_rad=sigma_rad/float(npart)
             velrad1=velrad1/float(npart)
             velrad2=velrad2/float(npart)
          end if

          mstar = mN + g_sigma*sigma_rad

          !********** Pauli blocking check, radial dependence: ****************************
          do i=1,6
             momentum(1:3)=0.
             momentum(0)=sqrt(mN**2+dot_product(momentum(1:3),momentum(1:3)))
             place(1:3)=0.
             if(i.eq.1) place(1)=float(k)*gridSpacing(1)
             if(i.eq.2) place(1)=-float(k)*gridSpacing(1)
             if(i.eq.3) place(2)=float(k)*gridSpacing(2)
             if(i.eq.4) place(2)=-float(k)*gridSpacing(2)
             if(i.eq.5) place(3)=float(k)*gridSpacing(3)
             if(i.eq.6) place(3)=-float(k)*gridSpacing(3)
             blockFlag=pauliBlocking(momentum,place,1,realparticles,probability(i))
          end do
          !********************************************************************************

          write(34,5) float(k)*gridSpacing(1), rhorad, rhorad_n, rhorad_p, velrad1, velrad2,&
               &probability   !, mstar

       end do


       !********** Pauli blocking check, momentum dependence: *****************************
       place(1:3)=0.
       do k=1,201
          do i=1,6
             momentum(1:3)=0.
             if(i.eq.1) momentum(1)=0.00240*float(k-1)
             if(i.eq.2) momentum(1)=-0.00240*float(k-1)
             if(i.eq.3) momentum(2)=0.00240*float(k-1)
             if(i.eq.4) momentum(2)=-0.00240*float(k-1)
             if(i.eq.5) momentum(3)=0.00240*float(k-1)
             if(i.eq.6) momentum(3)=-0.00240*float(k-1)
             momentum(0)=sqrt(mN**2+dot_product(momentum(1:3),momentum(1:3)))
             blockFlag=pauliBlocking(momentum,place,1,realparticles,probability(i))
          end do
          write(37,5) 0.00240*float(k-1),probability
       end do
       !********************************************************************************

       if(getRMF_flag()) then
          write(35,'(2A)') &
               & '# z, fm:     rhoz_bar:   rhoz_antibar:',&
               &'rhoz_Bar_points:  rhoz_AntiBar_points:  rhoz_BoundBar_points: rhoz_TargetBar_points:'
       else
          write(35,'(2A)') &
               & '# z, fm:     rhoz_bar:   rhoz_antibar:   rholrf:',&
               &'rhoz_Bar_points:  rhoz_AntiBar_points:  rhoz_BoundBar_points: rhoz_TargetBar_points:'
       end if

       Loop_over_zGrid1 : do k= -gridpoints(3),gridpoints(3)

          rhoz= densityField(0,0,k)%baryon(0)

          !proton & neutron densities
          rhoz_n=densityField(0,0,k)%neutron(0)
          rhoz_p=densityField(0,0,k)%proton(0)

          if(getRMF_flag()) then

             factor=ModificationFactor(nucleon,.true.)    ! Modification factor of the antinucleon coupling constants

             rhoz_bar=(rhoz+factor*totalDensity(0,0,k))/(1.+factor)  ! density of the baryons

             rhoz_antibar=(totalDensity(0,0,k)-rhoz)/(1.+factor)     ! density of the antibaryons

             pf=(1.5*pi**2*max(0.,rhoz_bar))**0.333333*hbarc

             pf_p=(3.*pi**2*max(0.,rhoz_p))**0.333333*hbarc
             pf_n=(3.*pi**2*max(0.,rhoz_n))**0.333333*hbarc


             mstar=mN+g_sigma*sigmaField(0,0,k)

             !           Ef=sqrt(pf**2+mstar**2)+g_omega*omegaField(0,0,k,0)

             Ef_p=sqrt(pf_p**2+mstar**2)+g_omega*omegaField(0,0,k,0)
             Ef_n=sqrt(pf_n**2+mstar**2)+g_omega*omegaField(0,0,k,0)

          else

             ! Remember: w/o RMF totalDensity is (baryon-antibaryon) density
             rhoz_bar=(rhoz+totalDensity(0,0,k))/2.         ! density of the baryons
             rhoz_antibar=(rhoz-totalDensity(0,0,k))/2.     ! density of the antibaryons

             rholrf= sqrt(max( densityField(0,0,k)%baryon(0)**2 &
                  -densityField(0,0,k)%baryon(1)**2 &
                  -densityField(0,0,k)%baryon(2)**2 &
                  -densityField(0,0,k)%baryon(3)**2, 0. ))

          end if

          rhoz_Bar_points=0.
          rhoz_AntiBar_points=0.
          rhoz_BoundBar_points=0.
          rhoz_TargetBar_points=0.
          !        ECoul = 0.0
          !        EcoulNorm = 0.0

          do i = 1,numensembles
             particle_loop2 : do j = 1,size(realParticles,dim=2)

                if (realParticles(i,j)%ID == EOV) exit particle_loop2
                if (realParticles(i,j)%ID == NOP) cycle particle_loop2

                if(realParticles(i,j)%ID >= pion) cycle particle_loop2   ! Count only (anti)baryons

                place(1:3)= realParticles(i,j)%position(1:3)

                if( abs(place(1)) <= 0.5*gridSpacing(1) .and. &
                     &abs(place(2)) <= 0.5*gridSpacing(2) .and. &
                     &abs(place(3)-float(k)*gridSpacing(3)) <= 0.5*gridSpacing(3) ) then

                   if( .not.realParticles(i,j)%antiparticle ) then
                      rhoz_Bar_points=rhoz_bar_points+1.
                      if(getRMF_flag()) then
                         call Particle4Momentum_RMF(realParticles(i,j),momentum)
                         energy=momentum(0)
                      else
                         energy=realParticles(i,j)%momentum(0)
                      end if
                      impuls(1:3)= realParticles(i,j)%momentum(1:3)
                      cpot = emfoca(place,impuls,realParticles(i,j)%charge,realParticles(i,j)%ID)
                      energy=energy+cpot
                      !for testing E_F=E_F(r)
                      !                   if (realParticles(i,j)%charge==1) then
                      !                      Ecoul = Ecoul + cpot
                      !                      EcoulNorm = EcoulNorm + 1.
                      !                   end if
                      ! write(*,*)'Id, energy:', realParticles(i,j)%ID, energy
                      if(energy-realParticles(i,j)%mass < 0) rhoz_BoundBar_points=rhoz_BoundBar_points+1.
                      if(realParticles(i,j)%event(1)==1) &
                           & rhoz_TargetBar_points=rhoz_TargetBar_points+1.
                   else
                      rhoz_antibar_points=rhoz_antibar_points+1.
                   end if

                end if

             end do particle_loop2
          end do

          if(getRMF_flag()) then
             !Coulomb contribution to the Fermi-Energy:
             place=(/0.,0.,float(k)*gridSpacing(3)/)
             impuls=(/0.,0.,pf/)
             Ecoul = emfoca(place,impuls,1,1)
          end if

          fnorm=1./float(numensembles)/gridSpacing(1)/gridSpacing(2)/gridSpacing(3)
          rhoz_Bar_points=rhoz_Bar_points*fnorm
          rhoz_AntiBar_points=rhoz_AntiBar_points*fnorm
          rhoz_BoundBar_points=rhoz_BoundBar_points*fnorm
          rhoz_TargetBar_points=rhoz_TargetBar_points*fnorm

          if(getRMF_flag()) then

             write(35,5)  float(k)*gridSpacing(3), &
                  & rhoz_bar, rhoz_antibar, &
                  & rhoz_Bar_points,rhoz_AntiBar_points,&
                  & rhoz_BoundBar_points,rhoz_TargetBar_points

             write(36,5) float(k)*gridSpacing(3), &
                  & pf_p, pf_n, rhoz_p, rhoz_n, &
                  & -g_sigma*sigmaField(0,0,k), &
                  & g_omega*omegaField(0,0,k,0), &
                  & g_rho*rhoField(0,0,k,0),ECoul, &
                  & (Ef_p+Ecoul+g_rho*rhoField(0,0,k,0)-mN), &
                  & (Ef_n-g_rho*rhoField(0,0,k,0)-mN)
          else

             write(35,5)  float(k)*gridSpacing(3),rhoz_bar,rhoz_antibar,rholrf,rhoz_Bar_points,&
                  &rhoz_AntiBar_points,rhoz_BoundBar_points,rhoz_TargetBar_points

          end if

       enddo Loop_over_zGrid1

    end if Final_Output

    if(flag_output) write(38,'(A)') '# z, fm:  x, fm:     rhoz_bar:   rhoz_antibar:'
    if(getRMF_flag()) factor=ModificationFactor(nucleon,.true.)
    rho_bar_max=-0.1
    rho_antibar_max=-0.1

    Loop_over_yGrid2 : do j= -gridpoints(2),gridpoints(2)
       Loop_over_zGrid2 : do k= -gridpoints(3),gridpoints(3)
          Loop_over_xGrid2 : do i= -gridpoints(1),gridpoints(1)

             rhoz= densityField(i,j,k)%baryon(0)

             if(getRMF_flag()) then
                rhoz_bar=(rhoz+factor*totalDensity(i,j,k))/(1.+factor)  ! density of the baryons
                rhoz_antibar=(totalDensity(i,j,k)-rhoz)/(1.+factor)     ! density of the antibaryons
             else
                rhoz_bar=(rhoz+totalDensity(i,j,k))/2.         ! density of the baryons
                rhoz_antibar=(rhoz-totalDensity(i,j,k))/2.     ! density of the antibaryons
             end if

             if(rhoz_bar.gt.rho_bar_max) rho_bar_max=rhoz_bar
             if(rhoz_antibar.gt.rho_antibar_max) rho_antibar_max=rhoz_antibar

             if(flag_output .and. j.eq.0) then
                write(38,5)  float(k)*gridSpacing(3),float(i)*gridSpacing(1),&
                     & rhoz_bar, rhoz_antibar
                if(i.eq.gridpoints(1))  write(38,*)
             end if

          enddo Loop_over_xGrid2
       enddo Loop_over_zGrid2
    enddo Loop_over_yGrid2

    write(39,5) time, rho_bar_max, rho_antibar_max

  end subroutine HeavyIon_evol


  !*** Determine whether j-th particle of i-th parallel ensemble is free or not
  logical function IsFree(realParticles,i,j,teilchen)

    use particleDefinition
    use densitymodule, only: Particle4Momentum_RMF
    use RMF, only: getRMF_flag
    use coulomb, only: emfoca

    type(particle), dimension(:,:), intent(in)  :: realparticles
    integer, intent (in) :: i, j
    type(particle), intent(in)  ::  teilchen  ! particle of interest

    integer, parameter :: imode=1   ! 1 - criterion according to dstmin
    ! 2 - criterion according to binding energy
    real, parameter :: dstmin = 3.  ! minimum distance for free particle (fm)
    real :: tmp,dist2,energy
    integer :: j1
    real, dimension(0:3) :: momentum
    real, dimension(1:3) :: place

    place=teilchen%position

    if(imode.eq.1) then

       tmp=100.
       do j1=1,size(realParticles,dim=2)
          if(j1.ne.j .or. teilchen%perturbative) then
             dist2=(realParticles(i,j1)%position(1)-place(1))**2 &
                  &+(realParticles(i,j1)%position(2)-place(2))**2 &
                  &+(realParticles(i,j1)%position(3)-place(3))**2
             if(dist2.lt.tmp) tmp=dist2
          end if
       end do
       if(tmp.gt.dstmin**2) then
          IsFree=.true.
       else
          IsFree=.false.
       end if

    else if(imode.eq.2) then

       if(getRMF_flag()) then
          call Particle4Momentum_RMF(teilchen,momentum)
          energy=momentum(0)
       else
          energy=teilchen%momentum(0)
       end if
       energy=energy+emfoca(place,(/0.,0.,0./),realParticles(i,j)%charge,realParticles(i,j)%ID)
       if(energy.gt.teilchen%mass) then
          IsFree=.true.
       else
          IsFree=.false.
       end if

    end if

  end function IsFree


end module HeavyIonAnalysis
