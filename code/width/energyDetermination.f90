!***************************************************************************
!****m* /energyCalc
! NAME
! module energyCalc
! PURPOSE
! Incorporates routines to evaluate the one-particle energy of a given
! particle and a routine for making a so called energy correction to a given
! final state.
!***************************************************************************
module energyCalc
  implicit none

  Private

  !*************************************************************************
  !****g* energyCalc/accuracy
  ! PURPOSE
  ! Determines accuracy of energy determination in the determination of
  ! the one particle energy. Units : GeV.
  ! SOURCE
  !
  real, parameter :: accuracy=0.0005
  !*************************************************************************


  PUBLIC :: energyDetermination
  PUBLIC :: energyCorrection
  PUBLIC :: updateEnergies


contains


  !*************************************************************************
  !****s* energyCalc/updateEnergies
  ! NAME
  ! subroutine updateEnergies(teilchen,check_in)
  ! PURPOSE
  ! Routine updates energies of the particles
  ! INPUTS
  ! * type(particle), dimension(:) :: teilchen  -- particle array
  ! * logical, OPTIONAL :: check_in -- flag for 'energyDetermination'
  ! OUTPUT
  ! * type(particle), dimension(:) :: teilchen  -- particle array
  !*************************************************************************
  subroutine updateEnergies(teilchen,check_in)
    use particleDefinition
    use output, only: DoPR
    type(particle),intent(inOut), dimension(:,:) :: teilchen
    logical,intent(in),optional :: check_in
    integer :: i,j
    logical :: check

    check = .false.
    if (present(check_in)) check = check_in

    if (DoPr(2)) write(*,*) 'Updating energies of Particles'

    Do i=1,Size(teilchen,dim=1)
       Do j=1,Size(teilchen,dim=2)
          if (teilchen(i,j)%ID < 0 ) exit
          if (teilchen(i,j)%ID == 0 ) cycle
          call energyDetermination(teilchen(i,j),check=check)
       end do
    end do

  end subroutine updateEnergies


  !*************************************************************************
  !****s* energyCalc/energyDetermination
  ! NAME
  ! subroutine energyDetermination(teilchen,betaToCF,countMax_in,error_out,check,warn,ForbidCoulomb)
  !
  ! PURPOSE
  ! This subroutine determines the one-particle energy E with a
  ! momentum-dependent scalar potential that is defined in the local
  ! rest frame.
  !
  ! NOTES
  ! * We keep p(1:3) and m constant.
  ! * To solve: p**2=E**2-p(1:3)**2=(m+s(E))**2
  ! * The scalar potential s is only defined in the LRF.
  ! * First we try to find the solution f(E)=0 via a
  !   the newton method,
  !     f(E)=m_eff**2-p**2=m_eff**2+p(1:3)**2-E**2
  !   with
  !     m_eff=m+s(E)
  !   by the iteration
  !     E(i+2)=E(i+1)-f(i+1)*(E(i+1)-E(i))/(f(i+1)-f(i)).
  ! * If the newton method fails, we try it again via a
  !   bisection method
  ! * At the end, the Coulomb potential is added.
  !
  ! INPUTS
  ! * type(particle) :: teilchen --  Particle to be considered.
  ! * real,dimension (1:3), OPTIONAL :: betaToCF -- Velocity of
  !   "Calculation Frame" in the frame where the particle's momentum is
  !    defined.
  ! * integer, OPTIONAL :: countMax_in  -- maximal number of cycles to find
  !   a solution (default=10000)
  ! * logical, OPTIONAL :: check -- check for energy conservation
  !   (Default: .false.)
  ! * logical, OPTIONAL :: warn -- if true, a warning message is issued,
  !   if the incomming momentum is too large (otherwise, this is skipped
  !   silently). (Default: .true.)
  ! * logical, OPTIONAL :: ForbidCoulomb -- if true, the Coulomb potential
  !   is NOT added to p0 (needed by 'propagate', where the coulomb force
  !   is treated seperately) (Default: .false.)
  ! OUTPUT
  ! * type(particle) :: teilchen
  ! * real, OPTIONAL    :: error_out    -- approximate error of E = sqrt(f(e))
  !*************************************************************************
  subroutine energyDetermination(teilchen, betaToCF, countMax_in, error_out, check, warn, ForbidCoulomb)
    use IdTable, only: rho,omegaMeson,phi,isMeson
    use particleDefinition
    use minkowski, only: abs4
    use offShellPotential, only: treatParticleOffShell, get_offshell_debug
    use output, only: WriteParticle_debug, WriteParticle
    use hist2Df90
    use callStack, only: traceback
    use dichteDefinition
    use densityModule, only: densityAt
    use coulomb, only: emfoca, getCoulombFlag

    type(particle), intent(inOut) :: teilchen
    real, optional, dimension (1:3), intent(in) :: betaToCF
    integer, optional, intent(in) :: countMax_in
    real, optional, intent(out)   :: error_out
    logical, optional, intent(in) :: check, warn, ForbidCoulomb

    integer, save :: failures=0, numtry=0  !,hugeCounter=0
    real :: mass,psqr,E0, E_min, Cpot
    real,dimension(1:3) :: f,E  ! f=f(E), need 1:3 for iteration
    logical :: CFBoostFlag      ! Determines whether it's necessary to boost to Calculation Frame
    integer :: counter,countmax
    integer,parameter :: countmax_default=100
    real,parameter :: maxError=accuracy**2   !wished accuracy of f(E) in GeV**2
    ! --> Energy is correct within (maxError)**(1/2)
    logical,parameter :: debugflag=.false.
    real,parameter :: delta_E = 0.01
    logical :: checkE,success,outOfBounds_offshell,warn_in,DoCoulomb
    type(dichte) :: density

    warn_in=.true.
    if(present(warn)) warn_in=warn

!     If(AbsMom(teilchen).gt.10000) then
!        !*********************************************************************
!        ! This is to explore some problem which seems to occur in  the
!        ! energyCorrection loops:
!        ! When no solution can be found, the regula falsi is producing a X to
!        ! scale the momentum which is way to high!!
!        ! This should be no problem to finished jobs. It just crashes the code
!        ! from time to time since we can not find a solution for 100 TEV
!        ! particles in this routine here.
!        ! Oliver
!        !*********************************************************************
!        if(warn_in) then
!           write(*,*)
!           write(*,*) 'WARNING: Particle with huge momentum in EnergyDetermination'
!           call errorMessage
!           hugeCounter=hugeCounter+1
!           write(*,*) 'It is the ',hugeCounter, 'th time'
!        end if
!        teilchen%momentum(0)=FreeEnergy(teilchen)
!        return
!        write(*,*)
!     end if

!!$    If(AbsMom(teilchen).gt.1000) then
!!$       write(*,*) "energyDetermination: Particle with large momentum: ",AbsMom(teilchen)
!!$       call WriteParticle(6,0,0,teilchen)
!!$       call traceback("calllist",-1)
!!$    endif


    countMax=countMax_default
    if(present(countMax_in))  countMax=countMax_in

    CFBoostFlag=.false.
    If (Present(betaToCF)) CFBoostFlag=.true.

    checkE=.false.
    if (present(check)) checkE=check

    DoCoulomb = .true.
    if (Present(ForbidCoulomb)) DoCoulomb = .not.ForbidCoulomb

    !----------------------------------------------

    E0=teilchen%momentum(0)

    mass=teilchen%mass
    psqr=Sum(teilchen%momentum(1:3)**2)
    if (getCoulombFlag()) then
      E_min = max(0.,sqrt(psqr)-0.2)  ! just a guess (potentials can be negative)
    else
      E_min = sqrt(psqr+0.001**2)
    end if
    density = densityAt(teilchen%position)

    if (isMeson(teilchen%ID) &
         .and. treatParticleOffShell(teilchen%ID,teilchen%OffShellParameter) &
         .and. checkE &
         .and. E0>sqrt(psqr)) then

       ! special treatment for offshell mesons, since f(E)=0 may have more than one solution:
       ! find a solution to f(E)=0, restricting E to E0 +- a few MeV

       if (get_offshell_debug()) then
          call MesonHistFailures(1)
       end if

       numtry=numtry+1

       success=bisection(E0-delta_E,E0+delta_E)

       if (.not. success) then
          failures=failures+1
          if (get_offshell_debug()) then !  .and. checkE
             teilchen%momentum(0)=E0
             call MesonHistFailures(2)
          end if
          write(*,*) "WARNING in EnergyDetermination: Impossible to determine energy for offshell meson!"
          write(*,*) "Number of Failures: ",failures,' (',float(failures)/float(numtry)*100., '%)'
          write(*,*)
          !call plot_f(teilchen)
          teilchen%ID = 0
          !call traceback('Error in energyDetermination',0)
       end if

    else

       ! default case

       outOFBounds_offShell=.false.
       success=newton()
       If (outOFBounds_offShell) success=vary_mass_offshell()
       if (.not. success) success=bisection()
       if (.not. success) then
          write(*,*) "WARNING in EnergyDetermination: Impossible to determine energy!"
          write(*,*)
          if(warn_in) then
             call plot_f(teilchen)
             call traceback('Error in energyDetermination')
          end if
       end if

    end if

    if (checkE .and. abs(teilchen%momentum(0)-E0)>delta_E) then
       write(*,*) "WARNING in EnergyDetermination: Energy not conserved!",teilchen%momentum(0)-E0,E0
       call WriteParticle_debug(teilchen)
       !call plot_f(teilchen)
       !stop
    end if


  contains

    !***********************************************************************
    subroutine MesonHistFailures(iTyp)
      integer, intent(in) :: iTyp

      logical, save :: firstTime=.true.
      type(histogram2D),save :: hist2D_rho,hist2D_omega,hist2D_phi

      select case (iTyp)
      case (1)
         if (firstTime) then
            call createHist2D(hist2D_rho,&
                 'energy violation as function of mass and momentum',&
                 (/0.,0./),(/2.,2./),(/0.02,0.02/))
            call createHist2D(hist2D_omega,&
                 'energy violation as function of mass and momentum',&
                 (/0.,0./),(/2.,2./),(/0.02,0.02/))
            call createHist2D(hist2D_phi,&
                 'energy violation as function of mass and momentum',&
                 (/0.,0./),(/2.,2./),(/0.02,0.02/))
            firstTime=.false.
         end if

         select case (teilchen%ID)
         case (rho)
            call AddHist2D(hist2D_rho,(/abs4(teilchen%momentum),absMom(teilchen)/),0.,1.)
         case (omegaMeson)
            call AddHist2D(hist2D_omega,(/abs4(teilchen%momentum),absMom(teilchen)/),0.,1.)
         case (phi)
            call AddHist2D(hist2D_phi,(/abs4(teilchen%momentum),absMom(teilchen)/),0.,1.)
         end select

      case (2)

         select case (teilchen%ID)
         case (rho)
            call AddHist2D(hist2D_rho,(/abs4(teilchen%momentum),absMom(teilchen)/),1.)
         case (omegaMeson)
            call AddHist2D(hist2D_omega,(/abs4(teilchen%momentum),absMom(teilchen)/),1.)
         case (phi)
            call AddHist2D(hist2D_phi,(/abs4(teilchen%momentum),absMom(teilchen)/),1.)
         end select

         if(mod(failures,100)==0) then
            call writeHist2D_Gnuplot(hist2D_rho,44,file='Energy2D_rho.dat')
            call writeHist2D_Gnuplot(hist2D_omega,45,file='Energy2D_omega.dat')
            call writeHist2D_Gnuplot(hist2D_phi,46,file='Energy2D_phi.dat')
         end if

      end select

    end subroutine MesonHistFailures

    !***********************************************************************

    ! Try maximizing p_0 such that equation is still fulfilled
    logical function vary_mass_offshell()
      use particleDefinition
      type(particle) :: part
      real, parameter :: delta_pNull=0.0001
      real :: error
      part=teilchen
      vary_mass_offshell=.true.
      Do
         part%momentum(0)=part%momentum(0)+delta_pNull
         error=m_eff_sqr(part,density)-(part%momentum(0)**2-psqr)
         if(abs(error).lt.maxError) then !  Successful!!!
            If (DebugFlag) Print *, 'SuccessFull with f(2)',f(2),E(2)
            If(Present(error_out)) error_out=sqrt(abs(f(2)))
            return
         else
            part%momentum(0)=part%momentum(0)-delta_pNull
            exit
         end if
         if(part%momentum(0).gt.10.) then
            ! FAILURE
            write(*,*) 'FAILURE in vary_mass_offshell'
            vary_mass_offshell=.false.
            exit
         end if
      End do

    end function vary_mass_offshell


    !***********************************************************************

    ! solve f(E)=0 by Newton Method
    logical function newton()

      newton=.true.

      if (DebugFlag) Print *, '**In EnergyDetermination (Newton)'

      !Make Vacuum Ansatz to evaluate starting values
      E(1)=SQRT(psqr+mass**2)
      teilchen%momentum(0)=E(1)
      f(1)=m_eff_sqr(teilchen,density)-(E(1)**2-psqr)

      if(abs(f(1)).lt.maxError) then !  Successful!!!
         If (DebugFlag) Print *, 'SuccesFull with f(1)',f(1)
         If(Present(error_out)) error_out=sqrt(abs(f(1)))
         return
      end if

      E(2)=sqrt(max(f(1)+E(1)**2,E_min**2))
      teilchen%momentum(0)=E(2)
      f(2)=m_eff_sqr(teilchen,density)-(E(2)**2-psqr)

      if(abs(f(2)).lt.maxError) then !  Successful!!!
         If (DebugFlag) Print *, 'SuccessFull with f(2)',f(2),E(2)
         If(Present(error_out)) error_out=sqrt(abs(f(2)))
         return
      end if

      If(debugFlag) Print * , 'Begin Iteration'

      ! Begin Iteration

      counter=0
      do
         if(Abs(f(2)-f(1)).eq.0) then
            if (DebugFlag) then
               write(*,*)'Error in EnergyDetermination (Newton). Derivative too small:', f(1:2),f(2)-f(1)
               call errorMessage
               write(*,*)'Counter:' , counter
            end if
            newton=.false.
            return
         End if

         E(3)=E(2)-f(2)*(E(2)-E(1))/(f(2)-f(1))
         teilchen%momentum(0)=E(3)
         f(3)=m_eff_sqr(teilchen,density)-(E(3)**2-psqr)

         If(debugFlag)  write(11,*) e(3),f(3)

         if(abs(f(3)).lt.maxError) then !Successful!!!
            If (DebugFlag) Print *, 'SuccesFull with f(3)',f(3)
            If (DebugFlag) Print *, 'counter=',counter
            If (Present(error_out)) error_out=sqrt(abs(f(3)))
            if (teilchen%momentum(0)<0.) newton=.false.
            return
         end if

         if(counter.gt.countmax) then
            if (DebugFlag) then
               write(*,*)'Error in Energydetermination (Newton):  counter.gt.countmax'
               call errorMessage
               call plot_f(teilchen)
            end if
            newton=.false.
            return
         end if
         counter=counter+1
         E(1)=E(2)
         E(2)=E(3)
         f(1)=f(2)
         f(2)=f(3)
      end do

    end function newton

    !***********************************************************************


    ! solve f(E)=0 by Bisection Method
    logical function bisection(E1,E2)

      real,intent(in),optional :: E1,E2 ! initial boundaries

      bisection=.true.

      if (Present(E1) .and. Present(E2)) then

         E(1)=max(E_min,E1)
         teilchen%momentum(0)=E(1)
         f(1)=m_eff_sqr(teilchen,density)-(E(1)**2-psqr)

         if(abs(f(1)).lt.maxError) then !  Successful!!!
            If(Present(error_out)) error_out=sqrt(abs(f(1)))
            return
         end if

         E(2)=E2
         teilchen%momentum(0)=E(2)
         f(2)=m_eff_sqr(teilchen,density)-(E(2)**2-psqr)

         if(abs(f(2)).lt.maxError) then !  Successful!!!
            If(Present(error_out)) error_out=sqrt(abs(f(2)))
            return
         end if

         counter=0
         do while (f(1)*f(2)>0. .and. counter<20)
            E(3)=(E(1)+E(2))/2.  ! half interval
            teilchen%momentum(0)=E(3)
            f(3)=m_eff_sqr(teilchen,density)-(E(3)**2-psqr)
            if(abs(f(3)).lt.maxError) then !  Successful!!!
               If(Present(error_out)) error_out=sqrt(abs(f(2)))
               return
            end if
            if (f(1)*f(3)<0.) then
               f(2)=f(3)
               E(2)=E(3)
            else if (f(2)*f(3)<0.) then
               f(1)=f(3)
               E(1)=E(3)
            else if (abs(f(1))<abs(f(2))) then
               f(2)=f(3)
               E(2)=E(3)
            else if (abs(f(2))<abs(f(1))) then
               f(1)=f(3)
               E(1)=E(3)
            else
               f(mod(counter,2)+1)=f(3)
               E(mod(counter,2)+1)=E(3)
            end if
            counter=counter+1
         end do

      else

         !if (E0==0. .and. debugFlag) write(*,*) 'EnergyDetermination (Bisection): initializing E!'

         ! set starting values: E(1) and E(2)
         E(1)=E_min
         teilchen%momentum(0)=E(1)
         f(1)=m_eff_sqr(teilchen,density)-(E(1)**2-psqr)

         if(abs(f(1)).lt.maxError) then !  Successful!!!
            If(Present(error_out)) error_out=sqrt(abs(f(1)))
            return
         end if

         E(2)=E(1)+0.001
         teilchen%momentum(0)=E(2)
         f(2)=m_eff_sqr(teilchen,density)-(E(2)**2-psqr)

         if(abs(f(2)).lt.maxError) then !  Successful!!!
            If(Present(error_out)) error_out=sqrt(abs(f(2)))
            return
         end if

         ! increase E(2)
         counter=0
         do while (f(1)*f(2)>0. .and. counter<20)
            E(2)=2*E(2)-E(1)  ! double interval
            If (DebugFlag) write(*,*) 'increasing E(2): ',E(2)
            teilchen%momentum(0)=E(2)
            f(2)=m_eff_sqr(teilchen,density)-(E(2)**2-psqr)
            if(abs(f(2)).lt.maxError) then !  Successful!!!
               If(Present(error_out)) error_out=sqrt(abs(f(2)))
               return
            end if
            counter=counter+1
         end do

      end if


      if (f(1)*f(2)>0) then
         write(*,*) "Error in EnergyDetermination (Bisection): Impossible to find starting values!"
         if (Present(E1) .and. Present(E2)) write(*,*) E1,E2
         call errorMessage
         bisection=.false.
         return
         !call plot_f(teilchen)
         !Stop 'Stop in EnergyDetermination (Bisection)'
      end if

      ! Begin Iteration
      !If (DebugFlag) write(*,*) "Bisection: starting iteration!"

      counter=0
      do
         if (f(1)*f(2)>0) then
            write(*,*)'Error in EnergyDetermination (Bisection):  Same Sign!'
            call errorMessage
            call plot_f(teilchen)
            Stop 'Stop in EnergyDetermination (Bisection)'
         end if

         E(3)=(E(1)+E(2))/2
         teilchen%momentum(0)=E(3)
         f(3)=m_eff_sqr(teilchen,density)-(E(3)**2-psqr)

         if(abs(f(3))<maxError .or. E(2)-E(1)<accuracy) then !Successful!!!
            If (DebugFlag) write(*,*) 'Bisection: Success!'
            If(Present(error_out)) error_out=sqrt(abs(f(3)))
            return
         end if

         if (counter>40) write(*,*) E,f

         if(counter.gt.countmax) then
            write(*,*)'Error in EnergyDetermination (Bisection):  counter.gt.countmax'
            call errorMessage
            bisection=.false.
            return
            !call plot_f(teilchen)
            !Stop 'Stop in EnergyDetermination (Bisection)'
         end if
         counter=counter+1

         if (f(1)*f(3)<0.) then
            If (DebugFlag) write(*,'(A)',ADVANCE='NO') 'L'
            E(2)=E(3)
            f(2)=f(3)
         else if (f(2)*f(3)<0.) then
            If (DebugFlag) write(*,'(A)',ADVANCE='NO') 'U'
            E(1)=E(3)
            f(1)=f(3)
         else
            write(*,*) 'Error in EnergyDetermination (Bisection): ZERO!'
            write(*,*) 'counter = ',counter
            call errorMessage
            call plot_f(teilchen)
            call traceback('Stop in EnergyDetermination (Bisection)')
         end if
      end do

    end function bisection

    !***********************************************************************

    ! This subroutine plots the function m_eff**2(teilchen,ener)-(Energy**2-Sum(p**2)) as a
    ! function of energy. Zoomes into the region where the function is zero.
    ! Useful for debugging!!
    subroutine plot_f(teilchen_in)
      use output

      type(particle),intent(in) :: teilchen_in
      type(particle) :: teilchen
      real :: dE
      real :: resu,energy,me2
      real :: down,ener!,old
      integer :: i!,j,ende
      integer :: max=10000

      teilchen=teilchen_in
      down=E_min
      dE=1./float(max)
      open(111,file='plot_f.dat')
      !do j=1,5 ! Loop which zooms into the region of the zero
      !ende=max
      do i=0,max ! Loop which steps through the energy
         energy=down+float(i)*dE
         teilchen%momentum(0)=energy
         me2 = m_eff_sqr(teilchen,density,ener)
         resu= me2 - (energy**2-psqr)
         !if(i.eq.0) old=resu
         write(111,'(4E12.4)') teilchen%momentum(0),resu,me2,ener
         !if(old*resu.lt.0) then
         !   ! Change of sign in the last dE step
         !   down=energy-dE
         !   ende=i+10 ! stop plotting after the next 10 steps
         !end if
         !old=resu
         !if(i.gt.ende) exit
      end do
      !dE=dE/float(max) ! Decrease step size during the zooming process
      !end do
      close(111)
    end subroutine plot_f

    !***********************************************************************

    subroutine errorMessage
      use output, only :  WriteParticle_debug

      write(*,'(a,3E20.7)')'SQRT( f(1:3) )       :', SQRT(Abs(f))
      write(*,'(a,3E20.7)')'f(1:3)               :', f
      write(*,'(a,3E20.7)')'E(1:3)               :', E
      write(*,'(a,i4)')    'Id of particle       :', teilchen%Id
      write(*,'(a,i4)')    'Charge of particle   :', teilchen%charge
      write(*,'(a,E12.4)') 'Mass of particle     :', teilchen%mass
      write(*,'(a,4E12.4)')'Momentum of particle :', teilchen%momentum
      write(*,'(a,E12.4)') 'Abs. mom. of particle:', absMom(teilchen)
      write(*,'(a,L2)')    'Perturbative         :', teilchen%perturbative
      write(*,'(a,E12.4)') 'Offshell parameter   :', teilchen%offshellParameter
      If(present(betaToCF)) then
         write(*,'(a,3E12.4)')'beta to Calculation frame :', betaToCF
      else
         write(*,*)'optional parameter "betaToCF" not used'
      end if

      call  WriteParticle_debug(teilchen)

    end subroutine errorMessage

    !***********************************************************************
    !****f* energyDetermination/m_eff_sqr
    ! NAME
    ! real function m_eff_sqr(teilchenIN,density,enerOut)
    ! PURPOSE
    ! Evaluate m_eff**2=(m+S)**2 of a given particle by boosting to the LRF-Frame,
    ! where S is the scalar potential.
    !***********************************************************************
    real function m_eff_sqr (teilchenIN, density, enerOut)
      use particleDefinition
      use potentialModule, only : potential_LRF, massDetermination
      use lorentzTrafo, only: lorentzCalcBeta, lorentz
      use offShellPotential, only : hamiltonFunc_offshell
      use minkowski, only: SP

      type(particle),intent(in):: teilchenIN
      type(dichte),intent(in) :: density
      type(particle) :: teilchen
      real, optional :: enerOut
      real, dimension(1:3) :: betaToLRF

      !Evaluate scalar potential in LRF:
      teilchen=teilchenIN
      ! (1.1) Boost to calculation frame
      if (CFBoostFlag) call lorentz(betaToCF,teilchen%momentum, 'energyDetermination,m_eff_sqr')
      ! (1.2) Boost from calculation frame to LRF
      if (density%baryon(0)>1E-8 .and. sum(density%baryon(1:3)**2)>1E-8) then
         betaToLRF = lorentzCalcBeta (density%baryon, 'energyDetermination,m_eff_sqr')
         call lorentz(betaToLRF,teilchen%momentum,'energyDetermination,m_eff_sqr')
      else
         betaToLRF=0.
      end if

      ! (1.3) Evaluate scalar potential V_S :(m+V_S)**2=p(0)**2-p(1:3)**2=p**2 =>
      if(treatParticleOffShell(teilchen%ID,teilchen%OffShellParameter)) then
         call massDetermination(teilchen,success=success,verbose=.false.)
         outOfBounds_offshell=.false.
         if(.not.success .or. SP(teilchen%momentum,teilchen%momentum)<0.) then
            ! BAD Solution !!
            teilchen%momentum(0)=999999999.
         else
            teilchen%momentum(0)=HamiltonFunc_offshell(teilchen,outOfBounds_offshell,.true.)
         end if
      else
         teilchen%momentum(0)=FreeEnergy(teilchen)+potential_LRF(teilchen,density,addCoulomb=.false.)
      end if

      if (DoCoulomb) then
         ! This Lorentz transformation to LRF is needed because Cpot from emfoca is defined in CF
         Cpot = emfoca(teilchen%position,teilchen%momentum(1:3),teilchen%charge,teilchen%ID) * sqrt(1.-sum(betaToLRF(1:3)**2))
         teilchen%momentum(0)=teilchen%momentum(0)+Cpot
      end if

      if(present(enerOut)) enerOUT=teilchen%momentum(0)

      m_eff_sqr=teilchen%momentum(0)**2-Dot_product(teilchen%momentum(1:3),teilchen%momentum(1:3))

    end function m_eff_sqr

  end subroutine energyDetermination



  !*************************************************************************
  !****s* energyCalc/energyCorrection
  ! NAME
  ! subroutine energyCorrection (srqts, betaToLRF, betaToCM, mediumAtCollision,
  ! finalState, successFlag, potentialFailure, verbose)
  !
  ! PURPOSE
  ! Gets the finalState particles with vacuum kinematics in the CM Frame.
  ! Also the real srts of the final state is handed over to the routine.
  ! Then this routine tries to solve energy and momentum conservation by the
  ! ansatz described in
  ! |html  <a href="../../Documentation_Extra/crossSections/Xsections/node25.html"> the cross section documentation. </a>.
  ! Therefore it evaluates the scaling factor "x" which is used to scale the
  ! CM-momenta of the final states. This is done by  a Regula-Falsi-Method.
  ! We search of a zero of the function "deltaSrts" which depends on this
  ! scaling factor x.
  !
  ! IMPORTANT :
  ! Output are the final state particles in the calculation frame!
  ! INPUTS
  ! * real, intent(in) :: srtS
  ! * real, intent(in), dimension(1:3) :: betaToLRF -- velocity of LRF - NOT USED right now!!!
  ! * real, intent(in), dimension(1:3) :: betaToCM  -- velocity of CM frame in calculation frame
  ! * type(medium), intent(in) :: mediumAtCollision
  ! * type(particle), dimension(:) , intent (INOUT)  :: finalState -- In CM frame
  ! * logical, OPTIONAL :: verbose -- flag print messages on/off (default: on)
  ! OUTPUT
  ! * type(particle), dimension(:) , intent (INOUT)  :: finalState -- In calculation frame
  ! * logical, intent(out) :: successFlag
  ! * logical, optional :: potentialFailure -- is returned .true. if it fails due to neglecting perturbative potentials
  ! NOTES
  ! If the iteration fails successFlag is set to .false. as output.
  ! Now is used also in the RMF mode (only function deltaSrts is modified).
  !*************************************************************************
  subroutine energyCorrection(srts, betaToLRF, betaToCM, mediumAtCollision, &
       & finalState, successFlag, potentialFailure, verbose)

    use mediumDefinition, only: medium
    use particleDefinition, only: particle,sqrts
    use LorentzTrafo, only: lorentz
    use random, only: rn
    use mesonPotentialModule, only: getNoPertPot_meson
    use baryonPotentialModule, only: getNoPertPot_baryon
    use RMF, only : getRMF_flag

    real, intent(in) :: srtS
    real, intent(in), dimension(1:3) :: betaToLRF, betaToCM
    type(medium), intent(in) :: mediumAtCollision
    type(particle), dimension(:), intent (inout)  :: finalState
    logical, intent(out) :: successFlag
    logical, optional,intent(out) :: potentialFailure
    logical, OPTIONAL,            intent(in)    :: verbose

    integer, parameter :: maxLoops=10  ! default maximal number of tries to find a solution for x

    ! local variables
    type(particle) , Allocatable, dimension(:):: teilchen
    integer :: i, counterEqual
    real, dimension (-1:1) :: f, x

    logical, parameter :: debug=.false.
    real, parameter :: accuracySiter=1E-10
    logical :: verb

    successFlag=.false.
    if(Present(potentialFailure)) potentialFailure=.false.

    verb = .true.
    if (present(verbose)) verb = verbose

    if(.not.getRMF_flag()) then
       if(getNoPertPot_meson().and.getNoPertPot_baryon()) then
          ! If perturbative potentials are switched off and finalState is perturbative,  then we can only initialize
          ! particles above the mass threshold.
          ! Using the following if-statement endless loops are prohibited and computing time shall be saved.
          if(finalState(1)%perturbative.and.(Sum(finalstate%mass).gt.srts)) then
             if(Present(potentialFailure)) potentialFailure=.true.
             return
          end if
       end if
    end if

    Allocate(teilchen(1:size(finalState,dim=1)))

    ! Initialize
    x(0)=1.
    f(0)=deltaSrts( x(0) )


    If (debug) write(*,*) 'f(0)=', f(0)
    If (debug) write(*,*) 'srst=', srts

    If (Abs(f(0)).lt. accuracySiter) then
       call setFinalState
       successFlag=.true.
       return
    end if

    x(1)=0.9 ! first wild guess for starting point
    f(1)=deltaSrts( x(1) )

    If (debug) write(*,*) 'f(1)=', f(1)

    If (Abs(f(1)).lt. accuracySiter) then
       call setFinalState
       successFlag=.true.
       return
    end if

    If(Abs(f(1)-f(0)).lt.1E-6) then
       ! Try to find new starting point for the regula-falsi-method
       counterEqual=0
       findX_loop : do
          x(1)=0.5+rn()
          f(1)=deltaSrts( x(1) )
          If(Abs(f(1)-f(0)).gt.1E-6) exit findX_loop
          If(counterEqual.gt.10) then
             if (verb) then
                write(*,*) ' Problem in energyCorrection . f(1) = f(0) :'
                write(*,*) f(0:1)
                write(*,*) ' x=',x(0:1)
                write(*,*) ' Hence iteration not possible'
             end if
             successFlag=.false.
             deAllocate(Teilchen)
             return
          end if
          counterEqual=counterEqual+1
       end do findX_loop
    end if


    ! Start Iteration
    Do i=1, maxLoops
       ! Store results of last loop
       x(-1:0)=x(0:1)
       f(-1:0)=f(0:1)
       ! Evaluate new guess
       x(1)=x(0)-f(0)*(x(0)-x(-1))/(f(0)-f(-1))
       If (debug) write(*,*) 'Loop', i, 'x=',x(1)
       If (x(1).lt.0) then
          x(1)=-x(1)/float(maxLoops)
          ! /float(maxLoops) because if there is a solution then probably close to 0
       end if
       f(1)=deltaSrts(x(1))
       If (debug) write(*,*) 'Loop', i, 'f(x)=',f(1)
       if(abs(f(1)).lt.accuracySiter) then
          call setFinalState
          successFlag=.true.
          return
       else if(f(1)-f(0).eq.0.) then
          if (verb) then
             write(*,*)' WARNING : In energycorrection f(1)-f(0).eq.0',f(0:1),x(0:1)
             write(*,*) 'Hence iteration not possible'
          end if
          successFlag=.false.
          deAllocate(Teilchen)
          return
       end if
    End do


    !    If(.not.successFlag) write(*,*) 'energyCorrection:  Energy correction failed'

  contains

    subroutine setFinalState
      integer :: i
      Do i=1,size(teilchen,dim=1)
         call lorentz(-betaToCM,teilchen(i)%momentum(0:3), 'energyDetermination,setFinalState')  ! Boost to Calculation frame
      End do
      finalState = teilchen
      If (debug) write(*,*) 'Setting finalState. Sqrts=',sqrtS(finalState(1),finalState(2))
      deAllocate(Teilchen)
    end subroutine setFinalState


    !***************************************
    !****f* energyCorrection/deltaSrts
    ! NAME
    ! real function deltaSrts(x)
    ! PURPOSE
    ! Evaluates sqrt(s)-Sum(Energies in CM frame). This value is the size of the energy non-conservation.
    ! The sum is depend on the scaling factor for the three momenta
    !***************************************
    real function deltaSrts(x)

      use densitymodule, only : energyDeterminationRMF

      real, intent(in) :: x
      integer :: i

      teilchen=finalState

      Do i=1, size( finalState, dim=1 )
         teilchen(i)%momentum(1:3)=x * teilchen(i)%momentum(1:3)
      End do

      If(debug) Print *, 'moms vorher=' ,  sum(teilchen%momentum(0)),sum(teilchen%momentum(1)), &
           & sum(teilchen%momentum(2)),sum(teilchen%momentum(3))

      Do i=1, size( finalState, dim=1 )
         IF(debug) Write(*,*) 'teilchen ',i,'vorher :',teilchen(i)%mass, teilchen(i)%momentum,teilchen(i)%perturbative
         if( .not.getRMF_flag() ) then
            call energyDetermination(teilchen(i),-betaToCM)
         else
            call energyDeterminationRMF(teilchen(i))
         end if
         IF(debug) Write(*,*) 'teilchen ',i,'nachher :',teilchen(i)%mass, teilchen(i)%momentum
      End do

      deltaSrts=srts-Sum(teilchen%momentum(0))

      If(debug) then
         write(*,*) 'moms nach energyDetermination=' ,  &
              & sum(teilchen%momentum(0)),sum(teilchen%momentum(1)), &
              & sum(teilchen%momentum(2)),sum(teilchen%momentum(3))
         write(*,*) 'real srst=', srts
         write(*,*) 'srts=', sqrtS(teilchen(1),teilchen(2))
         write(*,*) 'moms nachher=', &
              & sum(teilchen%momentum(0)),sum(teilchen%momentum(1)), &
              & sum(teilchen%momentum(2)),sum(teilchen%momentum(3))
         write(*,*) 'delta=',deltaSrts
         write(*,*) '-------------------------------------------------------------------------------------------------'
      end if
    end function deltaSrts


  end subroutine energyCorrection





end module energyCalc
