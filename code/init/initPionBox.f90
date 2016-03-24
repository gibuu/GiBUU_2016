!*********************************************************************
!****m* /initPionBox
! NAME
! module initPionBox
! PURPOSE
! Initializes pions for a box of pions
!*********************************************************************
module initPionBox

  implicit none

  Private

  !*******************************************************************
  !****g* initPionBox/thermalInit
  ! SOURCE
  !
  logical, save :: thermalInit = .false.
  ! PURPOSE
  ! flag how to initialize
  !*******************************************************************

  !*******************************************************************
  !****g* initPionBox/nDens
  ! SOURCE
  !
  real, save :: nDens = 1.0
  ! PURPOSE
  ! particle density [fm^-3]
  !*******************************************************************

  !*******************************************************************
  !****g* initPionBox/ChargeSelection
  ! SOURCE
  !
  integer, save :: ChargeSelection = 0
  ! PURPOSE
  ! define the type of the charge selection:
  ! * 0: only pi0
  ! * 1: 50% pi+, 50% pi-
  ! * 2: 33% for +,0,-
  !*******************************************************************

  !*******************************************************************
  !****g* initPionBox/pInit
  ! SOURCE
  !
  real, save :: pInit = 0.5
  ! PURPOSE
  ! initial momentum of particles [GeV/c]
  !*******************************************************************

  !*******************************************************************
  !****g* initPionBox/BoostZ
  ! SOURCE
  !
  real, save :: BoostZ = 0.0
  ! PURPOSE
  ! additional boost for all particles in z-direction
  !*******************************************************************

  !*******************************************************************
  !****g* initPionBox/Temp
  ! SOURCE
  real,dimension(101:122), save :: Temp = 0.0
  ! PURPOSE
  ! for thermal init: temperature of every meson species in GeV,
  ! if larger than 0. otherwise this species is not initialized
  !*******************************************************************

  !*******************************************************************
  !****g* initPionBox/Fugacity
  ! SOURCE
  real,dimension(101:122), save :: Fugacity = 1.0
  ! PURPOSE
  ! for thermal init: fugacity of every meson species.
  !*******************************************************************

  Public :: initializePionBox


  real,dimension(101:122), save :: densMes = 0

contains

  !*******************************************************************
  !****s* initPionBox/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Reads input out of jobcard. Namelist 'initBox'.
  !*******************************************************************
  subroutine initInput
    use output, only: Write_ReadingInput
    use callstack, only: traceBack
    use constants, only: mPi

    !*****************************************************************
    !****n* initPionBox/PionBox
    ! NAME
    ! NAMELIST PionBox
    ! PURPOSE
    ! Includes the input parameters:
    ! * thermalInit
    ! * nDens
    ! * ChargeSelection
    ! * pInit
    ! * BoostZ
    !*****************************************************************
    NAMELIST /PionBox/ thermalInit, &
         nDens,ChargeSelection,pInit,BoostZ, &
         Temp, Fugacity

    integer :: ios, iMes
    character(20), dimension(0:2), parameter :: cSel = (/ &
         'only pi0        ',&
         '50% pi+, 50% pi-',&
         '33% for +,0,-   ' /)
    real :: eDens,beta(1:2)

    call Write_ReadingInput('PionBox',0)
    rewind(5)
    read(5,nml=PionBox,iostat=ios)
    call Write_ReadingInput('PionBox',0,ios)

    if (thermalInit) then
       write(*,*) '### thermal init ###'

       do iMes=101,122
          if ((iMes>=114).and.(iMes<=121)) cycle
          if (Temp(iMes) >= 0.01) then
             densMes(iMes) = integrate(iMes, Temp(iMes)) * Fugacity(iMes)
          end if
       end do

       do iMes=101,122
          if (Temp(iMes) >= 0.01) then
             write(*,'("  ",i3,"  ",f5.3,"  ",f5.3," -> ",1P,e13.4)') iMes, &
                  Temp(iMes), Fugacity(iMes), densMes(iMes)
          end if
       end do

       write(*,*) 'total density of particles: ',sum(densMes)

    else
       write(*,*) '### non-thermal init ###'

       select case(ChargeSelection)
       case (0:2)
          write(*,*) 'charge selection: ',ChargeSelection,"= ",&
               cSel(ChargeSelection)
       case default
          write(*,*) 'charge selection: ',ChargeSelection,"= WRONG!"
          call traceBack('wrong input value')
       end select
       write(*,'(A,F8.3,A)') ' particle density: nDens =',nDens,' fm^(-3)'
       write(*,'(A,F8.3,A)') ' initial momentum: pInit =',pInit,' GeV/c'

       write(*,*)
       eDens = sqrt(mPi**2+pInit**2)*nDens

       ! the temperature is given by
       !   eDens/(m*nDens) == 3/betam + K1(betam)/K2(betam)
       ! As approx, we can give K1/K2 ~ betam/2
       ! This is then solved for beta:
       ! (The resulting T is off by 5-10 MeV)
       beta(1) = eDens + sqrt(eDens**2-6*(mPi*nDens)**2)
       beta(2) = eDens - sqrt(eDens**2-6*(mPi*nDens)**2)
       beta = beta/(mPi**2*nDens)
       write(*,'(A,F8.3,A)') ' energy density:   eDens =',eDens,' GeV fm^(-3)'
       write(*,'(A,F6.2,A)') ' Estimate temperature: T = ',1000/beta(2),' MeV'
       write(*,'(A,F6.2,A)') ' Estimate fugacity   : -- to be done-- '

       write(*,*)
       write(*,'(A,F8.3,A)') ' additional boost: betaZ =',BoostZ

    end if

    call Write_ReadingInput('PionBox',1)

  contains
    real function integrate(ID, T)

      use constants, only: pi
      use particleProperties, only: hadron
      use mesonWidth, only: fullWidthMeson

      integer, intent(in) :: ID
      real, intent(in) :: T

      integer :: im, nm
      real, parameter :: massMax = 5.0 ! upper boundary of integration
      real :: mass0, gamma0
      real :: mmin, mmax, m
      real :: ymin, ymax, dy, y
      real :: gamma, spectral, intfac, deg, sum
      real, parameter :: dy0 = pi/200.

      integrate = 0.0

      mass0  = hadron(id)%mass
      gamma0 = hadron(id)%width
      deg = (hadron(ID)%Spin*2+1)*(hadron(ID)%isoSpinTimes2+1)

      if (gamma0 < 1e-3) then
         if (massMax > mass0) then
            integrate = BoltzmannN(mass0, T)*deg
         end if
         return
      end if

      mmin = hadron(id)%minmass
      mmax = massMax
      if (mmax < mmin) return ! integral is zero

      ymax = 2*atan(2*(mmax-mass0) / gamma0)
      ymin = 2*atan(2*(mmin-mass0) / gamma0)

      nm = max(int((ymax-ymin)/dy0),1)
      dy  = (ymax-ymin)/float(nm)

      sum = 0.0
      do im=1,nm
         y = ymin+(float(im)-0.5)*dy
         m = 0.5*tan(0.5*y)*gamma0 + mass0
         m = min(max(m,mmin),mmax)
         gamma = fullWidthMeson(id, m)
         spectral = 2./pi * m**2 * gamma / ((m**2-mass0**2)**2+m**2*gamma**2)
         intfac = gamma0 / ((m-mass0)**2+gamma0**2/4.)
         sum = sum + spectral/intfac * BoltzmannN(m, T)
       end do
       integrate = sum * dy * deg

    end function integrate

    real function BoltzmannN(mass, T)
      use constants, only: pi
      use besselK, only: BesselK2

      real, intent(in) :: mass, T

      real :: betam, res
      real, parameter :: fak = 4*pi/(2*pi*0.197)**3

      betam  = mass/T
      res = BesselK2(betam)
      BoltzmannN = fak*res * mass**2 * T

      return
    end function BoltzmannN

  end subroutine initInput

  !*******************************************************************
  !****s* initPionBox/initializePionBox
  ! NAME
  ! subroutine initializePionBox(part)
  ! PURPOSE
  ! Initialize nucleons in a box
  !*******************************************************************
  subroutine initializePionBox(part)

    use particleDefinition
    use IdTable, only: pion
    use constants, only: mPi
    use densityModule, only: gridsize,get_densitySwitch
    use output, only: Write_InitStatus
    use callstack, only: traceBack
    use insertion, only: GarbageCollection
    use collisionNumbering, only: real_firstnumbering


    type(particle), dimension(:,:),intent(inOut) :: part

    integer :: dummy, iEns, iPart, nEns
    integer, dimension(0:3) :: nCharge

    call Write_InitStatus('box of pions',0)
    call initInput

    dummy = get_densitySwitch() ! force density module to be initialized
    nEns = size(part(:,1))

    write(*,'(A,3F9.3)') ' Gridsize   =',gridsize(:)
    write(*,*) 'Size of box= (8*Gridsize) = ',8.*gridsize(1)*gridsize(2)*gridsize(3)
    write(*,*) 'Number Ensembles,              nEns =', nEns
    write(*,*) 'Number Particles per Ensemble, nPart=', size(part(1,:))

    if (thermalInit) then
       call doThermal
    else
       call doNonThermal
    end if

    call GarbageCollection(part)
    call Write_InitStatus('box of pions',1)

  contains
    !*******************************************************************
    !****is* initializePionBox/doThermal
    ! NAME
    ! subroutine doThermal
    ! PURPOSE
    ! Do the thermal init
    !*******************************************************************
    subroutine doThermal

      use randomMaxwell, only: initMaxwell
      use AssignMassMC, only: AssignMass_1
      use mediumDefinition
      use particleProperties, only: hadron

      integer :: nPart0, nPart
      integer :: iMes
      real, dimension(0:3) :: pSum=0
      real :: mass
      real, dimension(0:3) :: momLRF = (/0,0,0,0/)

      nPart0 = 0
      pSum = 0

      do iMes=101,122
!         if ((iMes>=110).and.(iMes<=121)) cycle
!         if (densMes(iMes) <= 1e-8) cycle

         nPart = NINT(8.*gridsize(1)*gridsize(2)*gridsize(3)*densMes(iMes))
         write(*,'(" Number of ",i3," per ensemble:         ",i7,1P,e13.4)') &
              iMes,nPart,densMes(iMes)

         if (nPart < 1) cycle

         if (nPart0+nPart > size(part(1,:))) then
            call traceback('particle vector too small!')
         end if

         if (hadron(iMes)%width < 1e-03) then
            call initMaxwell(hadron(iMes)%mass, temp(iMes))
         end if

         do iEns=1,nEns

            select case(hadron(iMes)%isoSpinTimes2+1)
            case (1)
               call prepareChooseCharge1(nPart) ! only neutral
               
            case (2)
               if (hadron(iMes)%charm.ne.0) then
                  call traceback('charm not to be initialized!')
               end if
               if (hadron(iMes)%strangeness > 0) then
                  call prepareChooseCharge2plus(nPart) ! 0,+1
               else
                  call prepareChooseCharge2minus(nPart) ! -1,0
               end if
               
            case (3)
               call prepareChooseCharge3(nPart) ! -1,0,+1
               
            end select


            do iPart=nPart0+1,nPart0+nPart
               call setToDefault(part(iEns,iPart))
               call setNumber(part(iEns,iPart)) ! give it a unique number

               part(iEns,iPart)%event = real_firstnumbering()

               part(iEns,iPart)%ID = iMes
               if (hadron(iMes)%width < 1e-03) then
                  part(iEns,iPart)%mass = hadron(iMes)%mass
               else
                  call AssignMass_1(iMes, 5.0, momLRF,vacuum, mass)
                  part(iEns,iPart)%mass = mass
                  call initMaxwell(mass, temp(iMes)) ! re-init with new mass
               end if

               call ChooseCharge
               call ChoosePosition
               call ChooseMomentumTherm

               pSum = pSum + part(iEns,iPart)%momentum(0:3)

            end do ! iPart
         end do ! iEns

         nPart0 = nPart0+nPart
      end do ! iMes

      write(*,*)
      write(*,'(" Number of ALL per ensemble:         ",i7,1P,e13.4)') &
           nPart0,sum(densMes)
      write(*,'(" Size          per ensemble:         ",i7)') &
           size(part(1,:))

!      call traceback('not yet implemented!')
    end subroutine doThermal

    !*******************************************************************
    !****is* initializePionBox/doNonThermal
    ! NAME
    ! subroutine doNonThermal
    ! PURPOSE
    ! Do the non-thermal init
    !*******************************************************************
    subroutine doNonThermal
      use lorentzTrafo, only: lorentz

      integer :: nTot
      real, dimension(0:3) :: pSum=0
      real, dimension(1:3) :: betaVec=0

      nTot = NINT(8.*gridsize(1)*gridsize(2)*gridsize(3)*nDens)
      write(*,*) ' Number of pions per ensemble:        ', nTot

      if (nTot > size(part(1,:))) then
         call traceback('particle vector too small!')
      end if

      do iEns=1,nEns
         pSum = 0
         
         select case (ChargeSelection)
         case (0)
            call prepareChooseCharge1(nTot)
         case (1)
            call prepareChooseCharge2zero(nTot)
         case (2)
            call prepareChooseCharge3(nTot)
         end select

         do iPart=1,nTot
            call setToDefault(part(iEns,iPart))
            call setNumber(part(iEns,iPart)) ! give it a unique number

            part(iEns,iPart)%event = real_firstnumbering()
            part(iEns,iPart)%ID = pion
            part(iEns,iPart)%mass = mPi
            call ChooseCharge
            call ChoosePosition
            call ChooseMomentum

            pSum = pSum + part(iEns,iPart)%momentum(0:3)

         end do

!!$         ! correct for 'resting box':
!!$         pSum = pSum/nTot
!!$         do iPart=1,nTot
!!$            part(iEns,iPart)%momentum(1:3) = part(iEns,iPart)%momentum(1:3) - pSum
!!$            part(iEns,iPart)%momentum(0) = sqrt(mPi**2+sum(part(iEns,iPart)%momentum(1:3)**2))
!!$         end do

      end do


      if (BoostZ>0) then
         betaVec = (/ 0.0, 0.0, -BoostZ /)
         pSum = 0
         do iEns=1,nEns
            do iPart=1,nTot
               call lorentz(betaVec, part(iEns,iPart)%momentum )
               pSum = pSum + part(iEns,iPart)%momentum(0:3)
            end do
         end do
         pSum = pSum/(nEns*nTot)
         write(*,*) 'momentum/particle: ',pSum
      end if


    end subroutine doNonThermal

    !*******************************************************************
    !****is* initializePionBox/prepareChooseCharge1
    ! NAME
    ! subroutine prepareChooseCharge1(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-1,0,+1) = ( - , 100%, - )
    !*******************************************************************
    subroutine prepareChooseCharge1(nTot)
      integer, intent(in) :: nTot

      nCharge = (/ nTot, 0, nTot, 0 /)
      
    end subroutine prepareChooseCharge1

    !*******************************************************************
    !****is* initializePionBox/prepareChooseCharge2minus
    ! NAME
    ! subroutine prepareChooseCharge2minus(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-1,0,+1) = ( 50%, 50%, - )
    ! If the total number of particles is not even, the number of
    ! neutral particles is increased by one
    !*******************************************************************
    subroutine prepareChooseCharge2minus(nTot)
      integer, intent(in) :: nTot
      integer :: nTot2
      
      nTot2 = int(nTot/2)
      select case(nTot-2*nTot2)
      case (0)
         nCharge = (/ nTot, nTot2, nTot2, 0 /)
      case (1)
         nCharge = (/ nTot, nTot2, nTot2+1, 0 /)
      case DEFAULT
         write(*,*) nTot,nTot2
         call traceback('fishy nCharge')
      end select

    end subroutine prepareChooseCharge2minus
    
    !*******************************************************************
    !****is* initializePionBox/prepareChooseCharge2zero
    ! NAME
    ! subroutine prepareChooseCharge2zero(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-1,0,+1) = ( 50%, - , 50% )
    ! If the total number of particles is not even, the number of
    ! negative particles is increased by one, thus the symmetry is
    ! broken (but who cares?)
    !*******************************************************************
    subroutine prepareChooseCharge2zero(nTot)
      integer, intent(in) :: nTot
      integer :: nTot2
      
      nTot2 = int(nTot/2)
      select case(nTot-2*nTot2)
      case (0)
         nCharge = (/ nTot, nTot2, 0, nTot2/)
      case (1)
         ! this breaks the symmetry, because we have more pi- than pi+,
         ! but who cares? ;)
         nCharge = (/ nTot, nTot2+1, 0, nTot2 /)
      case DEFAULT
         write(*,*) nTot,nTot2
         call traceback('fishy nCharge')
      end select
      
    end subroutine prepareChooseCharge2zero

    !*******************************************************************
    !****is* initializePionBox/prepareChooseCharge2plus
    ! NAME
    ! subroutine prepareChooseCharge2plus(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-1,0,+1) = ( - , 50%, 50% )
    ! If the total number of particles is not even, the number of
    ! neutral particles is increased by one
    !*******************************************************************
    subroutine prepareChooseCharge2plus(nTot)
      integer, intent(in) :: nTot
      integer :: nTot2
      
      nTot2 = int(nTot/2)
      select case(nTot-2*nTot2)
      case (0)
         nCharge = (/ nTot, 0, nTot2, nTot2 /)
      case (1)
         nCharge = (/ nTot, 0, nTot2+1, nTot2 /)
      case DEFAULT
         write(*,*) nTot,nTot2
         call traceback('fishy nCharge')
      end select
      
    end subroutine prepareChooseCharge2plus

    !*******************************************************************
    !****is* initializePionBox/prepareChooseCharge3
    ! NAME
    ! subroutine prepareChooseCharge3(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-1,0,+1) = ( 33%, 33%, 33% )
    ! If the total number of particles is not dividable by 3, the
    ! number of neutral or the number of both charged states is
    ! increased by one in order to get the total number of particles,
    ! but always in charge neutral sum.
    !*******************************************************************
    subroutine prepareChooseCharge3(nTot)
      integer, intent(in) :: nTot
      integer :: nTot3
      nTot3 = int(nTot/3)
      select case(nTot-3*nTot3)
      case (0)
         nCharge = (/ nTot, nTot3, nTot3, nTot3 /)
      case (1)
         nCharge = (/ nTot, nTot3, nTot3+1, nTot3 /)
      case (2)
         nCharge = (/ nTot, nTot3+1, nTot3, nTot3+1 /)
      case DEFAULT
         write(*,*) nTot,nTot3
         call traceback('fishy nCharge')
      end select
    end subroutine prepareChooseCharge3

    !*******************************************************************
    !****is* initializePionBox/ChooseCharge
    ! NAME
    ! subroutine ChooseCharge
    ! PURPOSE
    ! Set the charge of the particle randomly according the numbers
    ! given in the array nCharge
    !*******************************************************************
    subroutine ChooseCharge
      use random, only: rn
      real :: r

      r = rn() * nCharge(0)
      if ( r < nCharge(1) ) then
         part(iEns,iPart)%charge = -1
         nCharge(1) = nCharge(1)-1
      else if ( r < nCharge(1)+nCharge(2) ) then
         part(iEns,iPart)%charge = 0
         nCharge(2) = nCharge(2)-1
      else
         part(iEns,iPart)%charge =  1
         nCharge(3) = nCharge(3)-1
      end if
      nCharge(0) = nCharge(0)-1
      
    end subroutine ChooseCharge
    
    !*******************************************************************
    !****is* initializePionBox/ChoosePosition
    ! NAME
    ! subroutine ChoosePosition
    ! PURPOSE
    ! Set the position of the particle randomly in the box
    !*******************************************************************
    subroutine ChoosePosition
      use random, only: rn

      part(iEns,iPart)%position(1)=(1.-2.*rn())*gridSize(1)
      part(iEns,iPart)%position(2)=(1.-2.*rn())*gridSize(2)
      part(iEns,iPart)%position(3)=(1.-2.*rn())*gridSize(3)

    end subroutine ChoosePosition
    
    !*******************************************************************
    !****is* initializePionBox/ChooseMomentum
    ! NAME
    ! subroutine ChooseMomentum
    ! PURPOSE
    ! Set the momentum of the particle within the non-thermal init.
    ! The direction is choosen isotropically, while the absolute value
    ! is fixed.
    !*******************************************************************
    subroutine ChooseMomentum
      use random, only: rnOmega

      part(iEns,iPart)%momentum(1:3) = rnOmega() * pInit
      part(iEns,iPart)%momentum(0) = sqrt(mPi**2+pInit**2)

    end subroutine ChooseMomentum

    !*******************************************************************
    !****is* initializePionBox/ChooseMomentumTherm
    ! NAME
    ! subroutine ChooseMomentumTherm
    ! PURPOSE
    ! Set the momentum of the particle within the thermal init.
    ! The direction is choosen isotropically. The absolute value of the
    ! momentum is choosen according a relativistic Maxwell distribution.
    !*******************************************************************
    subroutine ChooseMomentumTherm
      use random, only: rnOmega
      use randomMaxwell, only: rnMaxwell

      real :: p, m
      p = rnMaxwell()
      m = part(iEns,iPart)%mass
      part(iEns,iPart)%momentum(1:3) = rnOmega() * p
      part(iEns,iPart)%momentum(0) = sqrt(m**2+p**2)
    end subroutine ChooseMomentumTherm

  end subroutine initializePionBox

end module initPionBox
