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

  Public :: initializePionBox

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
    ! * nDens
    !*****************************************************************
    NAMELIST /PionBox/ nDens,ChargeSelection,pInit,BoostZ

    integer :: ios
    character(20), dimension(0:2), parameter :: cSel = (/ &
         'only pi0        ',&
         '50% pi+, 50% pi-',&
         '33% for +,0,-   ' /)
    real :: eDens,beta(1:2)

    call Write_ReadingInput('PionBox',0)
    rewind(5)
    read(5,nml=PionBox,iostat=ios)
    call Write_ReadingInput('PionBox',0,ios)

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

    call Write_ReadingInput('PionBox',1)

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
    use lorentzTrafo, only: lorentz
    
    type(particle), dimension(:,:),intent(inOut) :: part

    integer :: numberPions, dummy, iEns,iPart, nEns
    real, dimension(0:3) :: pSum=0
    real, dimension(1:3) :: betaVec=0
    integer, dimension(0:3) :: nCharge

    call Write_InitStatus('box of pions',0)
    call initInput

    dummy = get_densitySwitch() ! force density module to be initialized

    numberPions = NINT(8.*gridsize(1)*gridsize(2)*gridsize(3)*nDens)
    nEns = size(part(:,1))

    write(*,'(A,3F9.3)') '  Gridsize   =',gridsize(:)
    write(*,*) ' Size of box= (8*Gridsize) = ',8.*gridsize(1)*gridsize(2)*gridsize(3)
    write(*,*) ' Number Ensembles,              nEns =', nEns
    write(*,*) ' Number Particles per Ensemble, nPart=', size(part(1,:))
    write(*,*) ' Number of pions per ensemble:        ', numberPions

    if (numberPions > size(part(1,:))) then
       call traceback('particle vector too small!')
    end if

    do iEns=1,nEns
       pSum = 0
       call prepareChooseCharge
       do iPart=1,numberPions
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

       pSum = pSum/numberpions

       ! correct for 'resting box'

!!$       do iPart=1,numberPions
!!$          part(iEns,iPart)%momentum(1:3) = part(iEns,iPart)%momentum(1:3) - pSum
!!$          part(iEns,iPart)%momentum(0) = sqrt(mPi**2+sum(part(iEns,iPart)%momentum(1:3)**2))
!!$       end do

    end do

    

    if (BoostZ>0) then
       betaVec = (/ 0.0, 0.0, -BoostZ /)
       pSum = 0
       do iEns=1,nEns
          do iPart=1,numberPions
             call lorentz(betaVec, part(iEns,iPart)%momentum )
             pSum = pSum + part(iEns,iPart)%momentum(0:3)
          end do
       end do
       pSum = pSum/(nEns*numberpions)
       write(*,*) 'momentum/particle: ',pSum
    end if
    
    call GarbageCollection(part)
    call Write_InitStatus('box of pions',1)

  contains

    subroutine prepareChooseCharge

      nCharge(0) = numberPions

      select case (ChargeSelection)
      case (1)
         nCharge(1:3) = int(numberPions/2)
         nCharge(2) = 0
         select case(numberPions - sum(nCharge(1:3)))
         case (0)
            ! NOP
         case (1)
            ! this breaks the symmetry, because we have more pi- than pi+,
            ! but who cares? ;)
            nCharge(1) = nCharge(1)+1 ! pi- ++
         case DEFAULT
            write(*,*) nCharge
            call traceback('fishy nCharge (2)')
         end select

      case (2)
         nCharge(1:3) = int(numberPions/3)
         select case(numberPions - sum(nCharge(1:3)))
         case (0)
            ! NOP
         case (1)
            nCharge(2) = nCharge(2)+1 ! pi0 ++
         case (2)
            nCharge(1) = nCharge(1)+1 ! pi- ++
            nCharge(3) = nCharge(3)+1 ! pi+ ++
         case DEFAULT
            write(*,*) nCharge
            call traceback('fishy nCharge (3)')
         end select
      
      end select

    end subroutine prepareChooseCharge

    subroutine ChooseCharge
      use random, only: rn
      real :: r

      select case (ChargeSelection)
      case (0)

         part(iEns,iPart)%charge = 0

      case (1:2)

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

      case DEFAULT
         write(*,*) 'wrong ChargeSelection = ',ChargeSelection
         call traceback('correct input')
      end select

    end subroutine ChooseCharge

    subroutine ChoosePosition
      use random, only: rn

      part(iEns,iPart)%position(1)=(1.-2.*rn())*gridSize(1)
      part(iEns,iPart)%position(2)=(1.-2.*rn())*gridSize(2)
      part(iEns,iPart)%position(3)=(1.-2.*rn())*gridSize(3)

    end subroutine ChoosePosition

    subroutine ChooseMomentum
      use random, only: rnOmega
      
      part(iEns,iPart)%momentum(1:3) = rnOmega() * pInit
      part(iEns,iPart)%momentum(0) = sqrt(mPi**2+pInit**2)

    end subroutine ChooseMomentum

  end subroutine initializePionBox

end module initPionBox
