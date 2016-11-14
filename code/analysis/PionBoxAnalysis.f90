!*****************************************************************************
!****m* /PionBoxAnalysis
! NAME
! module BoxAnalysis
! PURPOSE
!*****************************************************************************
module PionBoxAnalysis

  use histf90
  use histMPf90
  use histMC

  implicit none

  Private

  PUBLIC :: DoPionBoxAnalysisTime

  !***************************************************************************
  !****t* PionBoxAnalysis/tTmunuNmu
  ! NAME
  ! type tTmunuNmu
  ! PURPOSE
  ! This stores Tmunu, Nmu and Jmu:
  ! * Tmunu = 1/(n_ens*V_box) sum_i (pmu pnu)/p0
  ! * Nmu = 1/(n_ens*V_box) sum_i pmu/p0
  ! * Jmu = 1/(n_ens*V_box) sum_i pmu/p0 * q_i
  !
  ! SOURCE
  !
  type tTmunuNmu
     real, dimension(0:9) :: Tmunu = 0. ! the Tmunu tensor
     real, dimension(0:3) :: Nmu   = 0. ! the particle flow
     real, dimension(0:3) :: Jmu   = 0. ! the electrical charge flow
  end type tTmunuNmu
  !***************************************************************************

  type(histogram), save :: hMassRho
  type(histogramMC), save :: hMCMomPion
  type(histogramMP), save :: hMP_pSet2, hMP_pSet4,hMP_ESet2, hMP_ESet4
  real, dimension(:,:), allocatable, save :: arrMultSet2, arrMultSet4
  type(tTmunuNmu), dimension(:), allocatable, save :: arrTmunuNmu
  type(tTmunuNmu), dimension(2), save :: arrTmunuNmu_hadr


  logical, save :: initFlag=.true.
  logical, parameter :: do_Tmunu_pirho=.false.

  !***************************************************************************
  !****g* PionBoxAnalysis/do_Tmunu
  ! SOURCE
  logical,save :: do_Tmunu=.false.
  ! PURPOSE
  ! Switch for Tmunu output. default: Only one file for all ensemble!
  ! you may change this with the flag perEnsemble_Tmunu
  !***************************************************************************

  !***************************************************************************
  !****g* PionBoxAnalysis/perEnsemble_Tmunu
  ! SOURCE
  logical,save :: perEnsemble_Tmunu=.false.
  ! PURPOSE
  ! Switch for Tmunu output. One file per ensemble!
  !
  ! NOTE
  ! this may slow down the execution dramatically, since huge output to the
  ! hard drive is induced.
  ! You may observe this, if e.g the cpu load drops permanently to 30%.
  ! Thus: switch it on, only if you want it!
  !***************************************************************************

  !***************************************************************************
  !****g* PionBoxAnalysis/do_P
  ! SOURCE
  logical,save :: do_P=.false.
  ! PURPOSE
  ! Switch for dN/p^2 dp output
  !***************************************************************************

  !***************************************************************************
  !****g* PionBoxAnalysis/do_velrel
  ! SOURCE
  logical,save :: do_velrel=.false.
  ! PURPOSE
  ! Switch for calculating velrel
  !***************************************************************************

contains
  !***************************************************************************
  !****s* PionBoxAnalysis/DoBoxAnalysisTime
  ! NAME
  ! subroutine DoPionBoxAnalysisTime(realPart,timestep)
  ! PURPOSE
  !***************************************************************************
  subroutine DoPionBoxAnalysisTime(realPart,timestep)

    use CallStack, only: TRACEBACK
    use densityModule, only: gridsize
    use history, only: history_getParents
    use output, only: Write_InitStatus, intToChar4, WriteParticleVector
    use particleDefinition

    type(particle),dimension(:,:),intent(in), target :: realPart
    integer, intent(in) :: timestep

    integer, save :: nHist
    real, save :: boxVol ! the volume of the box in fm^3

    integer, save :: nEns,nPart
    real, save :: mulfak

    integer :: i,j, iID, iCh
    real :: mom0,mom, mass
    type(particle), POINTER :: pPart
    integer :: parents(1:3)

    type(tTmunuNmu) :: TmunuNmu

    ! string constants may be broken over multiple continuation lines:
    character(len=*), parameter :: headTmunu = "# 1: timestep &
         &2: T00 3: T11 4: T22 5: T33 &
         &6: T01 7: T02 8: T03 &
         &9: T21 10: T31 11: T32 &
         &12: N0 13: N1 14: N2 15: N3 &
         &16: J0 17: J1 18: J2 19: J3"

    if (initFlag) then
       call Write_InitStatus('PionBoxAnalysis',0)
       initFlag=.false.

       call readInput

       if ((do_Tmunu) .and. (perEnsemble_Tmunu) .and. (nEns > 9999)) then
          call TRACEBACK("PionBoxAnalysis: Tmunu not prepared for more than 9999 ensembles.")
       end if

       nEns  = size(realPart,dim=1)
       nPart = size(realPart,dim=2)
       boxVol = 8.*gridsize(1)*gridsize(2)*gridsize(3)
       mulfak = 1.0/(nEns*boxVol)

       !----- mass distribution: -----
       call CreateHist(hMassRho, "mass(rho)", 0.0, 2.5, 0.01)
       call CreateHistMC(hMCMomPion, "momentum(pion)", 0.0, 2.5, 0.02, 6)
       hMCMomPion%yDesc(1:6) = (/ "original  ",  &
            "rho       ", "sigma     ", "other dec ", &
            "pi pi     ", "other coll" /)

       !----- multiplicities: -----

!       call CreateHistMP(hMP_pSet2, "dN/p^2 dp", 0.0, 2.5, 0.05, 2)
!       call CreateHistMP(hMP_pSet4, "dN/p^2 dp", 0.0, 2.5, 0.05, 4)
       call CreateHistMP(hMP_pSet2, "dN/p^2 dp", 0.0, 2.5, 0.02, 2)
       call CreateHistMP(hMP_pSet4, "dN/p^2 dp", 0.0, 2.5, 0.02, 4)
       call CreateHistMP(hMP_ESet2, "dN/pE dE", 0.0, 2.5, 0.02, 2)
       call CreateHistMP(hMP_ESet4, "dN/pE dE", 0.0, 2.5, 0.02, 4)

       nHist = Map2HistMP_getN(2)
       allocate( arrMultSet2(0:nHist, 2) )
       nHist = Map2HistMP_getN(4)
       allocate( arrMultSet4(0:nHist, 2) )

       open(123,file="PionBoxAnalysis_Mult_Set2.dat", status="unknown")
       call WriteHistMP_Names(hMP_pSet2,123)
       close(123)

       open(123,file="PionBoxAnalysis_Mult_Set4.dat", status="unknown")
       call WriteHistMP_Names(hMP_pSet4,123)
       close(123)

       !----- hydro tensors: -----
       if (do_Tmunu) then
          allocate( arrTmunuNmu(nEns) )

          open(123,file="PionBoxAnalysis_Tmunu.dat", status="unknown")
          write(123,'(A)') headTmunu
          close(123)
          if (do_Tmunu_pirho) then
             open(123,file="PionBoxAnalysis_Tmunu.pion.dat", status="unknown")
             write(123,'(A)') headTmunu
             close(123)
             open(123,file="PionBoxAnalysis_Tmunu.rho.dat", status="unknown")
             write(123,'(A)') headTmunu
             close(123)
          end if

          if (perEnsemble_Tmunu) then
             do i=1,nEns
                open(123,file="PionBoxAnalysis_Tmunu."//intToChar4(i)//".dat", status="unknown")
                write(123,'(A)') headTmunu
                close(123)
             end do
          end if
       end if

       call Write_InitStatus('PionBoxAnalysis',1)
    end if

    call ClearHistMP(hMP_pSet2)
    call ClearHistMP(hMP_pSet4)
    call ClearHistMP(hMP_ESet2)
    call ClearHistMP(hMP_ESet4)
    arrMultSet2 = 0.0
    arrMultSet4 = 0.0

    call ClearHistMC(hMCMomPion)
    call ClearHist(hMassRho)

    if (do_Tmunu) then
       arrTmunuNmu = TmunuNmu ! set all values to zero
       arrTmunuNmu_hadr = TmunuNmu ! set all values to zero
    end if

    !
    ! accumulate data:
    !

    do i=1,nEns
       do j=1,nPart
          pPart => realPart(i,j)
          if(pPart%Id <  0) exit
          if(pPart%Id <= 0) cycle

          mom = absMom(pPart)
          mom0 = pPart%momentum(0)

          select case(pPart%ID)
          case (101)
             parents = history_getParents(pPart%history)
             if (parents(2) == 0) then
                select case(parents(1))
                case (0)
                   iCh = 1
                case (103)
                   iCh = 2
                case (104)
                   iCh = 3
                case default
                   iCh = 4
                end select
             else
                if (parents(1)==101 .and. parents(2)==101) then
                   iCh = 5
                else
                   iCh = 6
                end if
             end if
             call AddHistMC(hMCMomPion, mom, iCh, 1.0/(mom**2))
          case (103)
             mass = sqrtS(pPart)
             call AddHist(hMassRho, mass, 1.0)
          end select

          call AddHistMP(hMP_pSet2, pPart, mom, 1.0/(mom**2), 1.0)
          call AddHistMP(hMP_pSet4, pPart, mom, 1.0/(mom**2), 1.0)

          call AddHistMP(hMP_ESet2, pPart, mom0, 1.0/(mom0*mom), 1.0)
          call AddHistMP(hMP_ESet4, pPart, mom0, 1.0/(mom0*mom), 1.0)

          arrMultSet2(0,1) = arrMultSet2(0,1) + 1.0
          arrMultSet4(0,1) = arrMultSet4(0,1) + 1.0

          iID = Map2HistMP_ID(pPart%ID,pPart%charge,pPart%antiparticle, 2)
          if (iID>0) then
             arrMultSet2(iID,1) = arrMultSet2(iID,1) + 1.0
          endif
          iID = Map2HistMP_ID(pPart%ID,pPart%charge,pPart%antiparticle, 4)
          if (iID>0) then
             arrMultSet4(iID,1) = arrMultSet4(iID,1) + 1.0
          endif

          ! fill Tmunu and Jmu:
          if (do_Tmunu) then
             call fillTmunu(TmunuNmu, pPart)
             if (do_Tmunu_pirho) then
                select case(pPart%ID)
                case (101)
                   call fillTmunu(arrTmunuNmu_hadr(1), pPart)
                case (103)
                   call fillTmunu(arrTmunuNmu_hadr(2), pPart)
                end select
             end if
             if (perEnsemble_Tmunu) call fillTmunu(arrTmunuNmu(i), pPart)
          end if

       end do
    end do

    !
    ! produce output:
    !

    if (do_P) then
       if (mod(timestep,5)==1) then
          call WriteHistMP(hMP_pSet2, file='p_Set2_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak, iColumn=1)
          call WriteHistMP(hMP_pSet4, file='p_Set4_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak, iColumn=1)
          call WriteHistMP(hMP_ESet2, file='E_Set2_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak, iColumn=1)
          call WriteHistMP(hMP_ESet4, file='E_Set4_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak, iColumn=1)

          call WriteHist(hMassRho, file='massRho_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak)
          !          call WriteParticleVector('parts_'//intTochar(timestep),realPart)
          !          call WriteHistMC(hMCMomPion, file='MomPion_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak)
       end if
    end if

    open(123,file="PionBoxAnalysis_Mult_Set2.dat",status="old",position='append')
    write(123,'(i11,1P,100E12.4,0P)') timestep, &
         & arrMultSet2(1:,1)*mulfak,arrMultSet2(0,1)*mulfak
    close(123)

    open(123,file="PionBoxAnalysis_Mult_Set4.dat",status="old",position='append')
    write(123,'(i11,1P,100E12.4,0P)') timestep, &
         & arrMultSet4(1:,1)*mulfak,arrMultSet4(0,1)*mulfak
    close(123)


    if (do_Tmunu) then

       open(123,file="PionBoxAnalysis_Tmunu.dat",status="old",position='append')
       write(123,'(i11,1P,100E14.6,0P)') timestep, &
            & TmunuNmu%Tmunu(:)*mulfak, &
            & TmunuNmu%Nmu(:)*mulfak, &
            & TmunuNmu%Jmu(:)*mulfak
       close(123)


       if (do_Tmunu_pirho) then
          open(123,file="PionBoxAnalysis_Tmunu.pion.dat",status="old",position='append')
          write(123,'(i11,1P,100E14.6,0P)') timestep, &
               & arrTmunuNmu_hadr(1)%Tmunu(:)*mulfak, &
               & arrTmunuNmu_hadr(1)%Nmu(:)*mulfak, &
               & arrTmunuNmu_hadr(1)%Jmu(:)*mulfak
          close(123)

          open(123,file="PionBoxAnalysis_Tmunu.rho.dat",status="old",position='append')
          write(123,'(i11,1P,100E14.6,0P)') timestep, &
               & arrTmunuNmu_hadr(2)%Tmunu(:)*mulfak, &
               & arrTmunuNmu_hadr(2)%Nmu(:)*mulfak, &
               & arrTmunuNmu_hadr(2)%Jmu(:)*mulfak
          close(123)
       end if

       if (perEnsemble_Tmunu) then
          ! since we print all information for every ensemble, we must not divide by nEns here
          do i=1,nEns
             open(123,file="PionBoxAnalysis_Tmunu."//intToChar4(i)//".dat", status="old",position='append')
             write(123,'(i11,1P,100E12.4,0P)') timestep, &
                  & arrTmunuNmu(i)%Tmunu(:)/boxVol, &
                  & arrTmunuNmu(i)%Nmu(:)/boxVol, &
                  & arrTmunuNmu(i)%Jmu(:)/boxVol
             close(123)
          end do
       end if
    end if

    if (do_velrel) then
       ! calculate average of velrel
       call CalcAverageVelRel(realPart,timestep,nEns,nPart)
    end if

  end subroutine DoPionBoxAnalysisTime

  !***************************************************************************
  !****s* PionBoxAnalysis/CalcAverageVelRel
  ! NAME
  ! subroutine CalcAverageVelRel(realPart,timestep,nEns,nPart)
  ! PURPOSE
  ! calculate the average v_rel of all particles with each other. Due to
  ! speed reasons, it may be a good idea to correlate only particles in the
  ! same ensemble, but this s only a approximation.
  !***************************************************************************
  subroutine CalcAverageVelRel(realPart,timestep,nEns,nPart)
    use particleDefinition
    use lorentzTrafo, only: eval_sigmaBoost

    type(particle),dimension(:,:),intent(in), target :: realPart
    integer, intent(in) :: timestep
    integer, intent(in) :: nEns,nPart

    integer :: iEns1,iEns2, iPart1,iPart2
    type(particle), POINTER :: pPart1, pPart2
    real :: sum0,sum1,velrel,m1,m2,s
    real :: ptot(0:3)
!    real :: vrel
!    real, dimension(1:3) :: vrel_vector


    ! due to speed reasons, I only correlate particles in the same ensemble

!    write(*,*) 'calculating velrel....'
    sum0 = 0.0
    sum1 = 0.0
    do iEns1=1,nEns
       do iPart1=1,nPart
          pPart1 => realPart(iEns1,iPart1)
          if(pPart1%Id <  0) exit
          if(pPart1%Id <= 0) cycle
          m1 = pPart1%mass**2

          iEns2 = iEns1
          do iPart2=iPart1+1,nPart

             pPart2 => realPart(iEns2,iPart2)
             if(pPart2%Id <  0) exit
             if(pPart2%Id <= 0) cycle
             m2 = pPart2%mass**2

             ptot = pPart1%momentum + pPart2%momentum
             s = ptot(0)**2-sum(ptot(1:3)**2)
             velrel = sqrt( max(0.0,(s-m1-m2)**2/4-m1*m2) )/( pPart1%momentum(0)*pPart2%momentum(0) )

!             vrel_vector=pPart1%velocity-pPart2%velocity
!             vrel = sqrt(Dot_product(vrel_vector,vrel_vector))

!             write(*,*) velrel,vrel,velrel/vrel,eval_sigmaBoost(pPart1%momentum,pPart2%momentum)

!             write(*,*) (velrel-vrel*eval_sigmaBoost(pPart1%momentum,pPart2%momentum))/velrel

             sum0 = sum0 + 1
             sum1 = sum1 + velrel

          end do
       end do
    end do

    write(*,'(A,f8.5,f15.0,i13)') 'Average velrel: ',sum1/sum0, sum0, timestep

  end subroutine CalcAverageVelRel

  !***************************************************************************
  !****s* PionBoxAnalysis/FillTmunu
  ! NAME
  ! subroutine FillTmunu(TmunuNmu, pPart)
  ! PURPOSE
  ! add the momentum of the given particle to the entries of Tmunu, Nmu, Jmu
  !***************************************************************************
  subroutine FillTmunu(TmunuNmu, pPart)
    use particleDefinition

    type(tTmunuNmu), intent(inOut) :: TmunuNmu
    type(particle), POINTER, intent(in) :: pPart

    real :: oneE
    oneE = 1.0/pPart%momentum(0)

    TmunuNmu%Tmunu(0) = TmunuNmu%Tmunu(0) + pPart%momentum(0) ! T00
    TmunuNmu%Tmunu(1:3) = TmunuNmu%Tmunu(1:3) + pPart%momentum(1:3)**2*oneE ! T11,T22,T33
    TmunuNmu%Tmunu(4:6) = TmunuNmu%Tmunu(4:6) + pPart%momentum(1:3) ! T01,T02,T03
    TmunuNmu%Tmunu(7)   = TmunuNmu%Tmunu(7)   + pPart%momentum(2)*pPart%momentum(1)*oneE ! T21
    TmunuNmu%Tmunu(8)   = TmunuNmu%Tmunu(8)   + pPart%momentum(3)*pPart%momentum(1)*oneE ! T31
    TmunuNmu%Tmunu(9)   = TmunuNmu%Tmunu(9)   + pPart%momentum(3)*pPart%momentum(2)*oneE ! T32

    TmunuNmu%Nmu(0:3) = TmunuNmu%Nmu(0:3) + pPart%momentum(0:3) * oneE
    if (pPart%charge .ne. 0) then
       TmunuNmu%Jmu(0:3) = TmunuNmu%Jmu(0:3) + pPart%momentum(0:3) * pPart%charge * oneE
    endif
  end subroutine FillTmunu

  !***************************************************************************
  !****s* PionBoxAnalysis/readInput
  ! NAME
  ! subroutine readInput
  !
  ! PURPOSE
  ! Reads input in jobcard out of namelist "pionBoxAnalysis"
  !***************************************************************************
  subroutine readInput

    use output, only: Write_ReadingInput

    !*************************************************************************
    !****n* PionBoxAnalysis/pionBoxAnalysis
    ! NAME
    ! NAMELIST pionBoxAnalysis
    ! PURPOSE
    ! Includes the switches:
    ! * do_Tmunu
    ! * perEnsemble_Tmunu
    ! * do_P
    ! * do_velrel
    !*************************************************************************
    NAMELIST /pionBoxAnalysis/ &
         do_Tmunu, do_P, do_velrel, perEnsemble_Tmunu
    integer :: ios

    call Write_ReadingInput('pionBoxAnalysis',0)
    rewind(5)
    read(5,nml=pionBoxAnalysis,IOSTAT=ios)
    call Write_ReadingInput('pionBoxAnalysis',0,ios)

    write(*,*) '  do Tmunu: ',do_Tmunu,'   perEnsemble: ',perEnsemble_Tmunu
    write(*,*) '  do P:     ',do_P
    write(*,*) '  do velrel:',do_velrel

    call Write_ReadingInput('pionBoxAnalysis',1)
  end subroutine readInput

end module PionBoxAnalysis
