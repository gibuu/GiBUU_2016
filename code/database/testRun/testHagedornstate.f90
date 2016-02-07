program testHagedornstate
  use inputGeneral
  use particleProperties, only: initParticleProperties, hadron, isCharmed
  use IdTable, only: isHadron, isBaryon
  use clebschGordan
  use output

  implicit none

!  integer, parameter :: qBmax=4, qSmax=12, qIx2max=12
!  integer, parameter :: qBmax=6, qSmax=18, qIx2max=18
  integer, parameter :: qBmax=7, qSmax=20, qIx2max=20

  real, dimension(-1:1,-3:3,0:3) :: Arr1 !B,S,I 
  real, dimension(-qBmax:qBmax,-qSmax:qSmax,0:qIx2max) :: ArrMinMass ! B,S,I 

  ! Database parameters:
  
  integer :: nPart=0
  real, dimension(1:400) :: ArrM, ArrW
  integer, dimension(1:400) :: arrID,arrB,arrS,arrI,arrJ
  character*18, dimension(1:400) :: ArrA
  
  ! Histogram parameters:
  
  real, parameter :: dM = 0.01
  integer, parameter :: nBin=1000
!  integer, parameter :: nBin=500
  real, dimension(1:nBin) :: HistAll,HistM,HistB,HistBSI,HistTheo

  real, dimension(-1:1,-3:3,0:3,1:nBin) :: ArrHadrons
  real, dimension(-qBmax:qBmax,-qSmax:qSmax,0:qIx2max,1:nBin,0:2) :: ArrBootstrap
  logical, dimension(-qBmax:qBmax,-qSmax:qSmax,0:qIx2max,0:2) :: ArrIsZero
  integer, dimension(-qBmax:qBmax,-qSmax:qSmax,0:qIx2max,0:2) :: ArrMinIM

!   real, dimension(nBin,nBin) :: arrPCM

  integer, parameter :: qB=0,qS=0,qIx2=0
!  integer, parameter :: qB=0,qS=0,qIx2=2
!  integer, parameter :: qB=0,qS=1,qIx2=1
!  integer, parameter :: qB=1,qS=-1,qIx2=0
!  integer, parameter :: qB=1,qS=-1,qIx2=2
!  integer, parameter :: qB=1,qS=0,qIx2=1
!  integer, parameter :: qB=1,qS=0,qIx2=3
!  integer, parameter :: qB=1,qS=-2,qIx2=1
!  integer, parameter :: qB=1,qS=-3,qIx2=0


!  call ReadFort
!  call DoFinalFolding
!  call PrintSomeGamma
!  stop

  call readInputGeneral
  call initParticleProperties

  call ReadFilePDG
  call PrintTable(12)

  call CalcMinMass
!  call CalcBootstrap
  call DoHamerFrautschi

!  call CalcTheoRho(qB,qS,qIx2)
!  call CalcExpRho(qB,qS,qIx2)

!  call PrintGamma(qB,qS,qIx2)
  

contains

  subroutine CalcTheoRho(qB,qS,qIx2)
    integer, intent(in) :: qB,qS,qIx2
    real :: m,mg
    integer :: iBin

    HistTheo = 0

    mg = funcMG(qB,qS,qIx2)

    do iBin=1,nBin
       m = iBin*dM
       HistTheo(iBin) = funcRho(m,mg)
    end do

  end subroutine CalcTheoRho


  subroutine CalcExpRho(qB,qS,qIx2)
    integer, intent(in) :: qB,qS,qIx2
    
    call ReadFilePDG
    call PrintTable(12)
    
    call FillHistMB
    call FillHistBSI(qB,qS,qIx2)
    call CalcTheoRho(qB,qS,qIx2)
    call WriteHist(13)
    call WriteHistBSI(15)
    call IntegrateHist
    call WriteHist(14)
    call WriteHistBSI(16)


    ArrW = 0.0
    call FillHistMB
    call FillHistBSI(qB,qS,qIx2)
    call CalcTheoRho(qB,qS,qIx2)
    call WriteHist(23)
    call WriteHistBSI(25)
    call IntegrateHist
    call WriteHist(24)
    call WriteHistBSI(26)

  end subroutine CalcExpRho


  subroutine IntegrateHist
    integer :: iBin
    do iBin=2,nBin
       HistM(iBin) = HistM(iBin) + HistM(iBin-1)
       HistB(iBin) = HistB(iBin) + HistB(iBin-1)
       HistAll(iBin) = HistAll(iBin) + HistAll(iBin-1)
       HistBSI(iBin) = HistBSI(iBin) + HistBSI(iBin-1)
       HistTheo(iBin) = HistTheo(iBin) + HistTheo(iBin-1)
    end do
    HistM = HistM*dM
    HistB = HistB*dM
    HistAll = HistAll*dM
    HistBSI = HistBSI*dM
    HistTheo = HistTheo*dM
    
  end subroutine IntegrateHist

  subroutine WriteHist(iFile)
    integer, intent(in) :: iFile
    integer :: iBin
    real :: m
    
    do iBin=1,nBin
       m = iBin*dM
       write(iFile,'(f12.3,1p,99e13.4,0p)') m,HistAll(iBin),HistM(iBin),HistB(iBin)
    end do
  end subroutine WriteHist

  subroutine WriteHistBSI(iFile)
    integer, intent(in) :: iFile
    integer :: iBin
    real :: m
    
    do iBin=1,nBin
       m = iBin*dM
       write(iFile,'(f12.3,1p,99e13.4,0p)') m,HistBSI(iBin),HistTheo(iBin)
    end do
  end subroutine WriteHistBSI

  subroutine FillHistMB
    integer :: iBin,iPart
    real :: m

    HistAll = 0
    HistM = 0
    HistB = 0

    do iBin=1,nBin
       m = iBin*dM
       do iPart=1,129
          HistM(iBin) = HistM(iBin) &
               + BWcontrib(m,m+dM, arrM(iPart),arrW(iPart)) &
               * (arrI(iPart)+1)*(arrJ(iPart)+1)
       end do
       do iPart=130,nPart
          HistB(iBin) = HistB(iBin) &
               + BWcontrib(m,m+dM, arrM(iPart),arrW(iPart)) &
               * (arrI(iPart)+1)*(arrJ(iPart)+1)
       end do
    end do
    HistAll = HistM+HistB
  end subroutine FillHistMB

  subroutine FillHistBSI(qB,qS,qIx2)
    integer, intent(in) :: qB,qS,qIx2

    integer :: iBin,iPart
    real :: m

    HistBSI = 0

    do iBin=1,nBin
       m = iBin*dM
       do iPart=1,nPart
          if (arrB(iPart)/=qB) cycle
          if (arrS(iPart)/=qS) cycle
          if (arrI(iPart)/=qIx2) cycle

          HistBSI(iBin) = HistBSI(iBin) &
               + BWcontrib(m,m+dM, arrM(iPart),arrW(iPart)) &
               * (arrI(iPart)+1)*(arrJ(iPart)+1)
       end do
    end do
  end subroutine FillHistBSI

  subroutine PrintTable(iFile)
    integer, intent(in) :: iFile
    integer :: i

    write(iFile,'(A)')"#Nr          KF    Mass   Width   B   S  2*I 2*J   Name"
    write(iFile,'(65("-"))')
    do i=1,nPart
       write(iFile,'(i3,i12,2f8.3,4i4,"    ",A)') &
            i,arrID(i),arrM(i),arrW(i),arrB(i),arrS(i),arrI(i),arrJ(i),arrA(i)
    end do
  end subroutine PrintTable


  subroutine ReadFilePDG
    character*200 BUF
    integer :: i
    character :: A
    integer,dimension(4) :: ID
    real :: val,err1,err2
    
    integer :: Q1,Q2,Q3,nJ,nB,nS,nI
!     integer :: defB=1, defS=-1, defI=2
    integer :: iC, iiC
    integer :: ios

    nPart = 0

    open(123,file="mass_width_2011.mcd", status="OLD")
    
    do i=1,33
       read (123,'(A)') BUF
    end do
    
    ! allPDG: 364

    do i=1,364
       read (123,'(BN, A1, 4I8, 1X, E15.0, 2(1X, E8.0), 1X, A21)',iostat=ios ) &
            & A,ID,val,err1,err2,BUF

       if (ios/=0) then
          write(*,*) 'ios=',ios,'. stop reading.'
          exit
       end if

       if (ID(1)==311) cycle ! 
       if (ID(1)==130) ID(1) = 311

100    continue ! REDO !!!

       call SplitKF(ID(1), Q1,Q2,Q3,nJ)
       if (nJ==0) cycle
       if (Q3==0) cycle
       if (Q2==0) cycle

       if (ID(1)==211) cycle ! pi+ treated differently
       if (ID(1)==321) cycle ! 
       if (ID(1)==323) cycle ! 
       if (ID(1)==325) cycle ! 
       if (ID(1)==327) cycle ! 
       if (ID(1)==411) cycle ! 
       if (ID(1)==413) cycle ! 
       if (ID(1)==415) cycle ! 
       if (ID(1)==521) cycle ! 
       if (ID(1)==523) cycle ! 
       if (ID(1)==525) cycle ! 
       if (ID(1)==2212) cycle ! p
       if (ID(1)==3222) cycle ! Sigma+
       if (ID(1)==3112) cycle ! Sigma-
       if (ID(1)==3224) cycle ! Sigma+
       if (ID(1)==3114) cycle ! Sigma-
       if (ID(1)==3312) cycle ! Xi-
       if (ID(1)==3314) cycle ! Xi-
       if (ID(1)==4222) cycle ! Sigma(c)(2455)     ++
       if (ID(1)==4212) cycle ! Sigma(c)(2455)      +
       if (ID(1)==4224) cycle ! Sigma(c)(2455)     ++
       if (ID(1)==4214) cycle ! Sigma(c)(2455)      +
       if (ID(1)==4232) cycle ! Xi(c)               +
       if (ID(1)==4322) cycle ! Xi(c)'              +
       if (ID(1)==4234) cycle ! Xi(c)               +
       if (ID(1)==4324) cycle ! Xi(c)'              +
       if (ID(1)==5222) cycle ! Sigma(b)            +
       if (ID(1)==5224) cycle ! Sigma(b)            +
       if (ID(1)==5114) cycle ! Sigma(b)            +
       if (ID(1)==5132) cycle ! Xi(b)               -

       select case(A)
       case ('M')
          nPart = nPart+1
          arrID(nPart) = ID(1)
!          write(*,*) nPart,ID(1)
       case ('W')
          arrW(nPart) = val
          cycle
       case DEFAULT
          cycle
       end select
       
!!$     call SplitKF(ID(1), Q1,Q2,Q3,nJ)
!!$     write(*,*) i,A,ID(1),nJ, Q1,Q2,Q3
       
       call SplitKF(ID(1), Q1,Q2,Q3,nJ)
       
       if (nJ==0) cycle
       if (Q3==0) cycle
       if (Q2==0) cycle
       
       nB=0
       if (Q1/=0) nB=1
       
       nS=0
       if (Q1==3) nS=nS-1
       if (Q1==-3) nS=nS+1
       if (Q2==3) nS=nS-1
       if (Q2==-3) nS=nS+1
       if (Q3==3) nS=nS-1
       if (Q3==-3) nS=nS+1
       
       call CalcIsospin(nB,nS,BUF, nI)
       
       arrM(nPart) = val
       arrB(nPart) = nB
       arrS(nPart) = nS
       arrI(nPart) = nI
       arrJ(nPart) = nJ-1 ! not 2J+1, but only 2J
       
       do iC=1,18
          if (BUF(iC:iC)==' ') then
             iiC = iC
             exit
          end if
       end do
       arrA(nPart) = BUF(1:iiC)
       
       
       if ((nB==0).and.(ID(1) > 0).and.(Q3/=-Q2)) then
          ID(1)=-ID(1)
          goto 100 ! REDO!!!
       end if

    end do
    close(123)
    
  end subroutine ReadFilePDG

  subroutine SplitKF(KF, Q1, Q2, Q3, nJ)
    integer, intent(in):: KF
    integer, intent(out):: Q1, Q2, Q3, nJ

    integer :: aKF

    aKF = mod(abs(KF), 10000)
    Q1 = aKF/1000
    Q2 = (aKF-1000*Q1)/100
    Q3 = (aKF-1000*Q1-100*Q2)/10
    if (Q1==0) Q2=-Q2 ! for mesons turn topmost quark

    nJ = mod(aKF,10)

    if (KF<0) then
       Q1 = -Q1
       Q2 = -Q2
       Q3=-Q3
    end if

  end subroutine SplitKF

  subroutine CalcIsospin(B,S,Name, I)
    integer, intent(in) :: B,S
    character*(*), intent(in) :: Name
    integer, intent(out) :: I

    I=-1 ! times 2 !

    select case(B)
    case(0)
       select case(Name(1:2))
       case ('et','om','ph', 'f(','h(','J/','ch','ps','Up')
          I = 0
       case ('K ','K*','K(', 'D ','D*','D(', 'B ','B*','B(')
          I = 1 ! i.e. 1/2
       case ('pi','rh', 'a(', 'b(')
          I = 2 ! i.e. 2/2=1
       case DEFAULT
          write(*,*) Name
          stop
       end select
    case(1)
       select case(Name(1:1))
       case ('n','p','N','X')
          I = 1 ! i.e. 1/2
       case ('D')
          I = 3 ! i.e. 3/2
       case ('L','O')
          I = 0
       case ('S')
          I = 2 ! i.e. 2/2=1
       case DEFAULT
          write(*,*) Name
          stop
       end select

    end select
  end subroutine CalcIsospin

  real function BWcontrib(m1,m2, m0,gamma)
    real, intent(in) :: m1,m2      ! lower/upper bin boundary
    real, intent(in) :: m0,gamma ! BW parameters

    real :: x1=0, x2=0
    real, parameter :: x0 = 1/(2*atan( 3.0 ))

    BWcontrib = 0.0

    if (gamma > 1e-10) then
       if (m2<m0-3.0*gamma) return
       if (m1>m0+3.0*gamma) return
       x1 = atan2(2*(m1-m0),gamma) 
       x2 = atan2(2*(m2-m0),gamma)
       BWcontrib = x0 * (x2-x1) / (m2-m1)
    else
       if ((m1<=m0).and.(m0<=m2)) then
          BWcontrib = 1.0 / (m2-m1)
       end if
    end if
  end function BWcontrib




  


  real function funcRho(m,mg)
    real, intent(in) :: m,mg
    real :: funcf, mArg

    real, parameter :: mr = 0.5

    real, parameter :: norm = 0.25*1.5 ! PalDan
    real, parameter :: a = 2.82 ! PalDan
    real, parameter :: TH = 0.170 ! PalDan
    real, parameter :: mth = 2.0 ! PalDan
!    real, parameter :: mth = 2.2 ! PalDan
    integer :: n = 1 ! PalDan

!!$    real, parameter :: norm = 10*0.0362 ! Beitel
!!$    real, parameter :: a = 1.03 ! Beitel
!!$    real, parameter :: TH = 0.180 ! Beitel
!!$    real, parameter :: mth = 1.8 ! Beitel
!!$    integer :: n = 1 ! Beitel

    funcf = 1.0/(1.0+(abs((m-mg)/mth))**n)
    mArg = m-mg*funcf
    funcRho = norm*exp(mArg/TH)/(mArg**2+mr**2)**a
  end function funcRho


  real function funcMG(B,S,Ix2)
    integer, intent(in) :: B,S,Ix2 ! ATTENTION: Ix2 = I*2

    real, parameter :: aQ=0.387
    real, parameter :: aS=0.459
    integer, parameter :: iParam = 3

    funcMG = 0
    select case(iParam)
    case (1)
       funcMG = 0
    case (2)
       if (abs(3*B+S) > Ix2) then
          funcMG = aQ*abs(3*B+S)+aS*abs(S)
       else
          funcMG = aQ*(Ix2)+aS*abs(S)
       end if
    case (3)
       funcMG = ArrMinMass(B,S,Ix2)
    case (4)
!       funcMG = 0.400 ! B=0, S=0, I=0
!       funcMG = 0.138 ! B=0, S=0, I=1
!       funcMG = 0.496 ! B=0, S=1, I=1/2
!       funcMG = 1.116 ! B=1, S=-1, I=0
!       funcMG = 1.189 ! B=1, S=-1, I=1
       funcMG = 0.938 ! B=1, S=0, I=1/2
!       funcMG = 1.076 ! B=1, S=0, I=3/2
!       funcMG = 1.254 ! B=1, S=-2, I=1/2
!       funcMG = 1.672 ! B=1, S=-3, I=0
!       funcMG = 1.200
    end select
  end function funcMG

  real function funcR2(m)
    real, intent(in) :: m

    real, parameter :: r0 = 1.0 ! fm
    real, parameter :: md = 1.0 ! GeV

    funcR2 = (r0*(m/md)**(1./3.))**2
  end function funcR2

  subroutine PrintGamma(qB,qS,qIx2)
    integer, intent(in) :: qB,qS,qIx2

!    integer :: qB = 0, qS = 0, qIx2 = 0
!    integer :: qB = 0, qS = -1, qIx2 = 1
!    integer :: qB = 1, qS = 0, qIx2 = 1
!    integer :: qB = 1, qS = 1, qIx2 = 2

    integer :: ID1,ID2
    real :: m1,m2
    integer :: qB1,qB2
    integer :: qS1,qS2
    integer :: qIx21,qIx22
    integer :: iM
    real :: m,mg  !, h
    real, dimension(3) :: gamma


    mg = funcMG(qB,qS,qIx2)

    do iM=50,1000
!    do iM=500,500
       m=iM*0.01

       gamma = 0

       !#### GAMMA 1 ####

       do ID1=1,121
          if (.not.(isHadron(ID1))) cycle
          if (isCharmed(ID1)) cycle

          m1 = hadron(ID1)%mass
          qS1 = hadron(ID1)%strangeness
          qIx21 = hadron(ID1)%isoSpinTimes2
          qB1 = 0
          if (isBaryon(ID1)) qB1 = 1
          
          do ID2=ID1,121
             if (.not.(isHadron(ID2))) cycle
             if (isCharmed(ID2)) cycle

             m2 = hadron(ID2)%mass
             qS2 = hadron(ID2)%strangeness
             qIx22 = hadron(ID2)%isoSpinTimes2
             qB2 = 0
             if (isBaryon(ID2)) qB2 = 1
    
             gamma(1) = gamma(1) + GammaI(m,qB,qS,qIx2, ID1,m1,qB1,qS1,qIx21, ID2,m2,qB2,qS2,qIx22)

             if (isBaryon(ID2)) then
                gamma(1) = gamma(1) + GammaI(m,qB,qS,qIx2, ID1,m1,qB1,qS1,qIx21, ID2,m2,-qB2,-qS2,qIx22)
             end if

          end do
       end do

       !#### GAMMA 2 ####

       do ID1=1,121
          if (.not.(isHadron(ID1))) cycle
          if (isCharmed(ID1)) cycle

          m1 = hadron(ID1)%mass
          qS1 = hadron(ID1)%strangeness
          qIx21 = hadron(ID1)%isoSpinTimes2
          qB1 = 0
          if (isBaryon(ID1)) qB1 = 1
          
          gamma(2) = gamma(2) + GammaII(m,qB,qS,qIx2,ID1,m1,qB1,qS1,qIx21)
          if (isBaryon(ID1)) then
             gamma(2) = gamma(2) + GammaII(m,qB,qS,qIx2,ID1,m1,-qB1,-qS1,qIx21)
          end if

       end do

       !#### GAMMA 3 ####

       gamma(3) = GammaIII(m,qB,qS,qIx2)

       write(*,'(f7.3,1p,3(e13.4))') m,gamma
       write(101,'(f7.3,1p,3(e13.4))') m,gamma
    end do


  end subroutine PrintGamma

  real function GammaIII(m,qB,qS,qIx2)
    use constants, only : twopi,hbarc

    real, intent(in) :: m
    integer, intent(in) :: qB,qS,qIx2
    real :: pp, mg
!     integer :: qIx2Z1, qIx2Z2
    real :: sum,sum1,sum2

    real :: m1,m1min,m1max,dm1,mg1
    real :: m2,m2min,m2max,dm2,mg2
    integer :: qB1,qS1,qIx21
    integer :: qB2,qS2,qIx22
    integer :: qIx22min, qIx22max
    integer :: im1,im2

    GammaIII = 0.0
    sum = 0
    do qB1=-7,7
       qB2 = qB - qB1
       do qS1 = -7,7
          qS2 = qS - qS1
          do qIx21=0,14
             qIx22min = abs(qIx2-qIx21)
             qIx22max = qIx2+qIx21
             
             mg1 = FuncMG(qB1,qS1,qIx21)
             m1min = max(2.0,mg1)
             m1max = m

             if (m1min>m1max) cycle
             dm1 = (m1max-m1min)/100


             do qIx22=qIx22min,qIx22max ! ??? ,2 ???

                ! summing over all qM1,qM2 is just the orthonormality condition
                ! thus: clebsch = 1, if "MayBranch()==.true."
                if (.not. MayBranch(qIx21,qIx22,qIx2)) cycle

                
                mg2 = FuncMG(qB2,qS2,qIx22)
                m2min = max(2.0,mg2)

                if (m1min+m2min > m) cycle

                sum1 = 0
                do im1=1,100
                   m1 = m1min+(float(im1)-0.5)*dm1
                   m2max = m - m1

                   if (m2min>m2max) exit
                   dm2 = (m2max-m2min)/100

                   sum2 = 0
                   do im2=1,100
                      m2 = m2min+(float(im2)-0.5)*dm2
                      pp = pstar2(m,m1,m2)
                      if (pp<0) exit
                      sum2 = sum2 + pp*funcRho(m2,mg2)
                   end do

                   sum1 = sum1 + sum2*dm2*funcRho(m1,mg1)

                end do ! im1
                sum = sum + sum1*dm1

!                write(*,'(4i5,1p,e12.4,0p,2f7.3)') qB1,qS1,qIx21,qIx22,sum1*dm1,mg1,mg2



             end do ! qIx22
          end do ! qIx21
       end do ! qS1
    end do ! qB1

    mg = funcMG(qB,qS,qIx2)
    GammaIII = sum * funcR2(m) &
         /(twopi*hbarc**2*funcRho(m,mg))
  end function GammaIII

  real function GammaII(m,qB,qS,qIx2, ID1,m1,qB1,qS1,qIx21)
    use constants, only : twopi,hbarc

    real, intent(in) :: m,m1
    integer, intent(in) :: qB,qS,qIx2,id1
    integer, intent(in) :: qB1,qS1,qIx21
    real :: pp, mg, mg2
!     integer :: qIx2Z1, qIx2Z2
    real :: sum,sum1

    real :: m2,m2min,m2max,dm2
    integer :: qB2,qS2,qIx22
    integer :: qIx22min, qIx22max
    integer :: im2

    GammaII = 0.0

    qB2 = qB - qB1
    qS2 = qS - qS1
    
    qIx22min = abs(qIx2-qIx21)
    qIx22max = qIx2+qIx21

    sum1=0
    do qIx22=qIx22min,qIx22max ! ??? ,2 ???

       ! summing over all qM1,qM2 is just the orthonormality condition
       ! thus: clebsch = 1, if "MayBranch()==.true."
       if (.not. MayBranch(qIx21,qIx22,qIx2)) cycle

       mg2 = FuncMG(qB2,qS2,qIx22)
       m2min = max(2.0,mg2)
       m2max = m - m1

       if (m2min>m2max) cycle
       dm2 = (m2max-m2min)/100
       sum = 0
       do im2=1,100
          m2 = m2min+(float(im2)-0.5)*dm2
          pp = pstar2(m,m1,m2)
          if (pp<0) exit
          sum = sum + pp*funcRho(m2,mg2)
       end do
       sum1 = sum1 + sum*dm2

    end do
    mg = funcMG(qB,qS,qIx2)
    GammaII = (hadron(ID1)%spin*2+1) * sum1 * funcR2(m) &
         /(twopi*hbarc**2*funcRho(m,mg))

  end function GammaII


  real function GammaI(m,qB,qS,qIx2, ID1,m1,qB1,qS1,qIx21, ID2,m2,qB2,qS2,qIx22)
    use constants, only : twopi,hbarc

    real, intent(in) :: m,m1,m2
    integer, intent(in) :: qB,qS,qIx2
    integer, intent(in) :: ID1,qB1,qS1,qIx21
    integer, intent(in) :: ID2,qB2,qS2,qIx22
    real :: pp, mg

    GammaI = 0.0

    if (qB1+qB2 /= qB) return
    if (qS1+qS2 /= qS) return

    pp = pstar2(m,m1,m2)
    if (pp<0) return

    ! summing over all qM1,qM2 is just the orthonormality condition
    ! thus: clebsch = 1, if "MayBranch()==.true."
    if (.not. MayBranch(qIx21,qIx22,qIx2)) return

    mg = funcMG(qB,qS,qIx2)

    GammaI = pp * funcR2(m) &
         *(hadron(ID1)%spin*2+1)*(hadron(ID2)%spin*2+1) &
         /(twopi*hbarc**2*funcRho(m,mg))
    
  end function GammaI

  real function pstar2(m,ma,mb)
    real, intent(in) :: m,ma,mb
    real :: h

    pstar2 = -99
    if (ma+mb > m) return

    h = (m**2-(ma+mb)**2)*(m**2-(ma-mb)**2)
    if (h > 0) then
       pstar2 = h/(2*m)**2
    end if
  end function pstar2


  logical function MayBranch(Iso1,Iso2,Iso)
    integer, intent(in) :: Iso1,Iso2,Iso

    integer :: IsoZ1, IsoZ2
    real :: clebsch

    MayBranch = .true.

    do IsoZ1=-Iso1,Iso1,2
       do IsoZ2=-Iso2,Iso2,2
          clebsch = clebschSquared(Iso1/2.,Iso2/2.,Iso/2.,IsoZ1/2.,IsoZ2/2.)
          if (clebsch.gt.0.0) return
       end do
    end do

    MayBranch = .false.
    return
  end function MayBranch


  subroutine CalcMinMass

    use inputGeneral
    use particleProperties
    use clebschGordan

    implicit none

    real, dimension(-qBmax:qBmax,-qSmax:qSmax,0:qIx2max) :: Arr2,Arr3 ! B,S,I 

    integer :: ID
    integer :: Iso1, Ib1, Is1
    integer :: Iso2, Ib2, Is2
    integer :: Iso, Ib, Is
    real :: m
    integer :: iIt
    logical :: DidReplace
    logical, parameter :: verbose = .false.

    
    !##### 1: fill the hadron array

    Arr1 = 99.9 ! dummy

    do ID=1,121
       if (.not.(isHadron(ID))) cycle
       if (isCharmed(ID)) cycle

       Iso = hadron(ID)%isoSpinTimes2
       Is = hadron(ID)%strangeness
       Ib = 0
       if (isBaryon(ID)) Ib = Ib+1

       m = hadron(ID)%minmass
       if (hadron(ID)%stability == 0) m = hadron(ID)%mass
!       m = hadron(ID)%mass

       if (m.lt.Arr1(Ib,Is,Iso)) Arr1(Ib,Is,Iso) = m

       if (Ib.gt.0) then
          if (m.lt.Arr1(-Ib,-Is,Iso)) Arr1(-Ib,-Is,Iso) = m
       end if

    end do

    Arr1(0,0,0) = 0.400 ! setting some special value

    !##### 2: Copy the arrays

    Arr3 = 99.9 ! dummy
    Arr2 = 99.9 ! dummy
    Arr2(-1:1,-3:3,0:3) = Arr1

    write(*,*) '##### Arr 1 #####'
    do Ib=-1,1
       write(*,*) 'Baryon = ',Ib
       do Is=-3,3
          write(*,'(i3," : ",4f8.3)') Is,Arr1(Ib,Is,0:3)
       end do
    end do

    !##### 3: Do the iteration

    do iIt=1,200

       if (verbose) write(*,*)
       if (verbose) write(*,*) '----- Iteration: ',iIt
       DidReplace = .false.

       do Ib1=-1,1
          do Is1=-3,3
             do Iso1 = 0,3

                if (Arr1(Ib1,Is1,Iso1).eq.99.9) cycle

                do Ib2=-qBmax,qBmax
                   Ib = Ib1+Ib2
                   if (Ib.lt.-qBmax) cycle
                   if (Ib.gt. qBmax) cycle

                   do Is2=-qSmax,qSmax
                      Is = Is1+Is2
                      if (Is.lt.-qSmax) cycle
                      if (Is.gt. qSmax) cycle

                      do Iso2 = 0,qIx2max
                
                         if (Arr2(Ib2,Is2,Iso2).eq.99.9) cycle

                         do Iso=abs(Iso1-Iso2),Iso1+Iso2
                            if (Iso.gt.qIx2max) cycle
                            if (.not.MayBranch(Iso1,Iso2,Iso)) cycle

                            m = Arr1(Ib1,Is1,Iso1)+Arr2(Ib2,Is2,Iso2)

                            if (m.lt.Arr3(Ib,Is,Iso)) then
                               if (Arr3(Ib,Is,Iso).eq.99.9) then
                                  if (verbose) write(*,'(i3,i3,i3," : new")') Ib,Is,Iso
                               else
                                  if (verbose) write(*,'(i3,i3,i3," : replace")') Ib,Is,Iso
                               end if

                               Arr3(Ib,Is,Iso) = m
                               DidReplace = .true.

                            end if

                         end do ! Iso
                      end do ! Iso2
                   end do ! Is2
                end do ! Ib2

             end do ! Iso1
          end do ! Is1
       end do ! Ib1



       Arr2 = Arr3

       if (DidReplace.and.verbose) then
          do Ib=-qBmax,qBmax
             write(*,'(16("#"))')
             write(*,'("## Baryon =",i2," ##")') Ib
             write(*,'(16("#"))')
             write (*,'(A)') "S=  !     I=0    I=1/2    I=1    I=3/2    I=2    I=5/2    I=3    I=7/2    I=4    I=9/2    I=5   I=11/2    I=6"
             write(*,'(110("-"))')
             
             do Is=-qSmax,qSmax
                write(*,'(i3," ! ",50f8.3)') Is,Arr3(Ib,Is,0:qIx2max)
             end do
             write(*,'(110("-"))')
          end do
       end if

       if (.not.DidReplace) exit

    end do ! iIt

    do Ib=-qBmax,qBmax
       write(*,'(16("#"))')
       write(*,'("## Baryon =",i2," ##")') Ib
       write(*,'(16("#"))')
       write (*,'(A)') "S=  !     I=0    I=1/2    I=1    I=3/2    I=2    I=5/2    I=3    I=7/2    I=4    I=9/2    I=5   I=11/2    I=6"
       write(*,'(110("-"))')
       
       do Is=-qSmax,qSmax
          write(*,'(i3," ! ",50f8.3)') Is,Arr3(Ib,Is,0:qIx2max)
       end do
       write(*,'(110("-"))')
    end do

    ArrMinMass = Arr3


  end subroutine CalcMinMass

  subroutine CalcBootstrap_FillArrays
    integer :: qB,qS,qIx2, im

    ArrBootstrap = 0.0
    ArrIsZero = .true.
    ArrMinIM = nBin+1

    do qB=-1,1
       do qS=-3,3
          do qIx2 = 0,3
             
             if (Arr1(qB,qS,qIx2).eq.99.9) cycle
             select case (qB)
             case (-1)
                call FillHistBSI(-qB,-qS,qIx2)
             case (1)
                call FillHistBSI(qB,qS,qIx2)
             case (0)
                call FillHistBSI(0,abs(qS),qIx2)
             end select

             ArrHadrons(qB,qS,qIx2,:) = HistBSI(:)
             ArrBootstrap(qB,qS,qIx2,:,0) = HistBSI(:)
             ArrIsZero(qB,qS,qIx2,0) = .false.

             do im=1,nbin
                if (HistBSI(im).gt.1e-3) then
                   ArrMinIM(qB,qS,qIx2,0) = im
                   exit
                end if
             end do
          end do
       end do
    end do

  end subroutine CalcBootstrap_FillArrays

  subroutine CalcBootstrap_WriteArrays(iIt,iArr)
    integer, intent(in) :: iIt, iArr

    integer :: qB,qS,qIx2, iM

    rewind(1000+iIt)
    rewind(2000+iIt)

    do qB=-qBmax,qBmax
       do qS=-qSmax,qSmax
          do qIx2 = 0,qIx2max
!             if (ArrIsZero(qB,qS,qIx2,iArr)) cycle
             if (ArrMinIM(qB,qS,qIx2,iArr)>nBin) cycle
             write(1000+iIt,'(3i5)') qB,qS,qIx2
             do iM=1,nBin
                write(2000+iIt,'(f12.3,1p,e13.4,0p,3i5)') &
                     iM*dM, ArrBootstrap(qB,qS,qIx2,iM,iArr),qB,qS,qIx2
             end do
          end do
       end do
    end do
  end subroutine CalcBootstrap_WriteArrays

  subroutine CalcBootstrap_DoConvolute(iArr,iArr0,iCount,iCountIso)

    integer, intent(in) :: iArr,iArr0
    integer, intent(inout) :: iCount, iCountIso

    integer :: qB,qS,qIx2, iM
    integer :: qB1,qS1,qIx21, iM1
    integer :: qB2,qS2,qIx22, iM2
    logical :: newIso 
    
    real :: m, m1, m2
    real :: r1, r2, Sum, h, pcm
    real, dimension(1:nBin) :: ArrIsoH

    do qB1=-qBmax,qBmax
       do qS1=-qSmax,qSmax
          do qIx21 = 0,qIx2max
             if (ArrIsZero(qB1,qS1,qIx21,iArr0)) cycle
             
             do qB2=-qBmax,qBmax
                   qB = qB1+qB2
                   if (qB.lt.-qBmax) cycle
                   if (qB.gt. qBmax) cycle

                   do qS2=-qSmax,qSmax
                      qS = qS1+qS2
                      if (qS.lt.-qSmax) cycle
                      if (qS.gt. qSmax) cycle

                      do qIx22 = 0,qIx2max

                         if (ArrIsZero(qB2,qS2,qIx22,3-iArr)) cycle

                         newIso = .true.
                         do qIx2=abs(qIx21-qIx22),qIx21+qIx22
                            if (qIx2.gt.qIx2max) cycle
                            if (.not.MayBranch(qIx21,qIx22,qIx2)) cycle

                            ! now we have the quantum numbers, the mass integrations may start...

                            write(*,'("Doing: ",i7,2("[",3i5,"]")," --> ",("[",3i5,"]"),20i5)') &
                                 iCount,qB1,qS1,qIx21,qB2,qS2,qIx22,qB,qS,qIx2
                            iCount = iCount + 1
                            ArrIsZero(qB,qS,qIx2,iArr) = .false.

                            if (NewIso) then
                               do im=1,nBin
                                  m = im*dM
                                  
                                  sum = 0
                                  
                                  do im1=1,im
                                     m1 = im1*dM
!                                     r1 = ArrHadrons(qB1,qS1,qIx21,im1)
                                     r1 = ArrBootstrap(qB1,qS1,qIx21,im1,iArr0)
                                     if (r1 < 1e-3) cycle
                                     
                                     do im2=1,im-im1-1
                                        m2 = im2*dM
                                        r2 = ArrBootstrap(qB2,qS2,qIx22,im2,3-iArr)
                                        if (r2 < 1e-3) cycle
                                        
!!$                                     E1 = dM*(im**2+im1**2-im2**2)/(2*im)
!!$                                     E2 = dM*(im**2-im1**2+im2**2)/(2*im)
!!$                                     h = ((im**2-im1**2-im2**2)**2-4*im1**2*im2**2)
!!$                                     if (h < 0) then 
!!$                                        write(*,*) '----------h<0: (Integer)'
!!$                                        write(*,*) m,m1,m2
!!$                                        write(*,*) im,im1,im2
!!$                                        write(*,*) h
!!$                                        write(*,*) im**2,im1**2,im2**2
!!$                                        write(*,*) (im**2-im1**2-im2**2),(im**2-im1**2-im2**2)**2
!!$                                        stop
!!$                                     end if
!!$                                     pcm = dM*sqrt(h)/(2*im)


!                                        E1 = (m**2+m1**2-m2**2)/(2*m)
!                                        E2 = (m**2-m1**2+m2**2)/(2*m)
                                        h = (m**2-m1**2-m2**2)**2-4*m1**2*m2**2
                                        if (h < 0) then 
                                           write(*,*) '----------h<0:'
                                           write(*,*) m,m1,m2
                                           write(*,*) im,im1,im2
                                           write(*,*) h
                                           cycle
                                        end if
                                        pcm = sqrt(h)/(2*m)

!                                        sum = sum + r1*r2*E1*E2*pcm/m
!                                        sum = sum + r1*r2*m1*m2*pcm/m
                                        sum = sum + r1*r2*m1*m2*pcm ! /m

                                     end do ! im2
                                  end do ! im1

                                  ArrIsoH(im) = sum/m

                               end do ! im
                               ArrIsoH = ArrIsoH * dM**2 * 6.24
                               ArrBootstrap(qB,qS,qIx2,:,iArr) = ArrBootstrap(qB,qS,qIx2,:,iArr) &
                                       + ArrIsoH(:)
                               NewIso = .false.
                               iCountIso = iCountIso + 1

                            else ! NewIso
                               ArrBootstrap(qB,qS,qIx2,:,iArr) = ArrBootstrap(qB,qS,qIx2,:,iArr) &
                                       + ArrIsoH(:)
                            end if

                         end do ! Iso
                      end do ! Iso2
                   end do ! Is2
                end do ! qB2

             end do ! qIx21
          end do ! qS1
       end do ! qB1


  end subroutine CalcBootstrap_DoConvolute

  subroutine CalcBootstrap
    integer :: iArr !, im, im1, im2
    integer :: iIt, iCount, iCountIso
    logical, parameter :: verbose = .true.

    !##### fill the array

    call CalcBootstrap_FillArrays

    iArr=1
    ArrBootstrap(:,:,:,:,iArr) = ArrBootstrap(:,:,:,:,0)
    ArrIsZero(:,:,:,iArr) = ArrIsZero(:,:,:,0)



    !##### 

    do iIt=1,200

       if (verbose) then
          write(*,*)
          write(*,*) '----- Iteration: ',iIt
          call CalcBootstrap_WriteArrays(iIt, iArr)
       end if
       call timeMeasurement(.true.)

       if (iIt==10) exit
       if (iIt==5) exit

       iArr = 3-iArr
       ArrBootstrap(:,:,:,:,iArr) = 0.0
       ArrIsZero(:,:,:,iArr) = .true.
       iCount = 1
       iCountIso = 1
       
       call CalcBootstrap_DoConvolute(iArr,0,iCount,iCountIso)
       write(123,*) 'iIt,iCount,iCountIso:',iIt,iCount,iCountIso
       call TimeMeasurement(iFile=123)

    end do ! iIt
  end subroutine CalcBootstrap

  subroutine CalcGamma(qB,qS,qIx2,mass, gamma,nDecay,ArrDecayGamma,ArrDecayPart, verbose)
    use constants, only : twopi,hbarc

    integer, intent(in) :: qB,qS,qIx2
    real, intent(in) :: mass
    real, intent(out) :: gamma
    integer, intent(out) :: nDecay
    real, dimension(:), intent(inout), OPTIONAL :: ArrDecayGamma
    integer, dimension(:,:,:), intent(inout), OPTIONAL :: ArrDecayPart ! nDecay x (qB,qS,qIx2) x 2
    logical, intent(in) :: verbose

    integer :: qB1,qS1,qIx21
    integer :: qB2,qS2,qIx22
    integer :: qIx22min, qIx22max
    integer :: im,im1,im2
    real :: fak,sum
    real :: r1,r2,m1,m2,pcm2

    nDecay = 0
    gamma = 0
    im = nint(mass/dM)
    if (im.gt.nBin) return

    fak = funcR2(mass)/(twopi*hbarc**2*ArrBootstrap(qB,qS,qIx2,im,1))
!    write(*,*) 'fak:',fak,ArrBootstrap(qB,qS,qIx2,im,1)

    do qB1=-qBmax,qBmax
       qB2 = qB-qB1
       if (qB2.lt. qB1) cycle ! avoid double counting
       if (qB2.gt. qBmax) cycle

       do qS1=-qSmax,qSmax
          qS2 = qS-qS1
          if (qS2.lt.-qSmax) cycle
          if (qS2.gt. qSmax) cycle
          if ((qB1==qB2).and.(qS2<qS1)) cycle ! avoid double counting

          do qIx21 = 0,qIx2max
             qIx22min = abs(qIx2-qIx21)
             qIx22max = min(qIx2+qIx21,qIx2max)
             
             do qIx22 = qIx22min,qIx22max

                if (.not.MayBranch(qIx21,qIx22,qIx2)) then
!                   write(*,*) 'not branch:',qIx21,qIx22,qIx2
                   cycle
                end if

                sum = 0
                do im1=1,im
                   do im2=1,im-im1
                      r1 = ArrBootstrap(qB1,qS1,qIx21,im1,1)
                      if (r1 < 1e-3) cycle
                      r2 = ArrBootstrap(qB2,qS2,qIx22,im2,1)
                      if (r2 < 1e-3) cycle
                      m1 = im1*dM
                      m2 = im2*dM
                      pcm2 = ((mass**2-m1**2-m2**2)**2-4*m1**2*m2**2)/(2*mass)**2
                      if (pcm2 < 0) cycle
                      sum = sum + r1*r2*pcm2
                      !                         write(*,*) m1,m2,r1,r2,pcm2
                   end do ! im2
                end do ! im1


                if (sum==0) cycle

                sum = sum*fak*dM**2 * 2 ! *2 because of double counting
                gamma = gamma+sum

                if (sum<1e-10) cycle 

                nDecay = nDecay+1
                if (present(ArrDecayGamma)) then
                   if (nDecay>size(ArrDecayGamma)) then
                      write(*,*) "Array too small:",nDecay,size(ArrDecayGamma)
                      stop
                   end if
                   ArrDecayGamma(nDecay) = gamma
                end if

                if (verbose) then
                   write(*,'(i4,2("[",3i5,"]"),1P,e13.5,0P)') nDecay,&
                        qB1,qS1,qIx21, &
                        qB2,qS2,qIx22,sum
                end if

                
             end do ! qIx22
          end do ! qIx21
       end do ! qS1
    end do ! qS2
    if (verbose) write(*,*) "===>",gamma

  end subroutine CalcGamma

  subroutine PrintGammaPlot(iFile,qB,qS,qIx2)
    integer, intent(in) :: iFile,qB,qS,qIx2

    integer :: iM
    real :: gamma
    integer :: nDecay

    write(*,*) 'Plot:',qB,qS,qIx2

    !    do iM = 1,nBin
    do iM = 1,nBin/2
       call CalcGamma(qB,qS,qIx2, iM*dM, gamma,nDecay,verbose=.false.)
       write(7000+iFile,*) iM*dM,gamma,nDecay
    end do


  end subroutine PrintGammaPlot

  subroutine ReadFort
    integer :: iIt,iArr,iM
    integer :: ios
    integer :: qB,qS,qIx2
    real :: x,y

    ArrIsZero = .true.
    ArrBootstrap = 0.0
    iArr = 1
    do iIt=1,20
       write(*,*) 'Reading ',iIt,':...'

       do
          read (1000+iIt,'(3i5)',iostat=ios) qB,qS,qIx2
          if (ios/=0) exit
          ArrIsZero(qB,qS,qIx2,iArr) = .false.

          do iM=1,nBin
             read (2000+iIt,*,iostat=ios) x,y,qB,qS,qIx2
             if (ios/=0) exit
             ArrBootstrap(qB,qS,qIx2,iM,iArr) = ArrBootstrap(qB,qS,qIx2,iM,iArr) + y
          end do

       end do

       !### copy iteration==1 into Hadrons array

       if (iIt==1) then
          do qB=-1,1
             do qS=-3,3
                do qIx2 = 0,3
                   ArrHadrons(qB,qS,qIx2,:) =  ArrBootstrap(qB,qS,qIx2,:,iArr)
                end do
             end do
          end do
       end if

    end do

    open(131,file='ArrBootstrap.dat',status='UNKNOWN',form='UNFORMATTED')
    rewind(131)
    write(131) qBmax, qSmax, qIx2max
    write(131) ArrIsZero
    write(131) ArrBootstrap
    close(131)

  end subroutine ReadFort


  subroutine PrintSomeGamma

    real :: gamma
    integer :: nDecay
    real, dimension(500) :: ArrDecayGamma
    integer, dimension(500,3,2) :: ArrDecayPart ! nDecay x (qB,qS,qIx2) x 2

    call CalcGamma(2,0,0, 2.6, gamma,nDecay,ArrDecayGamma,ArrDecayPart,.true.)
    write(*,*) '===='
    call CalcGamma(2,0,0, 4.1, gamma,nDecay,ArrDecayGamma,ArrDecayPart,.true.)
    write(*,*) '===='
    call CalcGamma(2,0,2, 2.6, gamma,nDecay,ArrDecayGamma,ArrDecayPart,.true.)
    write(*,*) '===='
    call CalcGamma(0,0,0, 2.6, gamma,nDecay,ArrDecayGamma,ArrDecayPart,.true.)

    call PrintGammaPlot(1, 2,0,0)
    call PrintGammaPlot(2, 2,0,2)

    call PrintGammaPlot(11, 3,0,1)
    call PrintGammaPlot(12, 3,0,3)
    call PrintGammaPlot(13, 3,0,5)
    
    call PrintGammaPlot(21, 0,0,0)
    call PrintGammaPlot(22, 0,0,2)
    call PrintGammaPlot(23, 0,1,1)
  end subroutine PrintSomeGamma

  subroutine DoFinalFolding

    integer :: iArr, iCount, iCountIso

    write(*,*) 'Doing a final convolutuion...'
    call timeMeasurement(.true.)
    iArr = 2
    ArrBootstrap(:,:,:,:,iArr) = 0.0
    ArrIsZero(:,:,:,iArr) = .true.
    iCount = 1
    iCountIso = 1
    call CalcBootstrap_DoConvolute(iArr,3-iArr,iCount,iCountIso)
    write(123,*) 'iIt,iCount,iCountIso:',999,iCount,iCountIso
    call TimeMeasurement(iFile=123)
    call CalcBootstrap_WriteArrays(999, iArr)
    
  end subroutine DoFinalFolding

  subroutine DoHamerFrautschi_DoConvolute(im)
    
    integer, intent(in) :: im

    integer :: qB,qS,qIx2
    integer :: qB1,qS1,qIx21, iM1
    integer :: qB2,qS2,qIx22, iM2
    logical :: newIso 
    integer :: iArr

    real :: mm1, mm2, mm
    real :: E1, E2
    real :: r1, r2, Sum, h, pcm
    real :: ArrIsoH
    integer :: iCount, iCountIso
    integer :: im1min, im2min

    iArr = 1
    iCount = 0
    iCountIso = 0

    call timeMeasurement(.true.)

    !    m = im*dM
    mm = (im*dM)**2

    do qB1=-qBmax,qBmax
       do qS1=-qSmax,qSmax
          do qIx21 = 0,qIx2max
             im1min = ArrMinIM(qB1,qS1,qIx21,iArr)
             if (im1min > im) cycle
!             if (ArrIsZero(qB1,qS1,qIx21,iArr)) cycle

             do qB2=-qBmax,qBmax
                   qB = qB1+qB2
                   if (qB.lt.-qBmax) cycle
                   if (qB.gt. qBmax) cycle

                   do qS2=-qSmax,qSmax
                      qS = qS1+qS2
                      if (qS.lt.-qSmax) cycle
                      if (qS.gt. qSmax) cycle

                      do qIx22 = 0,qIx2max
                         im2min = ArrMinIM(qB2,qS2,qIx22,iArr)
                         if (im2min > im) cycle
!                         if (ArrIsZero(qB2,qS2,qIx22,iArr)) cycle

                         newIso = .true.
                         do qIx2=abs(qIx21-qIx22),qIx21+qIx22
                            if (qIx2.gt.qIx2max) cycle
                            if (.not.MayBranch(qIx21,qIx22,qIx2)) cycle

                            ! now we have the quantum numbers, the mass integrations may start...

!                            write(*,'("Doing: ",i7,2("[",3i5,"]")," --> ",("[",3i5,"]"),20i5)') &
!                                 iCount,qB1,qS1,qIx21,qB2,qS2,qIx22,qB,qS,qIx2

                            if (NewIso) then
                               
                               sum = 0
                               
                               do im1=im1min,im-1
                                  r1 = ArrBootstrap(qB1,qS1,qIx21,im1,iArr)
                                  if (r1 < 1e-3) cycle
                                  mm1 = (im1*dM)**2
                                  
                                  do im2=im2min,im-im1-1
                                     r2 = ArrBootstrap(qB2,qS2,qIx22,im2,iArr)
                                     if (r2 < 1e-3) cycle
                                     mm2 = (im2*dM)**2
                                     
                                     h = (mm-mm1-mm2)**2-4*mm1*mm2
                                     if (h < 0) then 
                                        write(*,*) '----------h<0:'
                                        write(*,*) mm,mm1,mm2
                                        write(*,*) im,im1,im2
                                        write(*,*) h
                                        cycle
                                     end if
                                     pcm = sqrt(h*mm1*mm2) ! /(2*m)
                                     E1 = (mm+mm1-mm2) !/(2*m)
                                     E2 = (mm-mm1+mm2) !/(2*m)
                                     
!                                     sum = sum + r1*r2*m1*m2*pcm ! /m
!                                     sum = sum + r1*r2*pcm/(E1*E2) ! /m
                                     sum = sum + r1*r2*pcm ! /m

!!$                                     sum = sum + r1*r2*arrPCM(im1,im2)

                                  end do ! im2
                               end do ! im1

                               ArrIsoH = sum  * dM**2 * (1593.8/(16*3.14**3)) / (2*mm) * 4*3.14
!                               ArrIsoH = sum * (2*(im*dM))**2 * dM**2 *(547.9/(16*3.14**3)) ! see above: /m and /(2*m) 
!                               ArrIsoH = sum * (2*(im*dM)) * dM**2 *(1593.8/(16*3.14**3)) ! see above: /m and /(2*m) 
                               NewIso = .false.
                               if (ArrIsoH == 0.0) cycle
                               ArrBootstrap(qB,qS,qIx2,iM,iArr) = ArrBootstrap(qB,qS,qIx2,iM,iArr) &
                                       + ArrIsoH
                               iCountIso = iCountIso + 1

                            else ! NewIso
                               if (ArrIsoH == 0.0) cycle
                               ArrBootstrap(qB,qS,qIx2,iM,iArr) = ArrBootstrap(qB,qS,qIx2,iM,iArr) &
                                       + ArrIsoH
                            end if


                            iCount = iCount + 1
!                            ArrIsZero(qB,qS,qIx2,iArr) = .false.
                            
                            if (im < ArrMinIM(qB,qS,qIx2,iArr)) then
                               if ((ArrBootstrap(qB,qS,qIx2,iM,iArr) > 1e-3)) then
                                  ArrMinIM(qB,qS,qIx2,iArr) = im
                               end if
                            end if

                         end do ! qIx2
                      end do ! qIx22
                   end do ! qS2
                end do ! qB2

             end do ! qIx21
          end do ! qS1
       end do ! qB1

       write(123,*) 'iM,iCount,iCountIso:',iM,iCount,iCountIso
       call TimeMeasurement(iFile=123)

  end subroutine DoHamerFrautschi_DoConvolute

  subroutine DoHamerFrautschi_CalcGamma(im)
    use constants, only : twopi,hbarc

    integer, intent(in) :: im

    integer :: qB,qS,qIx2
    integer :: qB1,qS1,qIx21, iM1
    integer :: qB2,qS2,qIx22, iM2
    logical :: newIso 
    integer :: iArr

    real :: mm1, mm2, mm
    real :: r0, r1, r2, Sum, h
    real :: ArrIsoH
    integer :: iCount, iCountIso
    integer :: im1min, im2min
    real :: fak,fakSymmetry
    real :: E1, E2

    iArr = 1
    iCount = 0
    iCountIso = 0

    call timeMeasurement(.true.)

    fak = funcR2(im*dM)/(twopi*hbarc**2)

    !    m = im*dM
    mm = (im*dM)**2

    do qB1=-qBmax,qBmax
       do qS1=-qSmax,qSmax
          do qIx21 = 0,qIx2max
             im1min = ArrMinIM(qB1,qS1,qIx21,1)
             if (im1min > im) cycle

             do qB2=-qBmax,qBmax
                   qB = qB1+qB2
                   if (qB.lt.-qBmax) cycle
                   if (qB.gt. qBmax) cycle

                   do qS2=-qSmax,qSmax
                      qS = qS1+qS2
                      if (qS.lt.-qSmax) cycle
                      if (qS.gt. qSmax) cycle

                      do qIx22 = 0,qIx2max
                         im2min = ArrMinIM(qB2,qS2,qIx22,1)
                         if (im2min > im) cycle

                         newIso = .true.
                         do qIx2=abs(qIx21-qIx22),qIx21+qIx22
                            if (qIx2.gt.qIx2max) cycle
                            if (.not.MayBranch(qIx21,qIx22,qIx2)) cycle

                            ! now we have the quantum numbers, the mass integrations may start...

!                            write(*,'("Doing: ",i7,2("[",3i5,"]")," --> ",("[",3i5,"]"),20i5)') &
!                                 iCount,qB1,qS1,qIx21,qB2,qS2,qIx22,qB,qS,qIx2

                            r0 = ArrBootstrap(qB,qS,qIx2,iM,1)-ArrBootstrap(qB,qS,qIx2,iM,0)
                            if (r0<1e-3) cycle

                            if (NewIso) then

                               fakSymmetry = 1.0
                               if ((qB1==qB2).and.(qS1==qS2).and.(qIx21==qIx22)) then
                                  fakSymmetry = (2.0*(qIx21+1)**2-1.0)/(2.0*(qIx21+1)**2)
                               end if
                               
                               sum = 0
                               
                               do im1=im1min,im-1
                                  r1 = ArrBootstrap(qB1,qS1,qIx21,im1,1)
                                  if (r1 < 1e-3) cycle
                                  mm1 = (im1*dM)**2
                                  
                                  do im2=im2min,im-im1-1
                                     r2 = ArrBootstrap(qB2,qS2,qIx22,im2,1)
                                     if (r2 < 1e-3) cycle
                                     mm2 = (im2*dM)**2
                                     
                                     h = (mm-mm1-mm2)**2-4*mm1*mm2
                                     if (h < 0) then 
                                        write(*,*) '----------h<0:'
                                        write(*,*) mm,mm1,mm2
                                        write(*,*) im,im1,im2
                                        write(*,*) h
                                        cycle
                                     end if
                                     E1 = (mm+mm1-mm2) !/(2*m)
                                     E2 = (mm-mm1+mm2) !/(2*m)

!                                     sum = sum + r1*r2*h ! /m
                                     sum = sum + r1*r2*h *sqrt(mm1*mm2)/(E1*E2) ! /m

                                  end do ! im2
                               end do ! im1

                               fakSymmetry = 1.0
                               if ((qB1==qB2).and.(qS1==qS2).and.(qIx21==qIx22)) then
                                  fakSymmetry = (2.0*(qIx21+1)**2-1.0)/(2.0*(qIx21+1)**2)
                               end if

                               ArrIsoH = fak * sum * dM**2 / r0 * fakSymmetry

                               NewIso = .false.
                               if (ArrIsoH == 0.0) cycle
                               ArrBootstrap(qB,qS,qIx2,iM,2) = ArrBootstrap(qB,qS,qIx2,iM,2) &
                                       + ArrIsoH
                               iCountIso = iCountIso + 1

                            else ! NewIso
                               if (ArrIsoH == 0.0) cycle
                               ArrBootstrap(qB,qS,qIx2,iM,2) = ArrBootstrap(qB,qS,qIx2,iM,2) &
                                       + ArrIsoH
                            end if


                            iCount = iCount + 1
                            
                            if (im < ArrMinIM(qB,qS,qIx2,2)) then
                               if ((ArrBootstrap(qB,qS,qIx2,iM,2) > 1e-3)) then
                                  ArrMinIM(qB,qS,qIx2,2) = im
                               end if
                            end if

                         end do ! qIx2
                      end do ! qIx22
                   end do ! qS2
                end do ! qB2

             end do ! qIx21
          end do ! qS1
       end do ! qB1

       write(123,*) 'iM,iCount,iCountIso:',iM,iCount,iCountIso
       call TimeMeasurement(iFile=123)

     end subroutine DoHamerFrautschi_CalcGamma

  subroutine DoHamerFrautschi

    integer :: im


    call CalcBootstrap_FillArrays
    call CalcBootstrap_WriteArrays(2001,0)

    ! copy hadrons into Hagedorn array:
    ArrBootstrap(:,:,:,:,1) = ArrBootstrap(:,:,:,:,0)
    ArrIsZero(:,:,:,1) = ArrIsZero(:,:,:,0)
    ArrMinIM(:,:,:,1) = ArrMinIM(:,:,:,0)

    do im=1,nBin
       write(*,*) 'Doing iM=',iM

       call DoHamerFrautschi_DoConvolute(im)
       call DoHamerFrautschi_CalcGamma(im)

       if (mod(im,10)==0) then
          call CalcBootstrap_WriteArrays(2000,1)
          call CalcBootstrap_WriteArrays(2002,2)

          open(131,file='ArrBootstrap.dat',status='UNKNOWN',form='UNFORMATTED')
          rewind(131)
          write(131) qBmax, qSmax, qIx2max, nBin
          write(131) ArrBootstrap
          close(131)
       end if


    end do


  end subroutine DoHamerFrautschi


end program testHagedornstate
