!***************************************************************************
!****m* /barAntiBar
! PURPOSE
! Implements baryon+antibaryon -> X cross section
! NOTES
!***************************************************************************
module barAntiBar

  implicit none
  Private

  !*************************************************************************
  !****g* barAntiBar/fact_LambdaBar
  ! SOURCE
  !
  real, save :: fact_LambdaBar=1.
  ! PURPOSE
  ! Enhancement factor of pbar p -> Lambda LambdaBar cross section
  ! (for larger statistics)
  !*************************************************************************

  !*************************************************************************
  !****g* barAntiBar/fact_JPsi
  ! SOURCE
  !
  real, save :: fact_JPsi=1.
  ! PURPOSE
  ! Enhancement factor of pbar p -> J/Psi cross section (for larger statistics)
  !*************************************************************************

  !*************************************************************************
  !****g* barAntiBar/fact_JPsi_width
  ! SOURCE
  !
  real, save :: fact_JPsi_width=1.
  ! PURPOSE
  ! Enhancement factor of the J/Psi total width (for larger statistics)
  !*************************************************************************

  !*****************************************************************************
  !****g* barAntiBar/useAnni
  ! SOURCE
  !
  logical,save :: useAnni = .true.
  ! PURPOSE
  ! Flag whether to perform Baryon-Antibarion annihilation or not at all
  !*****************************************************************************


  logical, save :: initFlag=.true.

  Public :: sigmaBarAntiBar

contains

  !***************************************************************************
  !****s* barAntiBar/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads in namelist "barAntiBar_input"
  ! INPUTS
  ! * (none)
  ! OUTPUT
  ! * Initializes global module variables
  !***************************************************************************
  subroutine init
    use output

    integer :: ios

    !*************************************************************************
    !****n* barAntiBar/barAntiBar_input
    ! NAME
    ! NAMELIST barAntiBar_input
    ! PURPOSE
    ! Namelist which includes the input variables:
    ! * fact_LambdaBar
    ! * fact_JPsi
    ! * fact_JPsi_width
    ! * useAnni
    !*************************************************************************
    NAMELIST /barAntiBar_input/ fact_LambdaBar,fact_JPsi,fact_JPsi_width, &
         useAnni

    call Write_ReadingInput('barAntiBar_input',0)
    rewind(5)
    read(5,nml=barAntiBar_input,iostat=ios)
    call Write_ReadingInput('barAntiBar_input',0,ios)

    if(fact_LambdaBar.ne.1.) then
       write(*,*) ' ATTENTION: pbar p -> LambdaBar Lambda cross section'
       write(*,*) ' is rescaled by a factor of ', fact_LambdaBar
    end if

    if(fact_JPsi.ne.1.) then
       write(*,*) ' ATTENTION: pbar p -> J/Psi cross section'
       write(*,*) ' is rescaled by a factor of ', fact_JPsi
    end if

    if(fact_JPsi_width.ne.1.) then
       write(*,*) ' ATTENTION: J/Psi total width'
       write(*,*) ' is rescaled by a factor of ', fact_JPsi_width
    end if

    write(*,*) 'use Anni: ',useAnni

    call Write_ReadingInput('barAntiBar_input',1)

    initFlag = .false.
  end subroutine init


  !********************************************************************
  !****s* barAntiBar/sigmaBarAntiBar
  ! NAME
  ! subroutine sigmaBarAntiBar(srts,teilchenIN,mediumATcollision,sigTotal,
  ! sigElastic,sigChEx,sigAnnihilation,
  ! sigProduction,sigHyperon,sigLambdaBar,sigSigmaBar,sigXiBar,
  ! sigOmegaBar,sigJPsi)
  !
  ! PURPOSE
  ! Computes total, elastic and other cross sections for baryon+antibaryon
  ! collisions.
  ! INPUTS
  ! * real :: srts --- sqrt(s) of collision (GeV)
  ! * type(particle),dimension(1:2) :: teilchenIn --- colliding particles
  ! * type(medium) :: mediumATcollision  --- medium infos at collision point
  ! OUTPUT
  ! * real ::           sigTotal     --- total cross section (mb)
  ! * real ::           sigElastic   --- elastic cross section (mb)
  ! * real, optional :: sigChEx      ---
  !   charge exchange cross section (mb)
  ! * real, optional :: sigAnnihilation ---
  !   annihilation into mesons cross section (mb)
  ! * real, optional :: sigProduction   ---
  !   production Bbar+B -> B+Bbar+mesons cross section (mb)
  ! * real, optional :: sigHyperon      ---
  !   (anti)hyperon production Bbar+B -> Y+Ybar+mesons,
  !   B+Ybar+Kbar, Bbar+Y+K cross section (mb)
  ! * real, optional :: sigLambdaBar ---
  !   exclusive (anti)hyperon production:
  !   Bbar+B -> Lambda+Lambdabar
  ! * real, optional :: sigSigmaBar ---
  !   exclusive (anti)hyperon production :
  !   Bbar+B -> Lambda+Sigma0Bar, LambdaBar+Sigma0 (mb)
  ! * real, optional :: sigXiBar ---
  !   exclusive (anti)cascade production:
  !   Bbar+B -> Xi+XiBar
  ! * real, optional :: sigOmegaBar ---
  !   exclusive (anti)Omega (S=-3) production:
  !   Bbar+B -> Omega+OmegaBar
  ! * real, optional :: sigJPsi ---
  !   J/Psi production cross section (mb):
  !   Bbar+B -> J/Psi
  !********************************************************************
  subroutine sigmaBarAntiBar(srts,teilchenIN,mediumATcollision,sigTotal, &
                           & sigElastic,sigChEx,sigAnnihilation, &
                           & sigProduction,sigHyperon, &
                           & sigLambdaBar,sigSigmaBar,sigXiBar,sigOmegaBar, &
                           & sigJPsi)

    use mediumDefinition
    use particleDefinition
    use IdTable, only : nucleon,delta,Xi,JPsi
    use constants, only : rhoNull, mN, pi
    use twoBodyTools, only : pcm
    use ParticleProperties, only : hadron

    real, intent(in)   ::           srts
    type(particle),dimension(1:2), intent(in)    :: teilchenIn
    type(medium),                  intent(in)    :: mediumATcollision
    real, intent(out)  ::           sigTotal
    real, intent(out)  ::           sigElastic
    real, optional, intent(out)  :: sigChEx
    real, optional, intent(out)  :: sigAnnihilation
    real, optional, intent(out)  :: sigProduction
    real, optional, intent(out)  :: sigHyperon
    real, optional, intent(out)  :: sigLambdaBar
    real, optional, intent(out)  :: sigSigmaBar
    real, optional, intent(out)  :: sigXiBar
    real, optional, intent(out)  :: sigOmegaBar
    real, optional, intent(out)  :: sigJPsi

    real :: m1,m2,s,p,sigCEX,sigANN,sigPROD,sigY,u,exc,p12,Gamma_JPsi
    integer :: totCharge, totId, nDelta, I3_Delta

    ! momentum above which the PDG parameterization is working (GeV/c):
    real, parameter :: pCut=1.85
    real :: E,v_rel

    ! J/Psi total width (GeV) and branching ratio to pbar p (nbar n):
    real, parameter :: Gamma_JPsi0=92.9e-06, Br=2.2e-03

    ! J/Psi partial width to pbar p (nbar n):
    real, parameter :: Gamma_JPsi_pbarp=Gamma_JPsi0*Br

    ! If .true. then annihilation cross section is medium-modified according to
    ! E. Hernandez and E. Oset, Z. Phys. A 341, 201 (1992)
    logical, parameter :: densityDependence_Oset=.false.


    if(initFlag) call init

    sigCEX = 0.
    sigANN = 0.
    sigProd = 0.
    sigY = 0.

    s=srts**2
    m1=teilchenIn(1)%mass
    m2=teilchenIn(2)%mass
    E=(s-m1**2-m2**2)/(2.*m2)       ! Energy of 1-st particle in lab frame
    p=sqrt(E**2-m1**2)              ! Momentum of 1-st particle in lab frame
    v_rel=p/E
    p=v_rel*mN/sqrt(1-v_rel**2)  ! Momentum of antiproton for the same relative velocity

    p12=pcm(srts,m1,m2)   ! c.m. momentum of colliding particles

    !s0=4.*0.938**2
    !p=sqrt(max(s**2/s0-s,1e-06))           ! lab. momentum

    totCharge= sum(teilchenIn(1:2)%charge)
    totId= sum(teilchenIN(1:2)%ID)

    if(totId.eq.3) then
       if(teilchenIN(1)%ID.eq.2) then
          nDelta=1
       else
          nDelta=2
       end if
       if(teilchenIN(nDelta)%antiparticle) then
          I3_Delta=2*teilchenIN(nDelta)%charge+1  ! twice I3 projection by Gell-Mann - Nishidjima formula
       else
          I3_Delta=2*teilchenIN(nDelta)%charge-1
       end if
    end if

    if( totId.eq.2 .and. totCharge.eq.0 &
                                   & .or.&
       &totId.eq.3 .and. abs(totCharge).ne.2 ) then

       ! Parameterizations for pbar+p -> nbar+n from J. Cugnon and J. Vandermeulen,
       ! Annales de Physique (France) 14, 49 (1989);
       ! see also C.B. Dover et al., Prog. Part. Nucl. Phys. 29, 87 (1992):

       if(present(sigChEx) .or. p <= pCut) then
         if( p < 0.5 ) then
           if( p <= 0.1 ) then
             sigCEX= 0.
           else
             sigCEX= 10.9*(p-0.1)/p**1.6
           end if
         else
           sigCEX= 7.1/p**0.9       ! This is actually not good at p > 3 GeV/c
                                    ! (above experiment)
         end if
       end if

    else

       sigCEX= 0.

    end if


    if(present(sigAnnihilation) .or. p <= pCut) then
      if(p.lt.0.51) then
        if( totId.eq.2 .and. totCharge.ne.0 &
           &.and. p.lt.0.382 ) then
           ! Parameterization for low energy nbar+p annihilation cross section
           ! from T. Armstrong et al., PRD 36, 659 (1987):
           sigANN= 41.4 + 29./p
        else
           ! pbar+p data parameterization by A.L.:
           sigANN= 51.52/p**0.85 + 0.034/p**2.94
        end if
      else if(p.lt.6.34) then
        sigANN=88.8/p**0.4 - 24.2             ! A.L.
      else
        sigANN= 38./p**0.5 + 24./p**1.1         ! Cugnon
      end if
      if(densityDependence_Oset.and.mediumATcollision%useMedium) then
        u=mediumATcollision%density/rhoNull
        sigANN=sigANN*( 1. + 4.59*u + 10.6*u**2 + 12.8*u**3 )
      end if
    end if

    if(present(sigProduction) .or. p <= pCut) then
      if(totId <= 3) then
         if( p <= 0.793 ) then
           sigPROD= 0.
         else
           sigPROD= 30.*(p-0.793)**1.5/(2.+(p-0.793)**1.5)
         end if
      else
         sigPROD= 0.
      end if
    end if

    if(present(sigHyperon) .or. p <= pCut) then
      if(totId <= 3) then
         if( p <= 1.435 ) then
           sigY= 0.
         else
           sigY= 3.*(p-1.435)/(10.+(p-1.435))
         end if
      else
         sigY= 0.
      end if
    end if

    if (present(sigLambdaBar)) then

       if (totId.eq.2 .and. totCharge.eq.0 ) then

          if (s <= 4.982) then !Lambda+LambdaBar Threshold [GeV**2]
             sigLambdaBar = 0.0
          else if (s < 5.071) then
             exc = sqrt(s) - 2.232 !excess energy
             sigLambdaBar = 0.0357345*sqrt(exc)+9.00539*exc**1.5                !mb
          else
             sigLambdaBar = max(0.010,0.725872*(s/4.982-1.)**0.774485*(4.982/s)**3.35044)  !mb
          endif

       else

          sigLambdaBar = 0.0

       endif

    endif


    if (present(sigSigmaBar)) then

       if ( totId <= 3 .and. abs(totCharge) <= 1 ) then

          if (s <= 5.327) then !Lambda+Sigma0 Threshold [GeV**2]
             sigSigmaBar = 0.0
          else
             !sigSigmaBar is fitted to Lambda+Sigma0Bar+c.c.
             !sigSigmaBar = 2.43292*(s/5.327-1.)**1.2983*(5.327/s)**11.1481    !mb
             sigSigmaBar = 0.183665*(s/5.327-1.)**0.436713*(5.327/s)**1.85006    !mb
          endif

          if(totCharge /= 0) then   ! Multiplicative factors from isospin relations
             if(totId .eq. 2) then
                sigSigmaBar=2.*sigSigmaBar
             else if(totId .eq. 3) then
                if(abs(I3_Delta).eq.3) then
                    sigSigmaBar=1.5*sigSigmaBar
                else
                    sigSigmaBar=0.5*sigSigmaBar
                end if
             end if
          end if

       else

          sigSigmaBar = 0.0

       endif

    endif


    if (present(sigXiBar)) then
       if( totId <= 3 .and. abs(totCharge) <= 1 ) then
          if(srts <= 2.630) then        ! Xi+XiBar threshold [GeV]
             sigXiBar=0.
          else
             sigXiBar=0.004  !mb
          end if
          if(totId .eq. 3) then ! DeltaBar+N or Delta+Nbar collisions --> multiplicative factors from isospin relations
             if(totCharge == 0) then
                sigXiBar=0.5*sigXiBar
             else if(abs(I3_Delta).eq.1) then
                sigXiBar=0.25*sigXiBar
             else
                sigXiBar=0.75*sigXiBar
             end if
          end if
       else
          sigXiBar=0.
       end if
    end if

    if (PRESENT(sigOmegaBar)) then

       if (totId==2 .and. totCharge == 0) then
          if(srts <= 3.344) then        ! Omega+OmegaBar threshold [GeV]
             sigOmegaBar=0.0
          else
             sigOmegaBar=7.6e-05*(s/11.1823-1.)**1.7212*(11.1823/s)**6.5033  !mb
          end if
       else
          sigOmegaBar = 0.0
       endif

    endif


    if (present(sigJPsi)) then
       if( totId.eq.2 .and. totCharge.eq.0 ) then
          Gamma_Jpsi=Gamma_Jpsi0*fact_JPsi_width
          sigJPsi=3.*pi*0.389/p12**2*s*Gamma_JPsi_pbarp*Gamma_Jpsi &
                &/((s-hadron(JPsi)%mass**2)**2+s*Gamma_Jpsi**2)
       else
          sigJPsi=0.
       end if
    end if


    if( p < 2.03 ) then
       ! Parameterization from J. Cugnon and J. Vandermeulen,
       ! Annales de Physique (France) 14, 49 (1989);
       ! see also C.B. Dover et al., Prog. Part. Nucl. Phys. 29, 87 (1992):
       !       sigElastic= 42.3/p**0.54 + 4.3*exp(-(p-1.5)**2)

       sigElastic= 40.0/p**0.56 + 5.8*exp(-(p-1.85)**2)   ! A.L.

    else
       ! Use PDG parameterization:
       sigElastic=10.2+52.7*p**(-1.16)+0.125*log(p)**2-1.28*log(p)
    end if


    if (.not.useAnni) sigANN = 0.

    if( p <= pCut ) then

      sigTotal= sigElastic + sigCEX + sigANN + sigPROD + sigY

    else

      ! Parameterization of total cross sections from PDG,
      ! L. Montanet et al., PRD 50, 1173 (1994) (see p. 1335);
      ! see also T. Falter et al., PRC 70, 054609 (2004):

      sigTotal=38.4+77.6*p**(-0.64)+0.26*log(p)**2-1.2*log(p)

    end if

    if(present(sigChEx)) sigChEx= sigCEX
    if(present(sigAnnihilation)) sigAnnihilation= sigANN
    if(present(sigProduction)) sigProduction= sigPROD
    if(present(sigHyperon)) sigHyperon= sigY

    ! Enhancement factors (to get more statistics):
    sigLambdaBar=sigLambdaBar*fact_LambdaBar
    sigJPsi=sigJPsi*fact_JPsi

  end subroutine sigmaBarAntiBar

end module barAntiBar
