program test
  use inputGeneral

  implicit none

  write(*,*) '**************************************************'
  write(*,*) 'Initializing database'
  call readinputGeneral
  call init_Database
  write(*,*) '**************************************************'

  write(*,*) '**************************************************'
  write(*,*) 'Testing Phi Nuk -> '
  write(*,*) '**************************************************'

  call testPhiNuc


  contains


  subroutine testPhiNuc
    use IDTABLE
    use mediumDefinition
    use particleDefinition
    use particleProperties
    use phiNucleon
    use preEventDefinition
    implicit none
    real                                          :: srts                  ! sqrt(s) in the process
    type(particle),dimension(1:2)   :: teilchenIn        ! colliding particles
    type(medium)                      :: mediumATcollision    ! Medium informations at the position of the collision

    logical                         :: plotFlag          ! Switch on plotting of the  Xsections
    real,dimension(0:3)       :: momentumLRF        ! Total Momentum in LRF
    type(preEvent),dimension(1:3) :: teilchenOut     ! colliding particles
    real                            :: sigmaTot         ! total Xsection
    real                            :: sigmaElast      ! elastic Xsecti
    integer :: i,j,chargeMeson,chargeNuk, k
    real :: plab,ekin,mass
    real :: piN,phiN, pipiN,r_p13
    real :: dens
    integer :: numTries=1000

    NAMELIST /initTest/ chargeMeson, chargeNuk, dens
    write (*,*)
    write (*,*) '**Initializing testing parameter'
    write(*,*)  ' Reading input ....'
    rewind(5)
    read(5,nml=initTest)
    write(*,*) ' Set charge to ', chargeMeson,'.'
    write(*,*) ' Set nuk charge to ',  chargeNuk,'.'
    write(*,*) ' Set density to ',  dens,'.'

    mediumAtCollision%useMedium=.true.
    mediumAtCollision%densityProton=dens/2.
    mediumAtCollision%densityNeutron=dens/2.

    teilchenIN(1)%Id=phi
    teilchenIN(1)%charge=chargeMeson
    teilchenIN(1)%mass=meson(phi)%mass

    teilchenIN(2)%Id=nucleon
    teilchenIN(2)%charge=chargeNuk
    teilchenIN(2)%mass=0.938
    teilchenIN(2)%momentum=(/teilchenIN(2)%mass,0.,0.,0./)
    teilchenIN(2)%velocity=(/0.,0.,0./)

    Open(301,file='PhiNucleonXsections.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    mesonCharge=', chargeMeson
    Do i=1,100
      plab=i*0.02
      Do j=1,100
        mass=j*0.02
        teilchenIN(1)%mass=mass
        teilchenIN(1)%momentum=(/sqrt(mass**2+plab**2),plab,0.,0./)
        teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

       srts=SQRT((SQRT(mass**2+plab**2)+baryon(nucleon)%mass)**2-plab**2)
       momentumLRF=(/SQRT(mass**2+plab**2)+baryon(nucleon)%mass ,plab,0.,0./)
       ekin=SQRT(mass**2+plab**2)-mass
       call phiNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.2,.true.)

       write(301,'(5F12.3)') plab,mass,srts,sigmaTot,sigmaElast
      end do
      write(301,*)
      print *,i
    end do
    close(301)

    piN=0.
    phiN=0.
    pipiN=0.
    Do i=1,numTries
       plab=0.2
       srts=SQRT((SQRT(meson(phi)%mass**2+plab**2)+baryon(nucleon)%mass)**2-plab**2)
       teilchenIN(1)%momentum=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

       momentumLRF=(/SQRT(meson(phi)%mass**2+plab**2)+baryon(nucleon)%mass ,plab,0.,0./)
       call phiNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.3,.false.)
       If ((teilchenOut(1)%ID.eq.pion).and.(teilchenOut(2)%ID.eq.pion).and.(teilchenOut(3)%ID.eq.nucleon)) then
          pipiN= pipiN+sigmatot
       else If ((teilchenOut(1)%ID.eq.phi).and.(teilchenOut(2)%ID.eq.nucleon)) then
          phiN=phiN+sigmatot
       else If ((teilchenOut(1)%ID.eq.P13)) then
          R_p13=R_p13+sigmatot
       else If ((teilchenOut(1)%ID.eq.pion).and.(teilchenOut(2)%ID.eq.nucleon)) then
          piN= piN+sigmatot
       end if
    end do
    write(*,*) 'Xsections at plab= 0.2 GeV'
    write(*,*) 'Simga for P13. :', r_p13/float(numTries)
    write(*,*) 'Simga for piN. :', piN/float(numTries)
    write(*,*) 'Simga for phi nucleon prod. :', phiN/float(numTries)
    write(*,*) 'Simga for piPiN prod. :', pipiN/float(numTries)




    teilchenIN(1)%Id=phi
    teilchenIN(1)%charge=chargeMeson
    teilchenIN(1)%mass=meson(phi)%mass

    teilchenIN(2)%Id=nucleon
    teilchenIN(2)%charge=-chargeNuk
    teilchenIN(2)%antiparticle=.true.
    teilchenIN(2)%mass=0.938
    teilchenIN(2)%momentum=(/teilchenIN(2)%mass,0.,0.,0./)
    teilchenIN(2)%velocity=(/0.,0.,0./)

    Open(301,file='PhiNucleonXsections_anti.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    mesonCharge=', chargeMeson
    Do i=1,250
       plab=i*0.01

       teilchenIN(1)%momentum=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

       srts=SQRT((SQRT(meson(phi)%mass**2+plab**2)+baryon(nucleon)%mass)**2-plab**2)
       momentumLRF=(/SQRT(meson(phi)%mass**2+plab**2)+baryon(nucleon)%mass ,plab,0.,0./)
       ekin=SQRT(meson(phi)%mass**2+plab**2)-meson(phi)%mass
       call phiNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.3,.false.)

       write(301,'(5F8.3)') plab,sigmaTot,sigmaElast
    end do
    close(301)



  end subroutine testPhiNuc

end program test
