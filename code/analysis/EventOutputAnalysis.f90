!*******************************************************************************
!****m* /EventOutputAnalysis
! NAME
! module EventOutputAnalysis
!
! PURPOSE
! This module provides routines for writing the outgoing particle vector
! to output files (split into events).
!
! Currently the following formats are supported:
! * "Les Houches" event files
! * "OSCAR 2013" event files, compatible with SMASH output
! * "Shanghai 2014" even files, used for the code-comparison project at the
!   Transport2014 workshop in Shanghai
!
! For a description of the Les Houches format, please refer to:
! * http://arxiv.org/abs/hep-ph/0609017
! * https://gibuu.hepforge.org/trac/wiki/LesHouches
!
! For a description of the OSCAR 2013 format, see:
! * http://phy.duke.edu/~jeb65/oscar2013
!
! For a description of the Shanghai 2014 format, see:
! * http://www.physics.sjtu.edu.cn/hic2014/node/12
!
! INPUTS
! (none)
!
! NOTES
! These analysis routines are independent of the specific initialization
! and should work for all event types.
!*******************************************************************************
module EventOutputAnalysis

  implicit none
  private

  !*****************************************************************************
  !****g* EventOutputAnalysis/WritePerturbativeParticles
  ! SOURCE
  !
  logical, save :: WritePerturbativeParticles = .false.
  !
  ! PURPOSE
  ! Flag to write out the perturbative particle vector to an output file.
  ! The switch 'EventFormat' determines which format is used.
  !*****************************************************************************

  !*****************************************************************************
  !****g* EventOutputAnalysis/WriteRealParticles
  ! SOURCE
  !
  logical, save :: WriteRealParticles = .false.
  !
  ! PURPOSE
  ! Flag to write out the real particle vector to an output file.
  ! The switch 'EventFormat' determines which format is used.
  !*****************************************************************************

  !*****************************************************************************
  !****g* EventOutputAnalysis/EventFormat
  ! SOURCE
  !
  integer, save :: EventFormat = 1
  !
  ! PURPOSE
  ! This switch determines the format of the event output files. Possible values:
  ! * 1 = Les Houches format (default)
  ! * 2 = OSCAR 2013 format
  ! * 3 = Shanghai 2014 format
  !
  ! NOTES
  ! For Les Houches, the output will be written to files called
  ! EventOutput.Pert.lhe and EventOutput.Real.lhe.
  ! For OSCAR, the output files are called EventOutput.Pert.oscar and
  ! EventOutput.Real.oscar.
  ! For Shanghai, the output files are called EventOutput.Pert.dat and
  ! EventOutput.Real.dat.
  !*****************************************************************************

  !*****************************************************************************
  !****g* EventOutputAnalysis/Interval
  ! SOURCE
  !
  integer, save :: Interval = 0
  !
  ! PURPOSE
  ! Interval for event output, i.e. number of timesteps after which output is
  ! written. If zero, only final output at the end of the time evolution is
  ! produced.
  !*****************************************************************************


  logical, save :: init = .true.
  integer, save :: nCall = 0

  public :: DoEventOutput

contains

  !*****************************************************************************
  !****s* EventOutputAnalysis/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Read namelist 'EventOutput' from jobcard.
  !*****************************************************************************
  subroutine initInput
    use output, only: Write_ReadingInput

    !***************************************************************************
    !****n* EventOutputAnalysis/EventOutput
    ! NAME
    ! namelist /EventOutput/
    ! PURPOSE
    ! Namelist for EventOutput includes:
    ! * WritePerturbativeParticles
    ! * WriteRealParticles
    ! * EventFormat
    ! * Interval
    !***************************************************************************
    namelist /EventOutput/ WritePerturbativeParticles, WriteRealParticles, EventFormat, Interval

    integer :: ios

    call Write_ReadingInput('EventOutput',0)
    rewind(5)
    read(5,nml=EventOutput,iostat=ios)

    call Write_ReadingInput('EventOutput',0,ios)

    write(*,*) 'Event output of final particles (real,pert): ', &
               WriteRealParticles, WritePerturbativeParticles
    write(*,*) 'Event-output format: ', EventFormat
    write(*,*) 'Interval: ', Interval
    call Write_ReadingInput('EventOutput',1)

    init = .false.

  end subroutine initInput


  !*****************************************************************************
  !****s* EventOutputAnalysis/DoEventOutput
  ! NAME
  ! subroutine DoEventOutput(realPart, pertPart)
  ! PURPOSE
  ! Do the actual writing out, if desired (as indicated in namelist).
  !*****************************************************************************
  subroutine DoEventOutput(realPart, pertPart, timestep)
    use particleDefinition
    use inputGeneral, only: numTimeSteps

    type(particle), intent(in), dimension(:,:) :: realPart, pertPart
    integer, intent(in) :: timestep

    logical :: doOutput

    if (init) call initInput
    if (.not. WritePerturbativeParticles .and. .not. WriteRealParticles) return

    if (Interval>0) then
      ! output at fixed interval
      doOutput = (mod(timestep,Interval)==0)
    else
      ! only final output
      doOutput = (timestep>numTimeSteps)
    end if

    if (.not. doOutput) return

    nCall = nCall+1

    select case (EventFormat)
    case (1)
      call DoOutput_LesHouches(realPart, pertPart)
    case (2)
      call DoOutput_Oscar(realPart, pertPart)
    case (3)
      call DoOutput_Shanghai(realPart, pertPart, timestep)
    case default
      write(*,*) "Bad EventFormat in EventOutputAnalysis: ", EventFormat
      stop
    end select

  end subroutine DoEventOutput


  !*****************************************************************************
  !****s* EventOutputAnalysis/DoOutput_LesHouches
  ! NAME
  ! subroutine DoOutput_LesHouches(realPart, pertPart)
  ! PURPOSE
  ! Write out event files in Les Houches format.
  !*****************************************************************************
  subroutine DoOutput_LesHouches(realPart, pertPart)
    use EventOutput, only: LHOutputFile
    use particleDefinition

    type(particle), intent(in), dimension(:,:) :: realPart, pertPart

    type(LHOutputFile) :: LH_pert
    type(LHOutputFile) :: LH_real

    !***************************************************************************
    !****o* EventOutputAnalysis/EventOutput.Pert.lhe
    ! NAME
    ! file EventOutput.Pert.lhe
    ! PURPOSE
    ! Contains all perturbative particles of a given run in Les Hoches format.
    ! Can be enabled by the switch WritePerturbativeParticles.
    ! For documentation of the file format see https://gibuu.hepforge.org/trac/wiki/LesHouches.
    ! For each subsequent run a separate file will be produced:
    !  * EventOutput.Pert.00000001.lhe
    !  * EventOutput.Pert.00000002.lhe
    !  * etc
    !***************************************************************************
    if (WritePerturbativeParticles) then
      call LH_pert%open(.true., nCall)
      call LH_pert%write_pert(pertPart)
      call LH_pert%close()
    end if

    !***************************************************************************
    !****o* EventOutputAnalysis/EventOutput.Real.lhe
    ! NAME
    ! file EventOutput.Real.lhe
    ! PURPOSE
    ! Contains all real particles of a given run in Les Hoches format.
    ! Can be enabled by the switch WriteRealParticles.
    ! For documentation of the file format see https://gibuu.hepforge.org/trac/wiki/LesHouches.
    ! For each subsequent run a separate file will be produced:
    !  * EventOutput.Real.00000001.lhe
    !  * EventOutput.Real.00000002.lhe
    !  * etc
    !***************************************************************************
    if (WriteRealParticles) then
      call LH_real%open(.false., nCall)
      call LH_real%write_real(realPart)
      call LH_real%close()
    end if

  end subroutine DoOutput_LesHouches


  !*****************************************************************************
  !****s* EventOutputAnalysis/DoOutput_Oscar
  ! NAME
  ! subroutine DoOutput_Oscar(realPart, pertPart)
  ! PURPOSE
  ! Write out event files in OSCAR 2013 format.
  !*****************************************************************************
  subroutine DoOutput_Oscar(realPart, pertPart)
    use EventOutput, only: OscarOutputFile
    use particleDefinition

    type(particle), intent(in), dimension(:,:) :: realPart, pertPart

    type(OscarOutputFile), save :: Oscar_pert
    type(OscarOutputFile), save :: Oscar_real

    !***************************************************************************
    !****o* EventOutputAnalysis/EventOutput.Pert.oscar
    ! NAME
    ! file EventOutput.Pert.oscar
    ! PURPOSE
    ! Contains all perturbative particles of a given run in OSCAR 2003 format.
    ! Can be enabled by the switch WritePerturbativeParticles.
    ! The data from all subsequent runs will be written into a single file.
    !***************************************************************************
    if (WritePerturbativeParticles) then
      call Oscar_pert%open(.true., nCall)
      call Oscar_pert%write_pert(pertPart)
      call Oscar_pert%close()
    end if

    !***************************************************************************
    !****o* EventOutputAnalysis/EventOutput.Real.oscar
    ! NAME
    ! file EventOutput.Real.oscar
    ! PURPOSE
    ! Contains all real particles of a given run in OSCAR 2003 format.
    ! Can be enabled by the switch WriteRealParticles.
    ! The data from all subsequent runs will be written into a single file.
    !***************************************************************************
    if (WriteRealParticles) then
      call Oscar_real%open(.false., nCall)
      call Oscar_real%write_real(realPart)
      call Oscar_real%close()
    end if

  end subroutine DoOutput_Oscar


  !*****************************************************************************
  !****s* EventOutputAnalysis/DoOutput_Shanghai
  ! NAME
  ! subroutine DoOutput_Shanghai(realPart, pertPart)
  ! PURPOSE
  ! Write out event files in Shanghai 2014 format.
  !*****************************************************************************
  subroutine DoOutput_Shanghai(realPart, pertPart, timestep)
    use EventOutput, only: ShanghaiOutputFile
    use particleDefinition

    type(particle), intent(in), dimension(:,:) :: realPart, pertPart
    integer, intent(in) :: timestep

    type(ShanghaiOutputFile), save :: Shanghai_pert, Shanghai_real

    !***************************************************************************
    !****o* EventOutputAnalysis/EventOutput.Pert.dat
    ! NAME
    ! file EventOutput.Pert.dat
    ! PURPOSE
    ! Contains all perturbative particles of a given run in Shanghai format.
    ! Can be enabled by the switch WritePerturbativeParticles.
    ! The data from all subsequent runs will be written into a single file.
    !***************************************************************************
    if (WritePerturbativeParticles) then
      call Shanghai_pert%open(.true., timestep)
      call Shanghai_pert%write_pert(pertPart)
      call Shanghai_pert%close()
    end if

    !***************************************************************************
    !****o* EventOutputAnalysis/EventOutput.Real.dat
    ! NAME
    ! file EventOutput.Real.dat
    ! PURPOSE
    ! Contains all real particles of a given run in Shanghai format.
    ! Can be enabled by the switch WriteRealParticles.
    ! The data from all subsequent runs will be written into a single file.
    !***************************************************************************
    if (WriteRealParticles) then
      call Shanghai_real%open(.false., timestep)
      call Shanghai_real%write_real(realPart)
      call Shanghai_real%close()
    end if

  end subroutine DoOutput_Shanghai


end module EventOutputAnalysis
