#include "finclude/petscdef.h"
#include "parser_defaults.h"
!     GLOBAL constants and global parameters
      module global_params
      PetscScalar,parameter :: ZERO=0.0d0 
      PetscScalar,parameter :: PI=3.14159265358979d0 
      PetscScalar,parameter :: TWODPI=0.63661977236758138d0   ! 2/pi
      PetscScalar,parameter :: L=20.0d0 ! length of bar in test problems
      PetscScalar,parameter :: R=1.0d0  ! rad of sphere in test problems
      PetscScalar           :: DirichletPenalty  ! penalty term Dirichlet conditions
      PetscTruth ::  TAO_MONITOR     ! PETSC_TRUE ==> write parameters
      PetscTruth ::  OPTIMIZE_mesh   ! Logical variables to control mesh
      PetscTruth ::  DualEqPplusOne  !    optimization by GOALS algorithm
      PetscTruth ::  OwnElemNodeVert ! Control domain decomposition
      PetscTruth ::  Control_Task    ! true ==> this task is the Control Task

      integer ::  ITERNO=0   ! counter for optimization iterations
      ! MAXSTEPS:  maximum # of time steps
      ! GroupID:   compute group ID
      PetscInt :: MAXSTEPS, GroupID

      character(len=32)::  VISUALASE_FILE ! visualase file name
      PetscScalar :: VISUALASE_MAX_POW  ! visualase maximum power 
      PetscScalar :: VISUALASE_MAX_TEMP ! when this maximum temperature detected
                                        ! shut off the visualase
      character(len=32)::  PROFILEID ! used for file identification /tmp 
      PetscTruth  :: SKIPTMPCHK ! TRUE ==> additional error checking

      end module global_params
!----------------------------------------------------------------------
!   function name   - getpenalty    (last revision: dec, 2009)
!   purpose         - get penalty term
!   arguments       - none
!----------------------------------------------------------------------
      function getpenalty()
        use global_params , only: DirichletPenalty
        implicit none
        PetscScalar :: getpenalty 
        getpenalty = DirichletPenalty
      end function getpenalty
!----------------------------------------------------------------------
!
!   routine name       - opfil
!
!   purpose            - routine opens up files
!
!   usage              - call opfil
!
!----------------------------------------------------------------------
      subroutine opfil(Location,ProfileID)
 
      use global_params, only: GroupID

      implicit none
#include "finclude/petsc.h"
#include "cinout.blk"

      character(len=MAXLEN),intent(in) :: Location ! local location of data
      character(len=32),intent(in)  :: ProfileID 
      character(len=32)  :: qoiID , procid , idchargroup
!      
      integer :: irank,ierr,NTIMERS,NCOMMS,NIO,NDIFF,NGRAPHOUT,NGRAPHIN

      write(idchargroup,*) GroupID
      qoiID="qoi_"//trim(adjustl(idchargroup))

      call mpi_comm_rank(PETSC_COMM_WORLD,irank,ierr)
      write(procid,*) irank
      NIN  = 10
      open(unit=NIN,file='files/input',                                      &
           form='formatted',access='sequential',status='unknown')
!
      NOUT = 11
!  ...parallel version appends process number to output file
      open(unit=NOUT,file='files/'//trim(adjustl(qoiID))//'out_'//           &
                  trim(adjustl(procid))//'.o'//trim(adjustl(ProfileID)),     &
                  form='formatted',access='sequential',status='unknown')
      NDUMP=12
!
      NHIST=13
!  ...parallel version appends process number to output file
      open(unit=NHIST,file='files/'//trim(adjustl(qoiID))//                  &
                                   'history'//trim(adjustl(procid)),         &
                  form='formatted',access='sequential',status='unknown')
!  ...file 'files/control' is open from main/read_control...
      NWE = 14
      open(unit=NWE,file='files/input_new',form='formatted',                 &
           access='sequential',status='unknown')
!
      NPROFB = 15
      open(unit=NPROFB,file='files/'//trim(adjustl(qoiID))//                 &
                                   'profilb',form='formatted',               &
           access='sequential',status='unknown')
!
      KIN = 16
      open(unit=KIN,file=trim(Location)//'/input_compact',                   &
            form='formatted',access='sequential',status='unknown')
!
      end subroutine opfil
!----------------------------------------------------------------------
!     initialize visualase data structures
!----------------------------------------------------------------------
      subroutine initialize_visualase(lmonitor,CntlTask,Opt_mesh,Idgrp)

      use parse_ini
      use global_params, only : VISUALASE_FILE,VISUALASE_MAX_POW,         &
                                VISUALASE_MAX_TEMP , PROFILEID,           &
                                DirichletPenalty, GroupID, Control_Task,  &
                                OPTIMIZE_mesh, TAO_MONITOR  ,SKIPTMPCHK

      implicit none
#include "finclude/petsc.h"

      PetscTruth  ,intent(in):: lmonitor , CntlTask, Opt_mesh
      PetscInt    ,intent(in):: Idgrp
      PetscInt :: petintdum
      integer  :: intdum
      PetscScalar :: petdbldum
      double precision  :: dbldum
      character(len=MAXLEN) :: inifile   ! INI filename
      character(len=MAXLEN) :: powerfile ! text file containing power vs time
      character(len=MAXLEN) :: compfilelocation ! local location of data
      logical :: istat ! I/O status flag
      PetscErrorCode :: ierr 
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !  make sure that data type sizes are consistent
      !  petsc ensures that PetscScalar and PetscInt at same when passing
      !  variables between C and Fortran
      !  BUT when interfacing with hp3d could be possible at some point
      !  that sizeof(PetscScalar) .ne. size(double precision)
      !  double precision is implicitly used everywhere in hp3d
      if(sizeof(petdbldum).ne.sizeof(dbldum))then
         write(*,*) "double precision/PetscScalar data type error"
         call abort
      endif
      if(sizeof(petintdum).ne.sizeof(intdum))then
         write(*,*) "integer/PetscInt data type error"
         call abort
      endif

      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      GroupID       = Idgrp    !set group id
      Control_Task  = CntlTask !set control task id
      OPTIMIZE_mesh = Opt_mesh !set mesh optimization flag
      TAO_MONITOR   = lmonitor !set tao_monitor
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      inifile  ="files/control.ini"
      !Try to read the INI file
      call openIniFile(inifile, istat)
      !If unsuccessfully read,stop
      if(.not. istat) then
        write(*,*) 'could not open', inifile
        stop
      end if
      ! must get visualase data before initializing mrti data structures
      ! visualase data filename of visualase file
      call getIniString('visualase','VISUALASE_FILE',VISUALASE_FILE,'mlat.dat',istat)
      !  maximum visualase power
      call getIniReal('visualase','VISUALASE_MAX_POW',VISUALASE_MAX_POW , 15.0d0,istat)

      !  maximum temperature allowed in thermal images above this 
      !  temperature must shut off visualase
      call getIniReal('visualase','VISUALASE_MAX_TEMP',VISUALASE_MAX_TEMP , 350.0d0,istat)
      ! close ini file

      ! mesh file location
      call getIniString('compexec','compfilelocation', compfilelocation,'files/',istat)

      call getIniString('compexec','profileID',PROFILEID,'',istat)
      call opfil(compfilelocation,PROFILEID);

      ! additional error checking 
      call getIniPetscTruth('compexec','SKIPTMPCHK',SKIPTMPCHK)

      ! penalty term
      call getIniReal('dirichlet','Penalty',DirichletPenalty,1.0d10,istat)

      ! close file 
      call closeIniFile()
      if(Control_Task .eqv. PETSC_TRUE) then 
         call PetscPrintf(PETSC_COMM_SELF, &
                               'visualase: VISUALASE_FILE='//VISUALASE_FILE//char(10)//char(0),ierr)
         call printpetscscalar('visualase: VISUALASE_MAX_POW='//char(0) ,VISUALASE_MAX_POW )
         call printpetscscalar('visualase: VISUALASE_MAX_TEMP='//char(0),VISUALASE_MAX_TEMP)
         call printpetscscalar('visualase: DirichletPenalty='//char(0),DirichletPenalty)
         call printpetscint(   'visualase: SKIPTMPCHK        ='//char(0),SKIPTMPCHK        )
      endif
      end subroutine initialize_visualase
!----------------------------------------------------------------------
!     return visualase maximum temperature to c++ code
!----------------------------------------------------------------------
      function get_visualase_max_temp()
        use global_params , ONLY: VISUALASE_MAX_TEMP  
        implicit none
        PetscScalar :: get_visualase_max_temp  
        get_visualase_max_temp = VISUALASE_MAX_TEMP 
      end function get_visualase_max_temp  
!----------------------------------------------------------------------
!     reset the iteration count before each optimization solve
!----------------------------------------------------------------------
      subroutine reset_iterno
        use global_params , only: ITERNO
        implicit none
        ITERNO=0
      end subroutine reset_iterno
!----------------------------------------------------------------------
!   routine name       - cross_product
!   purpose            - routine evaluates cross product of two
!                        vectors in R^3
!   arguments :
!     in:              - a(3),b(3)
!     out:             - c(3)
!----------------------------------------------------------------------
        subroutine cross_product(a,b,c)
        implicit none
        PetscScalar ::  a(3),b(3),c(3)
        c(1) =   a(2)*b(3) - a(3)*b(2)
        c(2) = - a(1)*b(3) + a(3)*b(1)
        c(3) =   a(1)*b(2) - a(2)*b(1)
        end subroutine cross_product
