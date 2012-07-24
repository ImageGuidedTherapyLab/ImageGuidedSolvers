#include "finclude/petscdef.h"
#include "parser_defaults.h"
!----------------------------------------------------------------------
!   module  - parameters for plotting
!
!   latest revision    - Jun 06
!
!   purpose            - store parameters for plotting
!----------------------------------------------------------------------
      module plot_params
      implicit none
! the save attibute allows for the plot number to be save from one
! iteration to the next
      save
      !              break up hppatch  field width of plot file #
      !              by PLOTFACT**3     (must be same as mrti data for vis)
      !                        |          |
      !                        |          |
      !                        |          |
      !                       \|/        \|/
      integer, parameter :: PLOTFACT=1,NFLDPLOT=4
!  ... data for plotting
      integer, private :: iii,jjj
      double precision,dimension(9,0:9*PLOTFACT),parameter::              &
      xsi =  RESHAPE(                                                     &
      (/((dble(jjj)/dble(iii*PLOTFACT) ,iii=1,9),jjj=0,9*PLOTFACT)/),     &
                                                   (/9,9*PLOTFACT+1/) )
      PetscInt:: NELMTOTAL,   &  !sum of number of elements total
                 NUNIQUE          !number of unique nodes
      character(len=32)::  FEMFILE ! fem visualization file name
      character(len=MAXLEN):: COMPWRITEFEM, & ! where the fem vis files are written
                              COMPWRITEMRI    ! where the mri vis files are written
      character(len=MAXLEN)::  FEMVIS_SIZEFILE ! fem vis file size file name
      PetscTruth::       PLOTGMV, &  ! .true. ==> write out gmv files
                         PLOTRAW, &  ! .true. ==> write out raw files
                         PLOTAVS     ! .true. ==> write out avs files
      character(len=10),parameter:: CZEROS='0000000000'
      character(len=6):: BYTEID
      !module variables
      integer,parameter:: IDFIELDS=6 ! size integers below to define the data
                                     ! type used for MPI_Type_struct
      type  gmvvertex
      sequence 
          real    :: x,y,z
          integer :: nodetype   !  0 = vertex node
                                !  1 = non-vertex node
                                !  2 = node used for plotting
          integer :: nodenumber !  # of vertex/non-vertex node in 
                                !    hp3d data structures
                                !  -1 if node for plotting
          integer :: iii,jjj,kkk,mdle! how to find the node to plot variables
      end type
      type(gmvvertex), allocatable , dimension(:) :: NODES

      integer,allocatable,dimension(:,:) :: CONNECTIVITIES, & ! cell connectivites
                                            CONNECTEDELEMS    ! store the # of 
                                                              ! elements attached
                                                              ! to a particular
                                                              ! unique node
      ! max NUMBER OF elements attached to a plotting NODE
      integer , parameter :: MAXREPEAT=60
      ! store information about how many repeate nodes there are
      integer,allocatable,dimension(:) :: REPEATEDNODES
      ! array to write out variable data
      real,allocatable,dimension(:,:) :: NODVAR
      ! buffer to store exact solution for verification problems
      real,allocatable,dimension(:) :: EXACTVAR
      ! NUMBER OF NODE VARIABLES FOR AVS BINARY FILES
      !  this # should be kept 
      !         SMALL      SMALL      SMALL      SMALL 
      !                                 to keep the data transfer fast
      integer , parameter :: NUMNODVARSAVS=3
      ! buffer to write avs files
      character,allocatable,dimension(:) :: AVSBUFFER
      integer :: IPOS, AVSBUFFERSIZE ! buffer byte position and size
      PetscTruth ::   BYTESWAP     ! PETSC_TRUE ==> convert endian
      ! mpi derived data types are integers in fortran
#define MPI_Datatype integer
      MPI_Datatype :: MPI_HEXas,MPI_GMVVERTex,MPI_CELLINfo,MPI_NODVar
      end module plot_params
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine init_plot_params(DDDAS_COMM_WORLD)
        use parse_ini
        use plot_params

        implicit none
#include "finclude/petsc.h"
        PetscMPIInt , intent(in) :: DDDAS_COMM_WORLD
        character(len=MAXLEN),parameter::default_ini='files/control.ini'
        logical :: istat ! I/O status flag
        PetscInt :: irank
        PetscErrorCode :: ierr
        PetscTruth, external :: islittleendian



        call mpi_comm_rank(DDDAS_COMM_WORLD,irank,ierr)

        ! all processors Try to read the INI file
        call openIniFile(default_ini, istat)
        !If unsuccessfully read,stop
        if( .not. istat ) then
          write(*,*) 'could not open', default_ini
          call abort
        endif

        ! location of file writing
        call getIniString('output','compwritefem',COMPWRITEFEM,'femvis/',istat)
        call getIniString('output','compwritemri',COMPWRITEMRI,'mrivis/',istat)

        ! read in parameters for output files
        ! endianness
        call getIniPetscTruth('output','BYTESWAP',BYTESWAP)
        if(BYTESWAP .eqv. PETSC_TRUE) write(*,*) "Byte SWAPPING!!!!!!!!!!!"
        if(islittleendian() .eqv. PETSC_TRUE ) then
          BYTEID="litend" ; if(BYTESWAP .eqv. PETSC_TRUE) BYTEID="bigend"
        else
          BYTEID="bigend" ; if(BYTESWAP .eqv. PETSC_TRUE) BYTEID="litend"
        endif
        ! file write control
        call getIniString('output','FEM_FILENAME',FEMFILE,'fem',istat)
        ! gmv parameters
        call getIniPetscTruth('output','FEMGMV',PLOTGMV)
        ! avs parameters
        call getIniPetscTruth('output','FEMAVS',PLOTAVS)
        ! raw parameters
        call getIniPetscTruth('output','FEMRAW',PLOTRAW)
        
        ! fem vis size file 
        call getIniString('avs','femvis_sizefile',FEMVIS_SIZEFILE,'ucdvis.size',istat)
        ! close file 
        call closeIniFile()
      if(irank.eq.0) then 
         call printpetscint( 'init_plot_params: [output] BYTESWAP='//char(0), BYTESWAP )
         call PetscPrintf(PETSC_COMM_SELF, &
                             'init_plot_params: [output] FEMFILE='//trim(FEMFILE)//char(10)//char(0),ierr)
         call printpetscint( 'init_plot_params: [output] PLOTGMV='//char(0), PLOTGMV )
         call printpetscint( 'init_plot_params: [output] PLOTAVS='//char(0), PLOTAVS )
         call printpetscint( 'init_plot_params: [output] PLOTRAW='//char(0), PLOTRAW )
!         the full power data structure is written before to qoi_?power*.dat
      endif
      end subroutine init_plot_params
