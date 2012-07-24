// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// C include files 
#include <sys/types.h>
#include <sys/stat.h>

// Basic include files
#include "mesh_refinement.h"
#include "explicit_system.h"
#include "exodusII_io.h"
#include "exodusII_io_helper.h"
#include "numeric_vector.h"
#include "petsc_vector.h"

// Some (older) compilers do not offer full stream
// functionality, OStringStream works around this.
#include "o_string_stream.h"

// local class definition
#include "femInterface.h"
#include "tttkUtilities.h"

// declare static vars
bool                    FiniteElementInterface::m_instanceFlag = false;
FiniteElementInterface* FiniteElementInterface::m_single       = NULL;
Mesh*                   FiniteElementInterface::m_mesh         = NULL;
EquationSystems*        FiniteElementInterface::m_eqn_systems  = NULL;
LibMeshInit*            FiniteElementInterface::m_init         = NULL;
GetPot*                 FiniteElementInterface::m_controlfile  = NULL;
Imaging*                FiniteElementInterface::m_images       = NULL;
ExodusII_IO*            FiniteElementInterface::m_exodus_transient_file = NULL;

// profiling
std::vector<PetscInt>  FiniteElementInterface::logstages;
std::vector<PetscInt>  FiniteElementInterface::logevents;

// used for file/job ID
std::string FiniteElementInterface::profileID;


/** call private constructor if not initialized  */
#undef __FUNCT__
#define __FUNCT__ "FiniteElementInterface::getInstance"
FiniteElementInterface* FiniteElementInterface::getInstance()
{
    if(! m_instanceFlag)
    {
        m_single = new FiniteElementInterface();
        m_instanceFlag = true;
        return m_single;
    }
    else
    {
        return m_single;
    }
}

#undef __FUNCT__
#define __FUNCT__ "FiniteElementInterfaceGetInstance"
FiniteElementInterface* FiniteElementInterfaceGetInstance()
{
  return FiniteElementInterface::getInstance();
}


#undef __FUNCT__
#define __FUNCT__ "FiniteElementInterface::~FiniteElementInterface"
FiniteElementInterface::~FiniteElementInterface()
{ // clean up pointers

  std::cout << "cleaning up exodus " << std::endl;
  if(m_exodus_transient_file) delete m_exodus_transient_file;
  m_exodus_transient_file = NULL;

  std::cout << "cleaning up mesh " << std::endl;
  if(m_mesh)        delete m_mesh;
  m_mesh = NULL;

  std::cout << "cleaning up equation systems" << std::endl;
  if(m_eqn_systems) delete m_eqn_systems;
  m_eqn_systems = NULL;

  std::cout << "cleaning up control file " << std::endl;
  if(m_controlfile) delete m_controlfile;
  m_controlfile = NULL;

  std::cout << "not cleaning up libmesh due to seg faults " << std::endl;
  //if(m_init)        delete m_init;
  m_instanceFlag = false;
}

#undef __FUNCT__
#define __FUNCT__ "FiniteElementInterface::CleanUplibMesh"
PetscErrorCode FiniteElementInterface::CleanUplibMesh()
{ // clean up pointers
  PetscFunctionBegin; 
  std::cout << "cleaning up exodus " << std::endl;
  if(m_exodus_transient_file) delete m_exodus_transient_file;
  m_exodus_transient_file = NULL;

  std::cout << "cleaning up mesh " << std::endl;
  if(m_mesh)        delete m_mesh;
  m_mesh = NULL;

  std::cout << "cleaning up equation systems" << std::endl;
  if(m_eqn_systems) delete m_eqn_systems;
  m_eqn_systems = NULL;

  std::cout << "cleaning up control file " << std::endl;
  if(m_controlfile) delete m_controlfile;
  m_controlfile = NULL;

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- 
     Setup the majority of the variables
     for this class after the debugger is setup
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FiniteElementInterface::Setup"
PetscErrorCode FiniteElementInterface::Setup(GetPot  &controlfile)
{
  PetscFunctionBegin;

  // make sure the output directory exist, using the cross
  // platform tools: itksys::SystemTools. In this case we select to create
  // the directory if it does not exist yet.
  //
  // \index{itksys!SystemTools}
  // \index{itksys!MakeDirectory}
  // \index{SystemTools}
  // \index{SystemTools!MakeDirectory}
  // \index{MakeDirectory!SystemTools}
  // \index{MakeDirectory!itksys}
  itksys::SystemTools::MakeDirectory( "femvis" );

  // profile id
  profileID = controlfile("compexec/profileid","");

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- 
     Setup profiling data
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FiniteElementInterface::SetupProfile"
void FiniteElementInterface::SetupProfile()
{
 PetscFunctionBegin;

  const int NSTAGESIZE = 7;
  const int NEVENTSIZE = 33;
  logstages.resize(NSTAGESIZE,0);
  logevents.resize(NEVENTSIZE,0);

  //profiling
  /*profiling */
#if PETSC_VERSION_LESS_THAN(3,0,0)
  PetscLogStageRegister(&logstages[0] ,"Initialization"      );
  PetscLogStageRegister(&logstages[1] ,"Mesh Generation"     );
  PetscLogStageRegister(&logstages[2] ,"function evaluation" );
  PetscLogStageRegister(&logstages[3] ,"gradient evaluation" );
  PetscLogStageRegister(&logstages[4] ,"hessian  evaluation" );
  PetscLogStageRegister(&logstages[5] ,"Prediction    "      );
  PetscLogStageRegister(&logstages[6] ,"State Update  "      );

  PetscLogEventRegister(&logevents[0] ,"tao param xfer",PETSC_NULL);
  PetscLogEventRegister(&logevents[1] ,"assemble fnc  ",PETSC_NULL);
  PetscLogEventRegister(&logevents[2] ,"assemble jac  ",PETSC_NULL);
  PetscLogEventRegister(&logevents[3] ,"read disk     ",PETSC_NULL);
  PetscLogEventRegister(&logevents[4] ,"adj assemble  ",PETSC_NULL);
  PetscLogEventRegister(&logevents[5] ,"libMesh Solve ",PETSC_NULL);
  PetscLogEventRegister(&logevents[6] ,"gradient comp ",PETSC_NULL);
  PetscLogEventRegister(&logevents[7] ,"elemfnc setup ",PETSC_NULL);
  PetscLogEventRegister(&logevents[8] ,"elemfnc assble",PETSC_NULL);
  PetscLogEventRegister(&logevents[9] ,"elemfnc bndry ",PETSC_NULL);
  PetscLogEventRegister(&logevents[10],"elemjac setup ",PETSC_NULL);
  PetscLogEventRegister(&logevents[11],"elemjac assble",PETSC_NULL);
  PetscLogEventRegister(&logevents[12],"elemjac bndry ",PETSC_NULL);
  PetscLogEventRegister(&logevents[13],"elem diff     ",PETSC_NULL);
  PetscLogEventRegister(&logevents[14],"elem reac     ",PETSC_NULL);
  PetscLogEventRegister(&logevents[15],"data transfer ",PETSC_NULL);
  PetscLogEventRegister(&logevents[16],"data write    ",PETSC_NULL);
  PetscLogEventRegister(&logevents[17],"fnc evaluation",PETSC_NULL);
  PetscLogEventRegister(&logevents[18],"filters       ",PETSC_NULL);
  PetscLogEventRegister(&logevents[19],"field setup   ",PETSC_NULL);
  PetscLogEventRegister(&logevents[20],"init libMesh  ",PETSC_NULL);
  PetscLogEventRegister(&logevents[21],"setup libMesh ",PETSC_NULL);
  PetscLogEventRegister(&logevents[22],"init variable ",PETSC_NULL);
  PetscLogEventRegister(&logevents[23],"qoi eval      ",PETSC_NULL);
  PetscLogEventRegister(&logevents[24],"setup parallel",PETSC_NULL);
  PetscLogEventRegister(&logevents[25],"sort and merge",PETSC_NULL);
  PetscLogEventRegister(&logevents[26],"setup serial  ",PETSC_NULL);
  PetscLogEventRegister(&logevents[27],"eval mesh fcn ",PETSC_NULL);
  PetscLogEventRegister(&logevents[28],"find neighbour",PETSC_NULL);
  PetscLogEventRegister(&logevents[29],"set  neighbour",PETSC_NULL);
  PetscLogEventRegister(&logevents[30],"genorientation",PETSC_NULL);
  PetscLogEventRegister(&logevents[31],"update  gdof  ",PETSC_NULL);
  PetscLogEventRegister(&logevents[32],"set coordinate",PETSC_NULL);
#else
  PetscLogStageRegister("Initialization"      ,&logstages[0] );
  PetscLogStageRegister("Mesh Generation"     ,&logstages[1] );
  PetscLogStageRegister("function evaluation" ,&logstages[2] );
  PetscLogStageRegister("gradient evaluation" ,&logstages[3] );
  PetscLogStageRegister("hessian  evaluation" ,&logstages[4] );
  PetscLogStageRegister("Prediction    "      ,&logstages[5] );
  PetscLogStageRegister("State Update  "      ,&logstages[6] );

  PetscLogEventRegister("tao param xfer",PETSC_VIEWER_COOKIE,&logevents[0] );
  PetscLogEventRegister("assemble fnc  ",PETSC_VIEWER_COOKIE,&logevents[1] );
  PetscLogEventRegister("assemble jac  ",PETSC_VIEWER_COOKIE,&logevents[2] );
  PetscLogEventRegister("read disk     ",PETSC_VIEWER_COOKIE,&logevents[3] );
  PetscLogEventRegister("adj assemble  ",PETSC_VIEWER_COOKIE,&logevents[4] );
  PetscLogEventRegister("libMesh Solve ",PETSC_VIEWER_COOKIE,&logevents[5] );
  PetscLogEventRegister("gradient comp ",PETSC_VIEWER_COOKIE,&logevents[6] );
  PetscLogEventRegister("elemfnc setup ",PETSC_VIEWER_COOKIE,&logevents[7] );
  PetscLogEventRegister("elemfnc assble",PETSC_VIEWER_COOKIE,&logevents[8] );
  PetscLogEventRegister("elemfnc bndry ",PETSC_VIEWER_COOKIE,&logevents[9] );
  PetscLogEventRegister("elemjac setup ",PETSC_VIEWER_COOKIE,&logevents[10]);
  PetscLogEventRegister("elemjac assble",PETSC_VIEWER_COOKIE,&logevents[11]);
  PetscLogEventRegister("elemjac bndry ",PETSC_VIEWER_COOKIE,&logevents[12]);
  PetscLogEventRegister("elem diff     ",PETSC_VIEWER_COOKIE,&logevents[13]);
  PetscLogEventRegister("elem reac     ",PETSC_VIEWER_COOKIE,&logevents[14]);
  PetscLogEventRegister("data transfer ",PETSC_VIEWER_COOKIE,&logevents[15]);
  PetscLogEventRegister("data write    ",PETSC_VIEWER_COOKIE,&logevents[16]);
  PetscLogEventRegister("fnc evaluation",PETSC_VIEWER_COOKIE,&logevents[17]);
  PetscLogEventRegister("filters       ",PETSC_VIEWER_COOKIE,&logevents[18]);
  PetscLogEventRegister("field setup   ",PETSC_VIEWER_COOKIE,&logevents[19]);
  PetscLogEventRegister("init libMesh  ",PETSC_VIEWER_COOKIE,&logevents[20]);
  PetscLogEventRegister("setup libMesh ",PETSC_VIEWER_COOKIE,&logevents[21]);
  PetscLogEventRegister("init variable ",PETSC_VIEWER_COOKIE,&logevents[22]);
  PetscLogEventRegister("qoi eval      ",PETSC_VIEWER_COOKIE,&logevents[23]);
  PetscLogEventRegister("setup parallel",PETSC_VIEWER_COOKIE,&logevents[24]);
  PetscLogEventRegister("sort and merge",PETSC_VIEWER_COOKIE,&logevents[25]);
  PetscLogEventRegister("setup serial  ",PETSC_VIEWER_COOKIE,&logevents[26]);
  PetscLogEventRegister("eval mesh fcn ",PETSC_VIEWER_COOKIE,&logevents[27]);
  PetscLogEventRegister("find neighbour",PETSC_VIEWER_COOKIE,&logevents[28]);
  PetscLogEventRegister("set  neighbour",PETSC_VIEWER_COOKIE,&logevents[29]);
  PetscLogEventRegister("genorientation",PETSC_VIEWER_COOKIE,&logevents[30]);
  PetscLogEventRegister("update  gdof  ",PETSC_VIEWER_COOKIE,&logevents[31]);
  PetscLogEventRegister("set coordinate",PETSC_VIEWER_COOKIE,&logevents[32]);
#endif


 PetscFunctionReturnVoid();
}

/* -------------------------------------------------------------------
     initialize imaging
   ------------------------------------------------------------------- */
PetscErrorCode FiniteElementInterface::SetupImaging(int inputmethod)
{
  PetscFunctionBegin; 

  // pass controlfile 
  if(!this->m_controlfile)
    {
      std::cerr << "controlfile found" << std::endl; libmesh_error();
    }

  GetPot &controlfile = *this->m_controlfile;

  
  switch(static_cast<ImageAcquisitionType> (inputmethod))
   {
    case NO_IMAGING: // no imaging
      this->m_images = NULL;
      break;
    case DICOM2D: // 2d dicom series return snr and temperature as raw data
      this->m_images = new ImagingComplex(controlfile);
      std::cout << "conventional complex difference snr temperature maps..."
                << std::endl ;  
      break;
    case VOLUME3D:  // 3d Volume image
      this->m_images = new ImagingPreProcessed(controlfile);
      std::cout << "expecting preprocess temperature/snr maps..."
                << std::endl ;  
      break;
    case EIGHT_ECHO_CSI:   // CSI w/ 8 echos
      this->m_images = new ImagingChemicalShift< 8 >(controlfile);    
      break;
    case SIXTEEN_ECHO_CSI: // CSI w/ 16 echos
      this->m_images = new ImagingChemicalShift< 16 >(controlfile);
      break;
    default: 
      std::cout << "unknown inputmethod "<< inputmethod 
                << std::endl << std::flush; abort();
   }
 // get basic header info
 if(this->m_images) this->m_images->GetHeaderData
                      (static_cast<ImageAcquisitionType> (inputmethod));
 PetscFunctionReturn(0);
}
