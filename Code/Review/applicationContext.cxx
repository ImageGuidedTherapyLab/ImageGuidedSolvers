#include <vector>
#include "tao.h" // petsc solvers w/ Fortran Interface
#include "libmesh.h" // petsc solvers w/ Fortran Interface
#include "exodusII_io_helper.h" // petsc solvers w/ Fortran Interface
#include "getpot.h" // file parser
#include "equation_systems.h"
#include "nonlinear_implicit_system.h"
#include "mesh.h"
#include "diff_solver.h"
#include "euler_solver.h"
#include "steady_solver.h"
#include "petsc_diff_solver.h"

// local
#include "applicationContext.h"
#include "pennesVerification.h"
#include "pennesInverseSystem.h"
#include "KalmanFilter.h" 
#include "Utilities.h"
#include "parser_defaults.h"

// declare static vars
bool AppSolve::m_instanceFlag = false;
AppSolve* AppSolve::m_single = NULL;
std::vector<qoiBaseClass*> AppSolve::m_optimizerQOIS; 
TransientFEMSystem* AppSolve::m_pdeSolver = NULL;
PetscInt   AppSolve::Total_Nopt = 0; 
PetscInt   AppSolve::MAXIDEAL = -100000; // initialize maximum ideal step
PetscInt   AppSolve::MINIDEAL =  100000; // initialize minimum ideal step
PetscInt    AppSolve::m_NumQOI = 0; 
PetscInt    AppSolve::MRTI_nzero, 
            AppSolve::MRTI_ntime, // thermal image time bounds
            AppSolve::NEXACT, // thermal image time bounds
            AppSolve::MAXSTEPS,    // maximum number of time steps
            AppSolve::m_Refine,      // refinement level flag
            AppSolve::TIMELAG_IC=0, // time-lag for the initial condition
            AppSolve::m_ISTEPS_PER_IDEAL, // fem time steps per ideal timestep
            AppSolve::m_WriteInterval, // number of times to skip before plotting
            AppSolve::NfieldElem=1,  // number of field elements for a single param
            AppSolve::n_block=1;  // # of mesh blocks

PetscInt AppSolve::m_NoptSteps   = 0 ;
PetscInt AppSolve::m_IDEAL_NZERO = 0 ;
PetscInt AppSolve::m_NOFFSET     = 0 ;
PetscInt AppSolve::m_IDEAL_NTIME = 0 ;
PetscInt AppSolve::m_NUMIDEAL    = 0 ;

PetscInt AppSolve::controlfileID  = 0;
PetscInt AppSolve::outputfileID   = 0;
PetscInt AppSolve::m_checkpoint     = 0;
PetscInt AppSolve::m_checkpointiter = 0;

PetscScalar AppSolve::IDEAL_DT, // ideal deltat
            AppSolve::m_MaxTime,
            AppSolve::FEM_DT; // ideal deltat
PetscTruth  AppSolve::m_optimizeMesh          = PETSC_FALSE;
PetscTruth  AppSolve::m_COMPUTEFINALOBJECTIVE = PETSC_FALSE;
PetscTruth  AppSolve::m_IC_MRTI               = PETSC_FALSE;
PetscTruth  AppSolve::m_PLOTITER              = PETSC_FALSE;
PetscTruth  AppSolve::m_Control_Task          = PETSC_FALSE;

// used for restarts
PetscTruth  AppSolve::m_IDrestart             = PETSC_FALSE;
std::string AppSolve::restartFile;

// used for file/job ID
std::string AppSolve::profileID;
std::string AppSolve::localDisk;

PetscInt  AppSolve::GroupID      = 0, //store GroupID
          AppSolve::IDOPT        = 0, //initialize
          AppSolve::ISTEP        = 0; //initialize

// profiling
std::vector<PetscInt>  AppSolve::logstages;
std::vector<PetscInt>  AppSolve::logevents;
/* ------------------------------------------------------------------- 
     Default Constructor for AppSolve class
      minimal variables setup before debugger setup
   ------------------------------------------------------------------- */
// constructor
#undef __FUNCT__
#define __FUNCT__ "AppSolve::AppSolve"
/** minimal default initializations */
AppSolve::AppSolve()
{
}
/** call private constructor if not initialized  */
AppSolve* AppSolve::getInstance()
{
    if(! m_instanceFlag)
    {
        m_single = new AppSolve();
        m_instanceFlag = true;
        return m_single;
    }
    else
    {
        return m_single;
    }
}

void AppSolve::printSelf(std::ostream& os)
{
 PetscFunctionBegin; 
   
 std::vector<PetscScalar>::iterator TimeIter;
 // for(TimeIter=GenInfo.TIME.begin();TimeIter!=GenInfo.TIME.end(); TimeIter++)
 //      printf( "dddas: TIME[%d]=%f\n", 
 //               distance(GenInfo.TIME.begin(),TimeIter),*TimeIter);
  // echo data
  if( !libMesh::processor_id() )
   {
    os << "\nAppSolve: MINIDEAL        ="<<MINIDEAL
       << "\nAppSolve: optimizeMesh    ="<<m_optimizeMesh
       << "\nAppSolve: IDrestart       ="<<m_IDrestart  
       << "\nAppSolve: MAXIDEAL        ="<<MAXIDEAL
       << "\nAppSolve: MRTI_nzero      ="<<MRTI_nzero
       << "\nAppSolve: MRTI_ntime      ="<<MRTI_ntime
       << "\nAppSolve: NEXACT          ="<<NEXACT
       << "\nAppSolve: MAXSTEPS        ="<<MAXSTEPS
       << "\nAppSolve: Refine          ="<<m_Refine
       << "\nAppSolve: ISTEPS_PER_IDEAL="<<m_ISTEPS_PER_IDEAL
       << "\nAppSolve: WriteInterval   ="<<m_WriteInterval 
       << "\nAppSolve: IDEAL_DT        ="<<IDEAL_DT
       << "\nAppSolve: FEM_DT          ="<<FEM_DT
       << "\nAppSolve: profileID       ="<<profileID  
       << "\nAppSolve: localDisk       ="<<localDisk  
       << "\nAppSolve: n_block         ="<<n_block   
       << "\nAppSolve: NfieldElem      ="<<NfieldElem
       << "\nAppSolve: PLOTITER        ="<<m_PLOTITER      
       << "\nAppSolve: checkpoint      ="<<m_checkpoint    
       << "\nAppSolve: checkpointiter  ="<<m_checkpointiter     << std::endl;
    os << "AppSolve: NoptSteps         ="<<m_NoptSteps        << std::endl;
    os << "AppSolve: IDEAL_NZERO       ="<<m_IDEAL_NZERO      << std::endl;
    os << "AppSolve: NOFFSET           ="<<m_NOFFSET          << std::endl;
    os << "AppSolve: IDEAL_NTIME       ="<<m_IDEAL_NTIME      << std::endl;
    os << "AppSolve: NUMIDEAL          ="<<m_NUMIDEAL         << std::endl;
    os << "AppSolve: MaxTime           ="<<m_MaxTime          << std::endl;
    os << "AppSolve: FinalOpt_Ntime()  ="<<m_single->FinalOpt_Ntime()<< std::endl;
    os << "AppSolve: FinalOpt_Nzero()  ="<<m_single->FinalOpt_Nzero()<< std::endl;
   } 
 PetscFunctionReturnVoid(); 
}
/* ------------------------------------------------------------------- 
     Setup the majority of the variables
     for this class after the debugger is setup
   ------------------------------------------------------------------- */
PetscErrorCode AppSolve::Setup(GetPot  &controlfile)
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

  try
   {
    ExodusII_IO_Helper exio_helper;
    // Open the exodus file, if possible
    exio_helper.open(controlfile("compexec/meshdata","./mesh.e"));
    exio_helper.read_header();
    n_block = exio_helper.get_num_elem_blk();
    std::cout << "blocks " << n_block << std::endl ;
   }
  catch(...) // default to one
   {
     std::cout << "not an exodus file" << std::endl;
     n_block=1;
   }

  // minimum/maximum mrti time instance
  MRTI_nzero = controlfile("mrti/nzero",       0      );
  MRTI_ntime = controlfile("mrti/ntime",dfltmrti_ntime);
  // verification problem id
  NEXACT = controlfile("method/nexact",0 );

  /* Get timestep information about this Computational Group*/
  // number of optimization steps
  m_NoptSteps  = controlfile("timestep/noptsteps",1);
  // increment of IDEAL_NZERO at each optimization solve
  m_NOFFSET    = controlfile("timestep/noffset",dfltnoffset);  
  // lower bound of optimization window
  m_IDEAL_NZERO= controlfile("timestep/ideal_nzero_init",dfltideal_nzero);
  // increment of IDEAL_NTIME at each optimization solve
  m_NUMIDEAL   = controlfile("timestep/numideal",dfltnumideal);
  // upper bound of optimization window
  m_IDEAL_NTIME= controlfile("timestep/ideal_ntime_init",dfltideal_ntime);

  // domain decomposition control
  if(controlfile("dom_decomp/test_only",false)) NEXACT=-1; 

  // plot iteration and checkpoint control
  m_PLOTITER = controlfile("compexec/plotiter" ,  false) ? PETSC_TRUE:PETSC_FALSE;
  // default is to plot after five iterations 
  m_checkpointiter = controlfile("compexec/checkpointiter",  5);

  // restart from previous solution
  m_IDrestart= controlfile("compexec/restart" ,  false) ? PETSC_TRUE:PETSC_FALSE;

  // get consitutive data from disk
  restartFile = controlfile("compexec/restartfile","femvis/restartdata.e");

  // plot file write control
  m_ISTEPS_PER_IDEAL=controlfile("timestep/isteps_per_ideal",dfltistepsperideal);
  m_WriteInterval = controlfile("output/writeinterval",m_ISTEPS_PER_IDEAL);

  // time step per ideal time instance
  IDEAL_DT  = controlfile("timestep/ideal_dt", dfltideal_dt);

  // FEM time stepping...
  FEM_DT = IDEAL_DT/m_ISTEPS_PER_IDEAL;

  // default is no refinement for the mesh
  m_Refine= controlfile("compexec/meshrefine",0);

  // profile id
  profileID = controlfile("compexec/profileid","");

  // local disk storage 
  localDisk = controlfile("compexec/localdisk","/tmp");


  // error checking to enforce that start time is less than the finish time
  if(m_NOFFSET > m_NUMIDEAL) m_NOFFSET = m_NUMIDEAL;

  // error checking to enforce that start time is not more than the finish time
  if(m_IDEAL_NZERO > m_IDEAL_NTIME) m_IDEAL_NTIME=m_IDEAL_NZERO;

  /* when using the thermal images as the initial condition must be sure
          NOT to request and time instance of the thermal image that is
          out of the array bounds [0, MRTI_ntime] */
  if(m_IC_MRTI){
     if(m_IDEAL_NZERO + m_NOFFSET * (m_NoptSteps-1) 
                      - TIMELAG_IC > MRTI_ntime){
        if(!m_NOFFSET){
           printf("optvars.cpp: image initial condition input error \n");
           printf(" IDEAL_NZERO+NOFFSET*(NoptSteps-1)-TIMELAG_IC>MRTI_ntime\n");
           printf(" NOFFSET = %d\n",m_NOFFSET); abort();
        }
        //FLOOR property of int division ensures 
        //final time bounded by MRTI_ntime
        m_NoptSteps = ( MRTI_ntime - m_IDEAL_NZERO 
                                   + TIMELAG_IC) / m_NOFFSET + 1;
     }
  }

 m_MaxTime = m_single->setupTime(m_single->FinalOpt_Ntime(), 
                                 m_single->FinalOpt_Nzero(),1);

 PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- 
     Setup profiling data
   ------------------------------------------------------------------- */
void AppSolve::setupProfile()
{
 PetscFunctionBegin;

#include "profile_params.h"
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
     Setup basic time info
   ------------------------------------------------------------------- */
PetscScalar AppSolve::setupTime(const PetscInt FinalOptNtime,
                                const PetscInt IdealNzero,
                                const PetscInt MaxRefine)
{
 PetscFunctionBegin;

 /* get group-wise maximum and minimum of the last ideal step */
 MAXIDEAL=MAXIDEAL < FinalOptNtime ?  
                     FinalOptNtime : MAXIDEAL ;
 MINIDEAL=MINIDEAL > IdealNzero ?
                     IdealNzero : MINIDEAL ;

 //error check make sure that read in the lower/upper # thermal image that is needed
 if( MRTI_nzero > MINIDEAL ) MRTI_nzero = MINIDEAL;
 if( MRTI_ntime < MAXIDEAL ) MRTI_ntime = MAXIDEAL;

 // allocate memory all the way out to the final time instance of the 
 // optimization time interval ONLY
 MAXSTEPS= MaxRefine * m_ISTEPS_PER_IDEAL * MAXIDEAL;
 
 // fem max time this time is used to initialize the power as a 
 // function of time
 PetscScalar MaxTime = MAXIDEAL < MRTI_ntime ? 
                       IDEAL_DT * MRTI_ntime :
                       IDEAL_DT * MAXIDEAL   ;

 PetscFunctionReturn(MaxTime);
}
//----------------------------------------------------------------------
//   routine name       - indicateError
//   purpose            - indicate error for verification problems
//----------------------------------------------------------------------
void AppSolve::indicateError(const PetscScalar solnnorm,
                             const PetscScalar tol)
{
 PetscFunctionBegin; 
  if(solnnorm <  tol) 
    {
     std::cout<< "PROBLEM"               << std::endl
              << "        V"             << std::endl
              << "         E"            << std::endl
              << "          R"           << std::endl
              << "           I"          << std::endl
              << "            F"         << std::endl
              << "             I"        << std::endl
              << "              E"       << std::endl
              << "               D"      << std::endl
              << "PROBLEM"               << std::endl
              << "        V"             << std::endl
              << "         E"            << std::endl
              << "          R"           << std::endl
              << "           I"          << std::endl
              << "            F"         << std::endl
              << "             I"        << std::endl
              << "              E"       << std::endl
              << "               D"      << std::endl
              << "solnnorm = "<< solnnorm<< std::endl;
    }
   else if(solnnorm == tol) { /*do nothing*/ }
   else
    {
     std::cout <<  "PROBLEM  "                              << std::endl
               <<  "         N"                             << std::endl
               <<  "          O"                            << std::endl
               <<  "           T"                           << std::endl
               <<  "             VERIFIED!!!!!!!!!!!!!!!!!" << std::endl
               <<  "PROBLEM  "                              << std::endl
               <<  "         N"                             << std::endl
               <<  "          O"                            << std::endl
               <<  "           T"                           << std::endl
               <<  "             VERIFIED!!!!!!!!!!!!!!!!!" << std::endl
               <<  "PROBLEM  "                              << std::endl
               <<  "         N"                             << std::endl
               <<  "          O"                            << std::endl
               <<  "           T"                           << std::endl
               <<  "             VERIFIED!!!!!!!!!!!!!!!!!" << std::endl
               <<  "PROBLEM  "                              << std::endl
               <<  "         N"                             << std::endl
               <<  "          O"                            << std::endl
               <<  "           T"                           << std::endl
               <<  "             VERIFIED!!!!!!!!!!!!!!!!!" << std::endl
               <<  "solnnorm = " << solnnorm << std::endl ; 
    }
 PetscFunctionReturnVoid(); 
}

/**
 *  Setup the QOI 
 */
#undef __FUNCT__
#define __FUNCT__ "AppSolve::SetupQOI"
void AppSolve::SetupQOI( EquationSystems &eqn_systems )
{
  PetscFunctionBegin;

  GetPot &controlfile = *eqn_systems.parameters.get<GetPot*>("controlfile") ;

  /* check that AppSolve::SetupFEMSystem has been called */
  if( ! m_single->m_pdeSolver )
    {
     std::cout << "pdeSolver pointer not set "
               << "call AppSolve::SetupFEMSystem first "
               << std::endl << std::flush; libmesh_error();
    }
  /* instantiate class to handle the # of control parameters a pointer 
     to this class is needed by FormObjective and FormGradient 
     ***NOTE*** this class is instantiated within this
      scope for memory management                                   */
  const int DONTOPTIMIZE             = 0, // no optimization
            SPACETIMETEMP            = 1, // spacetime norm
            DDDASSPACETIMETEMP       = 2, // dddas spacetime norm
            VERIFICATIONQOI          = 3, // verification problem
            ONLYMRTIQOI              = 4, // only load MRTI data onto mesh
            SPACETIMEFDVERIFICATION  = 5,// verification problem w/ mrti data
            SPACETIMEHESSIANFROMDISK = 6,// build hessian from disk
            SPACETIMEHESSIANRECOMPUTE= 7,// recomp sensitiv to build hessian 
            SPARSEKALMAN             = 8,
            DIRECTKALMAN             = 9;
  // default is not to optimize
  m_NumQOI= controlfile("method/numqoi",1); // 
  PetscInt  qoiMethod= controlfile("method/qoi",DONTOPTIMIZE); // 

  // loop over the QOI
  for (PetscInt Ii = 0 ; Ii < m_NumQOI; Ii++)
   {
    qoiBaseClass *qoiOptimizer=NULL; 
    switch(qoiMethod)
     {
      case DONTOPTIMIZE:
        std::cout << "Forward SOLVE ONLY - NOT OPTIMIZING" << std::endl;
        qoiOptimizer = new noOptQOI< TransientFEMSystem > (m_single,controlfile,Ii);
        break;
      case SPACETIMETEMP:      // minimize space time norm of temperature data
        std::cout <<"Setting up ||u-u*||_{L_2([0,T],L_2(\\OMEGA))}"<<std::endl;
        qoiOptimizer = new spaceTimeQOI < PennesInverseSystem<PennesInverseModel> >(m_single,controlfile,Ii);
        break;
      case DDDASSPACETIMETEMP: 
        std::cout <<"Setting up DDDAS"<<std::endl;
        std::cout <<"Setting up ||u-u*||_{L_2([0,T],L_2(\\OMEGA))}"<<std::endl;
        //qoiOptimizer = new dddasQOI(m_single,controlfile,Ii);
        break;
      case VERIFICATIONQOI:   
        std::cout <<"Setting up Verification"<<std::endl;
        qoiOptimizer = new VerificationQOI< PennesInverseSystem < PennesVerification > >(m_single,controlfile,Ii);
        break;
      case ONLYMRTIQOI:   
        std::cout <<"Setting up MRTI data only"<<std::endl;
        qoiOptimizer = new onlyMRTIQOI< TransientFEMSystem >(m_single,controlfile,Ii);
        break;
      case SPACETIMEFDVERIFICATION:   
        std::cout <<"Setting up space Time FD Verification"<<std::endl;
        qoiOptimizer = new spaceTimeFDVerification< PennesInverseSystem < PennesVerification > >(m_single,controlfile,Ii);
        break;
      case SPACETIMEHESSIANFROMDISK:   
        std::cout <<"Setting up space Time Hessian from Disk "<<std::endl;
        //qoiOptimizer = new qoiSpaceTimeDiskRead(m_single,controlfile,Ii);
        break;
      case SPACETIMEHESSIANRECOMPUTE:   
        std::cout <<"Setting up space Time Hessian recompute sensitivity"<<std::endl;
        // qoiOptimizer = new qoiSpaceTimeRecompute(m_single,controlfile,Ii);
        break;
      case SPARSEKALMAN:   
        {
         // default to sparse kalman solver w/ maximum sparsity = F*F
         PetscInt maxSpread = 2; 
         PetscOptionsGetInt(PETSC_NULL,"-maxspread",&maxSpread,PETSC_NULL);
         if( maxSpread > 0 )
          {
           std::cout <<"Sparse Kalman..."<<std::endl;
           qoiOptimizer = new SparseKalmanFilter< TransientFEMSystem >(
                                             m_single,controlfile,Ii,maxSpread);
          }
         else if( maxSpread == 0 )
          {
           std::cout <<"Uncorrelated Kalman..."<<std::endl;
           qoiOptimizer = new UncorrelatedKalmanFilter< TransientFEMSystem >(
                                             m_single,controlfile,Ii,maxSpread);
          }
         else 
          {
           std::cout << "maxSpread >= 0 " 
                     << std::endl << std::flush; abort();
          }
        }
        break;
      case DIRECTKALMAN:   
        std::cout <<"Direct Kalman..."<<std::endl;
        qoiOptimizer = new DenseKalmanFilter< TransientFEMSystem >(m_single,controlfile,Ii);
        break;
      default:
        std::cout << "unknown qoi to optimize " << qoiMethod 
                  << std::endl << std::flush; abort();
     }

     // initialize any imaging data structures if any
     qoiOptimizer->init_imaging(qoiMethod,controlfile);

     // setup any auxillary systems
     qoiOptimizer->SetupAuxSystems(eqn_systems);

     // store this pointer for future access
     m_optimizerQOIS.push_back(qoiOptimizer); 
   }

 PetscFunctionReturnVoid(); 
}


/* ------------Setup the initial Mesh Domains ---------- */
#undef __FUNCT__
#define __FUNCT__ "InitMeshDomains"
PetscErrorCode InitMeshDomains(MeshBase &mesh,GetPot &controlfile)
{
  PetscFunctionBegin; 

  // Domain zero will always be a field domain
  // all elements with elem->subdomain_id() - 1 < n_fieldElemDomain will
  // be counted in the field optimization
  const unsigned int n_fieldDomains = controlfile("field/ndomain",1);
  libmesh_not_implemented();
  //libmesh_assert(n_fieldDomains>0 && 
  //               n_fieldDomains<= (unsigned int) user->get_num_elem_blk() );

  /* APPLY FIELD IN THE TRANSFORMED DOMAIN 
  ** Pre processes all the elements in the Serial mesh.
  */
  MeshBase::const_element_iterator       global_el   = mesh.elements_begin();
  const MeshBase::const_element_iterator global_el_end = mesh.elements_end();

  PetscInt n_elemFieldDomain = 0; // count # of elems in field domain

  for (  ; global_el != global_el_end ; ++global_el  )
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      Elem* elem = *global_el;

      // output different subdomains
      const unsigned int subdomain_id = elem->subdomain_id() - 1;
      if( !libMesh::processor_id() && subdomain_id) 
        std::cout <<"\nInitMeshDomains: element "<< elem->id() 
                  <<" has subdomain id "  << subdomain_id; 

      // count the elements in domain n_fieldDomains  
      if( elem->subdomain_id() - 1 < n_fieldDomains ) n_elemFieldDomain++;
    }

  // get data to setup up the possibly spatially varying fields
  PetscReal field_rad  = controlfile("field/field_rad",0.0);

  // error check should have no field elements if there are no field
  // optimization parameters 
  // FIXME maybe the error check should be pde dependent
  if( !controlfile("perfusion/w_0_vary",false)
                       &&
      !controlfile("thermal_conductivity/k_0_vary",false)
                       &&
      !controlfile("optical/mu_a_vary",false)
                       &&
      !controlfile("optical/mu_s_vary",false)
                       &&
      !controlfile("electric_conductivity/s_0_vary",false)
    ) field_rad = 0.0 ;

  // use laser position as default
  PetscReal laser_x_0  = controlfile("probe/x_0",0.0);
  PetscReal laser_y_0  = controlfile("probe/y_0",0.0);
  PetscReal laser_z_0  = controlfile("probe/z_0",0.0);
  Point field_center(controlfile("field/field_x_0",laser_x_0),
                     controlfile("field/field_y_0",laser_y_0),
                     controlfile("field/field_z_0",laser_z_0) );

  //PetscLogEventBegin(FiniteElementInterface::logevents[19],0,0,0,0); // field setup
    
  // count the local field elements on this proc
  PetscInt  localElemField = 0 ;
  const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
  for ( MeshBase::const_element_iterator  
         el = mesh.active_local_elements_begin(); el != el_end ; ++el  )
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      Elem* elem = *el;

      // counting assumes that a constant portion
      // of the field exists in entries [0,n_block]. if this assumption is 
      // incorrect it is caught following the loop over elements 
      if( (elem->centroid()-field_center).size() < field_rad 
                              && 
           elem->subdomain_id() - 1 < n_fieldDomains             )
        {
          localElemField++; 
        }
    }

  // determine how many field parameters are stored on all previous ranks
  std::vector<PetscInt> procElemIdStart( libMesh::n_processors() + 1, 0 );
  for(unsigned int  iii = libMesh::processor_id() + 1 ; 
                    iii < procElemIdStart.size() ; iii++ )
   {
      procElemIdStart[iii] = localElemField;
   }
  Parallel::sum(procElemIdStart);

  /* initialize counter such that first element in field optimization 
     begins with user->get_num_elem_blk()  */
  PetscInt icountfield = 0 ; 
  for ( MeshBase::const_element_iterator  
         el = mesh.active_local_elements_begin(); el != el_end ; ++el  )
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      Elem* elem = *el;

      // counting assumes that a constant portion
      // of the field exists in entries [0,n_block]. if this assumption is 
      // incorrect it is caught following the loop over elements 
      if( (elem->centroid()-field_center).size() < field_rad 
                              && 
           elem->subdomain_id() - 1 < n_fieldDomains         )
        {
          // TODO update to work w/ ExplicitSystem dof numbering
          //elem->_mat_data_field = user->get_num_elem_blk()   // begin at # of domains
          //                     // offset by other processors count of variables
          //                      +  procElemIdStart[ libMesh::processor_id() ] 
          //                     // offset by local counter
          //                      + icountfield;
          icountfield++; 
        }
    }

  // store the # of field elements for a single parameter
  //user->set_num_field_elem( procElemIdStart.back() );

  //// error check that that there is at least
  //// one element that belongs to the constant field section 
  //// TODO: seems that can relax this condition and if all elems in a domain 
  //// are used then the constant field section is allocated but unused
  //if( user->get_num_field_elem() >= n_elemFieldDomain )  
  //  {
  //   std::cout 
  //        << "NO ELEMENTS BELONG TO THE CONSTANT FIELD SECTION \n"
  //        << "in DOMAINS 0-"<< n_fieldDomains << std::endl
  //        << "This code assumes that at least two elements share a constant \n"
  //        << "part of the field in each domain, reduce field radius \n" 
  //        << "LEAVE THIS ALONE!!!!! \n" 
  //        << "the logic to solve the inverse problem with and without\n" 
  //        << "a constant field is unneccessarily complicated \n" 
  //        << "worst case scenario just have 2 elements share the constant field\n" ;
  //   std::cout  << std::flush; 
  //   abort();
  //  }

  if( !libMesh::processor_id() )
   { // echo parameters
     std::cout << "\nSetupModel: field_rad    ="<< field_rad 
          << "\nSetupModel: field_center ="<< field_center
          << std::endl;
   }
  CHKMEMQ; // check for memory corruption use -malloc_debug to enable
  //PetscLogEventEnd(FiniteElementInterface::logevents[19],0,0,0,0); // field setup

  PetscFunctionReturn(0);
}
