#ifndef __realtimeApplicationContext_h
#define __realtimeApplicationContext_h
// C++ include files that we need
#include <iostream>
#include <vector>

// boost includes
#include <boost/lexical_cast.hpp>

// petsc/tao  includes
#include "tao.h" 

// itk includes
#include <itksys/SystemTools.hxx>
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"

// libmesh
#include "libmesh.h"
#include "equation_systems.h"
#include "petsc_macro.h"
#include "getpot.h"

// libmesh
#include "transient_fem_system.h"

/**forward declaration */
class qoiBaseClass;
class ExactSolution;

/** default pixel type is double */
typedef PetscScalar                                              InputPixelType;
typedef float                                                   OutputPixelType;
const unsigned int          Dimension = 3;
typedef itk::Image< InputPixelType, Dimension >                  InputImageType;
typedef itk::ImageSeriesReader<    InputImageType >             InputReaderType;


/**
 *   this class structure is a singleton design pattern 
 *     http://www.codeproject.com/KB/cpp/singletonrvs.aspx?msg=242944#xx242944xx
 *   that contains general information meant to 
 *   coordinate a group of processor groups w/ different qoi.
 *   
 *   this is the main data structure that allows 
 *   communication between libMesh and TAO.
 *   This is essentially the brute force data structure
 *   written with the assumption that no matter how hard
 *   you try the Data Structures will NEVER be perfectly
 *   contained and modularize so a pointer the other data
 *   structures is stored here
 * 
 *   @todo {need to fully encapsulate this class... some members are still
 *          exposed b/c i'm lazy}
 */
class AppSolve 
{
public:
  /** calls private constructor 
    * minimal variables setup before debugger setup
    */
  static AppSolve* getInstance();
  ~AppSolve()
  {
      m_instanceFlag = false;
  }

  /** setup the majority of the variables after the debugger setup */
  PetscErrorCode Setup(GetPot &); 
  PetscScalar    setupTime(const PetscInt , const PetscInt ,
                           const PetscInt );
  void           setupProfile();
  /**
   * Read the controlfile and get the FEM System
   */
  TransientFEMSystem* SetupFEMSystem( EquationSystems &);

  /**
   * Read the controlfile and get the FEM System
   */
  void SetupQOI(EquationSystems &);

  /** echo params */
  void printSelf(std::ostream& );

  // basic functions
  PetscInt getMinFEMTimeStep(){return MRTI_nzero * m_ISTEPS_PER_IDEAL;} 
  PetscInt getMaxFEMTimeStep(){return MAXIDEAL   * m_ISTEPS_PER_IDEAL;} 
  PetscInt getFEMTimeStep(PetscInt step){return step*m_ISTEPS_PER_IDEAL;} 
  PetscScalar getTime(PetscInt istep){return istep* FEM_DT;} 
  bool modWrite(PetscInt istep){return istep%m_WriteInterval==0;} 
  /** IDEAL_NTIME at final optimization step  */
  PetscInt FinalOpt_Ntime() {return m_IDEAL_NTIME + m_NUMIDEAL*(m_NoptSteps-1);}
  /** IDEAL_NZERO at final optimization step  */
  PetscInt FinalOpt_Nzero() {return m_IDEAL_NZERO + m_NOFFSET *(m_NoptSteps-1);}
  /**
    * power time instance. floor fem time step 
    * by ISTEPS_PER_IDEAL to get power step
    */
  PetscInt get_id_power(){   return 
                               ISTEP/m_ISTEPS_PER_IDEAL;} 

  PetscTruth  OptimizeMesh(){         return m_optimizeMesh;}  
  PetscTruth  ComputeFinalObjective(){return m_COMPUTEFINALOBJECTIVE;}
  PetscTruth  PlotIteration(){        return m_PLOTITER;}
  PetscTruth  ControlTask(){          return m_Control_Task;}
  PetscTruth  Restart(){              return m_IDrestart;}
  PetscTruth  ImageInitialCondition(){return m_IC_MRTI ; }
  PetscInt    get_num_elem_blk(){     return n_block ; }
  PetscInt    get_num_field_elem(){   return NfieldElem ; }
  void        set_num_field_elem(PetscInt InputElem)
      { 
           NfieldElem = InputElem; 
           return ; 
      }
  void  UpdateIdealWindow()
   {
    m_IDEAL_NZERO += m_NOFFSET ;
    m_IDEAL_NTIME += m_NUMIDEAL;
    return ;
   }

  PetscScalar MaxTime(){          return m_MaxTime;}
  PetscScalar get_fem_dt()  {     return FEM_DT ; }
  PetscScalar get_ideal_dt(){     return IDEAL_DT; }
  PetscInt    NoptSteps(){        return m_NoptSteps;}
  PetscInt    IdealNzero(){       return m_IDEAL_NZERO;}
  PetscInt    Noffset(){          return m_NOFFSET;}
  PetscInt    IdealNtime(){       return m_IDEAL_NTIME;}
  PetscInt    NumIdeal(){         return m_NUMIDEAL;}
  PetscInt    IstepsPerIdeal(){   return m_ISTEPS_PER_IDEAL;}
  PetscInt get_max_ideal(){       return MAXIDEAL;} 
  PetscInt get_max_steps(){       return MAXSTEPS;} 
  PetscInt get_num_MRTI(){        return MRTI_ntime;} 
  PetscInt get_num_zero(){        return MRTI_nzero;} 
  PetscInt get_num_exact(){       return NEXACT;} 
  PetscInt WriteInterval(){       return m_WriteInterval;} 
  PetscInt Refine(){              return m_Refine;} 
  PetscInt NumQOI(){              return m_NumQOI;} 
  qoiBaseClass* qoiOptimizer()
    {   
     if(m_NumQOI > 1) 
       return m_optimizerQOIS[IDOPT];
     else if(m_NumQOI == 1) 
       return m_optimizerQOIS[0];
     else 
       return NULL;
    }
  TransientFEMSystem* pdeSolver(){return m_pdeSolver;    }
  EquationSystems & get_equation_systems() 
    {
     return m_pdeSolver->get_equation_systems();
    }
  MeshBase & get_mesh() 
    {
     return m_pdeSolver->get_equation_systems().get_mesh();
    }

  // checkpointing
  bool CheckPoint(PetscInt iter)
    { return iter >= m_checkpoint;}
  void UpdateCheckPoint(PetscInt iter)
    { m_checkpoint = iter + m_checkpointiter; return; }
  /**
    * pre FEMSystem current fem time step 
    *  @todo {remove this a public variable... but
	      it is embedded in the 
              interface of pdeBaseClass
             }
    */
  static PetscInt ISTEP; 

  // verification problems
  static void indicateError(const PetscScalar ,const PetscScalar );

  // For Profiling  
#if PETSC_VERSION_LESS_THAN(3,0,0)
  static std::vector<PetscInt>  logstages;
  static std::vector<PetscInt>  logevents;
#else
  static std::vector<PetscLogStage>  logstages;
  static std::vector<PetscLogEvent>  logevents;
#endif

  static PetscInt GroupID;                 ///< group ID

  /** current optimization step */
  static PetscInt IDOPT;

  // used for file/job ID and local disk location
  static std::string profileID,localDisk,restartFile;

  static PetscInt controlfileID,outputfileID ; // file ID's
private:
  /** private constructor (minimal variables setup before debugger setup) */
  AppSolve();
  /** private pointer to application context */
  static AppSolve *m_single;
  static bool m_instanceFlag;

  static TransientFEMSystem *m_pdeSolver;    ///< pointer to pde solver

  /** vector of pointers to qoi solvers */
  static std::vector<qoiBaseClass*> m_optimizerQOIS; 


  static PetscInt MRTI_nzero, MRTI_ntime, // thermal image time bounds
                  MAXSTEPS,    ///< maximum number of time steps
                  MINIDEAL,    ///< minimum number of IDEAL steps
                  MAXIDEAL,    ///< maximum number of IDEAL steps
                  m_ISTEPS_PER_IDEAL, ///< fem time steps per ideal timestep
                  m_checkpoint , m_checkpointiter, //checkpointing variables
                  m_Refine,      ///< refinement level flag
                  m_WriteInterval, ///< number of times to skip before plotting
                  NEXACT,      ///< exact solution
                  TIMELAG_IC , ///< time-lag for the initial condition
                  Total_Nopt,  ///< total number of optimization steps
                  NfieldElem,  ///< number of field elements for a single param
                  m_NumQOI,    ///< number of QOI's
                  n_block;     ///< # of mesh blocks

  /**  used in parameter studies to 
    *  compute the value of the objective 
    *  function over the entire time window 
    *  when the optimization is done only on a
    *  subset of the window
    */
  static PetscTruth  m_COMPUTEFINALOBJECTIVE; 
  static PetscTruth  m_optimizeMesh, 
                     m_PLOTITER, // plot file control
                     m_Control_Task, // PETSC_TRUE if control task
                     m_IDrestart,    // PETSC_TRUE if this is a restart 
                     m_IC_MRTI ; // IC_MRTI = PETSC_TRUE ==> get IC thermal image
  static PetscScalar IDEAL_DT, ///< # ideal deltat
                     FEM_DT;      ///< \f$ \Delta t \f$

  static PetscScalar m_MaxTime;  ///< max time
  static PetscInt m_NoptSteps,   ///< # of optimization steps
         m_IDEAL_NZERO, ///< lower bound of optimization window
         m_NOFFSET,     ///< increment of IDEAL_NZERO at each optimization solve
         m_IDEAL_NTIME, ///< upper bound of optimization window
         m_NUMIDEAL;    ///< increment of IDEAL_NTIME at each optimization solve

};
#endif // #define __realtimeApplicationContext_h
