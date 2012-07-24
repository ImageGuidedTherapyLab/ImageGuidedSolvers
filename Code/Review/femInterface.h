#ifndef __FiniteElementInterface_h
#define __FiniteElementInterface_h
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
 * c++ class meant to interface FEM routines with python
 * 
 *   @todo {need to fully encapsulate this class... some members are still
 *          exposed b/c i'm lazy}
 */
// libmesh includes
#include "mesh.h"
#include "equation_systems.h"
#include "getpot.h"

// forward declarations
class Imaging; 
namespace libMesh
{
class ExodusII_IO;
}


class FiniteElementInterface
{
public:

 /** calls private constructor 
   * minimal variables setup before debugger setup
   */
 static FiniteElementInterface* getInstance();
 ~FiniteElementInterface(); ///< destructor

 /** print local information */ 
 void printSelf (std::ostream&);
 void printSelf (){this->printSelf(std::cout);}

 /** setup libMesh */ 
 PetscErrorCode SetuplibMesh(MPI_Comm);

 /** clean up data for parameter Study */ 
 PetscErrorCode CleanUplibMesh();

 /** setup FEM mesh for libMesh */ 
 PetscErrorCode SetupStructuredGrid(int   , int   , int   ,
                                    double, double, 
                                    double, double, 
                                    double, double, 
                                    std::vector<int> * );
 
 /** setup FEM mesh for libMesh */ 
 PetscErrorCode SetupUnStructuredGrid(char *, int ,
                                    double, double, double, 
                                    double, double, double, 
                                    double, double, double, 
                                    double, double, double );

 /** add an explict system */ 
 PetscErrorCode AddExplicitSystem (char *,int,int);

 /** add a background system */ 
 PetscErrorCode AddBackgroundSystem (char *,int);

 /** add a Pennes system */ 
 PetscErrorCode AddPennesSystem (int,double); 

 /** update time step info */ 
 PetscErrorCode UpdateTransientSystemTimeStep(char *, int );
 PetscErrorCode StoreTransientSystemTimeStep( char *, int );
 PetscErrorCode StoreSystemTimeStep(          char *, int );

 /** setup ini file */ 
 PetscErrorCode SetupIniFile (char *);
 PetscErrorCode SetupCmdLineIni ();
 PetscErrorCode SetIniValue (char *,char *);

 /** parse control file */ 
 PetscErrorCode Setup( GetPot& ); 

 /** setup Imaging data structures */ 
 PetscErrorCode SetupImaging( int ); 
 PetscErrorCode SetupImaging( int, int, int,
	                      double , double , double , 
                              double , double , double  ) ; 

 /** setup Imaging to FEM data structures */ 
 PetscErrorCode ProjectImagingToFEMMesh (char *, Vec );

 /** copy PetscVec to python */ 
 PetscErrorCode GetSolutionVector(char *, Vec* );

 /** copy python to PetscVec */ 
 PetscErrorCode SetSolutionVector(char *, Vec );

 /** setup profiling data */ 
 void  SetupProfile();

 /** initialize Equation Systems */ 
 PetscErrorCode InitializeEquationSystems();

 /** Solve the System */ 
 PetscErrorCode SystemSolve(char *);

 /** write libMesh data structures to view in paraview */ 
 PetscErrorCode WriteEquationSystems( char *);
 PetscErrorCode WriteParameterSystems(char *);
 PetscErrorCode WriteTimeStep( char *, int , double );

 /** get reference to FEM Mesh data structure */ 
 MeshBase & get_mesh() { return *this->m_mesh; }

 /** get reference to equation systems data structure */ 
 EquationSystems & get_equation_systems() { return *this->m_eqn_systems; }

 /** get imaging pointer */
 Imaging *ImagingPointer(){return m_images;}

 // For Profiling  
 #if PETSC_VERSION_LESS_THAN(3,0,0)
 static std::vector<PetscInt>  logstages;
 static std::vector<PetscInt>  logevents;
 #else
 static std::vector<PetscLogStage>  logstages;
 static std::vector<PetscLogEvent>  logevents;
 #endif

 // used for job ID
 static std::string profileID;

private:

 // TODO: create separate python classes to interface each of these pointers
 static Mesh *m_mesh; ///< pointer for mesh object
 static ExodusII_IO *m_exodus_transient_file; ///< pointer for exodus output 
 static EquationSystems *m_eqn_systems; ///< pointer for equation systems object for state
 static LibMeshInit *m_init;///< pointer for equation systems object
 static GetPot *m_controlfile;///< pointer to parameter file
 static Imaging *m_images;///< pointer to imaging data structures 

 /** 
  * private constructor (minimal variables setup before debugger setup) 
  */
 FiniteElementInterface(){};
 /** private pointer to application context */
 static FiniteElementInterface *m_single;
 static bool m_instanceFlag;
};
// access from python
FiniteElementInterface* FiniteElementInterfaceGetInstance();
// Wrapper to evaluate exact solution
template< typename SystemModelType >
Number pdeExactSolution(const Point& p,
                        const Parameters& parameters,
                        const std::string& ,
                        const std::string& unknown_name)
{ // return the exact solution
  FiniteElementInterface *user = FiniteElementInterface::getInstance();
  SystemModelType &tmpSystem =  
                 user->get_equation_systems().get_system<SystemModelType>(unknown_name) ;
  return tmpSystem.exactSolution(p,parameters,unknown_name);
}  	 	 

#endif
