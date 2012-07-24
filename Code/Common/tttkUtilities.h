#ifndef __realtimetttkUtilities_h
#define __realtimetttkUtilities_h

#include "numeric_vector.h"
#include "petsc_vector.h"

// utility for printing
template<class T>
void printStdVector(std::ostream& os, const char* Name, std::vector<T> &x)
{
  for(unsigned int Ii = 0 ; Ii < x.size(); Ii ++)
      os << Name <<Ii<<"]="<< x[Ii] << std::endl;
}

PetscErrorCode MatDataInfo(Mat , const char [] );
PetscErrorCode MatDataSparseLossInfo(Mat , Mat ,const char []);
/**
 * print matrix norm info to determine the relative effect
 *   of matrix projection into the constrained sparsity pattern
 */
PetscErrorCode VecDataInfo(Vec , const char [] );
PetscErrorCode SetupUnStructuredGrid(Mesh *, char *, int ,
                                     double , double , double ,
                                     double , double , double ,
                                     double , double , double ,
                                     double , double , double  );
PetscErrorCode GenerateStructuredGrid( Mesh &,
                                    int , int , int ,
				    double , double , 
				    double , double , 
				    double , double , 
                                    std::vector<int> &); 
PetscErrorCode SetParameterIni( Parameters *, GetPot *);
PetscErrorCode PetscPopErrorHandlerAndDebug();
PetscErrorCode FlushStdCoutCerr();
PetscErrorCode StoreSystemTimeStep(System* , int );
PetscErrorCode AddStorageVectors(  System* , char * , int );
PetscErrorCode  GetSolutionVector(System& , Vec *   );
PetscErrorCode  SetSolutionVector(System& , Vec     );
PetscErrorCode CopySolutionVector(System& , System& );
System* GetSystem             ( char *, EquationSystems *);
System* AddExplicitSystem     ( char *, EquationSystems *);
System* AddBackgroundSystem   ( char *, EquationSystems *);
System* AddPennesDeltaPSystem ( char *, EquationSystems *,double);
System* AddPennesSDASystem    ( char *, EquationSystems *,double);
System* AddPennesRFSystem     ( char *, EquationSystems *,double);
unsigned int AddConstantMonomialVariable(System *, char *);
unsigned int AddFirstLagrangeVariable(   System *, char *);
PetscScalar WEIGHTEDL2Norm(   System*,char *, System*,char *, System*,char * );
Mat GETFactorMat(  Mat ,char *,char *);
PetscErrorCode PetscFEMSystemUpdateTimeStep(        System *, int);
PetscErrorCode PetscFEMSystemSetupInitialConditions(System *);
PetscErrorCode PetscFEMSystemGetSolnSubVector(      System *, int,Vec *);
PetscInt       PetscFEMSystemCreateNodeSetFromMask( System* ,double,int);
PetscErrorCode PennesSDASystemUpdateLaserPosition(System* ,
                      double,double,double,double,double,double );
PetscErrorCode PennesSDASystemUpdateLaserPower(System* , double,int);
#endif
