//prevent multiple inclusions of header file
#ifndef optimizationParameter_H
#define optimizationParameter_H

// petsc includes
#include "petsc_macro.h" 

//use Explicit system as the data storage...
#include "explicit_system.h"
#include "numeric_vector.h"

// arguement list of functers for Gradient and Hessian evaluation
#define OPTGAUSSARG const unsigned int &, const Point &
// typedefs to make function pointer declarations more readable
typedef PetscScalar (*OptGaussType)(OPTGAUSSARG) ;

class MeshBase; // forward declaration
class qoiBaseClass; // forward declaration
class PDEModelBaseClass;
typedef PetscScalar (qoiBaseClass::*dqoi_dmMemFn)(OPTGAUSSARG);
typedef PetscScalar (PDEModelBaseClass::*dpde_dmMemFn)(OPTGAUSSARG);
typedef PetscScalar (PDEModelBaseClass::*d2pde_dudmMemFn)(OPTGAUSSARG);

class optimizationParameter {

   friend class qoiBaseClass;

public:
 optimizationParameter( const char *, bool , dpde_dmMemFn, d2pde_dudmMemFn,
                        PetscScalar, PetscScalar, PetscScalar, PetscScalar);
 std::vector<PetscInt> dofs;
 bool spatial_field;   // if true this parameter may vary spatially
 bool time_vary;       // if true this parameter varies in time
 void getParam_lb(   PetscScalar &data,PetscInt &id){data=_value_lb[id];}
 void getParam_ub(   PetscScalar &data,PetscInt &id){data=_value_ub[id];}
 /* all function pointers for
    optimization should have the SAME argument lists given by OPTGAUSSARG */ 
 /* QOI, gradient, and hessian contributions of this parameter  
    these terms appear from the regularization terms 
 */ 
 PetscScalar qoi(PetscInt &,PetscScalar &)
  { libmesh_not_implemented(); return -1;
 // return  
 // std::exp( Scalefact/(_value_ub[id]-_value_lb[id])*(_value[id]-_value_ub[id]))
 //+std::exp(-Scalefact/(_value_ub[id]-_value_lb[id])*(_value[id]-_value_lb[id]));
  }
 /* QOI first derivates, NOTE that the integration rule DRAMATICALLY effects
    the time point used i.e.  k-1/2 vs k...  */ 
 dqoi_dmMemFn dqoi_dm;
 

 /* gradient, and hessian contributions of this parameter  
    these terms are inherent to the PDE 

    list of function pointers for dpde_dm
       error_estimate 
       dpde_dmu_a dpde_dmu_s dpde_dk_0 dpde_dw_0
       dwdw0 dwdwn dwdwi dwdwd dwdw2 dwdwni dwdwid
       dkdk0 dkdk1 dkdk2 dkdk3
       dqlaserdmu_a dqlaserdmu_s dqlaserdx dqlaserdy dqlaserdz dqlaserdpow
       dqmontedx dqmontedy dqmontedz dqmontedpow
 */ 
 dpde_dmMemFn dpde_dm;

 d2pde_dudmMemFn d2pde_du_dm, d2qoi_du_dm; /* d2pde_du_dk_0 d2pde_du_dw_0 */

 /* function to setup all function pointers and 
    verification suite data for one control variable*/
 virtual PetscErrorCode OneVarSetup(EquationSystems &) = 0 ;

 bool Optimize;

 /**
  * Arbitrary, user-specified name of the variable.
  */
 const std::string & name() const { return _name; }

 // upper and lower bounds
 PetscScalar&  lb(const int id){    return _value_lb.at(id); }
 PetscScalar&  ub(const int id){    return _value_ub.at(id); }
 PetscScalar&  dbeta(){ return _perturbation; }
 PetscScalar&  verif(){ return _verifparam; }

 std::string _name;   ///< name for id
private:
 std::vector< PetscScalar > _value_lb, _value_ub;
 PetscScalar _verifparam,   
             _perturbation;
} ;
class discreteParameter : public optimizationParameter 
{
public:
 discreteParameter (const char *,bool , dpde_dmMemFn,  d2pde_dudmMemFn,
                    PetscScalar, PetscScalar, PetscScalar, PetscScalar );
 PetscScalar&  operator[](const int id) { return _value[id]; }
 PetscScalar&  at(const int id) { return _value.at(id); }
 PetscInt size() const { return static_cast<PetscInt>(_value.size()); }
 void resize(unsigned int id,PetscScalar Xx){_value.clear();
                                             _value.resize(id,Xx);return;}
 void clear(){_value.clear();return;}
 void push_back(const PetscScalar Value){_value.push_back(Value);return;}
 virtual PetscErrorCode OneVarSetup(EquationSystems &)  ;
 void printStdVector(std::ostream& os, const char* Name)
 {
   for(unsigned int Ii = 0 ; Ii < _value.size(); Ii ++)
       os << Name <<Ii<<"]="<< _value[Ii] << std::endl;
 }
private:
 std::vector< PetscScalar > _value; /// local data depracated and trying to remove
} ;
class spatialParameter  : public optimizationParameter 
{
public:
 //overload operators
 PetscScalar  operator[] (const int id) { return m_data->current_solution(id); }
 void putParam(const PetscScalar &data,PetscInt &id){m_data->solution->set(id,data);}
 void getParam(      PetscScalar &data,PetscInt &id){data= m_data->current_solution(id); }

 // create a local copy of the global solution on all procs
 PetscErrorCode ScatterGlobalToAll();

 // get an entry for the local global vector
 PetscScalar GetGlobalSolution(const unsigned int);

 // constructor
 spatialParameter (const char *,EquationSystems &, bool , dpde_dmMemFn, d2pde_dudmMemFn,
                   PetscScalar, PetscScalar, PetscScalar, PetscScalar, bool);
 // clean up
 ~spatialParameter(); 
 // setup the variable
 virtual PetscErrorCode OneVarSetup(EquationSystems &)  ;
 void printStdVector(std::ostream& os, const char* Name)
 {
   os << Name  << std::endl;
   if(m_data)
    {
     os << Name << "current_local_solution"<< std::endl;
     m_data->current_local_solution->print(os); 
     os << Name << "solution"<< std::endl;
     m_data->solution->print(os); 
    }
   else
    {
     os << "not setup..."<< std::endl;
    }
   return; 
 }
protected:
 ExplicitSystem *m_data;
 VecScatter  m_globalScatter;
 Vec m_LocalGlobalSolution;
} ;


#endif
