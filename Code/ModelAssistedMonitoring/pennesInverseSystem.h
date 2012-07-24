#ifndef __pennes_inverse_system_h__
#define __pennes_inverse_system_h__

#include "thermal_therapy_system.h"  // interface to pennes solver 
#include "pennesSystem.h"  // interface to pennes solver 
/**@ingroup InverseSolver 
 * base clase for pde constitutitve data for an arbitrary number of variables.
 * contains basic methods for solution evaluation and boundary condition
 * evaluations assuming that systems with multiple variables have uncoupled
 * boundary conditions
 *
 *   bridge pattern implemented for qoi and pde heirarchies
 *
 *   http://www.parashift.com/c++-faq-lite/multiple-inheritance.html
 *
 *   for example see:
 *
 *             spaceTimeQOI::accumulateAdjointQOI
 *
 *   the pdeBaseClass pointer is used to call any simple operations on the
 *   solutions at the gauss points, ie
 *
 *   purpose            - this is the main data structure that allows 
 *                        communication between libMesh/PETSC and TAO.
 *                        this is used for the optimizations
 * 
 * Pennes Bioheat Equation with SDA laser is treated as the base class
 * for the variety of laser source models 
 * the default QOI is assumed to be the space time norm of the difference 
 * between an ideal solution. this will need to be overridden
 * in a derived class for a new qoi; FIXME: may or may not work out well.
 */ 
template< typename MathematicalModel  >
class PennesInverseSystem: public LITTSystem < MathematicalModel > 
{
public:

 //constructor
 PennesInverseSystem(EquationSystems& ,const std::string& , const unsigned int );

 // transfer C++ to F90 data structures 
 virtual void XferData(){};

 /** wrapper to old assembly */
 virtual bool element_time_derivative (bool , DiffContext &);

 /** pre FEMSystem adjoint assembly */
 virtual void assemble_adjoint(EquationSystems&);

 /** pre FEMSystem sensitivity load */
 virtual void sensitivityLoad( NumericVector<Number>& );

 /** regression */
 virtual void verifySolution( EquationSystems &eqnSystem )
    { this->m_MathModel.verifySolution(eqnSystem); return; }

 /** adjoint gradient */
 virtual PetscErrorCode FormGradient(TAO_APPLICATION ,Vec , Vec );

 /** hessian vector product */
 virtual PetscErrorCode hessianVectorProduct(Mat , Vec , Vec );

private:
  /** The type of the parent  */
  typedef LITTSystem< MathematicalModel > Parent;
  
  // boundary condition flag
  PetscInt dirichletID; 

};

void initial_sensitivity_matrix (   EquationSystems&, const std::string& );
void assemble_adjoint(  EquationSystems&, const std::string& );

#endif
