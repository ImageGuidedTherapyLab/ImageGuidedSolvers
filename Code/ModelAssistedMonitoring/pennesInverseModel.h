#ifndef __pennes_inverse_model_h__
#define __pennes_inverse_model_h__

#include "pennesModel.h"  // interface to pennes solver 

/* 
 * steady state coupled fluence solver
 */ 
class PennesInverseModel: public PennesStandardDiffusionApproximation 
{
public:
 // constructor
 PennesInverseModel(GetPot &, libMesh::MeshBase& );

 /** backward compatibility... no longer supported */
 virtual void accumulateResidual(unsigned int , 
                                 DiffContext &,TransientFEMSystem &) 
   { libmesh_error(); return;}

 /** backward compatibility... no longer supported */
 virtual void accumulateJacobian(unsigned int , 
                                 DiffContext &,TransientFEMSystem &) 
   { libmesh_error(); return;}

 //accumulate the adjoint matrix and rhs inherent to the pde
 virtual void accumulateAdjointPDE(const QGauss &, 
                                const unsigned int &, 
                                const std::vector<Real>&,
                                std::vector< SubMatrixVector > &,
                                std::vector< DenseSubVector<Number> > &,
                                TransientFEMSystem &,
                                TransientFEMSystem &) ;

 //accumulate the load for the sensitivity solve  useing previous sensitivity
 virtual void accumulatePreviousSensitivity(const QGauss &, 
                                 const unsigned int &, 
                                 const std::vector<Real>&,
                                 const std::vector<Point>&,
                                 std::vector< DenseSubVector<Number> > &,
                                 TransientFEMSystem &,
                                 TransientFEMSystem  &);

 // mass
 virtual PetscScalar getMass(){return rho*specific_heat;}

 virtual void SetupOptimizationVariables( libMesh::MeshBase &,
                 std::vector<optimizationParameter*> &);

 //mainly for fortran compatibility with old verification suite
 // to pass data buffer back and forth
 virtual void fillSolutionBuffer( std::vector<PetscScalar> &,const unsigned int);

 //default is to recompute jacobian matrices
 //virtual PetscTruth linearStatePDE(){  return PETSC_FALSE;} 
 virtual PetscTruth linearAdjointPDE(){return PETSC_FALSE;} 
   
 // sensitivity gradient for verification 
 virtual PetscScalar dqoidu_dudm( const unsigned int &);

protected:

  /* deriv WRT mu_a */
  inline Real dqlaserdmu_a(int domainId,const Point& qpoint, const int timeId)
  {
    Real mu_tr=mu_a[domainId]+mu_s[domainId]*(1.0e0-anfact);
    Real mu_eff=std::sqrt(3.0e0*mu_a[domainId]*mu_tr);
    Real source=0.0;
    for (unsigned int Ii = 0; Ii < m_volumeFraction.size(); Ii++ )
     {
      Point tmp(x_0[Ii],y_0[Ii],z_0[Ii]);
      Point dist = qpoint - tmp;
      source += 0.75*  Power[timeId] * m_volumeFraction[Ii] *
                 std::exp(-mu_eff*dist.size())/libMesh::pi/dist.size() * 
                 (mu_tr+mu_a[domainId]-mu_tr*mu_a[domainId]*
                  dist.size()*1.5*(mu_a[domainId]+mu_tr)/mu_eff);
     }
    return source;
  }

  /* deriv WRT mu_s */
  inline Real dqlaserdmu_s(int domainId,const Point& qpoint, const int timeId)
  {
    Real mu_tr=mu_a[domainId]+mu_s[domainId]*(1.0e0-anfact);
    Real mu_eff=std::sqrt(3.0e0*mu_a[domainId]*mu_tr);
    Real source=0.0;
    for (unsigned int Ii = 0; Ii < m_volumeFraction.size(); Ii++ )
     {
      Point tmp(x_0[Ii],y_0[Ii],z_0[Ii]);
      Point dist = qpoint - tmp;
      source += 0.75 * Power[timeId] * m_volumeFraction[Ii] *
                  std::exp(-mu_eff*dist.size())/libMesh::pi/dist.size()* 
                  mu_a[domainId]*(1.0-anfact)*
                  (1.0-1.5*mu_tr*dist.size()*mu_a[domainId]/mu_eff);
     }
    return source;
  }

  /* 2nd deriv WRT mu_a */
  inline Real d2qlaserd2mu_a(int domainId,const Point& qpoint, const int timeId)
  {
    /* brute force matlab symbolic toolbox
       d2qd2mu_a=simplify(diff(diff(qlaser,mu_a),mu_a));
       ccode(d2qd2mu_a) 
       TODO - if change, be sure to replace
               *  /mu_\([as]\)/mu_\1[domainId]/gc
               *  /Power/Power[timeId]/gc
               *  /exp/std::exp/
               *  /sqrt/std::sqrt/
               *  /pow/std::pow/
               *  /1\//1.0\//
               *  /dist/dist.size()/
     */

    Real t0,source=0.0;
    for (unsigned int Ii = 0; Ii < m_volumeFraction.size(); Ii++ )
     {
      Point tmp(x_0[Ii],y_0[Ii],z_0[Ii]);
      Point dist = qpoint - tmp;
      t0 = (Power[timeId]*std::exp(-std::sqrt(3.0)*dist.size()*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))))*(3.0/2.0))/(libMesh::pi*dist.size())+(Power[timeId]*dist.size()*std::exp(-std::sqrt(3.0)*dist.size()*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))))*std::pow(mu_a[domainId]*2.0-mu_s[domainId]*(anfact-1.0),2.0)*(9.0/1.6E1))/libMesh::pi-(std::sqrt(3.0)*Power[timeId]*mu_a[domainId]*std::exp(-std::sqrt(3.0)*dist.size()*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))))*1.0/std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)))*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*(3.0/4.0))/libMesh::pi-(std::sqrt(3.0)*Power[timeId]*std::exp(-std::sqrt(3.0)*dist.size()*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))))*1.0/std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)))*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*(mu_a[domainId]*2.0-mu_s[domainId]*(anfact-1.0))*(3.0/4.0))/libMesh::pi-(std::sqrt(3.0)*Power[timeId]*mu_a[domainId]*std::exp(-std::sqrt(3.0)*dist.size()*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))))*1.0/std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)))*(mu_a[domainId]*2.0-mu_s[domainId]*(anfact-1.0))*(3.0/4.0))/libMesh::pi+(std::sqrt(3.0)*Power[timeId]*mu_a[domainId]*std::exp(-std::sqrt(3.0)*dist.size()*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))))*1.0/std::pow(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)),3.0/2.0)*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*std::pow(mu_a[domainId]*2.0-mu_s[domainId]*(anfact-1.0),2.0)*(3.0/1.6E1))/libMesh::pi;
      source += t0 * m_volumeFraction[Ii] ;
     }
    return source;
  }

  /* 2nd deriv WRT mu_s */
  inline Real d2qlaserd2mu_s(int domainId,const Point& qpoint, const int timeId)
  {
    /* brute force matlab symbolic toolbox
       d2qd2mu_s=simplify(diff(diff(qlaser,mu_s),mu_s));
       ccode(d2qd2mu_s) 
       TODO - if change, be sure to replace
               *  /mu_\([as]\)/mu_\1[domainId]/gc
               *  /Power/Power[timeId]/gc
               *  /exp/std::exp/
               *  /sqrt/std::sqrt/
               *  /pow/std::pow/
               *  /1\//1.0\//
               *  /dist/dist.size()/
     */

    Real t0,source=0.0;
    for (unsigned int Ii = 0; Ii < m_volumeFraction.size(); Ii++ )
     {
      Point tmp(x_0[Ii],y_0[Ii],z_0[Ii]);
      Point dist = qpoint - tmp;
      t0 = (Power[timeId]*std::exp(-dist.size()*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*3.0))*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)))*1.0/std::pow(mu_a[domainId]-mu_s[domainId]*(anfact-1.0),2.0)*std::pow(anfact-1.0,2.0)*(std::sqrt(3.0)*(mu_a[domainId]*mu_a[domainId])-std::sqrt(3.0)*mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*4.0+std::sqrt(3.0)*mu_a[domainId]*mu_s[domainId]+dist.size()*mu_a[domainId]*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)))*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*3.0-std::sqrt(3.0)*anfact*mu_a[domainId]*mu_s[domainId])*(3.0/1.6E1))/libMesh::pi;
      source += t0 * m_volumeFraction[Ii] ;
     }
    return source;
  }

  /* mixed deriv WRT mu_a and mu_s */
  inline Real d2qlaserdmu_a_dmu_s(int domainId,const Point& qpoint, const int timeId)
  {
       
    /* brute force matlab symbolic toolbox
       d2qd2mu_a_mu_s=simplify(diff(diff(qlaser,mu_a),mu_s));
       ccode(d2qd2mu_a_mu_s) 
       TODO - if change, be sure to replace
               *  /mu_\([as]\)/mu_\1[domainId]/gc
               *  /Power/Power[timeId]/gc
               *  /exp/std::exp/
               *  /sqrt/std::sqrt/
               *  /pow/std::pow/
               *  /1\//1.0\//
               *  /dist/dist.size()/
     */

    Real t0,source=0.0;
    for (unsigned int Ii = 0; Ii < m_volumeFraction.size(); Ii++ )
     {
      Point tmp(x_0[Ii],y_0[Ii],z_0[Ii]);
      Point dist = qpoint - tmp;
      t0 = (Power[timeId]*1.0/(mu_a[domainId]*mu_a[domainId])*std::exp(-dist.size()*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*3.0))*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)))*1.0/std::pow(mu_a[domainId]-mu_s[domainId]*(anfact-1.0),2.0)*(anfact-1.0)*(mu_a[domainId]*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)))*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*4.0+std::sqrt(3.0)*dist.size()*(mu_a[domainId]*mu_a[domainId]*mu_a[domainId]*mu_a[domainId])*2.0+(dist.size()*dist.size())*(mu_a[domainId]*mu_a[domainId]*mu_a[domainId])*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)))*(mu_a[domainId]*6.0-mu_s[domainId]*(anfact-1.0)*6.0)-std::sqrt(3.0)*dist.size()*(mu_a[domainId]*mu_a[domainId]*mu_a[domainId])*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*1.0E1+std::sqrt(3.0)*dist.size()*(mu_a[domainId]*mu_a[domainId]*mu_a[domainId])*mu_s[domainId]*3.0+std::sqrt(3.0)*dist.size()*(mu_a[domainId]*mu_a[domainId])*(mu_s[domainId]*mu_s[domainId])-std::sqrt(3.0)*anfact*dist.size()*(mu_a[domainId]*mu_a[domainId])*(mu_s[domainId]*mu_s[domainId])*2.0+(dist.size()*dist.size())*(mu_a[domainId]*mu_a[domainId])*mu_s[domainId]*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)))*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*3.0+std::sqrt(3.0)*(anfact*anfact)*dist.size()*(mu_a[domainId]*mu_a[domainId])*(mu_s[domainId]*mu_s[domainId])-std::sqrt(3.0)*dist.size()*(mu_a[domainId]*mu_a[domainId])*mu_s[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*6.0-std::sqrt(3.0)*anfact*dist.size()*(mu_a[domainId]*mu_a[domainId]*mu_a[domainId])*mu_s[domainId]*3.0+std::sqrt(3.0)*anfact*dist.size()*(mu_a[domainId]*mu_a[domainId])*mu_s[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*6.0-anfact*(dist.size()*dist.size())*(mu_a[domainId]*mu_a[domainId])*mu_s[domainId]*std::sqrt(mu_a[domainId]*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0)))*(mu_a[domainId]-mu_s[domainId]*(anfact-1.0))*3.0)*(-3.0/1.6E1))/(libMesh::pi*dist.size());
      source += t0 * m_volumeFraction[Ii] ;
     }
    return source;
  }

  // optimization functions
  virtual PetscScalar dpde_derr(         OPTGAUSSARG );
  virtual PetscScalar dpde_dk_0(         OPTGAUSSARG );
  virtual PetscScalar dpde_dw_0(         OPTGAUSSARG );
  virtual PetscScalar dpde_dmu_s(        OPTGAUSSARG );
  virtual PetscScalar dpde_dmu_a(        OPTGAUSSARG );
  virtual PetscScalar dpde_dprobeTemp(   OPTGAUSSARG );
  virtual PetscScalar d2pde_du_dk_0(     OPTGAUSSARG );
  virtual PetscScalar d2pde_du_dw_0(     OPTGAUSSARG );
  virtual PetscScalar d2pde_d2mu_a(      OPTGAUSSARG );
  virtual PetscScalar d2pde_d2mu_s(      OPTGAUSSARG );
  virtual PetscScalar d2pde_dmu_a_dmu_s( OPTGAUSSARG );

};
/* 
 * Direct coupled fluence solver
 */ 
class PennesFluencePDE : public PennesInverseModel
{
public:

 //constructor
 PennesFluencePDE(GetPot &,libMesh::MeshBase&);

 //virtual void accumulateResidual(const QGauss &, 
 //                                const unsigned int &, 
 //                                const std::vector<Real>&,
 //                                const std::vector<Point>&,
 //                                std::vector< DenseSubVector<Number> > &,
 //                                TransientFEMSystem &) ;
 //virtual void accumulateJacobian(const QGauss &,
 //                                const unsigned int &, 
 //                                const std::vector<Real>& ,
 //                                std::vector< SubMatrixVector > &,
 //                               TransientFEMSystem &) ;

 //print data 
 virtual void printSelf(std::ostream& os=std::cout) ; 

 // overwrite exactsolution
 virtual Number exactSolution(const Point& p,
                              const Parameters& parameters,
                              const std::string& var);   	 	 

 //get initial condition data
 Real getInitFluence(const unsigned int domainID){return initialFluence[domainID];}

 // verify the solution
 virtual void verifySolution( EquationSystems &eqnSystem );

protected:

  // initial fluence
  std::vector< Real > initialFluence;
 
  // laser parameters
  std::vector< Real > speedLight,alpha;
};
/* 
 * steady state coupled fluence solver
 */ 
class PennesFluenceSS : public PennesFluencePDE 
{
public:

 //constructor
 PennesFluenceSS(GetPot &,libMesh::MeshBase&);

 //virtual void accumulateResidual(const QGauss &, 
 //                                const unsigned int &, 
 //                                const std::vector<Real>&,
 //                                const std::vector<Point>&,
 //                                std::vector< DenseSubVector<Number> > &,
 //                                TransientFEMSystem &) ;
 //virtual void accumulateJacobian(const QGauss &,
 //                                const unsigned int &, 
 //                                const std::vector<Real>& ,
 //                                std::vector< SubMatrixVector > &,
 //                               TransientFEMSystem &) ;

 // overwrite exactsolution
 virtual Number exactSolution(const Point& p,
                              const Parameters& parameters,
                              const std::string& var);   	 	 
protected:
 
private:
};

#endif
