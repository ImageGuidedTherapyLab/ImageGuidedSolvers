// Copyright 2008 Google Inc.
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Author: wan@google.com (Zhanyong Wan)
// Author: vladl@google.com (Vlad Losev)

// This provides interface PrimeTable that determines whether a number is a
// prime and determines a next prime number. This interface is used
// in Google Test samples demonstrating use of parameterized tests.

#ifndef GTEST_PENNES_REGRESSION_H_
#define GTEST_PENNES_REGRESSION_H_

/**@ingroup TreatmentPlanning 
 * Developing a manufactured solution is crucial to debugging the output
 * of the Kalman Filter. A simple 1D problem would be sufficient.
 * \f{figure}[h]
 * \centering
 * \[
 *   u(x,t,\omega) = q (x-L) ^2 t + w(x,\omega)
 * \qquad
 *   w(x_i,\omega) \sim N(0,Q) \qquad \forall x_i
 * \]
 * \begin{picture}(140,32)
 *   \multiput(0,10)(30,0){3}{\line(0,1){10}}
 *   \put(0,15){\line(1,0){120}}
 *   \put(2,5){$x_0$}
 *   \put(32,5){$x_1$}
 *   \put(62,5){$x_2$}
 *   \put(82,6){$\dots$}
 *   \put(120,10){\line(0,1){10}}
 *   \put(120,5){$x_n$}
 * \end{picture}
 * \caption{1D manufactured solution}
 * \f}
 * 
 * where
 *   \f[
 *     \vec{u}(t) = \left[\begin{split} 
 *                                       u(x_0,t,\omega) \\
 *                                       u(x_1,t,\omega) \\
 *                                       u(x_2,t,\omega) \\
 *                                               .       \\
 *                                               .       \\
 *                                               .       \\
 *                                       u(x_n,t,\omega) \\
 *                                   \end{split} \right]
 * \qquad
 *     f_{w(x_0),w(x_1),...,w(x_n)}(\xi_0,\xi_1,...,\xi_n)  = 
 *       \left[ (2\pi)^n  |P|^{1/2}\right]^{-1}
 *       \exp\left\{ -1/2 (\xi - 0)^T P^{-1}(\xi - 0)\right\}
 *    \f]
 * the statistics for the manufactueed solution are straight forward
 * \f[
 * \hat{u} (t) \equiv E[\vec{u}(t)] = 
 *                  \left[\begin{split} 
 *                                       q(x_0-L)^2 t \\
 *                                       q(x_1-L)^2 t \\
 *                                       q(x_2-L)^2 t \\
 *                                               .    \\
 *                                               .    \\
 *                                               .    \\
 *                                       q(x_n-L)^2 t \\
 *                                   \end{split} \right]
 * \qquad
 *  E[(\vec{u}(t)-\hat{u} (t))(\vec{u}(t)-\hat{u} (t) )^T] =  Q
 * \f]
 * Given a measurement model
 * \f[
 *   z(x,t,\omega) = q (x-L) ^2 t + v(x,\omega) \qquad v(x,\omega) \sim N(0,R)
 * \f]
 * $Q = \sigma I - \Phi \Phi^T $ and $R = \gamma I$ are
 * \textit{specifically chosen} such that the state covariance = $I$ at every
 * time step~\footnote{This is an quick check for simple bugs.  Is a more
 * detailed analytic form of the covariance values possible?}. 
 * \f[
 *     P_i^- = \Phi P_{i-1}^- \Phi^T + Q_i
 *           = \Phi I \Phi^T + \sigma I -  \Phi \Phi^T 
 *           =  \sigma I 
 * \qquad
 *     K_i = P_i^- \left( P_i^- + R_i \right)^{-1}
 *         = \sigma I \left( \sigma I + \gamma I \right)^{-1}
 *         = \frac{\sigma}{\sigma +\gamma} I 
 * \f]
 * \f[
 *  P(t^+) =  \sigma I -  \frac{\sigma}{\sigma +\gamma} I  \sigma I 
 *         =   \frac{\gamma \sigma}{\sigma +\gamma} I 
 * \qquad
 * \text{ let } \gamma = \frac{\sigma}{\sigma-1}
 * \f]
 * The \textit{analytical} form of the Kalman Filter may be used as a simple
 * verification of the dense linear algebra routines
 * \f[
 *  x(t^+) =  q (x-L) ^2 t 
 * \qquad
 *  P(t^+) =  I
 * \f]
 */
class VerifyPennesConstantSourceTerm : public PennesBioheatModel 
{
public:
  // constructor
  VerifyPennesConstantSourceTerm( GetPot &controlfile,EquationSystems &es) : 
                PennesBioheatModel(       controlfile,                 es)
    {
      PetscFunctionBegin; 
      PetscFunctionReturnVoid(); 
    }
  
  /** inline functions should be in the header see @ref PDEModelInlineStrategy */
  virtual PetscScalar PennesSource( const unsigned int &, 
                                    const Real&, const Real&, const Real&,
                                    const Point &, const int istep ) 
    { return this->Power[istep] ; }
  virtual PetscScalar dPennesSourcedu( const unsigned int &, const Real&,
                                    const Real&,const Real&, const Real&,
                                    const Point &, const int ) 
                                   {return 0.0;}

  /** set to full damage to test nonlinear solve */
  virtual PetscScalar getInitialDamage(unsigned int,unsigned int,
                                       const Point&,const Parameters& )
       { // for use in regression testing
         return 1.e10;
       }
 
  /** Print the Constitutive Data */
  virtual void printSelf(std::ostream &os) 
    {
      PetscFunctionBegin; 
      // print base class info
      this->PennesBioheatModel::printSelf(os); 
      //os << "VerifyPennesConstantSourceTerm.m_q0  =" <<   m_q0 << std::endl;
      PetscFunctionReturnVoid(); 
    }

};

/**
 * Linear version for testing
 */
class VerifyLinearPennesConstantSourceTerm : public VerifyPennesConstantSourceTerm 
{
public:
  // constructor
  VerifyLinearPennesConstantSourceTerm( GetPot &controlfile,EquationSystems &es) : 
        VerifyPennesConstantSourceTerm(         controlfile,                 es)
    {
      PetscFunctionBegin; 
      PetscFunctionReturnVoid(); 
    }
  
  /** inline functions should be in the header see @ref PDEModelInlineStrategy */
  virtual Real perfusion(    const unsigned int &domainId, const Real &, const Real &) 
    { 
      Real value = w_0[domainId];
      return value ; // faster to multiply
    }
  virtual Real dperfusiondu( const unsigned int &        , const Real &, const Real &, const Real &) 
    {return 0.0;}


};
class VerifyPennesExponentialSourceTerm: public PennesDeltaP1 
{
public:
  // constructor
  VerifyPennesExponentialSourceTerm( GetPot &controlfile,EquationSystems &es) : 
                PennesDeltaP1 (       controlfile,                 es)
    {
      PetscFunctionBegin; 
      PetscFunctionReturnVoid(); 
    }
  
  /**  special case neuman BC */
  virtual PetscScalar residualNeumannBC(const unsigned int i_var,
                                        const Real  &,
                                        int , const int ,
                                        const Point &p ,
                                        const Parameters& parameters)
    {
      PetscScalar time_current = parameters.get<PetscScalar>("time")  ;
      PetscScalar time_previous = parameters.get<PetscScalar>("timePrev");
      PetscScalar mu_a_healthy = parameters.get<PetscScalar>("mu_a_healthy") ; 
      PetscScalar k_0_healthy  = parameters.get<PetscScalar>("k_0_healthy" ) ; 
      PetscScalar xmin = parameters.get<PetscScalar>("xmin") ; 
      PetscScalar xmax = parameters.get<PetscScalar>("xmax") ; 
      const unsigned int u_var= parameters.get<unsigned int>("u_var") ; 
      if( i_var == u_var )
        {
        PetscScalar flux =  0.5* (  
          std::exp( - mu_a_healthy * p(0) )
          * ( std::exp( k_0_healthy * mu_a_healthy * mu_a_healthy * time_current )-1.0)
                     + 
          std::exp( - mu_a_healthy * p(0) )
          * ( std::exp( k_0_healthy * mu_a_healthy * mu_a_healthy * time_previous )-1.0)
                    );
        // have to explicitly account for the normal sign direction
        if( std::abs(p(0) - xmin) < 1.e-4 ) {flux = flux * -1.0;}
        else if( std::abs( p(0) - xmax ) < 1.e-4) {flux = flux * 1.0;}
        else { std::cout << p(0) << " " << xmin << " " << xmax << std::endl << std::flush; libmesh_error();}
        return  flux;
        }
      else
        return  m_NeumannFlux[i_var]; 
    }

  /** Print the Constitutive Data */
  virtual void printSelf(std::ostream &os) 
    {
      PetscFunctionBegin; 
      // print base class info
      this->PennesDeltaP1::printSelf(os); 
      os << "VerifyPennesExponentialSourceTerm.ResidualBC[0] =" <<  ResidualBC[0] << std::endl;
      os << "VerifyPennesExponentialSourceTerm.ResidualBC[1] =" <<  ResidualBC[1] << std::endl;
      os << "VerifyPennesExponentialSourceTerm.ResidualBC[2] =" <<  ResidualBC[2] << std::endl;
      PetscFunctionReturnVoid(); 
    }

};
class VerifyPennesPlanarSourceTerm: public PennesDeltaP1 
{
public:
  // constructor
  VerifyPennesPlanarSourceTerm( GetPot &controlfile,EquationSystems &es) : 
                PennesDeltaP1 (       controlfile,                 es)
    {
      PetscFunctionBegin; 

      // compute needed coeffs and store 
      Real mu_a      = es.parameters.get<PetscScalar>("mu_a_healthy") ; 
      Real mu_s      = es.parameters.get<PetscScalar>("mu_s_healthy") ; 
      Real mu_s_star = mu_s*(1.0e0-m_ForwardScatterFraction);
      Real mu_t_star = mu_a + mu_s_star ;
      Real mu_tr     = mu_a + mu_s*(1.0e0-anfact);
      Real mu_eff    = std::sqrt(3.0e0*mu_a*mu_tr);
      m_alpha     = 3.0e0 * mu_s_star * (mu_t_star + m_ScatteringAsymmetry * mu_a ) 
                       / (mu_eff*mu_eff - mu_t_star * mu_t_star ); 

      Real h = 2.0/3.0 / mu_tr;
      m_beta = - ( m_alpha * (1.0 + m_ReflectionFactor * h * mu_t_star) 
                     +3.0 * m_ReflectionFactor * h * m_ScatteringAsymmetry * mu_s_star  )
                  / (1.0 + m_ReflectionFactor * h * mu_eff) ; 

      // store
      es.parameters.set<PetscScalar>("mu_t_star") = mu_t_star ;
      es.parameters.set<PetscScalar>("mu_eff")    = mu_eff    ;
      es.parameters.set<PetscScalar>("alpha")     = m_alpha     ;
      es.parameters.set<PetscScalar>("beta")      = m_beta      ;

      PetscFunctionReturnVoid(); 
    }
  
  /** infinite  planar source at z=0 */
  virtual Real ExternalIrradiance( int domainId, const Real&, const Real&, const Point& qpoint, const int timeId)
  { 
     Real mu_t_star=mu_a_0[domainId]
                   +mu_s_0[domainId]*(1.0e0-m_ForwardScatterFraction);
     return Power[timeId] * (1.0 - m_SpecularReflectance) 
                          * std::exp( - mu_t_star * qpoint(0) ) ;
  }

  /**  special case neuman BC */
  virtual PetscScalar residualNeumannBC(const unsigned int i_var,
                                        const Real  &,
                                        int , const int ,
                                        const Point &p ,
                                        const Parameters& parameters)
    {
      Real mu_a      = parameters.get<PetscScalar>("mu_a_healthy") ; 
      Real mu_s      = parameters.get<PetscScalar>("mu_s_healthy") ; 
      const unsigned int z_var= parameters.get<unsigned int>("z_var") ; 
      Real mu_s_star = mu_s*(1.0e0-m_ForwardScatterFraction);
      Real mu_t_star = mu_a + mu_s_star ;
      Real mu_tr     = mu_a + mu_s*(1.0e0-anfact);
      Real mu_eff    = std::sqrt(3.0e0*mu_a*mu_tr);
      if( i_var == z_var )
        {
        PetscScalar flux = - mu_t_star * m_alpha * std::exp( - mu_t_star * p(0) )   
                           - mu_eff    * m_beta  * std::exp( - mu_eff    * p(0) ) ;
        return -flux; // MINUS Sign FUN! due to IBP for diffusion term
        }
      else
        return  m_NeumannFlux[i_var]; 
    }

  /** Print the Constitutive Data */
  virtual void printSelf(std::ostream &os) 
    {
      PetscFunctionBegin; 
      // print base class info
      this->PennesDeltaP1::printSelf(os); 
      os << "VerifyPennesPlanarSourceTerm.m_alpha =" <<  m_alpha << std::endl;
      os << "VerifyPennesPlanarSourceTerm.m_beta  =" <<  m_beta  << std::endl;
      PetscFunctionReturnVoid(); 
    }
private:
   Real m_alpha ,m_beta    ;
};
#endif  // GTEST_PENNES_REGRESSION_H_
