#ifndef __PennesModel_h
#define __PennesModel_h
#include "pdeBaseClass.h"


/**@ingroup TreatmentPlanning 
 * The PennesBioheatModel class is intended to be a default
 * class for simulating MRgLITT. It is intended to hold constitutive data ONLY
 * and is intended for use with either finite difference or finite element
 * solvers. The default is a linear Pennes Model and
 * arrhenius thermal dose.  Derived classes should add different forms of
 * constitutive nonlinearities and constitutive relations that depend on
 * damage.
 * 
 * @section PennesModelDerivation The Pennes Model
 * @htmlonly
 *   see pdf for in depth derivation of the Pennes Model
 * @endhtmlonly
 * 
 * @latexonly
 * This section focuses on the complete derivation of the Pennes bioheat 
 * transfer equation modified to include a laser source term and heterogeneous
 * nonlinear tissue properties. 
 * Consider a material body $\Omega \subset \mathbb{R}^3$. The
 * time span of interested is denoted $(0,T]$ and the space-time
 * domain is denoted $\Omega_T \equiv \Omega \times (0,T] $.
 * The biological domain $\Omega$ is 
 * assumed to be a Lipschitz domain to provide meaning to the trace
 * of the functions defined in the variational problem. 
 * The integral equations for conservation of
 * energy are now invoked under the following assumptions:
 * 
 * 
 * \begin{itemize}
 *   \item[-] no mass flux across boundary, $\partial \Omega$ 
 *   \item[-] no motion, deformation, or applied forces
 *   \item[-] muscle, fat blood, etc. is represented as a single homogenized tissue
 *            that may have spatially varying properties.
 *   \item[-] the presence of the blood within the homogenized tissue acts as a 
 * volumetric heating term at every point within the body $\Omega$
 *   \item[-] laser heating acts as a spatially and time varying heat source term
 * \end{itemize}
 * Note that the conservation of mass and momentum are identically 
 * satisfied under the first two assumptions.
 * In words, the conservation of energy states that the change of energy within 
 * an arbitrary subvolume of the 
 * body, $\hat{\omega} \subset \Omega$, is equal to the energy lost across
 *  the boundary, 
 * $\partial \hat{\omega}$, plus energy generated within the body,
 * Figure~\ref{derivebioheat}.
 * 
 * \begin{equation}\label{pennesintegral} 
 *   \frac{d}{dt} \int_{\hat{\omega} } \rho e \; dx  =
 *   - \int_{\partial \hat{\omega} } \textbf{q} \cdot \textbf{n} \;dA
 *   + \int_{\hat{\omega}} Q_{blood} + Q_{laser} \; dx  
 * \end{equation}
 * where $Q_{blood}$ is the volumetric heating due to blood perfusion
 * in the tissue and $Q_{laser}$ is the volumetric heating of the laser
 * source term. The density is denoted $\rho$, $e$ is the internal energy, and 
 * $\textbf{q}$ is the heat flux.
 * 
 * Integrals in \eqn{pennesintegral} may be simplified as follows.
 * \[
 *  \Omega \text{ fixed,} \qquad \rho = \rho(\textbf{x}) \qquad \Rightarrow
 * \frac{d}{dt} \int_{\hat{\omega} } \rho e \; dx  =
 * \int_{\hat{\omega} } \rho \frac{\partial e}{\partial t} \; dx  
 * \]
 * \[
 *  \text{Divergence theorem} \qquad \Rightarrow
 *    \int_{\partial \hat{\omega} } \textbf{q} \cdot \textbf{n} \;dA =
 * \int_{\hat{\omega} } \nabla \cdot \textbf{q} \; dx  
 * \]
 * Since $\hat{\omega}$ is arbitrary, this implies that the local 
 * form of the principle of conservation of energy is
 * \[
 * \rho \frac{\partial e}{\partial t} =
 *  -\nabla \cdot \textbf{q} 
 *  + Q_{blood} + Q_{laser}  \qquad \text{a.e. in } \Omega
 * \]
 * For the constitutive equations take
 * \begin{itemize}
 *   \item[$\bullet$] $e = c_p \; u(\textbf{x},t)$, where $c$ is the specific heat 
 *   $ \left[\frac{J}{kg \; K}\right] $ and
 * $u(\textbf{x},t)$ denotes the temperature of the homogenized tissue
 *     at point $\textbf{x} \in \bar{\Omega}$ at time $t$.
 *     Since there is no deformation the material may be considered
 *     incompressible and the specific heat at constant pressure is the same as
 *     the specific heat at constant volume. The constitutive equation follows 
 *     from the definition of the  specific heat
 *        $c_p = \left. \frac{\partial e}{ \partial u} \right|_p $
 * 
 *   \item[$\bullet$] $ q = -k(u,\textbf{x}) \nabla u $, 
 *    where $k(u,\textbf{x})$ is the scalar coefficient
 *  of thermal conductivity.
 * Assume a non-linear form for,
 * $k \left[\frac{J}{s \cdot m \cdot K}\right]$, 
 * as shown in Figure~\ref{coeffinvert}.
 * \[ 
 *  k(u,\textbf{x})   =  k_0(\textbf{x}) 
 *                    + k_1 \; \text{atan} ( \tilde{k}_2 ( u - \hat{k}_3))  
 * \] 
 * where $k_0(\textbf{x}) \left[\frac{J}{s \cdot m \cdot K}\right]$,
 * $k_1 \left[\frac{J}{s \cdot m \cdot K}\right] $,
 * $\tilde{k}_2  \left[\frac{1}{K}\right] $,
 * $\hat{k}_3  \left[K\right] \in \mathbb{R}$.
 * $k_0(\textbf{x})$ is allowed to vary over
 * the spatial dimension to capture the biological tissue heterogeneity.
 * The media is thus assumed to be thermally isotropic (referring to heat flux)
 * but spatially heterogeneous.
 * 
 * \end{itemize}
 * For the volumetric heating terms:
 * \begin{itemize}
 *   \item[$\bullet$] 
 *  The scattering and absorption of photons is the
 *  fundamental mechanism of heat transfer for the laser source term.
 *   The following form of the laser source term may be derived from transport
 * theory using the diffusion approximation~\cite{Welch}.
 * \[ \begin{split}
 * Q_{laser}(\textbf{x},t)  & =
 * 3 P(t) \mu_a \mu_{tr} \frac{\exp(-\mu_{eff} \|\textbf{x} -\textbf{x}_0\|) }{4\pi \|\textbf{x}-\textbf{x}_0\|} \\ \\
 *   \mu_{tr}  = \mu_a +  &\mu_s (1-g) \qquad
 *   \mu_{eff} = \sqrt{3 \mu_a \mu_{tr}}
 * \end{split} \]
 *   $P(t) $ is the laser power as a function of time, $\mu_a$ and $\mu_s $ 
 * are laser coefficients related to laser wavelength and give probability of
 *  absorption and scattering of photons, respectively. The anisotropic factor is 
 * denoted $g$ and $\textbf{x}_0$ denotes the position of laser photon source.
 *   \item[$\bullet$] Assume, as Pennes did in his original work in 1948, that
 * arterial blood acts as a heat source and blood is isothermal until it 
 * reaches the capillaries with a diameter of $\approx$ 8 $\mu$m.
 * Figure~\ref{tempdistribution} shows a schematic of the blood temperature 
 * distribution in the circulatory system. As shown, blood is equilibriated with 
 * surrounding tissue by the time it reaches arterioles with diameters 
 * $\approx$ 60 $\mu$m, but we will use the assumption of Pennes irregardless.
 *   \item[] 
 * The volumetric heating of the tissue by perfused blood is driven by
 * mass flow of blood to tissue and the temperature difference between 
 * arterial temperature, $u_a \left[ K \right]$, and local homogenized 
 * tissue temperature, $u \left[ K \right]$.
 * \[
 *   Q_{blood} =  - \omega c_{blood} (u - u_a )
 * \]
 * where $\omega \left[\frac{kg}{s \; m^3}\right]$ is the perfusion coefficient. 
 * Perfusion has units of mass flow (of blood to tissue) per unit volume of
 * tissue.  Empirical evidence suggests that
 *  $\omega \left[\frac{kg}{s \; m^3 }\right]$, 
 * is a smooth, monotone increasing bounded function of
 * temperature, Figure~\ref{coeffinvert}.
 * \[
 *  \omega(u,\textbf{x})   =  \omega_0(\textbf{x}) 
 *             + \omega_1 \; \text{atan} ( \tilde{\omega}_2 ( u - \hat{\omega}_3))
 * \]
 * where $\omega_0 \left[\frac{kg}{s \; m^3 }\right]$,
 * $\omega_1 \left[\frac{kg}{s \; m^3}\right] $,
 * $\tilde{\omega}_2  \left[\frac{1}{K}\right] $,
 * $\hat{\omega}_3  \left[K\right] \in \mathbb{R}$.
 * Note that $\omega_0(\textbf{x})$ is allowed to vary over
 * the spatial dimension as the blood perfusion within the necrotic core of
 * a cancerous tumor or the blood perfusion within a damaged tissue
 * is expected to be significantly lower than the 
 * surrounding healthy tissue. The specific heat of blood is denoted
 * $c_{blood} \left[\frac{J}{kg \cdot K}\right]$ 
 * \end{itemize}
 * 
 * Let $\beta$ denote an array of all the bioheat transfer model parameters,
 * \[ \boxed{
 * \beta \equiv ( k_0(\textbf{x}), k_1, \tilde{k}_3, \hat{k}_3,
 *           \omega_0(\textbf{x}), \omega_1, \tilde{\omega}_3, \hat{\omega}_3,
 *            P(t),\mu_a,\mu_s,\textbf{x}_0 )
 * } \]
 * Combining all terms, the problem statement is as follows:
 * \begin{itemize}
 * \item[]
 * \begin{minipage}{\linewidth}
 * Given a set of bioheat transfer model parameters, $\beta$, 
 * \[ \boxed{
 * \beta \equiv ( k_0(\textbf{x}), k_1, \tilde{k}_3, \hat{k}_3,
 *           \omega_0(\textbf{x}), \omega_1, \tilde{\omega}_3, \hat{\omega}_3,
 *            P(t),\mu_a,\mu_s,\textbf{x}_0 )
 * } \]
 *  find the spatially and temporally varying temperature field
 *    $ u(\beta,\textbf{x},t)$  such that
 * \begin{equation} \label{classicalpennes} \begin{split}
 * & \rho  c_p \frac{\partial u}{\partial t}
 *  -\nabla \cdot (k(u,\textbf{x},\beta) \nabla u) \\
 * &
 *  + \omega(u,\textbf{x},\beta) c_{blood} (u - u_a )
 *  = Q_{laser}(\beta,\textbf{x},t) \quad \text{in } \Omega
 * \end{split} \end{equation}
 * given the Cauchy and Neumann boundary conditions
 *   \[ -  k(u,\textbf{x},\beta) \nabla u \cdot \textbf{n} = h (u - u_{\infty})
 *         \qquad \text{on } \partial \Omega_C
 *   \]
 *   \[ -  k(u,\textbf{x},\beta) \nabla u \cdot \textbf{n} = \mathcal{G}
 *         \qquad \text{on } \partial \Omega_N
 *   \]
 * and the initial condition
 *  \[
 *    u(\textbf{x},0) = u^0 \qquad \text{in } \Omega
 *  \]
 * \end{minipage}
 * \end{itemize}
 * The domains of the Cauchy and Neumann boundary conditions satisfy
 * \[
 * \partial \Omega = \bar{\partial \Omega_C \cup \partial \Omega_N} 
 * \qquad \qquad
 * \partial \Omega_C \cap \partial \Omega_N = \emptyset 
 * \]
 * On the Cauchy boundary, $h$ is the coefficient of cooling and $u_\infty$
 * is the ambient temperature. The prescribed heat flux on the Neumann boundary
 * is denoted $\mathcal{G}$.
 * Notice the explicit dependence of thermal conductivity,
 * $k(u,\textbf{x},\beta)$, perfusion, $\omega(u,\textbf{x},\beta)$, and
 * laser source, $Q_{laser}(\beta,\textbf{x},t)$, on $\beta$.
 * A summary of the  nomenclature used throughout the remainder
 * of this dissertation is given in Table~\ref{nomenclature}.
 * \begin{table}[H] \centering
 * \caption{Nomenclature.} \label{nomenclature}
 * \begin{tabular}{lll}
 *    & &\\
 *    $u$         =  temperature & 
 *    $h$         =  coefficient of cooling&  
 *    $c_p$       =  tissue specific heat \\
 *    $\omega$    =  blood perfusion  &
 *    $k$         =  thermal conductivity&  
 *    $c_{blood}$ =  blood specific heat \\
 *    $\rho$      =  tissue density &
 *    $u_\infty$  =  ambient temperature&
 *    $Q_{laser}$ =  laser source term\\
 * \end{tabular}
 * \end{table}
 * 
 * 
 * Typical values of the parameters appearing in the Pennes model obtained by 
 * Rylander~\cite{Rylander2005}  are given in 
 * Tables~\ref{textbooklaser}~and~\ref{textbookthermal}.
 * Later, these parameters will be determined from 
 * real-time thermal MRTI data using calibration methods based 
 * on inverse analysis.
 * 
 * \begin{table}[h]
 * \caption{Model Coefficients (From Rylander~\cite{Rylander2005})}
 * \label{textbooklaser}
 * \centering
 * \begin{tabular}{|l|c||l|c|} \hline
 * & & & \\ 
 * Laser Wavelength &  810 nm $\qquad$ &
 * $\rho_{tissue}$  &  1045 $ \frac{kg}{m^3}$ \\
 * & & & \\ \hline
 * & & & \\ 
 * $\mu_s$          &  14.74cm   $\qquad$&
 * $\rho_{blood}$   &  1058 $ \frac{kg}{m^3}$ \\ 
 * & & & \\ \hline
 * & & & \\ 
 * $\mu_a$          &  0.046 cm $\qquad$&
 * $c_p^{tissue}$   &  3600 $ \frac{J}{kg \cdot K}$ \\
 * & & & \\ \hline
 * & & & \\ 
 * $u_a$            &  310 K  $\qquad$&
 * $c_p^{blood}$    &  3840 $ \frac{J}{kg \cdot K}$ \\ 
 * & & & \\ \hline
 * \end{tabular}
 * \end{table}
 * 
 *  
 * \begin{table}[h]
 * \caption{Model Coefficients}
 * \label{textbookthermal}
 * \centering
 * \begin{tabular}{|l|c|l|c|} \hline
 * \multicolumn{2}{|c|}{ } & \multicolumn{2}{c|}{ } \\
 * \multicolumn{2}{|c|}{$k(u) = k_0 + k_1 \cdot \text{atan}(k_2(u-k_3)) $} & 
 * \multicolumn{2}{c|}{$\omega(u) = \omega_0 + \omega_1 \cdot \text{atan}(\omega_2(u-\omega_3)) $} 
 *         \\
 * \multicolumn{2}{|c|}{ } & \multicolumn{2}{c|}{ } \\ \hline
 * $k_0$            &   $ 0.6489 \left[\frac{J}{s \cdot m \cdot K}\right]  \qquad$ &
 * $\omega_0$       &   $ 0.6267  \left[\frac{kg}{s \; m^3 }\right]       $ \\
 * & & & \\ \hline
 * $k_1$            &   $ 0.0427   \left[\frac{J}{s \cdot m \cdot K}\right]   \qquad$ &
 * $\omega_1$       &   $ -0.137    \left[\frac{kg}{s \; m^3}\right]         $ \\
 * & & & \\ \hline
 * $k_2$            &   $ 0.02529  \left[\frac{1}{K}\right]  \qquad$ &
 * $\omega_2$       &   $ 2.35893    \left[\frac{1}{K}\right]        $ \\
 * & & & \\ \hline
 * $k_3$            &   $ 315.314  \left[K\right] \qquad$ &
 * $\omega_3$       &   $ 314.262  \left[K\right]      $ \\
 * & & & \\ \hline
 * \end{tabular}
 * \end{table}
 * 
 * @endlatexonly
 *
 * @section ArrheniusDamageEquation Arrhenius Thermal Dose
 * 
 * \f[
 *  \Omega = \int_0^T A \exp \left( \frac{-E_a}{R \; u(t)}  \right) \; dt
 * \f]
 *
 */
class PennesBioheatModel : public PDEModelBaseClass 
{
public:
  // constructor
  PennesBioheatModel (GetPot &, EquationSystems & );

  /** Print the Constitutive Data */
  virtual void printSelf(std::ostream& );

  /** set power data info */
  void GetPowerData(GetPot &);

  /** Setup auxillary data structures needed for solve */
  virtual void SetupState(EquationSystems &es)
   {
     w_0.OneVarSetup(es);
     k_0.OneVarSetup(es);
   }

  /** scatter global solution locally */
  virtual void ScatterParametersLocally()
   {
     w_0.ScatterGlobalToAll();
     k_0.ScatterGlobalToAll();
   }
  /**
   * From jstafford@mdanderson.org Mon Jun 28 11:06:54 2010
   * Date: Mon, 28 Jun 2010 11:06:49 -0500
   * 
   * So, many authors use a "perfusion model" which at least shuts down with
   * damage. Have you implemented, or investigated if it has an impact on the
   * results?
   * 
   *  Perfusion constitutive equation and its solution derivative 
   *  Assume the nonlinearity is wrt to damage \in [0,1] of the form
   *     w   = w_0 +  damage(  u) ( w_1 - w_0 )
   *    dwdu =       DdamageDu(u) ( w_1 - w_0 )
   *
   *     damage(  u) =  arrhenius   /         (std::log(2.0) + arrhenius )
   *    DdamageDu(u) = DarrheniusDu / std::pow(std::log(2.0) + arrhenius ,2)
   *  inline functions should be in the header see @ref PDEModelInlineStrategy
   */
  virtual Real  perfusion(  const unsigned int &, const Real &, const Real &);
  virtual Real dperfusiondu(const unsigned int &, const Real &, const Real &, const Real &); 

  /** 
   *  Should be inlined reaction term to simplicity and backward compatibility
   *  inline functions should be in the header see @ref PDEModelInlineStrategy
   */
  Real  PennesReactionTerm(  const unsigned int &, const Real &, const Real &) ;
  Real dPennesReactionTermdu(const unsigned int &, const Real &, const Real &, const Real &) ;

  /* Thermal conductivity constitutive equation and its solution derivative 
   * Assume the nonlinearity is linear wrt to temperature
   */
  Real  ThermalConductivity(  const unsigned int &, const Real &, const Real &) ;
  Real dThermalConductivitydu(const unsigned int &, const Real &, const Real &, const Real &) ;
  
  /** initial condition
   *  return body  temp on subdomain 0 1 
   *  return probe temp on subdomain 2 3 
   */
  PetscScalar getInitialTemperature(unsigned int domainId )
     { 
       return ( domainId < m_ApplicatorDomain ) ? m_bodyTemp : m_probeTemp; 
     }
  virtual PetscScalar getInitialTemperature(unsigned int domainId, 
                                            unsigned int,
                                            const Point&  , const Parameters& )
     { 
       return this->getInitialTemperature(domainId);
     }

  /** get power */
  PetscScalar getPower(unsigned int idpow) {return this->Power[idpow];}

  /** set power */
  PetscErrorCode setPower(unsigned int idpow,PetscScalar NewPower)
    {
      this->Power[idpow]=NewPower;
      this->Power.printStdVector( std::cout, "PennesModel:  Power[");
      return 0;
    }

  Gradient BulkFluidFlow(unsigned int idDomain)
         {return m_BulkFluidFlow[idDomain] ; }

  Gradient DiffusionDirection(unsigned int idDomain)
         {return m_DiffusionDirection[idDomain] ; }

  void ThermalDose(Vec damage, PetscScalar deltat)
     {
       PetscFunctionBegin;
       // input should be temperature in degC
       // convert to Kelvin
       VecShift(damage,m_BaseTemperature);
       // 1/T 
       VecReciprocal(damage);
       // - E_a/ (RT) 
       VecScale(damage, - m_FrequencyFactor/m_GasConstant );
       // exp ( - E_a/ (RT)  ) 
       VecExp(damage);
       // A exp ( - E_a/ (RT)  ) \Delta t
       VecScale(damage, m_ActivationEnergy * deltat );
       PetscFunctionReturnVoid();
     }

  void DoseDerivative(Vec derivative, Vec workVec)
     {
       PetscFunctionBegin;
       // input should be damage and temperature: 
       //     derivative = A exp ( - E_a/ (RT)  ) \Delta t
       //     workVec    = temperature
       // 
       // convert to Kelvin
       VecShift(workVec,m_BaseTemperature);
       //  T^2
       VecPointwiseMult(workVec,workVec,workVec);
       //  1/T^2
       VecReciprocal(workVec);
       //  E_a/ (RT^2) 
       VecScale(workVec,  m_FrequencyFactor/m_GasConstant );
       // A E_a/ (RT^2) exp ( - E_a/ (RT)  ) \Delta t
       VecPointwiseMult(derivative,derivative,workVec);
       PetscFunctionReturnVoid();
     }
protected:
  Real rho,            ///<  \f$ \rho   \f$ 
       specific_heat;  ///<  \f$ c_p    \f$ 

  // convection term
  std::vector<RealGradient> m_BulkFluidFlow;

  // diffusion direction 
  std::vector<RealGradient> m_DiffusionDirection;
  // perfusion parameters
  Real bloodspecificheat,  ///<  \f$ c_{blood}   \f$ 
       u_artery;           ///<  \f$ u_a   \f$ 

  // perfusion and thermal conductivity
  spatialParameter k_0,  ///<  \f$ k_0    \f$ 
                   w_0;  ///<  \f$ \omega_0   \f$ 
  Real w_1,         ///<  \f$ \omega_1 \f$ 
       k_1,              ///<  \f$ k_1    \f$ 
       m_ActivationEnergy,      ///<  \f$ A      \f$ 
       m_FrequencyFactor ,      ///<  \f$ E_a    \f$ 
       m_GasConstant     ,      ///<  \f$ R      \f$ 
       m_bodyTemp        ,      ///<  initial temperature
       m_probeTemp       ,      ///<  probe temperature
       m_BaseTemperature ;      ///<  conversion to Kelvin

  // the data structure is setup as follows and ASSUMES equispace time 
  //     distances, IDEAL_DT   \forall i
  //  
  //                                           *NOTE* closed at beginning
  // time = 0    ---------                        BUT open   at end
  //                 |                                |
  //                 |                               \|/
  //              Power[0]    power between time = [0,1)  is Power[0]
  //                 |            
  //                 |
  // time = 1    ---------
  //                 |
  //                 |          
  //              Power[1]    power between time = [1,2)  is Power[1]   
  //                 |         
  //                 |
  // time = 2    ---------
  //                 |
  //                 |        
  //              Power[2]    power between time = [2,3)  is Power[2]
  //                 |       
  //                 |
  // time = 3    ---------
  //           .
  //           .
  //           .
  //           .

  discreteParameter Power;
 
private: 
  /** identify applicator domain */
  unsigned int m_ApplicatorDomain  ;
  /** The type of the parent  */
  typedef PDEModelBaseClass Parent;
};

/**
 * base class for laser treatment 
 */
class PennesStandardDiffusionApproximation : public PennesBioheatModel 
{
public:
  // constructor
  PennesStandardDiffusionApproximation(GetPot &, EquationSystems & );

  /** Print the Constitutive Data */
  virtual void printSelf(std::ostream& );

  /** Setup auxillary data structures needed for solve */
  virtual void SetupState(EquationSystems &es)
   {
     Parent::SetupState(es);
     mu_a_0.OneVarSetup(es);
     mu_s_0.OneVarSetup(es);
   }

  /** scatter global solution locally */
  virtual void ScatterParametersLocally()
   {
     Parent::ScatterParametersLocally();
     mu_a_0.ScatterGlobalToAll();
     mu_s_0.ScatterGlobalToAll();
   }
  /**
   * return the source given the position coordinates and the time step
   */
  virtual PetscScalar PennesSource( const unsigned int &domainId,
                     const Real &u_theta, const Real &damage, 
                     const Real &, const Point &qpoint, const int istep ) 
   {
      return this->absorption(domainId,u_theta,damage) *
        this->interstitialFluence(domainId,u_theta,damage,qpoint,istep);
   }

  virtual PetscScalar dPennesSourcedu( const unsigned int &,const Real&,
   const Real &, const Real &, const Real &, const Point &, const int ) 
   {
      // FIXME optical non linearities not coded
      // TODO  code optical non linearities
      return 0.0;
   }

  /* ODA laser source */
  virtual Real interstitialFluence(int ,const Real &, const Real &, const Point& , const int );

  /** nodal dirichlet data */
  virtual bool dirichletNodalData( const Point &);

  /**
   *  absorption and scatter constitutive equation and its solution derivative 
   *  Assume the nonlinearity is wrt to damage \in [0,1] of the form
   *     mu_a   = mu_a_0 + damage(  u) ( mu_a_1 - mu_a_0   )
   *    dmu_adu =         DdamageDu(u) ( mu_a_1 - mu_a_0   )
   *
   *     mu_s   = mu_s_0 + damage(  u) ( mu_s_1 - mu_s_0   )
   *    dmu_sdu =         DdamageDu(u) ( mu_s_1 - mu_s_0   )
   *
   *     damage(  u) =  arrhenius   /         (std::log(2.0) + arrhenius )
   *    DdamageDu(u) = DarrheniusDu / std::pow(std::log(2.0) + arrhenius ,2)
   *  inline functions should be in the header see @ref PDEModelInlineStrategy
   */
  Real absorption(   const unsigned int &, const Real &, const Real &);
  Real dabsorptiondu(const unsigned int &, const Real &, const Real &, const Real &); 
  Real scattering(   const unsigned int &, const Real &, const Real &);
  Real dscatteringdu(const unsigned int &, const Real &, const Real &, const Real &); 

  /** setup the laser domain */
  PetscErrorCode UpdateLaserPosition(PetscScalar,PetscScalar,PetscScalar,
                                     PetscScalar,PetscScalar,PetscScalar);
protected:

  PetscInt m_ProbeDomain;

  // laser parameters
  spatialParameter  mu_a_0, ///<  \f$\mu_{a,0}\f$ 
                    mu_s_0; ///<  \f$\mu_{s,0}\f$ 

  Real anfact; ///<  \f$ g   \f$ 
  Real mu_a_1,         ///<  \f$ \mu_{a,1} \f$ 
       mu_s_1;         ///<  \f$ \mu_{s,1} \f$ 
  discreteParameter x_0, y_0, z_0;
  // volume fractions for power source
  std::vector< PetscScalar > m_volumeFraction;
  // store associated elements 
  std::vector< Elem* > m_elementProbe;
 
  /** identify dirichlet domain */
  PetscScalar m_laserTip[3]; 
  PetscScalar m_unitVec[3];
private: 
  /** identify dirichlet domain */
  PetscScalar m_diffusingradius;
  PetscScalar m_diffusinglength;
  PetscInt    m_nelemtip;

  /** The type of the parent  */
  typedef PennesBioheatModel Parent;
};
/**@ingroup TreatmentPlanning 
 * @htmlonly
 *   see pdf for in depth derivation of the RF coupled Pennes Model
 * @endhtmlonly
 * 
 * @latexonly
 * The PennesVoltage class is intended for for simulating RFA with an
 * applied voltage. It is intended to hold constitutive data ONLY and is
 * intended for use with either finite difference or finite element solvers. 
 * The equations governing the resistive heating that provides the thermal
 * source for the bioheat equation~\eqn{bioheatmodel} may be obtained from
 * Maxwell's equation~\cite{demkowicz2006cha} using Faraday's
 * law~\eqn{FaradayLaw}, Ampere's law~\eqn{AmpereLaw}, and the continuity
 * of free charge~\eqn{ChargeContinuity}.
 * \begin{equation} \label{FaradayLaw}
 *  \nabla \times \vec{E} = -\frac{\partial}{\partial t} \vec{B}
 * \end{equation}
 * \begin{equation} \label{AmpereLaw}
 *  \nabla \times \vec{B} = 
 *         \mu \vec{J} + \mu \epsilon \frac{\partial}{\partial t} \vec{E}
 * \end{equation}
 * \begin{equation} \label{ChargeContinuity}
 *     \nabla \cdot \vec{J} + \frac{\partial}{\partial t} \rho_{charge} = 0 
 * \end{equation}
 * Here, the medium is assumed linear and isotropic with homogeneous
 * permittivity, $\epsilon$, and permeability, $\mu$.
 * The charge density is denoted $\rho_{charge}$.
 * The electric and magnetic fields are denoted
 * $\vec{E}$ and $\vec{B}$, respectively. 
 * The current density, $J$, is assumed related to the electric field
 * through Ohm's law, $\vec{J}=\sigma \vec{E}$, where $\sigma$ denotes
 * the electrical conductivity.
 * In the quasistatic case, the electric field is approximated as
 * irrotational and Maxwell's system decouples
 * \[ 
 *  \nabla \times \vec{E}   = 0 
 *  \qquad \Rightarrow \qquad
 *  \nabla \cdot \sigma \vec{E}   = 0 
 * \]
 * For an arbitrary control volume, the Poynting vector defines the amount
 * of electromagnetic power, $Q_{em}$, which crosses the surface
 * \[
 *  Q_{em} = \int_{\text{surface}} 
 *                 \left( \vec{E} \times \vec{B} \right) dA
 * \]
 * From the conservation of energy for an arbitrary control volume, the
 * mechanical energy and electromagnetic energy are coupled~\cite{pao1967fd}
 * \[
 *  \rho  c_p \frac{\partial u}{\partial t} + 
 *            \frac{\partial u_{em}}{\partial t} + 
 *   -\nabla  \cdot k \nabla u + \omega c_{blood} (u - u_a )
 *  = -Q_{em}(x,t) \\
 * \]
 * Here, the temperature is denoted $u$,
 * $\rho$  is the density of the continuum, which is liver tissue in
 * our case, and  $c_{blood}$ is the specific heat of blood. 
 * The perfusion coefficient, $\omega$, and thermal conductivity, $k$,
 * are assumed linear and
 * spatially homogeneous model coefficients. 
 * Combining the conservation of energy with Poynting's 
 * theorem~\cite{jackson1999cer} yields
 * \[
 *  \rho  c_p \frac{\partial u}{\partial t} + 
 *            \frac{\partial u_{em}}{\partial t} + 
 *   -\nabla  \cdot k \nabla u + \omega c_{blood} (u - u_a )
 *  = \vec{E} \cdot \vec{J} + \frac{\partial u_{em}}{\partial t} 
 * \]
 * From Poynting's theorem, the energy density of the
 * electromagnetic field does not change
 * \[
 * \begin{split}
 *  \frac{\partial u_{EM}}{\partial t}  
 *     & = 
 *         - \nabla \cdot \frac{1}{\mu}\left( \vec{E} \times \vec{B} \right) 
 *         - \vec{J} \cdot \vec{E}
 *  \\
 *     &  = 
 *         - \frac{1}{\mu}\left( \vec{B} \cdot \nabla \times \vec{E} 
 *                             - \vec{E} \cdot \nabla \times \vec{B} \right) 
 *         - \vec{J} \cdot \vec{E}
 *  \\
 *     &  = 
 *           \frac{1}{\mu} \vec{E} \cdot \left( \mu J \right) 
 *         - \vec{J} \cdot \vec{E}
 *  \\
 *     &  =  0 
 * \end{split}
 * \]
 * hence, the energy input into the tissue is
 * dissipated through the resistive heating.
 * \[
 * Q_{RF}(x,t) = \vec{E} \cdot \vec{J} 
 *             = \sigma \vec{E} \cdot \vec{E} 
 * \]
 * The relationship of the voltage, $\phi$, as the scalar potential for the
 * electric field, $\nabla \phi = \vec{E}$, leads to the time dependent system
 * of coupled equations governing the temperature and voltage.
 * \begin{equation} \label{bioheatmodel}
 * \left.
 * \begin{split}
 *  \rho  c_p \frac{\partial u}{\partial t}
 *  & -\nabla  \cdot k \nabla u + \omega c_{blood} (u - u_a )
 *  = Q_{RF}(x,t) \\
 *  & Q_{RF}(x,t) = \sigma |\nabla \phi|^2, \qquad 
 *         \nabla \cdot \sigma \nabla \phi = 0
 * \end{split}
 * \right\}
 * \end{equation}
 * Dimensionally, the units of $Q_{RF}$ are consistent with power density.
 * \[
 * \frac{S}{m} \frac{V^2}{m^2} 
 *   = \frac{\frac{A}{V}}{m} \frac{V^2}{m^2}  
 *   = \frac{A V}{m^3} 
 *   = \frac{ \frac{C}{s} \frac{J}{C} }{m^3} 
 *   = \frac{W}{m^3} 
 * \]
 * Neumann and Dirichlet boundary conditions are used on the surface of
 * the biological domain of interest, $\partial \Omega$, for the 
 * the temperature and voltage, respectively. 
 * \[
 * \left.
 * \begin{split}
 *  - k(u)      \nabla   u  \cdot \textbf{n}  & = 0 
 * \\
 *   \phi = 0
 * \end{split} 
 * \right\}
 * \qquad \text{on } \partial \Omega
 * \]
 * A summary of the 
 * constitutive data used is presented in Table~\ref{modeldata}.
 * 
 * \begin{table}[h]
 * \caption{Constitutive Data~\cite{Handbook05,stauffer2003paa,mcnichols2004mtb}}\label{modeldata}
 * \centering
 * \begin{tabular}{|c|c|c|c|c|c|c|c|c|} \hline
 * $\sigma$  &  $k$   &  $\rho$                &  $u_a$       &
 * $c_{blood}$                  &  $c_p$  & $A$ & $E_a$  & $R$  \\ \hline
 *  0.69 $\frac{S}{m}$  &  0.5  $\frac{W}{m K}$ &  1045 $\frac{kg}{m^3}$ &  308  $K$    &  3840 $ \frac{J}{kg \cdot K}$ &  3600 $ \frac{J}{kg \cdot K}$ 
 *   &  3.1e98 $ \frac{1}{s}$ &  6.28e5 $ \frac{J}{mol}$ &8.314 $\frac{J}{kg \; mol}$ \\ \hline
 * \end{tabular}
 * \end{table}
 *
 * @endlatexonly
 * 
 */
class PennesVoltage : public PennesBioheatModel 
{
public:
  
  /** Constructor. */
  PennesVoltage(GetPot &, EquationSystems & );

  /** Print the Constitutive Data */
  virtual void printSelf(std::ostream& );

  /** Setup auxillary data structures needed for solve */
  virtual void SetupState(EquationSystems &es)
   {
     Parent::SetupState(es);
     m_ElectricConductivity_0.OneVarSetup(es);
   }
  /**  return resistive heating source */
  virtual PetscScalar PennesSource( const unsigned int &domainId,
                     const Real &u_theta, const Real &damage, 
                     const Real &electricPower, const Point &, const int) 
   {
      return this->ElectricConductivity(domainId,u_theta,damage) *
             electricPower; 
   }

  virtual PetscScalar dPennesSourcedu( const unsigned int &domainId,
          const Real &u_theta, const Real &damage, const Real &DdamageDu, 
          const Real &electricPower, const Point &, const int) 
   {
      return this->dElectricConductivitydu(domainId,u_theta,damage,DdamageDu) *
             electricPower; 
   }

  /**  return spatially dependent electrical conductivity at a given point */
  virtual PetscScalar ElectricConductivity(const unsigned int &, 
                                           const PetscScalar &,  
                                           const PetscScalar &) ; 

  /**  return derivative of electrical  conductivity WRT temperature */
  virtual PetscScalar dElectricConductivitydu(const unsigned int &, 
                                              const PetscScalar &, 
                                              const PetscScalar &, 
                                              const PetscScalar &) ; 

  /** initial condition
   *  return body  temp on subdomain 0 1 
   *  return probe temp on subdomain 2 3 
   */
  virtual PetscScalar getInitialVoltage(unsigned int domainId, 
                                        unsigned int , 
                                        const Point&  , const Parameters& )
     { 
       return m_InitialVoltage[domainId];
     }
protected:
  std::vector< PetscScalar > m_InitialVoltage;
  spatialParameter m_ElectricConductivity_0; ///<  \f$ \sigma_0    \f$ 
  PetscScalar      m_ElectricConductivity_1; ///<  \f$ \sigma_1    \f$ 
private:
  /** The type of the parent  */
  typedef PennesBioheatModel Parent;
};
/**
 * \f$ \Delta P_1  \f$ model
 *
 */ 
class PennesDeltaP1 : public PennesStandardDiffusionApproximation
{
public:

 //constructor
 PennesDeltaP1(GetPot &,EquationSystems &);

 //print data 
 virtual void printSelf(std::ostream& os=std::cout) ; 

 //generate the attenuation data 
 PetscScalar ComputeRadialLoss( const Point &);
 PetscScalar ComputeTotalAttenuation( const Elem* , const Elem* );

 /** primary fluence */
 virtual PetscScalar getExternalIrradiance(    unsigned int,unsigned int, 
                                               const Point&,const Parameters&);
 virtual PetscScalar getExternalFlux_X(        unsigned int,unsigned int, 
                                               const Point&,const Parameters&);
 virtual PetscScalar getExternalFlux_Y(        unsigned int,unsigned int, 
                                               const Point&,const Parameters&);
 virtual PetscScalar getExternalFlux_Z(        unsigned int,unsigned int, 
                                               const Point&,const Parameters&);
 virtual PetscScalar getInterstitialIrradiance(unsigned int,unsigned int,
                                               const Point&,const Parameters&);
 virtual PetscScalar getInterstitialFlux_X(    unsigned int,unsigned int, 
                                               const Point&,const Parameters&);
 virtual PetscScalar getInterstitialFlux_Y(    unsigned int,unsigned int, 
                                               const Point&,const Parameters&);
 virtual PetscScalar getInterstitialFlux_Z(    unsigned int,unsigned int, 
                                               const Point&,const Parameters&);
 /** \f$ \mu_a  \phi \f$ is the delivered fluence */
 PetscScalar ScatteredSource( const unsigned int &domainId, 
                        const Real &u_theta, const Real &damage,
                        const Real &ScatteredFluence ) 
 {
    return this->absorption(domainId,u_theta,damage) 
           * ScatteredFluence ;
 }
 virtual PetscScalar PennesSource( const unsigned int &domainId, 
                        const Real &u_theta, const Real &damage,
                        const Real &PrimaryFluence,
                        const Point &, const int istep) 
 {
    return this->absorption(domainId,u_theta,damage) *Power[istep] * PrimaryFluence; 
 }

 Real FluenceReactionTerm( const unsigned int &domainId,const Real &u_theta, const Real &damage, const Real &fluence) 
 {
   Real mu_tr=this->absorption(domainId,u_theta,damage)
             +this->scattering(domainId,u_theta,damage)*(1.0e0-anfact);
   Real mu_eff=std::sqrt(3.0e0 
                        *this->absorption(domainId,u_theta,damage)*mu_tr);
   return mu_eff * mu_eff * fluence; 
 }
 Gradient FluxSource( int domainId,const Real &u_theta, const Real &damage,const Gradient & flux, const int timeId)
 {
   Real mu_s_star =  this->scattering(domainId,u_theta,damage)*(1.0e0-m_ForwardScatterFraction);
   return 3.0 * m_ScatteringAsymmetry * mu_s_star * Power[timeId] *flux ; 
 }

 Real IrradianceSource( int domainId,const Real &u_theta, const Real &damage,const Real& irradianceSource , const int timeId)
 {
   Real mu_s_star =  this->scattering(domainId,u_theta,damage)*(1.0e0-m_ForwardScatterFraction);
   Real mu_tr=this->absorption(domainId,u_theta,damage)
             +this->scattering(domainId,u_theta,damage)*(1.0e0-anfact);
   return 3.0 * mu_tr *  mu_s_star * irradianceSource ; 
 }

 virtual Real ExternalIrradiance( int domainId,const Real &u_theta, const Real &damage,const Point& qpoint, const int timeId)
 { // TODO extend to volumetric source
    // external beam for NOW
    PetscScalar locposRefTip[3]=
                      {qpoint(0) - m_laserTip[0],
                       qpoint(1) - m_laserTip[1],
                       qpoint(2) - m_laserTip[2]};
    PetscScalar axialComponent =
                     locposRefTip[0]* m_unitVec[0] +
                     locposRefTip[1]* m_unitVec[1] +
                     locposRefTip[2]* m_unitVec[2] ;
    PetscScalar radialVec[3]=
        {locposRefTip[0] - axialComponent * m_unitVec[0],
         locposRefTip[1] - axialComponent * m_unitVec[1],
         locposRefTip[2] - axialComponent * m_unitVec[2]};
    // check if satisfy the radius criteria
    PetscScalar locRadius;
    locRadius = std::sqrt( std::pow(radialVec[0],2)   +
                           std::pow(radialVec[1],2)   +
                           std::pow(radialVec[2],2)   
                         );
    Real mu_t_star=this->absorption(domainId,u_theta,damage)
                  +this->scattering(domainId,u_theta,damage)
                  *(1.0e0-m_ForwardScatterFraction);
    Real PowerAttenuation = 0.0;
    if( axialComponent < m_AgarLength ) 
       PowerAttenuation = std::exp( - mu_t_star * axialComponent ) ;
    else
       PowerAttenuation = std::exp( - m_Agar_mu_t_star * m_AgarLength ) *  
                     std::exp( - mu_t_star * ( axialComponent - m_AgarLength )  )   ;
    return 2.0*Power[timeId]/libMesh::pi/m_GuassBeamRadius/m_GuassBeamRadius 
                         * (1.0 - m_SpecularReflectance) 
                         * PowerAttenuation 
                         * std::exp( -2.0 * locRadius * locRadius /
                                       m_GuassBeamRadius / m_GuassBeamRadius ) ;
 }
 /**
  * Cauchy boundary data for each variable, assumes INWARD Pointing normal for
  * minus sign FUN!!!
  * \f[
  *   - \nabla \phi \cdot \hat{n}  = \frac{1}{Ah} (\phi +  A h 3 g^*  \mu_s^* E(r,z))
  * \f]
  */
 virtual PetscScalar residualFluenceBC(const unsigned int i_var, const Real &fluence,
                                       int domainId, const int timeId,
                             const Point &qpoint, const Parameters& parameters )
 {
   //FIXME design flaw in BC, assume linear for now
   Real u_theta = m_bodyTemp, damage = 0.0;
   Real mu_tr=this->absorption(domainId,u_theta,damage)
             +this->scattering(domainId,u_theta,damage)*(1.0e0-anfact);
   Real h = 2.0/3.0 / mu_tr;
   //
   PetscScalar boundaryValue = 0.0;
   const unsigned int z_var= parameters.get<unsigned int>("z_var") ; 
   if(i_var == z_var ) boundaryValue  = 1.0/m_ReflectionFactor/h * 
                       (fluence + m_ReflectionFactor * h * 3.0 * m_ScatteringAsymmetry 
                                                     * scattering(domainId,u_theta,damage)*(1.0e0-m_ForwardScatterFraction) 
                                                     * this->ExternalIrradiance(domainId,u_theta,damage,qpoint,timeId) ) ;
   return boundaryValue ;
 }
 virtual PetscScalar jacobianFluenceBC(const unsigned int i_var, int domainId, 
                         const Point &, const Parameters& parameters)
 {
   //FIXME design flaw in BC, assume linear for now
   Real u_theta = m_bodyTemp, damage = 0.0;
   Real mu_tr=this->absorption(domainId,u_theta,damage)
             +this->scattering(domainId,u_theta,damage)*(1.0e0-anfact);
   Real h = 2.0/3.0 / mu_tr;
   //
   PetscScalar boundaryValue = 0.0;
   const unsigned int z_var= parameters.get<unsigned int>("z_var") ; 
   if(i_var == z_var ) boundaryValue  = 1.0/m_ReflectionFactor/h ;
   return boundaryValue ;
 }
 virtual PetscScalar residualCauchyBC(const unsigned int i_var, const Real &varValue,
                                       int domainId, const int ,
            const Point &, const Parameters& parameters)
 {
   //FIXME design flaw in BC, assume linear for now
   Real u_theta = m_bodyTemp, damage = 0.0;
   Real mu_tr=this->absorption(domainId,u_theta,damage)
             +this->scattering(domainId,u_theta,damage)*(1.0e0-anfact);
   Real h = 2.0/3.0 / mu_tr;
   //
   PetscScalar boundaryValue = m_newton_coeff[i_var]*(varValue-m_u_infty[i_var]) ; 
   const unsigned int z_var= parameters.get<unsigned int>("z_var") ; 
   if(i_var == z_var ) boundaryValue  = 1.0/m_ReflectionFactor/h * varValue;
   return boundaryValue ;
 }
 virtual PetscScalar jacobianCauchyBC(const unsigned int i_var, int domainId, 
                     const Point &, const Parameters& parameters)
 {
   //FIXME design flaw in BC, assume linear for now
   Real u_theta = m_bodyTemp, damage = 0.0;
   Real mu_tr=this->absorption(domainId,u_theta,damage)
             +this->scattering(domainId,u_theta,damage)*(1.0e0-anfact);
   Real h = 2.0/3.0 / mu_tr;
   //
   PetscScalar boundaryValue = m_newton_coeff[i_var]; 
   const unsigned int z_var= parameters.get<unsigned int>("z_var") ; 
   if(i_var == z_var ) boundaryValue  = 1.0/m_ReflectionFactor/h ;
   return boundaryValue ;
 }
protected:
 /*
  * Delta P1 laser source
  * Simulates the effect of the deltaP1 laser source
  */
 virtual Real interstitialFluence(int ,const Real &, const Real &, const Point& , const int )
 {
   //Point tmp(x_0[0],y_0[0],z_0[0]);
   //Point dist = qpoint - tmp;
   //Real mu_ss=scattering(domainId,damage)*(1.0e0-anfact)*(1.0e0-m_f);
   //Real mu_tr=absorption(domainId,damage)+scattering(domainId,damage)*(1.0e0-anfact);
   //Real mu_trs=absorption(domainId,damage)+scattering(domainId,damage)*(1.0e0-anfact)*(1.0e0-f);
   //Real mu_eff=std::sqrt(3.0e0*absorption(domainId,damage)*mu_tr);
   //Real alpha=(3.0e0*mu_ss*(mu_trs-gs*absorption(domainId,damage)))/(mu_eff*mu_eff-mu_trs*mu_trs);
   //Real beta=(-alpha*(1+2.54e0*(6.66e-1)*mu_tr*mu_trs)-3*2.54e0*(6.66e-1)*mu_tr*gs*mu_ss)/(1+2.54e0*(6.66e-1)*mu_tr*mu_eff);
   //Real fluence=std::exp(-mu_trs*dist.size())+
   //     (alpha*std::exp(-mu_trs*dist.size())+beta*std::exp(-mu_trs*dist.size()));
   //return Power[timeId]*absorption(domainId,damage)*1.13682e5*fluence; // The 1.06157e5 adjusts for the probe area.
   libmesh_not_implemented(); return 0.0;
 }
 
  Real m_ScatteringAsymmetry,    ///< scattering asymmetry, \f$ g^* \f$ 
       m_ForwardScatterFraction, ///< forward scatter fraction, \f$ f \f$ 
       m_ReflectionFactor,       ///< \f$ A = -0.13755 n^3 +4.33900 n^2 -4.90366 n +1.68960; \f$ 
       m_SpecularReflectance;    ///<                           \f$ R_s \f$ 

private:
  Real m_GuassBeamRadius ;       ///<                           \f$ r_0 \f$ 
  // need to explicitly offset the attenuation of the agar region
  Real m_Agar_mu_t_star , ///< total attenuation in agar region
       m_AgarLength     ; ///< length of agar region

  /** The type of the parent  */
  typedef PennesStandardDiffusionApproximation Parent;
};
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
Real PennesBioheatModel::PennesReactionTerm(
                             const unsigned int &domainId, 
                             const Real &u_theta, 
                             const Real &damage) 
{
  return this->bloodspecificheat * 
         this->perfusion(domainId,u_theta,damage) * (u_theta-this->u_artery); 
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
Real PennesBioheatModel::dPennesReactionTermdu(
                             const unsigned int &domainId, 
                             const Real &u_theta, 
                             const Real  &damage, 
                             const Real &DdamageDu) 
{
  return this->bloodspecificheat * ( this->perfusion(domainId,u_theta,damage) 
                                   +
         this->dperfusiondu(domainId,u_theta,damage,DdamageDu)*(u_theta-this->u_artery) 
                                   ) ; 
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline
Real PennesBioheatModel::perfusion( 
                           const unsigned int &domainId, const Real &, const Real &damage) 
{ 
  Real thermal_dose = damage/(std::log(2.0)+damage);
  Real value = w_0[domainId]+thermal_dose*(w_1 - w_0[domainId]); 
  return value ; // faster to multiply
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline
Real PennesBioheatModel::dperfusiondu(
                           const unsigned int &domainId, const Real &, const Real &damage, const Real &DdamageDu) 
{ 
  Real dthermal_dosedu = DdamageDu * std::log(2.0)/(std::log(2.0)+damage)/(std::log(2.0)+damage);
  Real value = dthermal_dosedu * (w_1 - w_0[domainId]); 
  return value ; // faster to multiply
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
Real PennesBioheatModel::ThermalConductivity(
                           const unsigned int &domainId, const Real &u_theta,  const Real &) 
{
  return k_0[domainId]+u_theta*k_1; 
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
Real PennesBioheatModel::dThermalConductivitydu(
                           const unsigned int &, const Real &, const Real &, const Real &) 
{ 
  return k_1; 
}

// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
Real PennesStandardDiffusionApproximation::absorption( 
                           const unsigned int &domainId, const Real &, const Real &damage) 
{ 
  Real thermal_dose = damage/(std::log(2.0)+damage);
  Real value = mu_a_0[domainId]+thermal_dose*(mu_a_1 - mu_a_0[domainId]); 
  return value ; // faster to multiply
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
Real PennesStandardDiffusionApproximation::dabsorptiondu(
                           const unsigned int &domainId, const Real &, const Real &damage, const Real &DdamageDu) 
{ 
  Real dthermal_dosedu = DdamageDu * std::log(2.0)/(std::log(2.0)+damage)/(std::log(2.0)+damage);
  Real value = dthermal_dosedu * (mu_a_1 - mu_a_0[domainId]); 
  return value ; // faster to multiply
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
Real PennesStandardDiffusionApproximation::scattering( 
                           const unsigned int &domainId, const Real &, const Real &damage) 
{ 
  Real thermal_dose = damage/(std::log(2.0)+damage);
  Real value = mu_s_0[domainId]+thermal_dose *(mu_s_1 - mu_s_0[domainId]); 
  return value ; // faster to multiply
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
Real PennesStandardDiffusionApproximation::dscatteringdu(
                           const unsigned int &domainId, const Real &, const Real &damage, const Real &DdamageDu) 
{ 
  Real dthermal_dosedu = DdamageDu * std::log(2.0)/(std::log(2.0)+damage)/(std::log(2.0)+damage);
  Real value = dthermal_dosedu * (mu_s_1 - mu_s_0[domainId]); 
  return value ; // faster to multiply
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
Real PennesStandardDiffusionApproximation::interstitialFluence(
    int domainId,const Real &u_theta, const Real &damage,const Point& qpoint, const int timeId)
                           
{
  Real mu_tr=this->absorption(domainId,u_theta,damage)
            +this->scattering(domainId,u_theta,damage)
            *(1.0e0-anfact);
  Real mu_eff=std::sqrt(3.0e0 *
             this->absorption(domainId,u_theta,damage)*mu_tr);
  Real source=0.0;
  if(!this->dirichletNodalData(qpoint))
    for (unsigned int Ii = 0; Ii < m_volumeFraction.size(); Ii++ )
     {
      Point tmp(x_0[Ii],y_0[Ii],z_0[Ii]);
      Point dist = qpoint - tmp;
      source += 0.75 * Power[timeId] * m_volumeFraction[Ii] * 
              mu_tr * std::exp(-mu_eff*dist.size())/(libMesh::pi*dist.size());
     }
  return source;
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
PetscScalar PennesVoltage::ElectricConductivity( 
            const unsigned int &domainId, const PetscScalar & temperature, const PetscScalar &  ) 
{
   Real value = 
        ( m_ElectricConductivity_0[domainId] + m_ElectricConductivity_1 * temperature ) ; 
   return value ; // faster to multiply
}
// inline functions should be in the header see @ref PDEModelInlineStrategy
inline 
PetscScalar PennesVoltage::dElectricConductivitydu( 
            const unsigned int &, const PetscScalar &, const PetscScalar &, const PetscScalar & ) 
{
   Real value = m_ElectricConductivity_1 ;
   return value ; // faster to multiply
}
#endif
