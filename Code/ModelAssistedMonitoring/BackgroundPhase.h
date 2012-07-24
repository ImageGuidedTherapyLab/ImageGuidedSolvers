/**
 * Laplace Solution for Background phase from MagnetoStatics.
 * The resulting boundary value problem is of the form of an elliptic
 * PDE with pure Neumann boundary conditions (KSPSetNullSpace in PETSc).
 * \f[
 * \boxed{
 * \begin{split}
 * \text{Find } & \psi_{obj}:
 * \\
 * & \Delta \psi_{obj} = -H_0 \frac{\partial \chi}{\partial z}
 *                                      \qquad \text{in } \Omega
 * \\
 * & \nabla  \psi_{obj}
 *  =
 * \begin{bmatrix}
 *  0 \\
 *  0 \\
 * \left( 1 - \delta - \frac{2}{3} \chi \right)
 *                 \left( 1+ \chi \right) \mu
 *  H_0
 *   - \frac{\varphi}{\gamma \; TE}
 *   \\
 * \end{bmatrix}
 *                                      \qquad \text{on } \partial \Omega
 * \end{split}
 * }
 * \f]
 * 
 * @htmlonly
 *   see pdf for in depth derivation of governing equations 
 * @endhtmlonly
 * 
 * @latexonly
As water heats up the weak hydrogen bonds that holds the water together
stretches. As the bonds are stretched the protons are pulled in closer
to the Oxygen molecules. The proximity of the proton to the Oxygen
causes a shielding effect that \textit{linearly} 
changes the resonance frequency with distance from the Oxygen.
A typical model is that the phase of the MR signal for a gradient echo imaging sequence
is proportional to the local magnetic field seen at the hydrogen nucleus along the direction 
of the static $B_0$ field.
\[
   \varphi =  \gamma \; TE \; B_{nuc,z} \qquad \text{(scalar valued)}
\Rightarrow
   \delta \varphi =  \gamma \; TE \; \delta B_{nuc,z} 
\]
MR temperature imaging measures temperature change through a model
of the temperature dependent magnetic field
change, $ \delta B_{nuc,z} (\delta T) $.
The  model of the field dependence on temperature may be obtain from
a quasi magneto-statics.

The measured phase angle, $\varphi$, of a gradient echo imaging sequence
can be taken as proportional to the axial component of the
magnetic field, $B_{nuc,z}$, seen at the signal
generating nucleus~\cite{salomir2003fast}.
\[
   \varphi =  \gamma \; TE \; B_{nuc,z} \qquad \text{(scalar valued)}
\]
Where $\gamma$ is the gyromagnetic ratio and $TE$ is the echo time of the
sequence. 

\paragraph{Maxwell Equations}
Maxwell equation'
(Faraday's law~\eqn{FaradayLaw}, Ampere's law~\eqn{AmpereLaw}, Gauss's
law~\eqn{GaussLawElectric} \&\eqn{GaussLawMagnet})~\cite{demkowicz2006cha},
and the continuity of free charge~\eqn{ChargeContinuity} may be used to
obtain a model for phase.
\begin{equation} \label{FaradayLaw}
 \nabla \times \vec{E} = -\frac{\partial}{\partial t} \vec{B}
\end{equation}
\begin{equation} \label{AmpereLaw}
 \nabla \times \vec{B} = 
        \mu \vec{J} + \mu \epsilon \frac{\partial}{\partial t} \vec{E}
\end{equation}
\begin{equation} \label{GaussLawElectric}
    \nabla \cdot \left( \epsilon \vec{E} \right)  =  \rho_{charge}
\end{equation}
\begin{equation} \label{GaussLawMagnet}
    \nabla \cdot \vec{B} = 0 
\end{equation}
\begin{equation} \label{ChargeContinuity}
    \nabla \cdot \vec{J} + \frac{\partial}{\partial t} \rho_{charge} = 0 
\end{equation}
Here, the medium is assumed linear and isotropic with homogeneous
permittivity, $\epsilon$, and permeability, $\mu = \mu_0 (1+\chi)$.
The charge density is denoted $\rho_{charge}$.
The electric and magnetic fields are denoted
$\vec{E}$ and $\vec{B}$, respectively. 
The current density is denoted $J$.
\paragraph{Steady State} In steady state, the system decouples
\[\begin{split}
 \nabla \times \vec{E} =   0 
\\
 \nabla \times \vec{B} = 
        \mu \vec{J} 
\\
    \nabla \cdot \left( \epsilon \vec{E} \right)  =  \rho_{charge}
\\
    \nabla \cdot \vec{B} = 0 
\\
    \nabla \cdot \vec{J}  = 0 
\end{split} \]
\paragraph{Magnetic Scalar Potential}
Assuming a magnetic scalar potential 
\[
    \vec{B} =  -  \nabla  \psi
\]
The current density, $J$, must necessarily vanish because the
identity that the curl of the gradient equals zero.
\[
  \vec{J}  = \frac{1}{\mu} \nabla \times \vec{B}   
           = \frac{1}{\mu} \nabla \times \left( - \nabla  \psi \right) = 0 
\]
Gauss law for magnetism \eqn{GaussLawMagnet}
 requires that the magnetic field is solenoidal
and leads to a laplace equation for the scalar potential, $\psi$
\[
    \nabla \cdot \vec{B} = 0 \qquad \Rightarrow \qquad - \Delta  \psi = 0 
\]


\paragraph{Relating Phase to Magnetic Field}
Under typical MR operating conditions the magnetic field in the transverse
plane may be considered second order effects to the magnetic field
along the axial dimension 
\[
\vec{B}_{nuc} \approx 
\begin{bmatrix}
 0 \\
 0 \\
   B_{nuc,z} \\
\end{bmatrix}
   \qquad \text{(this will be used later for the BC)}
\]
The magnetic field at the nucleus $\vec{B}_{nuc}$
is related to the macroscopic magnetic field, $\vec{B}_{mac}$,
and flux, $\vec{H}$, 
through the susceptibility, $\chi$, and chemical shift $\sigma$.
\[
\vec{B}_{nuc} = \left( 1 - \sigma - \frac{2}{3} \chi \right) 
                \vec{B}_{mac}
\qquad
\vec{B}_{mac} = \left( 1+ \chi \right) \mu \vec{H}
\qquad
\vec{B}_{nuc} = \left( 1 - \sigma - \frac{2}{3} \chi \right) 
                \left( 1+ \chi \right) \mu \vec{H}
\qquad
\vec{B} = \mu \vec{H} 
\]
As no current is flowing ($\vec{J}=0$) through an object 
in a typical MR image (ignore eddy currents),
the magnetic flux may be derived from a potential
$\vec{H} =  -  \nabla  \psi$. Further, for an MR acquisition, 
the potential may be considered the superpostion  of the main magnetic
field, $\psi_0$, the inhomogenieties do to improperty shimming, $\psi_{in}$,
and a demagnetizing field due to the object within the field of view
$\psi_{obj}$.
\[
\psi = 
 \psi_0
+
 \psi_{in}
+
 \psi_{obj}
\qquad
\vec{H} =  -  \nabla  \psi_0
           -  \nabla  \psi_{in}
           -  \nabla  \psi_{obj}
\qquad
\psi_0 = -H_0 z
\Rightarrow
\vec{H}_0 =  -  \nabla  \psi_0 = 
\begin{bmatrix}
 0 \\
 0 \\
 H_0 \\
\end{bmatrix}
\]
Gauss law appled to the macroscopic magnetic field $\vec{B}_{mac}$ yeilds
\[
\nabla \cdot \vec{B}_{mac} =  0 
= 
\nabla \cdot \mu (1+\chi) \vec{H} 
\qquad
\Rightarrow
\qquad
\nabla \chi \cdot 
\nabla \psi 
+
(1+ \chi )
\Delta \psi = 0 
\]
Under standard MR conventions~\cite{salomir2003fast} 
Neglecting second order terms
($h_{obj}\chi \approx 0$) yeilds  
a PDE for the magnetic potential of the object field
that is analogous to a electric field driven by a susceptibility source.
\begin{equation} \label{DemagnetizingFieldEqn}
\Delta \psi_{obj} = -H_0 \frac{\partial \chi}{\partial z}
\end{equation}

\paragraph{Temperature Measurement}

Under standard MR conventions~\cite{salomir2003fast}, our
model of the B-field at the nucleus of the signal generating Hydrogen atom
is given as the superposition of the static field and the demagnitizing
field
\[
   B_{nuc} \approx B_0 \left( 1 -  \sigma(T) + \frac{1}{3} \chi(T) \right) 
                   + \mu \nabla \psi_{obj}
\]
Temperature dependent changes, $\delta T$ in the B-field are induced by the temperature
dependent change in the shielding constant, 
$\alpha_\sigma = \frac{ \delta \sigma(T)}{\delta T } \left[ \frac{ppm}{^oC}\right]$, and susceptibility 
$\alpha_\chi = \frac{ \delta \chi(T)}{\delta T }\left[ \frac{ppm}{^oC}\right]$.
Note that units of ppm = 1e-6 are \textbf{dimensionless}.
Also notice that a potentially additional spatially 
dependent change in the B-field and measure phase may be induce through
the demagnetizing field through changes in the susceptibility
\eqn{DemagnetizingFieldEqn}.
\[
 \delta  B_{nuc} \approx B_0 \left(  - \delta \sigma(T) + \frac{1}{3}\delta \chi(T) \right) 
                   + \mu \delta \left( \nabla \psi_{obj} \right)
\]
In the case that we ignore susceptibility induce changes in the B-field, we
arrive at the standard equation for temperature induced phase change.
\[
   \delta \varphi =   \gamma \; TE \; \delta B_{nuc,z} 
                  = - \gamma \; TE \; B_0  \alpha_\sigma \delta T 
\qquad
\Rightarrow
\qquad
   \delta T  =  \frac{-\delta \varphi}{ \gamma \; TE \; B_0  \alpha_\sigma}
\]
In general, the measured phase change is perturbed by susceptibility 
effects.
\[
   \delta \varphi 
                  = - \gamma \; TE \; B_0  \alpha_\sigma \delta T 
                   + \gamma \; TE \; B_0 \frac{1}{3} \alpha_\chi \delta T 
                   + \gamma \; TE \; \mu \delta \left( \nabla \psi_{obj} \right)
\]
One approach to bound the influence of the susceptibility effects is to take
the Fourier Transform of the governing equation of the demagnetizing 
field~\eqn{DemagnetizingFieldEqn}. This allows an analytic representation 
of the resulting local B-field~\cite{salomir2003fast}.
\[
\begin{split}
& B_{nuc,z}= B_0 
\left\{ 
 1 - \sigma 
   + FT^{-1}\left[ \left(\frac{1}{3}- \frac{k_z^2}{k^2}\right) FT(\chi ) \right]
\right\}
\\
& \Rightarrow
\delta B_{nuc,z}= B_0 
\left\{ 
   - \alpha_\sigma \delta T
   + \alpha_\chi FT^{-1}\left[ \left(\frac{1}{3}- \frac{k_z^2}{k^2}\right) FT(\delta T) \right]
\right\}
\end{split}
\]
Alternatively, PDE theory~\cite{evans-partial} may be used to bound the
elliptic solution to the demagnetizing field~\eqn{DemagnetizingFieldEqn}. 
\[
   \| \psi_{obj}\|_{H^2} \leq C \left( 
                        \| H_0 \frac{\partial \chi}{\partial z} \|_{L^2} + \| \psi_{obj} \|_{L^2}
                                \right)
\]
However the constant, $C$, would have to be determined for an exact estimate.

The Fourier approach agrees with the experimentally observed spatially dependent susceptibility
artifacts~\cite{peters2000heat} that depend on the temperature gradients and
derivatives.
\[
\delta B_{nuc,z}= B_0 
\left\{ 
   - \alpha_\sigma \delta T
   + \alpha_\chi FT^{-1}\left[ \left(\frac{FT(\delta T)}{3}-
\frac{FT\left( \frac{\partial^2 \delta T}{\partial z^2} \right )}{k^2}\right)  \right]
\right\}
= B_0 
\left\{ 
   - \alpha_\sigma \delta T
   + \alpha_\chi  f ( \delta T, \nabla \delta T, \Delta \delta T)
\right\}
\]
Using the susceptibility coefficient from Peter~\cite{peters2000heat},
$\alpha_\chi \approx  0.0027  \left[\frac{10^{-6}}{^oC}\right]$,
TE=20ms, and  $\bar{\gamma}B_0 =\frac{\gamma}{2\pi}\gamma B_0$ =  63.87 MHz,
a temperature change of $\delta T = 1^oC$ would cause a phase angle change
of (considering temperature change only, no gradient dependence)
\[
\begin{split}
   \delta \varphi 
                  &= - \frac{1}{3} \gamma \; TE \; B_0  \alpha_\chi \delta T 
                  = - \frac{1}{3} 2 \pi \bar{\gamma} B_0 \; TE   \alpha_\chi \delta T 
\\
                  &= - \frac{1}{3} 2 \pi \left(63.87 MHz\right) \; \left(20 ms \right)\; 
                     \left(0.0027 \left[\frac{10^{-6}}{^oC}\right]\right)  \;
                     \left( 1 ^oC\right)
% 2 * pi * 63.87 * 0.0027 * 0.020 * 1 /3.
                  = 0.0072 \; rad = 0.414^o
\end{split}
\]
As a comparison for $\alpha_\sigma \approx  0.0097  \left[\frac{10^{-6}}{^oC}\right]$,
\[
\begin{split}
   \delta \varphi 
              &   = - \gamma \; TE \; B_0  \alpha_\sigma \delta T 
                  = - 2 \pi \bar{\gamma} B_0 \; TE   \alpha_\sigma \delta T 
\\
               &  = - 2 \pi \left(63.87 MHz\right) \; \left(20 ms \right)\; 
                     \left(0.0097 \left[\frac{10^{-6}}{^oC}\right]\right)  \;
                     \left( 1 ^oC\right)
% 2 * pi * 63.87 * 0.0097 * 0.020 * 1
                  = 0.0779 \; rad = 4.46^o
\end{split}
\]
The error in the termperature measurement when not accounting for susceptibility effects
may be seen from considering a phase change measurement $\delta \varphi$, with
a  true value of the temperature change, $\delta T^*$ and a temperature change when not
accounting for the susceptibility, $\delta T$
\[
\begin{split}
 - \gamma \; & TE \; B_0  \alpha_\sigma \delta T 
   = 
     \delta \varphi
   = - \gamma \; TE \; B_0  \alpha_\sigma \delta T^*
     + \gamma \; TE \; B_0 \frac{1}{3} \alpha_\chi \delta T^*
     + \gamma \; TE \; f  \left(\nabla \delta T^* , \Delta \delta T^*  \right)
\\
 & \Rightarrow
  \delta T^* - \delta T
   = 
     \frac{\alpha_\chi}{   \alpha_\sigma}
       \left(
             \frac{\delta T^*}{ 3 }   +  f\left(\nabla \delta T^* , \Delta \delta T^*  \right)
       \right)
\end{split}
\]

\paragraph{Background Correction using Phase on the Boundary}

Ignoring shimming inhomogenieties,  $\psi_{in} \approx 0$, and
second order terms, $\chi^2 \approx \sigma \chi \approx \sigma \nabla
\psi_{obj} \approx \chi \nabla \psi_{obj} \approx 0 $, boundary
conditions may be obtain from phase data.
\[
\vec{B}_{nuc} = 
\begin{bmatrix}
 0 \\
 0 \\
 \frac{\varphi}{\gamma \; TE}  \\
\end{bmatrix}
=
\left( 1 - \sigma - \frac{2}{3} \chi \right) 
                \left( 1+ \chi \right) \mu 
\left(  
\begin{bmatrix}
 0 \\
 0 \\
 H_0 \\
\end{bmatrix}
        -  \nabla  \psi_{obj} \right)
\approx 
\mu 
\left( 1 - \sigma + \frac{1}{3} \chi \right) 
\left(  
\begin{bmatrix}
 0 \\
 0 \\
 H_0 \\
\end{bmatrix}
        - \nabla  \psi_{obj} \right)
\]
The resulting boundary value problem is of the form of an elliptic
PDE with pure Neumann boundary conditions (KSPSetNullSpace in PETSc).
\[ 
\boxed{
\begin{split}
\text{Find } & \psi_{obj}:
\\
& \Delta \psi_{obj} = -H_0 \frac{\partial \chi}{\partial z} 
                                     \qquad \text{in } \Omega
\\
& \nabla  \psi_{obj}
 = 
\begin{bmatrix}
 0 \\
 0 \\
 H_0
  - \frac{\varphi}{\gamma \; TE
      \left( 1 - \sigma - \frac{2}{3} \chi \right) 
                \left( 1+ \chi \right) \mu }  
  \\
\end{bmatrix}
                                     \qquad \text{on } \partial \Omega
\end{split} 
}
\]
Recall, once the magnetic potential is known the corresponding phase
may be obtained as
\[
\begin{split}
 \varphi  & =  \gamma \; TE \; B_{nuc,z} 
=
\gamma \; TE
\left( 1 - \sigma - \frac{2}{3} \chi \right) 
                \left( 1+ \chi \right) \mu 
\left(  
\begin{bmatrix}
 0 \\
 0 \\
 H_0 \\
\end{bmatrix}
        -  \nabla  \psi_{obj} \right) \cdot \vec{k}
\\
& \approx
\gamma \; TE
\mu 
\left(  
\left( 1 - \sigma + \frac{1}{3} \chi \right)  
\begin{bmatrix}
 0 \\
 0 \\
 H_0 \\
\end{bmatrix}
        -  \nabla  \psi_{obj} \right) \cdot \vec{k}
\end{split}
\]

\paragraph{Solving for Phase Difference}
As in standard complex phase difference based imaging,
The resulting boundary value problem may be more appropriate
in terms of the phase change $\delta \varphi$ between two time instances.
\[ 
\boxed{
\begin{split}
\text{Find } & \delta \psi_{obj}:
\\
& \Delta \delta \psi_{obj} = -H_0 
       \delta \left(\frac{\partial \chi}{\partial z} \right)
                                     \qquad \text{in } \Omega
\\
& \nabla \delta  \psi_{obj}
 = 
\begin{bmatrix}
 0 \\
 0 \\
  - \frac{\delta \varphi}{\gamma \; TE
\left( 1 - \sigma - \frac{2}{3} \chi \right) 
                \left( 1+ \chi \right) \mu 
}  
  \\
\end{bmatrix}
                                     \qquad \text{on } \partial \Omega
\\
& \delta \left( \cdot \right) \equiv \text{ change in quantity }
                     \left(  \cdot \right)
       \text{ between time instance one and time instance two } 
\end{split} 
}
\]
The corresponding phase difference and induced temperature change
is obtained as
\[
 \delta  \varphi 
=
\gamma \; TE
\left( 1 - \sigma - \frac{2}{3} \chi \right) 
                \left( 1+ \chi \right) \mu 
\left(  -  \nabla  \delta \psi_{obj} \right) \cdot \vec{k}
\qquad
\Delta u = \frac{ \delta \varphi } { 2 \pi \alpha \cdot \gamma B_0 \cdot \text{TE}}
         = \frac{ 
\left( 1 - \sigma - \frac{2}{3} \chi \right) 
                \left( 1+ \chi \right) \mu } { 2 \pi \alpha \cdot B_0 }
\left(  -  \nabla  \delta \psi_{obj} \right) \cdot \vec{k}
\]

 * @endlatexonly
 */
#ifndef __BackgroundPhase_h
#define __BackgroundPhase_h
#include "petsc_fem_system.h"
class BackgroundPhaseModel : public PDEModelBaseClass
{
public:
  // Constructor
  BackgroundPhaseModel(GetPot &, EquationSystems &) ; 

  /**  overwrite neuman BC */
  virtual PetscScalar residualNeumannBC(const unsigned int i_var,
                                        const Real &phase,
                                        int , const int ,
                                        const Point & , const Parameters&)
    {
      libmesh_not_implemented(); // TODO  verify
      m_NeumannFlux[i_var] = - phase / (  // MINUS SIGN
       m_EchoTime*1.0e-3 * m_GyroMagRatio*1.0e6 * //MHz cancel ppm ?
       ( 1.0 - m_ChemicalShift - 2.0/3.0*m_Susceptibility[0] ) * 
           ( 1.0 + m_Susceptibility[0] ) * m_ElectricPermeability  
                                     );
      return  m_NeumannFlux[i_var]; 
    }
  PetscScalar load(){return m_StaticFieldB0 * m_dSusceptibilitydz; }
  /** Print the Constitutive Data */
  virtual void printSelf(std::ostream& );
private:
  
  spatialParameter m_Susceptibility;
  PetscScalar m_EchoTime,          ///<  \f$TE \f$  [milli seconds]
              m_StaticFieldB0,     ///<  [Tesla]
              m_dSusceptibilitydz, ///<  [ppm]
              m_ChemicalShift ,    ///<  [ppm]
              m_ElectricPermeability , ///< [T*m/A]
              m_GyroMagRatio;     ///<  \f$\gamma \f$ [MHz/T]
  
};

/**@ingroup ModelAssistedMonitoring 
 * @htmlonly
 *   see pdf for in depth derivation of the variational form
 * @endhtmlonly
 * 
 * @latexonly
 * Our backround phase model is of the form
 * \[ \begin{split}
 *  -\nabla \cdot (a \nabla u ) &  = f \text{ in }  \Omega 
 *  \\
 *  a \nabla u \cdot n  &   = g \text{ on } \partial \Omega
 *  \end{split} 
 * \]
 * While classical solutions using fundamental solutions, greens functions,
 * Laplace/Fourier transform, etc exists for a variety of PDE's of interest,
 * classical solutions that are differentiable up to the order of the strong
 * form of an PDE are \textit{in general} difficult to
 * find~\cite{evans1998partial}. Not only does (1) the solution need to be
 * constructed using some intuition but (2) it must also be shown that the
 * constructed solution satisfies the PDE and possesses the desired regularity
 * and smoothness, ie that derivatives of the constructed solution are indeed
 * well defined (typically a couple page proof).  As mathematicians considered
 * solutions to more general PDE's, the typical strategy became to separate the
 * existance proof from the smoothness proof. The key idea is that, using a
 * functional analysis approach, it should be easier to find a unique solution
 * within an appropriate class of generalized weak solutions in which little is
 * assumed as to the regularity of the solution. Given a well-posed solution
 * the details of the regularity may be revisited.  An excellent example as to
 * the power of this method is given by PDE's in which a classical solution
 * does not exist but a solution to the PDE exists in a broader class of
 * generalized functions.
 * A key difference between the functional analysis approach is that it does
 * not provide explicit formulas (like the classical approach) but rather
 * establishes that (1) indeed a solution does exist and (2) under given
 * regularity of the data (model parameters, forcing function) possess a
 * certain degree of smoothness.  
 *
 * Our FEM method mimics existence theory techniques using Galerkin expansions of the solution
 * about a particular basis and
 * our system of equations that are input into the computer are derived from a
 * variational form of our PDE.  Using the Gauss' theorem, 
 * \[ \begin{split}
 * -\int_{\Omega} \nabla \cdot ( a \nabla u) v \; dx 
 *  & =
 *  \int_{\Omega} a \nabla u \cdot \nabla  v \; dx
 * -\int_{\partial \Omega}  a \nabla u \cdot \textbf{n} v  \; dA
 *  \\
 *  & =
 *  \int_{\Omega} a \nabla u \cdot \nabla  v \; dx
 * -\int_{\partial \Omega}  g v  \; dA
 * \end{split} 
 * \]
 * and applying the boundary conditions, the structure
 * of the variational problem is as follows:
 * \begin{equation} \label{PhaseVariationalForm}
 *   \int_{\Omega} a  \nabla u  \cdot \nabla  v \; dx = 
 *   \int_{\Omega}  f  v \; dx + 
 *   \int_{\partial \Omega} g \; v  \; dA
 *  \end{equation}
 *
 * The following Galerkin representation of the temperature field and
 * test function is assumed.
 * 
 * \[
 * u(\mathbf{x},t) = \sum_{j=1}^{N_{dof}}
 *                     \alpha_j \phi_j(\textbf{x})
 * \qquad \qquad
 * v(\mathbf{x},t) = \sum_{i=1}^{N_{dof}}
 *                        v_i \phi_i(\textbf{x})
 * \]
 * where $N_{dof}$ is the number
 * of Galerkin coefficients, and $\phi_i$'s are the finite element shape functions
 * of polynomial order p=1,2,3...
 * Substituting our Galerkin discretization into our variational
 * form~ \eqn{PhaseVariationalForm} results in the following
 * system of equations
 * \[ \boxed{ \begin{split}
 *   \text{find } & \vec{u}=(u_1,u_2,...) \text{ s.t. } \\
 *    & {\bf K} \vec{u} = \vec{f}
 * \end{split} }  \]
 * The Jacobian of this system of equations is
 *  \[
 *    {\bf K}  = 
 * \begin{bmatrix}
 *  . & . & . \\
 *  . & \int_{\Omega} a \nabla \phi_i \nabla \phi_j & . \\
 *  . & . & . \\
 * \end{bmatrix}
 *  \qquad
 *  \vec{f} = 
 * \begin{bmatrix}
 *   . \\
 *  \int_{\Omega} f \phi_i dx 
 *   +  \int_{\partial \Omega} g \; \phi_i  \; dA \\
 *   . \\
 * \end{bmatrix}
 *  \]
 * Notice that this system of equations is ill-posed due to the pure Nuemann
 * boundary condition. Given a solution $u$, the solution plus and arbitrary
 * constant, $u+C$, is also a solution.
 * @endlatexonly
 */
template< typename MathematicalModel >
class BackgroundPhaseSystem : public PetscFEMSystem
{
public:
  // Constructor
  BackgroundPhaseSystem(EquationSystems& , 
                        const std::string& ,
			const unsigned int ) ;

  /** BVP no IC needed... do nothing */
  virtual void user_initialization (){}

  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context(DiffContext &context);

  // Constraint parts
  virtual bool element_constraint (bool request_jacobian,
                                   DiffContext& context);

  // update dirichlet bc before assembly
  virtual void assembly(bool , bool );

  // Indices for each variable;
  unsigned int b_var;

private:
  /** The type of the parent. */
  typedef PetscFEMSystem Parent;
  
  /** this class holds all constitutive data.
    * should be an instaniation of 
    *     - BackgroundPhaseModel
    */
  MathematicalModel m_MathModel; 

};
#endif
