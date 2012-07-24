#ifndef __ThermalTherapy_h
#define __ThermalTherapy_h
/**@ingroup TreatmentPlanning 
 * This class is responsible for assembling the system dynamics matrix and load
 * vector for the LITT simulation. The actual bioheat model is templated so a
 * particular source term and the constitutive models may be inlined at compile
 * time. For the linear case, just assemble System Matrix then RHS only once.
 * Then SCALE the RHS each time for the power
 * 
 * @section PennesVariationalForm Variational Form of the Pennes Model
 *
 * @htmlonly
 *   see pdf for in depth derivation of variational form of Pennes Model
 * @endhtmlonly
 * 
 * @latexonly
 * It is infeasible to provide
 * proof that a solution exists in the classical setting developed in
 * \eqn{classicalpennes}.
 *
 * A variational method is used to weaken the regularity requirements on the
 * solution and search for the solution in a broader space of functions.
 * In order to do this, the sense in which the equality in 
 * \eqn{classicalpennes} holds must be reinterpreted.
 * 
 * A space-time formulation is developed 
 * in anticipation of the dual/adjoint formulation used in the optimization
 * process. As will be seen, the space of test functions will be the same
 * space of the adjoint solution,
 * and a variational formulation is sought in which the space
 * of the trial and test functions are identical. 
 * Proceeding formally, one multiplies the Pennes equation by a
 * \textit{smooth space-time}
 * test function, $v(\textbf{x},t)$, and integrates over~$\Omega$ and $(0,T)$.
 * Using the Gauss' theorem, 
 * \[
 *  \int_{\Omega}-\nabla \cdot (k(u,\textbf{x},\beta) \nabla u) v \; dx =
 *  \int_{\Omega} k(u,\textbf{x},\beta) \nabla u \cdot \nabla  v \; dx
 * -\int_{\partial \Omega}  k(u,\textbf{x},\beta) \nabla u \cdot \textbf{n} v  \; dx
 * \]
 * and applying the boundary conditions, the structure
 * of the variational problem is as follows:
 * \begin{equation} \label{smoothvariational} \begin{split}
 *    \int_0^T \int_{\Omega} \rho  c_p \frac{\partial u}{\partial t} \;v \; dx
 * & 
 *   +\int_0^T \int_{\Omega} \omega(u,\textbf{x},\beta) c_{blood} (u - u_a) \; v \; dx 
 * \\
 * & 
 *   +\int_0^T \int_{\Omega} k(u,\textbf{x},\beta)  \nabla u  \cdot \nabla  v \; dx
 *   +\int_0^T \int_{\partial \Omega_C} h ( u - u_{\infty} ) \; v  \; dx
 * \\
 *  &  =\int_0^T \int_{\Omega} Q_{laser}(\beta,\textbf{x},t) v  \; dx
 *   -\int_0^T \int_{\partial \Omega_N} \mathcal{G} \; v  \; dA
 * \end{split} \end{equation}
 * Sufficient regularity has been tacitly assumed to
 * justify the above computations. 
 * To provide a precise mathematical meaning to the variational problem,
 * the proper setting in which each of the integrals in
 * \eqn{smoothvariational} holds must be determined.
 * 
 * Given the assumptions on
 * the thermal conductivity $k(u,\textbf{x},\beta)$ 
 * and the perfusion $\omega(u,\textbf{x},\beta)$,  the corresponding
 * integrals involving these terms are well defined provided that
 * $u,v \in H^1(\Omega)$ for a.e. $t \in (0,T]$. 
 * Here, $H^1(\Omega)$ is the classical Sobolev space of the equivalence
 * class of  $L_2(\Omega)$ functions whose distributional derivatives also 
 * belong to the equivalence class of $L_2(\Omega)$ functions.
 * \[
 *   H^1(\Omega) \equiv 
 *    \left\{ 
 *      f \in L_2(\Omega) :
 *       \frac{\partial f}{\partial x_i} \in L_2(\Omega)     
 *                       \text{ for } i = 1,2,3
 *    \right\} 
 * \]
 * More precisely,
 * interpreting $u,v$ as mappings from the time domain to the appropriate
 * function space
 * \[
 * u(\textbf{x},t) \equiv (u(t))(\textbf{x}): (0,T) \rightarrow H^1(\Omega) 
 * \qquad
 * \qquad
 * v(\textbf{x},t) \equiv (v(t))(\textbf{x}): (0,T) \rightarrow H^1(\Omega) 
 * \]
 * the integrals involving the thermal conductivity and perfusion in 
 * \eqn{smoothvariational} are well defined.
 * The appropriate space containing such mappings is defined below.
 * \begin{defn}
 * \[\boxed{ \begin{split}
 * \mathcal{V} \equiv L_2\left(0,T ; H^1(\Omega)\right) 
 * \end{split} } \]
 * where $ u \in L_2\left(0,T;H^1(\Omega)\right) \Rightarrow $
 * \begin{itemize}
 * \item
 * $u$ is a mapping from the time domain to the spatial domain,
 * \[
 *     u:(0,T) \rightarrow H^1(\Omega)
 * \]  
 * \item
 * the $H^1$ norm of $u$ is $L_2$ integrable in time,
 * \[
 *    \int_0^T \|u(\cdot,t)\|_{H^1}^2 dt < \infty
 * \]
 * \end{itemize}
 * \end{defn}
 * 
 * A precise form of \eqn{smoothvariational} also requires that the 
 * laser source term, $Q_{laser}$, and the time derivative
 * of the solution, $\frac{\partial u}{\partial t}$, exist in the
 * dual space of the test functions.
 * \[
 * Q_{laser}, \frac{\partial u}{\partial t}
 *   \in \mathcal{V}' = \left(L_2(0,T ; H^1(\Omega))\right)' 
 *                     = L_2\left(0,T ; (H^1(\Omega))' \right) 
 * \] 
 * To be precise,
 * \begin{equation} 
 * \label{dualcharacterization},
 * \begin{split}
 *   h \in & \left(L_2(0,T ; H^1(\Omega))\right)'   \Rightarrow 
 *    \exists g \in L_2(0,T ; (H^1(\Omega))'   ): \\
 * & <h,y> = \int_0^T <g(t),y(t)>_{(H^1(\Omega))' \times H^1(\Omega)} \; dt 
 *                             \qquad  y \in L_2(0,T ; H^1(\Omega))
 * \end{split} \end{equation}
 * and the derivative $\frac{\partial u}{\partial t}$ is understood in 
 * the sense of distributions.
 * 
 * Finally, the boundary terms in \eqn{smoothvariational} are well defined
 * provided  $\mathcal{G} \in (H^{\frac{1}{2}}(\partial \Omega))'$ 
 * and the values of the solution and test function
 * are interpreted in the trace sense.
 * \[  \begin{split}
 *  \left| \int_{\partial \Omega}  \gamma u \; \gamma v  \; dx \right|  
 * \leq  \|\gamma u \|_{L_2(\partial \Omega)}  \|\gamma v \|_{L_2(\partial \Omega)}
 * & \leq  \| \gamma u \|_{H^{\frac{1}{2}}(\partial \Omega)}   \|\gamma v \|_{H^{\frac{1}{2}}(\partial \Omega)}  
 * \\
 * & \leq   C \|u\|_{H^1(\Omega)}   \|v\|_{H^1(\Omega)}  
 * \end{split} \]
 * where 
 * $\gamma : H^1(\Omega) \rightarrow H^{\frac{1}{2}}(\partial \Omega)$
 * is the trace operator.
 * 
 * The variational problem may now be stated.
 * \begin{equation} \label{showalterweakform}
 * \boxed{
 * \begin{split}
 *   \text{ Given} & \text{ a set of model parameters, }  
 *                              \beta  \text{, find}  \\
 *  u(\textbf{x},&t) \in \mathcal{V} 
 *     \text{ with } \frac{\partial u}{\partial t} \in \mathcal{V}' : \\
 *   & \rho c_p \frac{\partial u}{\partial t}(\textbf{x},t) + \mathcal{A}(\beta)u(\textbf{x},t) 
 *                       =  \mathcal{F}(\beta) \in \mathcal{V}' \\
 *   &   \qquad  u(\textbf{x},0) = u^0
 * \end{split}
 * }
 * \end{equation}
 * where $\mathcal{A}: \mathcal{V} \rightarrow \mathcal{V'}$ 
 * and denoting the duality 
 * pairing by $<\cdot,\cdot>_{\mathcal{V}'\times \mathcal{V}}$,
 * \begin{equation} \label{variationaldualityleft}
 * \begin{split}
 *  \left< \mathcal{A}(\beta)u(\textbf{x},t) , v \right>_{\mathcal{V}'\times \mathcal{V}} = 
 *    \int_0^T \int_{\Omega}  &
 * \left[
 *   k(u,\textbf{x},\beta)  \nabla u   \cdot \nabla  v 
 *  + \omega(u,\textbf{x},\beta)  c_{blood} u  \; v \; \right] dx dt\\
 * &+ \int_0^T \int_{\partial \Omega_C} h u \; v  \; dA dt
 * \end{split}
 * \end{equation}
 * and
 * \begin{equation} \label{variationaldualitysource}
 * \begin{split}
 *  \left< \mathcal{F}(\beta) , v \right>_{V'\times V} = & \int_0^T \int_{\Omega} 
 *                          Q_{laser}(\beta,\textbf{x},t) v \; dx dt
 * + \int_0^T \int_{\partial \Omega_C} h u_{\infty} \; v  \; dA dt \\
 * & + \int_0^T \int_{\Omega} 
 *                 \omega(u,\textbf{x},\beta)  c_{blood}  u_a \; v \;  dx dt
 *  - \int_0^T \int_{\partial \Omega_N} \mathcal{G} \; v  \; dA
 * \end{split} \end{equation}
 * As is typical for a function defined in $L_2(0,T)$,
 * we have the right to speak of a function $L_2\left(0,T;H^1(\Omega)\right)$
 * only for a.e. $t \in (0,T]$.
 * However, important calculus results~\cite{Showalter97} 
 * on the space $L_2\left(0,T;H^1(\Omega)\right)$,
 * with time derivatives in $L_2\left(0,T;(H^1(\Omega)\right)')$ 
 * justify the right to speak of the solution point-wise in time so that
 * the initial condition $u(x,0)=u^0$ is well-defined.
 * 
 * 
 * The theory on parabolic equations as found in Showalter~\cite{Showalter97}
 * provides conditions under which 
 * \eqn{showalterweakform}
 * may be shown to have a solution. In particular,
 * if the operator $\mathcal{A}: \mathcal{V} \to \mathcal{V}'$ is
 * \begin{itemize}
 * \item bounded, i.e. $\exists$ $C$ s.t.
 * \[ 
 * \| \mathcal{A}u\|_{\mathcal{V}'} \leq  C \; \| u\|_{\mathcal{V}} 
 *    \qquad \forall u \in  \mathcal{V}
 * \]
 * \item coercive: $\exists$ $\alpha > 0$ s.t.
 * \[ 
 *  \left< \mathcal{A}v , v \right>_{\mathcal{V}'\times \mathcal{V}} 
 *   \geq
 *    \alpha   \; \| u \|^2_{\mathcal{V}} 
 *    \qquad \forall v \in  \mathcal{V}
 * \]
 * \item type $\mathcal{M}$:
 * \[ 
 * u_n \rightharpoonup u ,  
 * A u_n \rightharpoonup f ,
 * \text{ lim sup}A(u_n)u_n \leq f(u)
 * \qquad  
 * \Rightarrow
 * \qquad  
 * A u = f  
 * \]
 * \end{itemize}
 * then \eqn{showalterweakform} has a solution~\cite{Showalter97}.
 * Coercivity may be shown from the realistic 
 * assumptions that thermal conductivity, perfusion, and coefficient
 * of cooling are bounded away from zero.
 * \[
 * \begin{split}
 *  \left< \mathcal{A}(v),v \right>_{\mathcal{V}'\times \mathcal{V}} & = 
 *    \int_0^T \int_{\Omega}  
 * \left[
 *   k(v,\textbf{x},\beta)  \nabla v   \cdot \nabla  v 
 *  + \omega(v,\textbf{x},\beta)  c_{blood} v  \; v \; \right] dx dt
 * \\
 * & \hspace{2in} + \int_0^T \int_{\partial \Omega_C} h v \; v  \; dA dt
 * \\
 *  & \geq
 *    \int_0^T k_* \| \nabla v\|^2_{L_2(\Omega)} \; dt
 *    + \int_0^T \omega_* c_{blood} \| v \|^2_{L_2(\Omega)} \; dt
 * \\
 *  & \geq
 *     \text{min}(k_*,\omega_*) \| v\|^2_{\mathcal{V}}
 * \end{split}
 * \]
 * where positivity of the boundary integral was used and
 * $k_* $ and $\omega_*$ denote bounds on the thermal conductivity and perfusion.
 * \[
 * 0 < k_* < k(u,\textbf{x},\beta)
 *   \qquad
 * 0 < \omega_* < \omega(u,\textbf{x},\beta)
 * \]
 * Upper bounds on the thermal conductivity and perfusion 
 * \[
 *   k(u,\textbf{x},\beta) < k^* <  \infty
 *   \qquad
 *  \omega(u,\textbf{x},\beta)  < \omega^* < \infty
 * \]
 * may be used to show boundedness. Using H\"older's inequality and rearranging
 * \[ 
 * \forall u,v
 * \qquad 
 *  \frac{ \left| 
 *     \left< \mathcal{A}(u),v \right>_{\mathcal{V}'\times \mathcal{V}}  \right| }
 *      {\|v\|_{\mathcal{V}}}
 *  \leq 
 *          \text{max}(k^*,\omega^*) \|u\|_{\mathcal{V}}
 * \]
 * Thus, $ \text{max}(k^*,\omega^*) \|u\|_{\mathcal{V}}$ is an upper bound of
 * $ \frac{ \left| 
 *     \left< \mathcal{A}(u),v \right>_{\mathcal{V}'\times \mathcal{V}}  \right| }
 *      {\|v\|_{\mathcal{V}}} $ for all $v$. 
 * The supremum is the least upper bound which implies that
 * \[ 
 *   \sup_{v \in \mathcal{V}, v \neq 0}
 *  \frac{ \left| 
 *     \left< \mathcal{A}(u),v \right>_{\mathcal{V}'\times \mathcal{V}}  \right| }
 *      {\|v\|_{\mathcal{V}}}
 *  =
 *   \| \mathcal{A}u\|_{\mathcal{V}'} 
 *  \leq 
 *          \text{max}(k^*,\omega^*) \|u\|_{\mathcal{V}}
 * \qquad 
 * \forall u
 * \]
 * 
 * The Lipschitz continuity of the constitutive equations
 * and a generalization of the Rellich-Kondrachov theorem
 * to the space $\mathcal{V}$
 * \begin{equation} \label{rellichkondrachov}
 *   u_n \overset{ L_2\left(0,T ; H^1(\Omega)\right) }{\rightharpoonup} u 
 * \qquad \Rightarrow \qquad
 *   u_n \overset{ L_2\left(0,T ; L_2(\Omega)\right) }{\rightarrow} u 
 * \end{equation}
 * may be used to show that the operator
 * $\mathcal{A}: \mathcal{V} \to \mathcal{V}'$ is weakly closed
 * \[
 *  u_n \rightharpoonup u , \mathcal{A} u_n \rightharpoonup f
 *  \qquad
 *   \Rightarrow
 *  \qquad
 *  \mathcal{A}u = f 
 * \]
 * and, hence, is of type  $\mathcal{M}$. 
 * Weak convergence of space-time functions in $\mathcal{V}$ is a straight
 * forward generalization of the classical theory on the domain $\Omega$.
 * \[
 *   u_n \overset{\mathcal{V}}{\rightharpoonup} u 
 *    \Longleftrightarrow
 * \left|
 * \int_0^T \int_{\Omega}    
 *    \left(  (u-u_n)  \; v \; + \nabla (u-u_n)  \cdot \nabla  v  \right)
 *  dx dt
 * \right| \rightarrow 0  \qquad  \forall v \in \mathcal{V}
 * \]
 * From the definition of weak convergence and 
 * using the characterization of $\mathcal{V}'$ \eqn{dualcharacterization},
 * $\mathcal{A} u_n \rightharpoonup f$ is understood as
 * \begin{equation} \label{operatorconvergence} \left| \begin{split}
 * & \int_0^T \int_{\Omega}  
 * \left[
 * \begin{split}
 *   &  k(u_n,\textbf{x},\beta)  \nabla u_n   \cdot \nabla  v 
 * \\
 *  & +  \omega(u_n,\textbf{x},\beta)  c_{blood} u_n   \; v \; 
 * \end{split} \right] dx dt
 *   \\
 *  & + \int_0^T \int_{\partial \Omega_C} h u_n \; v  \; dA dt 
 * -
 * \int_0^T \int_{\Omega}    
 *    \left(  f  \; v \; + \nabla f  \cdot \nabla  v  \right)
 *  dx dt
 * \end{split} 
 * \right|
 * \rightarrow 0 \qquad \forall v
 * \end{equation}
 * An element of the weakly converging sequence, $u_n$,
 *  may be chosen such that $\mathcal{A} u - f$ is as small as desired.
 * \begin{equation}  \label{typemproof1}
 * \left| \begin{split}
 * & \int_0^T \int_{\Omega}  
 * \left[
 * \begin{split}
 *   &  k(u,\textbf{x},\beta)  \nabla u   \cdot \nabla  v 
 * \\
 *  & +  \omega(u,\textbf{x},\beta)  c_{blood} u   \; v \; 
 * \end{split} \right] dx dt
 *   \\
 *  & \hspace{1.0in} + \int_0^T \int_{\partial \Omega_C} h u \; v  \; dA dt 
 * \end{split} 
 * -
 * \int_0^T \int_{\Omega}    
 *    \left(  f  \; v \; + \nabla f  \cdot \nabla  v  \right)
 *  dx dt
 * \right|
 * \end{equation}
 * \[
 * =
 * \left| \begin{split}
 * & \int_0^T \int_{\Omega}  
 * \left[
 * \begin{split}
 *   & ( k(u,\textbf{x},\beta)  \nabla u  \pm k(u_n,\textbf{x},\beta)  \nabla u_n )
 *                   \cdot \nabla  v 
 * \\
 *  & + c_{blood}( \omega(u,\textbf{x},\beta)   u \pm 
 *                                 \omega(u_n,\textbf{x},\beta)   u_n)  \; v \; 
 * \end{split} \right] dx dt
 *   \\
 *  & \hspace{1.0in} + \int_0^T \int_{\partial \Omega_C} h (u \pm u_n) \; v 
 *                                                                    \; dA dt 
 * \\
 * & \hspace{1.0in}-
 * \int_0^T \int_{\Omega}    
 *    \left(  f  \; v \; + \nabla f  \cdot \nabla  v  \right)
 *  dx dt
 * \end{split} 
 * \right|
 * \]
 * Using the definition of the weak convergence of the
 *  operator~\eqn{operatorconvergence} and a trace theorem and
 * Rellich-Kondrachov~\eqn{rellichkondrachov} for the boundary term, 
 * \eqn{typemproof1} reduces to
 * \[
 * \leq \epsilon + 
 * \left|
 * \int_0^T \int_{\Omega}  
 * \left[
 * \begin{split}
 *   & ( k(u,\textbf{x},\beta)  \nabla u  - k(u,\textbf{x},\beta)  \nabla u_n )
 *                   \cdot \nabla  v 
 * \\
 *  & + ( k(u,\textbf{x},\beta) \nabla u_n - k(u_n,\textbf{x},\beta)  \nabla u_n )
 *                   \cdot \nabla  v 
 * \\
 *  & + c_{blood}( \omega(u,\textbf{x},\beta)   u - 
 *                                 \omega(u,\textbf{x},\beta)   u_n)  \; v \; 
 * \\
 *  & + c_{blood}( \omega(u,\textbf{x},\beta)   u_n - 
 *                                 \omega(u_n,\textbf{x},\beta)   u_n)  \; v \; 
 * \end{split} \right] dx dt
 * \right|
 * \]
 * and finally using that the constitutive equations are Lipschitz
 * and that a weakly convergent sequence is bounded
 * \[ \begin{split}
 * \leq  & \tilde{\epsilon} + 
 * \left|
 * \int_0^T \int_{\Omega}  
 * \left[
 * \begin{split}
 *  & ( k(u,\textbf{x},\beta)  - k(u_n,\textbf{x},\beta) ) \nabla u_n 
 *                   \cdot \nabla  v 
 * \\
 *  & + c_{blood}( \omega(u,\textbf{x},\beta)  - \omega(u_n,\textbf{x},\beta)  )
 *                                                             u_n  \; v \; 
 * \end{split} \right] dx dt
 * \right|
 * \\
 * \leq  & \tilde{\epsilon} + C 
 *      \left( \begin{split}
 *            &\| k(u,\textbf{x},\beta)  - k(u_n,\textbf{x},\beta)\|_{L_2\left(0,T ; L_2(\Omega)\right)}
 *           \\ &+
 *            \| \omega(u,\textbf{x},\beta)  - \omega(u_n,\textbf{x},\beta)   \|_{L_2\left(0,T ; L_2(\Omega)\right)}
 *      \end{split} \right)
 *  \\
 * \leq & \tilde{\epsilon} +  \tilde{C} \|u-u_n\|_{L_2\left(0,T ; L_2(\Omega)\right)}
 * \end{split} \]
 * the strong convergence in $L_2\left(0,T ; L_2(\Omega)\right) $
 * implies the result.  Additional assumptions
 * on the constitutive equations may be required to prove monotinicity
 * of the operator $\mathcal{A}: \mathcal{V} \to \mathcal{V}'$ 
 * and hence uniqueness.
 * 
 * From the above conditions, we may begin to identify a
 * well-defined parameter space
 * over which a solution is expected to exist.
 * The parameter space used in the remainder of this dissertation is
 * defined below.
 * \begin{equation} \label{parameterspace}
 *     \mathbb{P} = 
 *         \left\{ \begin{split}
 * \beta \in
 * & L_\infty(\Omega) \times \mathbb{R}^3  \times 
 *   L_\infty(\Omega) \times \mathbb{R}^3 \times 
 *   L_\infty([0,T] ) \times \mathbb{R}^5   : \\
 * & 0 < k_* < k_0(\textbf{x}) + k_1 \text{atan}(k_2(u-k_3)) < k^* <  \infty \\
 * & 0 < \omega_* < 
 *    \omega_0(\textbf{x}) + \omega_1 \text{atan}(\omega_2(u-\omega_3)) 
 *   < \omega^* < \infty \\
 *         \end{split} \right\}
 * \end{equation}
 * 
 * The sense in which a solution
 * satisfies \eqn{showalterweakform} is not strictly equivalent to the 
 * solution of the original PDE \eqn{classicalpennes}.
 * For a solution to the variational problem \eqn{showalterweakform},
 * classical arguments using the Lebesgue or Fourier Lemma may be used to recover
 * the PDE in the sense of distributions, i.e.
 * \[
 *  \left< \text{ Pennes PDE } , v \right>_{\mathcal{D}' \times \mathcal{D} } 
 *      =  0  \qquad  \forall v \in \mathcal{D}
 * \]
 * however the remaining boundary terms
 * \[
 *  \int_0^T \int_{\partial \Omega}  
 *     k(u,\textbf{x},\beta) \nabla u \cdot \textbf{n} \; v  \; dA dt 
 *  = 
 *  \int_0^T \int_{\partial \Omega_N} \mathcal{G} \; v  \; dA dt 
 *  + 
 *  \int_0^T \int_{\partial \Omega_C} h (u-u_{\infty}) \; v  \; dx dt 
 * \]
 * may not be used to recover the original boundary conditions of
 * the PDE because the trace of the normal derivative of an $H^1(\Omega)$
 * function may not be well-defined. Additional regularity on 
 * the solution is required
 * to recover the boundary conditions of the PDE \eqn{classicalpennes}.
 *
 *
 * @endlatexonly
 *
 * 
 * @section PennesAdjointGradient Adjoint Gradient
 * 
 * @htmlonly
 *   see pdf for in depth derivation of adjoint gradient
 * @endhtmlonly
 * 
 * @latexonly
 * 
 * The Adjoint Gradient of the quantity of interest is constructed from 
 * the derivative
 * of the discretized equations with respect to a single model
 * variable.
 * The chain rule is used to compute the gradient of the quantity of interest 
 * for the optimization. The initial condition does not
 * depend on the model parameters, $\frac{\partial u_0}{\partial \beta_i}=0$.
 * \[ 
 * \frac{\partial}{\partial \beta_i} Q(u(\beta,\textbf{x},t), \beta ) 
 *    = \sum_{k=1}^{N_{step}}  \frac{\partial Q}{\partial u_k}
 *                             \frac{\partial u_k}{\partial \beta_i}
 * \]
 * The derivative of the discretized state equations~\eqn{diseqn} with respect
 *  to a single model variable yields the following.
 * \begin{align*}
 * \frac{\partial C(u,\beta,v)}{\partial \beta_i} & = 
 * \Delta t_k \; \int_{\Omega}  \frac{\rho c_p}{\Delta t_k} 
 *   \left( \frac{\partial u_k}{\partial \beta} -
 *          \frac{\partial u_{k-1}}{\partial \beta} \right )  v_k \;dx \\
 * & + \Delta t_k \int_{\Omega}  
 *          \frac{\partial k}{\partial u}(u_{\kmhalf},\textbf{x},\beta)
 *            \frac{1}{2}   \left[ \frac{\partial u_{k-1}}{\partial \beta}  +
 *                          \frac{\partial u_k}{\partial \beta}  \right] 
 *                                                          \nabla u_{\kmhalf} 
 * \cdot \nabla  v_k \; dx \\
 * & + \Delta t_k \int_{\Omega}  
 *     k(u_{\kmhalf},\textbf{x},\beta) \frac{1}{2} 
 *     \left[    \nabla \frac{\partial u_{k-1}  }{\partial \beta}
 *             + \nabla \frac{\partial u_{k}  }{\partial \beta}
 *     \right] \cdot \nabla  v_k \; dx \\
 * \intertext{}
 * &  +\Delta t_k \int_{\Omega} c_{blood} 
 *    \frac{\partial \omega}{\partial u}(u_{\kmhalf},\textbf{x},\beta) 
 *    \frac{1}{2}  \left[    \frac{\partial u_{k-1}}{\partial \beta}
 *             +  \frac{\partial u_k  }{  \partial \beta} \right] 
 *                          \left(  u_{\kmhalf} - u_a   \right)
 *    v_k dx \\
 * &  +\Delta t_k \int_{\Omega} c_{blood} 
 *    \omega(u_{\kmhalf},\textbf{x},\beta) \frac{1}{2} 
 *    \left[  \frac{\partial u_{k-1}  }{\partial \beta} + 
 *            \frac{\partial u_{k}  }{\partial \beta} \right] 
 *            v_k \; dx \\
 * & +\Delta t_k \int_{\partial \Omega_C}  \frac{h}{2} 
 *          \left[ \frac{\partial u_{k-1}  }{\partial \beta} + 
 *                \frac{\partial u_k  }{\partial \beta} \right]\; v_k \; dA \\
 * & +\Delta t_k \int_{\Omega}  
 *    \frac{\partial k}{\partial \beta}(u_{\kmhalf},\textbf{x},\beta) \nabla u_{\kmhalf}
 *                    \cdot \nabla  v_k \; dx \\
 * &  +\Delta t_k \int_{\Omega} c_{blood}  
 *   \frac{\partial \omega}{\partial \beta}(u_{\kmhalf},\textbf{x},\beta) 
 *             \left(u_{\kmhalf}  - u_a \right)  v_k \; dx  \\
 * &  -\Delta t_k \; \int_{\Omega} 
 *  \frac{\partial Q_{laser}}{\partial \beta}(\beta,\textbf{x},t_k) v_k \; dx 
 * = 0 \qquad k  = 1,2, ..., N_{step} 
 * \end{align*}
 * Solving for the adjoint variable, $p_k$, such that
 * \[ \begin{split}
 * \Delta t_k \; \int_{\Omega}& \frac{\rho c_p}{\Delta t_k} 
 *               \hat{u} p_k 
 *  +\frac{1}{2} \frac{\partial k}{\partial u}(u_{\kmhalf},\textbf{x},\beta) 
 *         \hat{u} \nabla u_{\kmhalf}\cdot \nabla  p_k  \; dx
 * \\
 *  + &  \Delta t_k \; \int_{\Omega}
 *  \frac{1}{2} \frac{\partial \omega}{\partial u}
 *             (u_{\kmhalf},\textbf{x},\beta) \hat{u}
 *                         \left( u_{\kmhalf} - u_a \right) p_k  \; dx
 * \\
 * & + \Delta t_k \int_{\Omega}  
 *    \frac{1}{2} k(u_{\kmhalf},\textbf{x},\beta) \nabla \hat{u}
 * \cdot \nabla  p_k \; dx 
 *  + \frac{1}{2} \omega(u_{\kmhalf},\textbf{x},\beta) \hat{u}   p_k \; dx \\
 * & +\Delta t_k \int_{\partial \Omega_C}  \frac{h}{2} 
 *                \hat{u}  p_k \; dA
 * = \int_{\Omega}  \delta Q_k \hat{u} \; dx 
 *                   \qquad \forall \hat{u} , \qquad k = N_{step}
 * \end{split} \]
 * and
 * \[ \begin{split}
 * \Delta t_k \; \int_{\Omega}& \frac{\rho c_p}{\Delta t_k} 
 *               \hat{u} p_k 
 *  +\frac{1}{2} \frac{\partial k}{\partial u}(u_{\kmhalf},\textbf{x},\beta) 
 *         \hat{u} \nabla u_{\kmhalf}\cdot \nabla  p_k  \;
 * \\
 * & + \Delta t_k \; \int_{\Omega}
 *  \frac{1}{2} \frac{\partial \omega}{\partial u}
 *              (u_{\kmhalf},\textbf{x},\beta) \hat{u}
 *                         \left( u_{\kmhalf} - u_a \right) p_k  \; dx
 * \\
 * & + \Delta t_k \int_{\Omega}  
 *    \frac{1}{2} k(u_{\kmhalf},\textbf{x},\beta) \nabla \hat{u}
 * \cdot \nabla  p_k \; dx 
 *  + \frac{1}{2} \omega(u_{\kmhalf},\textbf{x},\beta) \hat{u}   p_k \; dx \\
 * & +\Delta t_k \int_{\partial \Omega_C}  \frac{h}{2} 
 *                \hat{u}  p_k \; dA 
 *   = \int_{\Omega}  \delta Q_k \hat{u} \; dx  \\
 * \end{split} \]
 * \[ \begin{split}
 * &\; -\left(
 * \begin{split}
 *   -\Delta t_{k+1} \; \int_{\Omega}& \frac{\rho c_p}{\Delta t_{k+1}} 
 *                  \hat{u} p_{k+1} 
 *     +\frac{1}{2} \frac{\partial k}{\partial u}(u_{\kphalf},\beta) 
 *            \hat{u} \nabla u_{\kphalf}\cdot \nabla  p_{k+1}  \\
 *    & + \Delta t_{k+1} \; \int_{\Omega} \frac{1}{2} 
 *                \frac{\partial \omega}{\partial u}(u_{\kphalf},\beta) \hat{u}
 *                            \left( u_{\kphalf} - u_a \right) p_{k+1}  \; dx
 *    \\
 *    & + \Delta t_{k+1} \int_{\Omega}  
 *       \frac{1}{2} k(u_{\kphalf},\beta) \nabla \hat{u}
 *    \cdot \nabla  p_{k+1} \; dx  \\
 *    & + \Delta t_{k+1} \int_{\Omega}  
 *       \frac{1}{2} \omega(u_{\kphalf},\beta) \hat{u}   p_{k+1} \; dx 
 *    \\
 *    &  +\Delta t_{k+1} \int_{\partial \Omega_C}  \frac{h}{2} 
 *                   \hat{u}  p_{k+1} \; dA \\
 * \end{split} \right ) \\
 * &\hspace{2.0in} \forall \hat{u} , \qquad k = N_{step}-1,N_{step}-2,...,1
 * \end{split} \]
 * implies that the numerical gradient of the quantity of interest may be 
 * computed as follows.
 * \[ 
 *   \frac{\partial Q(u,\beta)}{\partial \beta}    =  \sum_{k=1}^{N_{step}}
 * \left(
 * \begin{split}
 *  -\Delta t_k \int_{\Omega}  &
 *    \frac{\partial k}{\partial \beta}(u_{\kmhalf},\textbf{x},\beta) \nabla u_{\kmhalf}  
 *                               \cdot \nabla  p_k \; dx \\
 *   -\Delta t_k \int_{\Omega} & c_{blood}  
 *   \frac{\partial \omega}{\partial \beta}(u_{\kmhalf},\textbf{x},\beta) 
 *                                \left( u_{\kmhalf} - u_a \right)  p_k \; dx  \\
 *  + \Delta t_k \int_{\Omega} & 
 *   \frac{\partial Q_{laser}}{\partial \beta}(\beta,\textbf{x},t_k) p_k \; dx 
 * \end{split} 
 * \right)
 * \]
 * When the temperature field, $u(\textbf{x},t)$, and the adjoint variable,
 * $p(\textbf{x},t)$, are known, the explicit form of the first variation of 
 * any particular quantity of interest is as follows.
 * \[ \begin{split}
 * \frac{\delta Q}{\delta k_0} =
 *  -\int_0^T \int_{\Omega} &
 *                               \nabla u \cdot  \nabla p \; \hat{k}_0 \; dx dt \\
 * \frac{\delta Q}{\delta k_1} =
 *  -\int_0^T \int_{\Omega} &
 *          atan(k_2 \; (u-k_3))
 *                               \nabla u \cdot  \nabla p \; \hat{k}_1 \; dx dt \\
 * \frac{\delta Q}{\delta k_2} =
 *  -\int_0^T \int_{\Omega} &
 *          \frac{k_1 \; (u-k_3)}{1+k_2^2 \; (u-k_3)^2}
 *                               \nabla u \cdot  \nabla p \; \hat{k}_2 \; dx dt \\
 * \frac{\delta Q}{\delta k_3} =
 *  -\int_0^T \int_{\Omega} &
 *          \frac{-k_1 \; k_2}{1+k_2^2 \; (u-k_3)^2}
 *                               \nabla u \cdot  \nabla p \; \hat{k}_3 \; dx dt \\
 * \end{split} \]
 * \[ \begin{split}
 * \frac{\delta Q}{\delta \omega_0} =
 *  -\int_0^T \int_{\Omega} &
 *                        ( u - u_a )  p \; \hat{\omega}_0(\textbf{x}) \; dx dt 
 *  + \int_{\Omega} \omega_0(\textbf{x}) \; \hat{\omega}_0(\textbf{x}) \; dx \\
 * \frac{\delta Q}{\delta \omega_1} =
 *  -\int_0^T \int_{\Omega} &
 *          atan(\omega_2 (u-\omega_3))
 *                                 ( u - u_a )  p \; \hat{\omega}_1 \; dx dt \\
 * \frac{\delta Q}{\delta \omega_2} =
 *  -\int_0^T \int_{\Omega} &
 *          \frac{\omega_1 \; (u-\omega_3)}{1+\omega_2^2 (u-\omega_3)^2}
 *                                 ( u - u_a )  p \; \hat{\omega}_2 \; dx dt \\
 * \frac{\delta Q}{\delta \omega_3} =
 *  -\int_0^T \int_{\Omega} &
 *          \frac{-\omega_1 \; \omega_2} {1+\omega_2^2  (u-\omega_3)^2}
 *                                 ( u - u_a )  p \; \hat{\omega}_3 \; dx dt \\
 * \end{split} \]
 * {\small
 * \[ \begin{split}
 * \frac{\delta Q}{\delta P} =
 *   \int_0^T \int_{\Omega} &
 *       \frac{3\;\mu_a\;\mu_{tr} \exp(-\mu_{eff} \|\textbf{x}-\textbf{x}_0\|)}
 *            {4 \; \pi \;  \| \textbf{x} - \textbf{x}_0 \|} 
 *                                                  \; p \; \hat{P}(t) dx dt \\
 * \frac{\delta Q}{\delta \mu_a} =
 *   \int_0^T \int_{\Omega} &  \left(  
 *       \begin{split}
 *           3 & P(t) \; 
 *            \frac{\exp(-\mu_{eff} \|\textbf{x}-\textbf{x}_0\|)}
 *                 {4 \pi\| \textbf{x}-\textbf{x}_0\|} \;  \cdot 
 *          \\
 *           &  \left( 
 *          \mu_{tr} +\mu_a -\mu_{tr} \; \mu_a  \|\textbf{x}-\textbf{x}_0\|
 *                      \frac{3(\mu_a+\mu_{tr})}{2\mu_{eff}} 
 *            \right)
 *                                                  \; p \; \hat{\mu}_a 
 *       \end{split}
 *             \right) dx dt \\
 * \frac{\delta Q}{\delta \mu_s} =
 *   \int_0^T \int_{\Omega} & \left(
 *       \begin{split}
 *         3 & P(t) \; 
 *            \frac{\exp(-\mu_{eff} \| \textbf{x} - \textbf{x}_0 \|)}
 *                 {4 \; \pi \;  \| \textbf{x} - \textbf{x}_0 \|} \; \cdot
 *             \\
 *          &  \left(\mu_a -\mu_{tr} \; \mu_a  \|\textbf{x}-\textbf{x}_0\|
 *                      \frac{3(\mu_a(1-\gamma))}{2\mu_{eff}} \right)
 *                                                  \; p \; \hat{\mu}_s 
 *       \end{split}
 *               \right) dx dt \\
 * \frac{\delta Q}{\delta x_0} =
 *   \int_0^T \int_{\Omega} &
 *       3 P(t) \mu_a \mu_{tr}(\mu_{eff}\|\textbf{x}-\textbf{x}_0\| + 1) 
 *     \frac{ \exp(-\mu_{eff} \| \textbf{x} - \textbf{x}_0 \|)}
 *            {4 \; \pi  \| \textbf{x} - \textbf{x}_0 \|^3} (x-x_0) 
 *                                                   p  \hat{x}_0 dx dt \\
 * \frac{\delta Q}{\delta y_0} =
 *   \int_0^T \int_{\Omega} &
 *       3 P(t) \mu_a \mu_{tr}(\mu_{eff}\|\textbf{x}-\textbf{x}_0\| + 1) 
 *     \frac{ \exp(-\mu_{eff} \| \textbf{x} - \textbf{x}_0 \|)}
 *            {4 \; \pi  \| \textbf{x} - \textbf{x}_0 \|^3} (y-y_0) 
 *                                                   p  \hat{y}_0 dx dt \\
 * \frac{\delta Q}{\delta z_0} =
 *   \int_0^T \int_{\Omega} &
 *       3 P(t) \mu_a \mu_{tr}(\mu_{eff}\|\textbf{x}-\textbf{x}_0\| + 1) 
 *     \frac{ \exp(-\mu_{eff} \| \textbf{x} - \textbf{x}_0 \|)}
 *            {4 \; \pi  \| \textbf{x} - \textbf{x}_0 \|^3} (z-z_0) 
 *                                                   p \hat{z}_0 dx dt \\
 * \end{split} \]
 * }
 *  * @endlatexonly
 */
template< typename MathematicalModel  >
class LITTSystem : public ThermalTherapySystem < MathematicalModel  >
{

public:
  // Constructor
  LITTSystem(EquationSystems& ,const std::string& , const unsigned int );

  /** System initialization */
  virtual void init_data ();

  // Context initialization
  virtual void init_context(DiffContext &context);

  /** Dirichlet BC */
  virtual void SetupDirichlet(libMesh::MeshBase& );
  virtual void ApplyDirichlet (); 
  
  // Element residual and jacobian calculations
  // Time dependent parts (Doxygen at definition below)
  virtual bool element_time_derivative (bool , DiffContext& );

  // Indices for temperature, damage, damage derivative; respectively
  unsigned int u_var, a_var, b_var;

private:
  /** The type of the parent  */
  typedef ThermalTherapySystem< MathematicalModel > Parent;
  // index maps for damage 
  std::vector<unsigned int> temp_indices, damage_indices, deriv_indices ;
};

/**@ingroup TreatmentPlanning 
 * This class is responsible for assembling the system dynamics matrix and load
 * vector for the RFA simulation. The actual bioheat model is templated so a
 * particular source term and the constitutive models may be inlined at compile
 * time. 
 * Reimplements the element_time_derivative for the coupled pennes RF Solve
 */
template< typename MathematicalModel  >
class RFASystem : public LITTSystem < MathematicalModel > 
{

public:
  // Constructor
  RFASystem(EquationSystems& ,const std::string& , const unsigned int );

  /** System initialization */
  virtual void init_data ();

  /** Dirichlet BC */
  virtual void ApplyDirichlet (); 
  
  // Context initialization
  virtual void init_context(DiffContext &context)
   {
    Parent::init_context(context); // call LITTSystem
    PetscFEMContext &c = libmesh_cast_ref<PetscFEMContext&>(context);
    // use forward euler voltage
    c.SetThetaValue(this->z_var,1.0); 
   }

  /** Assemble System Dynamics Matrix and Load Vector */
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext& context);

  // Indices for voltage;  
  unsigned int z_var;

private:
  /** The type of the parent  */
  typedef LITTSystem < MathematicalModel > Parent;

  int m_ElectrodeNodeSet;
};
/**@ingroup TreatmentPlanning 
 * This class is responsible for assembling the system dynamics matrix and load
 * vector for the coupled bioheat RHTE simulation. The actual bioheat model is templated so a
 * particular source term and the constitutive models may be inlined at compile
 * time. 
 * Reimplements the element_time_derivative for the coupled pennes fluence Solve
 */
template< typename MathematicalModel  >
class RHTESystem : public LITTSystem < MathematicalModel > 
{

public:
  // Constructor
  RHTESystem(EquationSystems& ,const std::string& , const unsigned int );

  /** System initialization */
  virtual void init_data ();

  /** scatter global parameters to all locally */
  virtual void ScatterGlobalToAll()
   {
   }

  /** Assemble System Dynamics Matrix and Load Vector */
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext& context);

  // setup dirichlet
  virtual void SetupDirichlet(libMesh::MeshBase& );

  // Context initialization
  virtual void init_context(DiffContext &context)
   {
    Parent::init_context(context); // call LITTSystem
    PetscFEMContext &c = libmesh_cast_ref<PetscFEMContext&>(context);
    // use forward euler voltage
    c.SetThetaValue(this->z_var,1.0); 
   }

  unsigned int  z_var; // Indices for scattered fluence  
  unsigned int  e_var; // Indices for primary   fluence  
  unsigned int fx_var; // Indices for x-direction flux  
  unsigned int fy_var; // Indices for y-direction flux  
  unsigned int fz_var; // Indices for z-direction flux  

private:
  /** The type of the parent  */
  typedef LITTSystem < MathematicalModel > Parent;

  PetscTruth m_ExternalFiber;
};

#endif
