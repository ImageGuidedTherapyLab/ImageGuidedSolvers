#ifndef __realtimeimagingKalman_h
#define __realtimeimagingKalman_h
#include "Imaging.h"
/**@ingroup ModelAssistedMonitoring 
 * @author  David Fuentes
 *
 *  \f{itemize}
 *    \item[] Kalman, R. (1960).  {A new approach to linear filtering and prediction problems}. {\em Journal of basic Engineering}, 82(1):35--45.
 *  \f}
 * 
 * The Kalman filter follows Bayesian principles and provides the
 * fundamental framework of uncertainty quantification. 
 * It is not really a filter in the conventional sense; it
 * is a state estimator that can extract a signal (system state) from a
 * series of incomplete, corrupt,  and noisy measurements. 
 * Given initial probability distribution functions to describe a known
 * state, the Kalman algorithm provides a rigorous methodology and
 * mathematical framework to completely characterize the state
 * conditioned on all previous data, \f$ f_{x(t_i)|Z(t_i)}(\xi|Z_i) \f$. Then the 
 * "optimal" choice of the state can be chosen as the conditional expected value.
 * \f[
 *  E\{x(t_i)|Z(t_i)\}  =  \int_{-\infty}^\infty \xi f_{\xi|Z(t_i)}(\xi|Z_i) d\xi
 * \f]
 * For the particular case of Gaussian distribution, the mean and median
 * coincide. Further the condition expectation may be interprited as the
 * uncertainty estimate
 * \f[
 *  E\{(x(t_i)-E\{x(t_i)|Z(t_i)\}) (x(t_i)-E\{x(t_i)|Z(t_i)\})^T|Z(t_i)\}  
 *        =  \int_{-\infty}^\infty (\xi -  E\{x(t_i)|Z(t_i)\} ) 
 *                                 (\xi -  E\{x(t_i)|Z(t_i)\} )^T
 *                                                  f_{\xi|Z(t_i)}(\xi|Z_i) d\xi
 * \f]
 * Kalman theory is the culmanation of earlier work in filtering and
 * stochastic differential equations by
 * Wiener, Kolmogorov, Bode, Shannon, Pugachev and others with a modern
 * state space control theory approach.  The Kalman mathematical
 * technique was used for control of complex dynamic systems such as
 * spacecraft (Apollo), manufacturing processes, aircraft, or ships.
 *
 * @section KalmanSimpleExample Simple Example of Main Ideas of Kalman Filter
 *
 * A simple example will highlight the main ideas of the Kalman Filter,
 * consider the following position example.
 * \f{itemize}
 *   \item  At time $t_1$ you measure your location to be $z_1$ with 
 * some initial uncertainty, $\sigma_{z_1}$
 *   \item  Using this information we can establish the conditional probability
 *          $f(x|z_1)$ of the position, $x(t_1)$, 
 *          conditioned on the measurement, $z_1$
 * \f}
 * 
 * \image html  maybeckFig1_5.jpg "Maybeck Figure 1.5" width=5cm
 * \image latex maybeckFig1_5.jpg "Maybeck Figure 1.5" width=5cm
 *
 * \f{itemize}
 *   \item  At time $t_2\approx t_1$, trained navigator measures the
 *          location to be $z_2$ with  uncertainty, $\sigma_{z_2}<\sigma_{z_1}$
 *   \item  How do you combine both data into the probabilistic representation?
 * \f}
 * \f{itemize}
 * \item[] under appropriate assumptions, conditional pdf based on the two
 * measurements, $f(x|z_1,z_2)$, retains Gaussian form
 * \f}
 *
 * \image html  maybeckFig1_6.jpg "Maybeck Figure 1.6" width=5cm
 * \image latex maybeckFig1_6.jpg "Maybeck Figure 1.6" width=5cm
 *
 * \f{itemize}
 *   \item[]
 *  \[
 *     \mu = \left( \frac{\sigma^2_{z_2}}{\sigma^2_{z_1} + \sigma^2_{z_2}} \right) 
 *            z_1
 *         + \left( \frac{\sigma^2_{z_1}}{\sigma^2_{z_1} + \sigma^2_{z_2}} \right)
 *            z_2
 *  \]
 *  { \scriptsize
 *        $\sigma_{z_1} > \sigma_{z_2} \Rightarrow$  greater confidence in $z_2$
 *        }
 *   \item[] uncertainty decreases by combining the two pieces of information
 *  \[ 
 *         \frac{1}{\sigma^2} = 
 *         \frac{1}{\sigma^2_{z_1}} + \frac{1}{\sigma^2_{z_2}}
 *  \]
 *   \item[]
 * \f}
 * 
 * \f{itemize}
 *   \item Best estimate is
 * \[
 * \begin{split}
 *  \hat{x}(t_2) & = \mu 
 *      = \left( \frac{\sigma^2_{z_2}}{\sigma^2_{z_1} + \sigma^2_{z_2}} \right) 
 *           z_1
 *      + \left( \frac{\sigma^2_{z_1}}{\sigma^2_{z_1} + \sigma^2_{z_2}} \right)
 *           z_2
 * \\
 *    & = z_1
 *      + \left( \frac{\sigma^2_{z_1}}{\sigma^2_{z_1} + \sigma^2_{z_2}} \right)
 *            (z_2 - z_1)
 *      = {\hat{x}(t_1) + K(t_2) \left[  z_2 - \hat{x}(t_1) \right]}
 * \end{split}
 * \]
 *  \begin{itemize}
 *   \item $\hat{x}(t_1) = z_1$ is the best predicted values {before}
 *         the second measurement is made
 *   \item plus a correction term 
 *         $K(t_2)= \left( \frac{\sigma^2_{z_1}}{\sigma^2_{z_1} + \sigma^2_{z_2}} \right) $
 *         weighting difference between $z_2$ and best prediction $\hat{x}(t_1)$
 *  \end{itemize}
 *   \item Similarly, uncertainty update  in language of Kalman filter
 *         \[
 *             \sigma^2_x(t_2)=\sigma^2_x(t_1) - K(t_2) \sigma^2_x(t_1)
 *         \]
 * 
 *  \begin{itemize}
 *   \item $\hat{x}(t_1) = z_1$ is the best predicted values {before}
 *         the second measurement is made
 *   \item plus a correction term 
 *         $K(t_2)= \left( \frac{\sigma^2_{z_1}}{\sigma^2_{z_1} + \sigma^2_{z_2}} \right) $
 *         weighting difference between $z_2$ and best prediction $\hat{x}(t_1)$
 *  \end{itemize}
 *   \item Similarly, uncertainty update  in language of Kalman filter
 *         \[
 *             \sigma^2_x(t_2)=\sigma^2_x(t_1) - K(t_2) \sigma^2_x(t_1)
 *         \]
 * \f}
 * 
 * 
 * \f{itemize}
 *   \item The previous problem is referred to as the 
 *         static estimation problem. Want to incorporate dynamics
 *   \item Suppose you travel some time before taking another measurement
 *          \[
 *             \frac{dx}{dt}=u+w
 *          \]
 *     \begin{itemize}
 *         \item  Very simple model, $u$=velocity, and absorb any further
 *                higher order complexities of the modeling detail into the
 *                noise, $w$. For tractability assume Gaussian
 *                noise of zero mean and variance $\sigma_w^2$
 *     \end{itemize}
 *   \item Most importantly, the model provides a mechanism to
 *          propagate the uncertainty
 * \f}
 * 
 * \image html  maybeckFig1_7.jpg "Maybeck Figure 1.7" width=5cm
 * \image latex maybeckFig1_7.jpg "Maybeck Figure 1.7" width=5cm
 * \f{itemize}
 * \item[] at time $t_3^-$ right before the measurement at $t_3$, pdf moves
 * according to nominal dynamics model, and spreads out in time due to the
 * constant addition of uncertainty
 * \[
 *  \hat{x}(t_3^-) =  \hat{x}(t_2) + u (t_3 - t_2)
 *  \qquad
 *  \sigma_x^2(t_3^-) =  \sigma_x^2(t_2) + \sigma_w^2(t_3-t_2)
 * \]
 * \f}
 * \f{itemize}
 *   \item again take measurement $z_3$, $\sigma_{z_3}$, similar to static
 * case $(t_1,t_2) \leftrightarrow (t_3^-,t_3)$
 * \[ \begin{split} 
 *  \hat{x}(t_3) &  = \hat{x}(t_3^-) + K(t_3) \left[  z_3 - \hat{x}(t_3^-) \right]
 * \\
 *  \sigma_x^2(t_3)& = \sigma_x^2(t_3^-) - K(t_3) \sigma_x^2(t_3^-) 
 * \\
 * K(t_3) & = \left( 
 *            \frac{\sigma^2_x(t_3^-)}{\sigma^2_x(t_3^-) + \sigma^2_{z_3}} 
 *          \right) 
 * \end{split}
 * \]
 *   \item low confidence in measurement, use model
 *     \[
 *          \sigma^2_{z_3} \rightarrow \infty
 *           \Rightarrow 
 *          K(t_3) \rightarrow 0 
 *           \Rightarrow 
 *         \hat{x}(t_3)   = \hat{x}(t_3^-) 
 *     \]
 *   \item low confidence in model, use measurement
 *     \[
 *            \sigma^2_w        \rightarrow \infty,
 *            \sigma^2_x(t_3^-) \rightarrow \infty
 *           \Rightarrow 
 *          K(t_3) \rightarrow 1 
 *           \Rightarrow 
 *         \hat{x}(t_3)   = z_3
 *     \]
 * \f}
 *
 * @section KalmanDerivation Derivation of Kalman Filter
 * 
 * Derivation of the Kalman filter begins with a
 * continous time dynamics/discrete time measurement model,
 * Figure~\ref{continDynamicDiscreteMeasBlock}. 
 * The governing system pde model 
 *
 * Ignoring the mathematical technicalities of why
 * \f[
 *    \int_{t_0}^t w(\tau) dt \tau
 * \f]
 * is not well defined
 *   \f{itemize}
 *    \item dirac delta for correlation~\cite{kloeden1999numerical} pg.76
 *           \[
 *                E\{x(t_i)x(t_j)\} = \Psi(t_i) \delta(t_i-t_j)
 *           \]
 *    \item white noise should be thought of as a "generalized stochastic process"
 *          similiar to how dirac delta is a 
 *          "generalized function"~\cite{kloeden1999numerical} pg.44
 *    \item wiener process of unbounded variation~\cite{kloeden1999numerical} pg.76
 *   \f}
 * we will begin with 
 * \f[
 *   dx(t) =  F(t) x(t) dt + G(t) d\beta(t)
 * \f]
 * as our governing differential equation.
 * 
 * may be obtained from a method of lines
 * approach on the variational form corrupted by some additive noise. Most
 * importantly, the linear stochastic differential model provides a
 * mechanism to propagate the uncertainty. The
 * same discretization may be used to obtain the measurement model.
 * \f[ \begin{split}
 *   dx(t) &  = F(t) x(t) dt + B(t) u(t) dt + G(t) d\beta(t) 
 * \\
 *    z(t_i)  & = H(t_i) x(t_i) +  v(t_i)
 * \end{split}
 * \f]
 * where \f$H\f$ is the measurement matrix, \f$v\f$ is additive white noise
 * \f[
 * E\{v(t_i)\} = 0 
 * \qquad
 * E\{v(t_i)v^T(t_i)\} = 
 * \left\{
 * \begin{split} 
 * R(t_i)  \qquad t_i =    t_j \\
 * 0       \qquad t_i \neq t_j \\
 * \end{split} 
 * \right.
 * \f]
 * The characterization of the solution to the continuous dynamics model is
 * provided by the solution the linear Stochastic Differential
 * Equation~\eqn{sdeSoln}.  Derivation of the Kalman filter exploits the
 * linearity of the problem and allows the assumed intial Gaussian
 * distributions to remain Gaussian throughout the time evolution of the
 * problem.  The nice algebraic properties when manipulating linear
 * combinations and products of Gaussian properties that maintains the
 * Gaussian structyure is a theme that is continually exploited throughout
 * the derivation of the Kalman filter.  Zero mean noise implies that the
 * mean of the measurement process is given by
 * \f[
 *    m_z(t_i) = E\{z(t_i)\} 
 *             = E\{H(t_i)x(t_i)+v(t_i)\} 
 *             = H(t_i) m_x(t_i)
 * \f]
 * and the autocorrelation is given by
 * \f[
 *   \Psi_{zz}(t_j,t_i) = E\{z(t_j) z^T(t_i)\} 
 * \left\{
 *  \begin{split}
 *   H(t_j) E\{x(t_j) x^T(t_i)\} H^T(t_i)          & \qquad t_j \neq t_i
 *   \\
 *   H(t_j) E\{x(t_j) x^T(t_i)\} H^T(t_i) + R(t_i) & \qquad t_j = t_i
 *  \end{split}
 * \right.
 * \f]
 * and the covariance
 * \f[
 *   P_{zz}(t_j,t_i) = 
 * \left\{
 *  \begin{split}
 *   H(t_j) P(t_j,t_i) H^T(t_i)          & \qquad t_j \neq t_i
 *   \\
 *   H(t_j) P(t_j,t_i) H^T(t_i) + R(t_i) & \qquad t_j = t_i
 *  \end{split}
 * \right.
 * \f]
 * 
 * \image html  continDynamicDiscreteMeasBlock.jpg "Continuous-time dynamics coupled with discrete time measurement model" width=5cm
 * \image latex continDynamicDiscreteMeasBlock.jpg "Continuous-time dynamics coupled with discrete time measurement model" width=5cm
 * 
 * Using the solution of the linear SDE~\eqn{sdeSoln},
 * the state propagation from time \f$t_{i-1}\f$ to \f$t_i\f$ is given by
 * \f[
 * x(t_i) = \Phi(t_i,t_{i-1})x(t_{i-1}) + w_d(t_{i-1})
 * \f]
 * We are interested in the density of the state conditioned on the
 * measurement data \textit{immediately before} the latest measurement, 
 * \f$f_{x(t_i)|Z(t_{i-1})}\f$. To show that this is Gaussian,
 * \f$x(t_i)\f$ is a linear combination of 
 * \f$ x(t_{i-1})\f$ and \f$w_d \f$, \f$f_{x(t_i)|Z(t_{i-1})}\f$ may be shown to be
 * Gaussian if the joint pdf \f$f_{x(t_{i-1}),w_d(t_{i-1}),Z(t_{i-1})}\f$ is
 * Gaussian which is equivalent to saying that  
 * \f$b=[x(t_{i-1}),w_d(t_{i-1}),Z(t_{i-1})]\f$ is a gaussian variable, (ie
 * find the mean and covariance). This can be shown easily because
 * independent Gaussian random variables are jointly gaussian and linear
 * transforms of a gaussian is gaussian~(\cite{Maybeck79} Sec 3.11).
 * Using the charaterization of the conditional density ,\f$f_{x(t_i)|Z(t_{i-1})}\f$,
 * the expected value and covariance may be computed to characterize the state.
 * 
 * Finally, given the latest set of measurement data, the 
 * the conditional density immediately before the measurement may be updated 
 * according to Bayes rule as
 * \f{equation}\label{pdfUpdate}
 * f_{x(t_i)|Z(t_i)} 
 * = \frac{ f_{x(t_i)|Z(t_i)} } { f_{Z(t_i)} }
 * = ...
 * = \frac{ f_{z(t_i)|x(t_i),Z(t_{i-1})} f_{x(t_i)|Z(t_{i-1})} }
 *        { f_{z(t_i)|Z(t_{i-1})} }
 * \f}
 * As discussed, the linearity of the PDE model may be used to show that
 * the conditional density before the update, \f$ f_{x(t_i)|Z(t_{i-1})} \f$, is
 * Gaussian.  The linearity of the measurement model and  the Gaussian
 * assumptions on the noise terms may be used to shown that the remaining
 * terms in \eqn{pdfUpdate} are indeed Gaussian. Thus, working through the
 * algebraic manipulations of the product of two Gaussians divided by another
 * Gaussian distribution, the density conditioned on all known information
 * may be explicitly written as
 * \f[
 *     f_{x(t_i)|Z(t_i)}(\xi|Z_i)  = 
 *       \left[ (2\pi)^{n/2} |P(t_i^+)|^{1/2}\right]^{-1}
 *       \exp\left\{ -1/2 (\xi - \hat{x}(t_i^+) )^T 
 *                        P(t_i^+)^{-1}
 *                        (\xi - \hat{x}(t_i^+) )\right\}
 * \f]
 * where
 * \f[ \begin{split}
 * K(t_i)   & = P(t_i^-)H^T(t_i)[H(t_i)P(t_i^-)H^T(t_i) + R(t_i)]^{-1} \\
 * \hat{x}(t_i^+)  & = \hat{x}(t_i^-)  + K(t_i)[z_i - H(t_i)\hat{x}(t_i^-)] \\   
 * P(t_i^+) & = P(t_i^-)  - K(t_i)H(t_i)P(t_i^-)                       \\
 * \end{split} \f]
 * The goals of uncertainty quantification have been attained: the
 * expected value is the "optimal" estimate of the state and covariance
 * provides a measure of the uncertainty.
 *
 * @section KalmanSummary Summary of Kalman Algorithm
 *
 * Main data structure for Kalman filter based off of Maybeck 
 *  \f{itemize}
 *    \item[] Maybeck, P.~S. (1979).  \newblock {\em Stochastic models, estimation, and control}, volume 141 of {\em Mathematics in Science and Engineering}.
 *  \f}
 * 
 * The abstract base class assumes a stochastic system dynamics model of the
 * form
 *  \f[
 *     dx(t) =  F(t) x(t) dt + B(t) u(t) dt + d\beta(t)
 *  \f]
 * 
 * The system dynamics model of the time evolution of the state vector, \f$x\f$,
 * is assumed of the form of a linear stochastic differential equation.
 * \f[
 *   dx(t) = F(t) x(t) dt + B(t) u(t) dt + G(t) d\beta(t)
 * \f]
 * Where \f$F\f$, \f$B\f$, and \f$G\f$ are the \textit{known} the system
 * dynamics matrix, deterministic input matrix, and noise input matrix,
 * respectively.  The system inputs, modeling error, discretization error,
 * and noise that is beyond our control is represented as 
 * a vector valued Brownian motion process, \f$\beta\f$, with statistics
 * \f[
 *  E\{\beta(t)\} = 0
 * \qquad \qquad
 *  E\{ [ \beta(t)- \beta(t') ] [ \beta(t)- \beta(t') ]^T \} 
 *                                              = \int_{t'}^t Q(\tau) d\tau
 * \f]
 * Where \f$Q(t)\f$, is the \textit{symmetric positive definite} diffusion 
 * matrix and is a measure of the rate of
 * divergence of the mean squared value from its initial value.
 * Measurments, \f$z(t_i)\f$, are taken at discrete time points and are assumed
 * to be a linear transformation of the state vector, \f$x(t_i)\f$,
 * corrupted by white Gaussian noise,
 * \f$v(t_i)\f$
 * \f[
 *    z(t_i)   = H(t_i) x(t_i) +  v(t_i) \qquad i=1,2,3,...
 * \f]
 * where \f$H\f$ is the measurement matrix, and the statistics of the white
 * Gaussian noise \f$v\f$ is 
 * \f[
 * E\{v(t_i)\} = 0 
 * \qquad
 * E\{v(t_i)v^T(t_i)\} = 
 * \left\{
 * \begin{split} 
 * R(t_i)  \qquad t_i =    t_j \\
 * 0       \qquad t_i \neq t_j \\
 * \end{split} 
 * \right.
 * \f]
 * Initial conditions are given by
 * \f[
 *  \hat{x}(t_0)  = E\{x(t_0)\} = \hat{x}_0
 * \qquad \qquad
 *  P(t_0)  = E\{[x(t_0)-\hat{x}_0][x(t_0)-\hat{x}_0]^T\} = P_0
 * \f]
 * The Kalman algorithm begins by 
 * first propagating the state from \f$t_{i-1}\f$ to \f$t_i\f$
 * @anchor KalmanPrediction
 * \f{equation} \label{KalmanPrediction}
 * \boxed{
 * \begin{split}
 *  \hat{x}(t^-_i) & = \Phi(t_i,t_{i-1}) \hat{x}(t^+_{i-1}) 
 *                 + \int_{t_{i-1}}^{t_i} \Phi(t_i,\tau) B(\tau) u(\tau) d\tau 
 * \\
 *  P(t^-_i) & = \Phi(t_i,t_{i-1}) P(t^+_{i-1}) \Phi^T(t_i,t_{i-1})
 *                 + \int_{t_{i-1}}^{t_i} 
 *                 \Phi(t_i,\tau) G(\tau) Q(\tau) G^T(\tau) \Phi^T(t_i,\tau) d\tau 
 * \end{split}
 * }
 * \f}
 * Predictions are updated based on measurements at time \f$t_i\f$.
 * @anchor KalmanUpdate
 * \f{equation} \label{KalmanUpdate}
 * \boxed{
 * \begin{split}
 * K(t_i)   & = P(t_i^-)H^T(t_i)[H(t_i)P(t_i^-)H^T(t_i) + R(t_i)]^{-1} \\
 * \hat{x}(t_i^+)  & = \hat{x}(t_i^-)  + K(t_i)[z_i - H(t_i)\hat{x}(t_i^-)] \\   
 * P(t_i^+) & = P(t_i^-)  - K(t_i)H(t_i)P(t_i^-)                       \\
 * \end{split} 
 * }
 * \f}
 *
 * @section ForwardEulerKalman Forward Euler Kalman Model Propagation
 *
 * For a time invariant system dynamics matrix, \f$ F \f$, 
 * the state transition matrix is represented 
 * by a matrix exponential  (\cite{Maybeck79} Sec 2.3.11) and 
 * a Taylor
 * series approximation of the matrix exponential may be used to reduce the 
 * analytical solution 
 * \eqn{KalmanPrediction}
 * to a familar forward Euler time discretization.
 * \f{equation}\label{TaylorSeriesStateTransition}
 *    \Phi(t,t_0) = \Phi(t - t_0) = e^{F(t-t_0)} 
 *                = I + F (t-t_0) + \frac{1}{2!}F^2(t-t_0)^2  +  ...
 *                \approx I + F (t-t_0)  \text{ (small time step)} 
 * \f}
 * Substituting this approximation into the expected value of the state
 * over the time increment, \f$ \Delta t = t_i - t_{i-1} \f$
 * \f[
 *  \hat{x}(t^-_i) = \Phi(t_i,t_{i-1}) \hat{x}(t^+_{i-1}) 
 *                 + \int_{t_{i-1}}^{t_i} \Phi(t_i,\tau) B(\tau) u(\tau) d\tau 
 *    \qquad \Rightarrow \qquad
 *  \hat{x}(t^-_i) = \left( I + F \Delta t \right)  \hat{x}(t^+_{i-1}) 
 *    + \int_{t_{i-1}}^{t_i} \left( I + F (t_i-\tau) \right) B(\tau) u(\tau) d\tau 
 * \f]
 * and approximating the integral with a low order "left point" rule
 * \f[
 * \begin{split}
 *    \int_{t_{i-1}}^{t_i} \left( I + F (t_i-\tau) \right) B(\tau) u(\tau) d\tau 
 *  & \approx
 *    \left( I + F \Delta t \right) B u(t_{i-1}) \Delta t
 *  =
 *     B u(t_{i-1}) \Delta t+ F  B u(t_{i-1}) \Delta t^2
 * \\
 *  & \approx
 *     B u(t_{i-1}) \Delta t 
 *    \quad \text{ ignoring higher order terms}
 * \end{split}
 * \f]
 * yields the forward Euler approximation
 * \f[\boxed{
 *  \hat{x}(t^-_i) = \left( I + F \Delta t \right)  \hat{x}(t^+_{i-1})  
 *                 + B u(t_{i-1}) \Delta t 
 * }
 * \f]
 * Similarly for the covariance update, assuming \f$ G=I \f$
 * \f[
 *  P(t^-_i) \approx 
 * \left( I + F \Delta t \right)
 *  P(t^+_{i-1})
 * \left( I + F \Delta t \right)^T
 *      + \int_{t_{i-1}}^{t_i} 
 *        \left( I + F (t_i-\tau) \right) Q(\tau) \left( I + F (t_i-\tau)\right)^T 
 *        d\tau 
 * \f]
 * Again using the left point rule for the integral yields
 * \f{equation} \label{ForwardEulerCovariancePrediction}
 * \boxed{
 * \begin{split}
 *     P(t^-_i) &  \approx 
 *     ( I + F \Delta t )
 *     P(t^+_{i-1}) 
 *     ( I + F \Delta t )^T 
 *      + 
 *     ( I + F \Delta t )
 *     Q(t_{i-1}) 
 *     ( I + F \Delta t )^T \; \Delta t
 * \\
 *      &    \approx
 *     P(t^+_{i-1}) 
 *      + 
 *     \left(F P(t^+_{i-1}) + P(t^+_{i-1}) F^T + Q(t_{i-1}) \right) \; \Delta t
 *     \quad \text{ ignoring higher order terms}
 * \\
 *      &    =
 *     P(t^+_{i-1}) 
 *      + 
 *     \left(F P(t^+_{i-1}) + (F P(t^+_{i-1}) )^T + Q(t_{i-1}) \right) \; \Delta t
 * \end{split}
 * }
 * \f}
 *
 * @section CrankNicolsonKalman Crank Nicolson Kalman Model Propagation
 *
 * The following identity for the state transition matrix is needed to 
 * derive the Crank Nicolson time stepping.
 * \f[\begin{split}
 * \Phi(t,t_0) \Phi(t_0,t)  = I  
 *  \qquad & \Rightarrow \qquad 
 * \Phi^{-1}(t,t_0) = \Phi(t_0,t)
 * \\
 * \Phi(t_3,t_1) & = \Phi(t_3,t_2) \Phi(t_2,t_1)
 * \end{split} \f]
 * A backwards time stepping scheme for the solution is needed for the
 * approximation at the half time step.
 * \f[
 *    \Phi(t_0,t) \left( 
 *      x(t) = \Phi(t,t_0) x_0  
 *           + \int_{t_0}^{t} \Phi(t,\tau) B(\tau) u(\tau) d\tau 
 *                \right) 
 *     \Rightarrow  
 * \hspace{-2in}
 * \begin{split}
 *    x_0 &  = \Phi(t_0,t) x(t) 
 *           - \int_{t_0}^{t} \Phi(t_0,t) \Phi(t,\tau) B(\tau) u(\tau) d\tau 
 *     \\
 *        &  = \Phi(t_0,t) x(t) 
 *           + \int_{t}^{t_0} \Phi(t_0,\tau) B(\tau) u(\tau) d\tau 
 * \end{split} 
 * \f]
 * The strategy of the Crank Nicholson derivation is to propagate the 
 * solution forward and backward over the interval, \f$ \Delta t = t_i - t_{i-1} \f$,
 * beginning with the half time point 
 * \f$ \hat{x}(t^-_{i-\frac{1}{2}}) \equiv  \frac{ \hat{x}(t^-_i) + \hat{x}(t^+_{i-1}) }{2} \f$  
 * \f[\begin{split}
 *  \hat{x}(t^-_i) & = \Phi(t_i,t_{i-\frac{1}{2}}) \hat{x}(t^-_{i-\frac{1}{2}}) 
 *          + \int_{t_{i-\frac{1}{2}}}^{t_i} \Phi(t_i,\tau) B(\tau) u(\tau) d\tau 
 * \\
 *  \hat{x}(t^+_{i-1}) & = \Phi(t_{i-1},t_{i-\frac{1}{2}}) \hat{x}(t^-_{i-\frac{1}{2}}) 
 *          + \int^{t_{i-i}}_{t_{i-\frac{1}{2}}} \Phi(t_{i-1},\tau) B(\tau) u(\tau) d\tau
 * \end{split}
 * \f]
 * using \eqn{TaylorSeriesStateTransition} and the corresponding integration
 * rule at the half time step \f$ t_{i-\frac{1}{2}} \f$ to approximate the integral
 * \f{equation}\label{ForwardTaylorSeries}
 *  \hat{x}(t^-_i)  = 
 *       \left(I + F \frac{\Delta t}{2} \right) 
 *              \hat{x}(t^-_{i-\frac{1}{2}}) 
 *          + \int_{t_{i-\frac{1}{2}}}^{t_i} 
 *               \left(I + F (t_i-\tau)  \right) 
 *             B(\tau) u(\tau) d\tau 
 *                  \approx 
 *       \left(I + F \frac{\Delta t}{2} \right) 
 *       \left( \hat{x}(t^-_{i-\frac{1}{2}}) 
 *               + 
 *             B(t_{i-\frac{1}{2}}) u(t_{i-\frac{1}{2}}) \frac{\Delta t}{2}
 *       \right) 
 * \f}
 * \f{equation} \label{BackwardTaylorSeries}
 *  \hat{x}(t^+_{i-1})  = 
 *    \left(I - F \frac{\Delta t}{2} \right) 
 *              \hat{x}(t^-_{i-\frac{1}{2}}) 
 *          + \int^{t_{i-i}}_{t_{i-\frac{1}{2}}} 
 *               \left(I + F (t_{i-1}-\tau)  \right) 
 *             B(\tau) u(\tau) d\tau
 *                  \approx 
 *       \left(I - F \frac{\Delta t}{2} \right) 
 *       \left( \hat{x}(t^-_{i-\frac{1}{2}}) 
 *               -  
 *             B(t_{i-\frac{1}{2}}) u(t_{i-\frac{1}{2}}) \frac{\Delta t}{2}
 *       \right) 
 * \f}
 * Notice that the additional minus sign is due to integrands. 
 * Subtracting \eqn{BackwardTaylorSeries} from \eqn{ForwardTaylorSeries}
 * and dropping higher order terms yields 
 * the Crank Nicholson scheme.
 * \f[
 * \boxed{
 *  \hat{x}(t^-_i) - \hat{x}(t^+_{i-1})  = 
 *   F \Delta t \hat{x}(t^-_{i-\frac{1}{2}}) 
 *             + 
 *             B(t_{i-\frac{1}{2}}) u(t_{i-\frac{1}{2}}) \Delta t
 * }
 * \f]
 * Similarly for the covariance update at the half time point, 
 * \f$ P(t^-_{i-\frac{1}{2}})\equiv\frac{P(t^-_i)  + P(t^+_{i-1})}{2} \f$,
 * assuming \f$ G=I \f$ 
 * \f[
 *     P(t^-_i)   \approx 
 *     \left( I + F \frac{\Delta t}{2} \right)
 *     \left( P(t^-_{i-\frac{1}{2}})  + Q(t_{i-\frac{1}{2}}) \; \frac{\Delta t}{2} \right)
 *     \left( I + F \frac{\Delta t}{2} \right)^T 
 * \f]
 * \f[
 *     P(t^+_{i-1})   \approx 
 *     \left( I - F \frac{\Delta t}{2} \right)
 *     \left( P(t^-_{i-\frac{1}{2}})  - Q(t_{i-\frac{1}{2}}) \; \frac{\Delta t}{2} \right)
 *     \left( I - F \frac{\Delta t}{2} \right)^T 
 * \f]
 * Subtracting  and ignoring higher order terms yields 
 * \f[ \boxed{
 *     P(t^-_i)  - P(t^+_{i-1})  
 *        =
 *     \left(F P(t^-_{i-\frac{1}{2}}) + P(t^-_{i-\frac{1}{2}}) F^T + Q(t_{i-\frac{1}{2}}) \right) \; \Delta t
 *        =
 *     \left(F P(t^-_{i-\frac{1}{2}}) + (F P(t^-_{i-\frac{1}{2}}) )^T + Q(t_{i-\frac{1}{2}}) \right) \; \Delta t
 * }
 * \f]
 */
class KalmanFilter 
{

public:
  /** constructor */
  KalmanFilter(EquationSystems *);

  /** destructor */
  ~KalmanFilter();

  /** close structured grid infrastructure */
  PetscErrorCode FinalizeDA();

  void DebugOn(); // enable debugging

  // echo data
  virtual void printSelf(std::ostream& os=std::cout) ; 

  //setup infrastructure
  PetscInt CreateIdentityMeasurementMap(PetscInt);
  PetscInt CreateMeasurementMapFromImaging(char*,PetscInt);
  PetscErrorCode Setup(PetscInt);
  PetscErrorCode CreateROIAverageMeasurementMap(PetscInt,PetscInt,PetscInt,
                                                PetscInt,PetscInt,PetscInt,PetscInt, 
                                                PetscScalar,PetscScalar,PetscScalar,
                                                PetscScalar,PetscScalar,PetscScalar);
  PetscErrorCode CreateROINodeSetMeasurementMap(PetscInt);

  //extract image data into mx1 petsc Vec & mxm petsc Mat
  PetscErrorCode ExtractMeasurementData(Vec,char*);

  // extract variance data for plotting
  PetscErrorCode ExtractVarianceForPlotting(char*);

  // extract covariance data for plotting
  PetscErrorCode ExtractCoVarianceForPlotting(char*,PetscInt);

  // debugging routine to show effect of measurement matrix on 
  // uniform roi vector projecting back to full state
  PetscErrorCode MeasurementMatrixTranspose(char*,PetscScalar);
  PetscErrorCode ProjectMeasurementMatrix(Vec,Vec *);

  /* Evaluate initial conditions for state and covariance */
  virtual PetscErrorCode InitializeStateCov() = 0;

  /** @ref KalmanPrediction "Prediction of state"  */
  virtual PetscErrorCode StatePredict( int );

  /** @ref KalmanPrediction "Prediction of state"  */
  virtual PetscErrorCode CovariancePredict();

  /**
   * Prediction of state 
   * 
   * Begin Template Classes... Note. Am using workaround #1  using the this
   * pointer
   * 
   * http://www.parashift.com/c++-faq-lite/templates.html#faq-35.19
   * 
   * [35.19] Why am I getting errors when my template-derived-class uses a
   * member it inherits from its template-base-class?
   * 
   * Perhaps surprisingly, the following code is not valid C++, even though some
   * compilers accept it:
   * 
   * \verbatim 
   *  template<typename T>
   *  class B {
   *  public:
   *    void f() { }  <- member of class B<T>
   *  };
   *  
   *  template<typename T>
   *  class D : public B<T> {
   *  public:
   *    void g()
   *    {
   *      f();  <- bad (even though some compilers erroneously (temporarily?) accept
   * it)
   *    }
   *  };
   * \endverbatim 
   * 
   * This might hurt your head; better if you sit down.
   * 
   * Within D<T>::g(), the name f does not depend on template parameter T, so f
   * is known as a nondependent name. On the other hand, B<T> is dependent on
   * template parameter T so B<T> is called a dependent name.
   * 
   * Here's the rule: the compiler does not look in dependent base classes
   * (like B<T>) when looking up nondependent names (like f).
   * 
   * This doesn't mean that inheritance doesn't work. Class D<int> is still
   * derived from class B<int>, the compiler still lets you implicitly do the
   * is-a conversions (e.g., D<int>* to B<int>*), dynamic binding still works
   * when virtual functions are invoked, etc. But there is an issue about how
   * names are looked up.
   * 
   * Workarounds:
   *     - Change the call from f() to this->f(). Since this is always
   *       implicitly dependent in a template, this->f is dependent and the
   *       lookup is therefore deferred until the template is actually
   *       instantiated, at which point all base classes are considered.  
   *     - Insert using B<T>::f; just prior to calling f().  
   *     - Change the call from f() to B<T>::f(). Note however that this might
   *       not give you what you want if f() is virtual, since it inhibits the
   *       virtual dispatch mechanism.
   */
  PetscErrorCode StateUpdate(char *);

  // solve for the auxillary measurement vector
  virtual PetscErrorCode MeasurementSolve(char*,Vec) = 0 ; 

  /**
   * @ref KalmanUpdate "update of covariance from measurement"  
   * 
   * MR thermal images measure the state variable directly and the 
   * measurement matrix is taken as \f$H=I\f$. The  corresponding
   * @ref KalmanUpdate  "Kalman update" reduces to 
   * \f[
   * \boxed{
   * \begin{split}
   * K(t_i)   & = P(t_i^-)[P(t_i^-) + R(t_i)]^{-1} 
   * \hat{x}(t_i^+)  & = \hat{x}(t_i^-)  + K(t_i)[z_i - \hat{x}(t_i^-)] 
   * P(t_i^+) & = P(t_i^-)  - K(t_i)P(t_i^-)                       
   * \end{split} 
   * }
   * \f]
   *  @todo {KSP MatMatSolve is perfectly parallel
   *         implement OpenMP scheduling for the group of small solves}
   */
  PetscErrorCode CovarianceUpdate(); 

  // invert the covariance sum
  virtual PetscErrorCode InvertSumCov(Mat&) = 0; 

  // structured grid infrastructure
  Vec globalVec; 

  // control over direct solver
  char solvertype[PETSC_MAX_PATH_LEN],
       ordertype[PETSC_MAX_PATH_LEN];
 
protected:
  Mat           CovP,   ///< covariance of the state vector estimate
                CovR,   ///< measurement noise covariance 
                CovQ,   ///< model error covariance 
                eyeState,   ///< template nxn identity matrix
                eyeMeasure, ///< template mxm identity matrix
                MeasurementMat,  ///< measurement matrix
                SumCov, ///< data structure used to hold P + R and factor
                tmpMatDenState; ///< dense matrix used as temporary storage in state update
  // m_CovJacobian p^i = m_CovRHSMat p^{i-1} + m_MassMatrix q
  Mat         m_CovJacobian, m_CovRHSMat, m_MassMatrix;
  PetscScalar modelcov,  ///< model covariance
              statecov;  ///< initial state covariance
  PetscScalar m_fill,    ///< fill for matrix multiply 
              m_zerotol; ///< tolerance for zero on coversion from dense to sparse

  // member functions to initialize covariance
  PetscErrorCode InitializeSparseUncorrCov(Mat&,Mat&,PetscScalar);
  PetscErrorCode InitializeDenseUncorrCov( Mat&,Mat&,PetscScalar);

  // member functions to convert whatever is store in the dense temporary
  // storage to the sparse format
  PetscErrorCode GetSparseMatrixFromDense(Mat &,Mat &);

  // data structures for plotting
  Vec          imageVec, 
               measurementVec, ///< mx1 measurement data
               covMaxVec;      ///< mx1 measurement data
  VecScatter   gather;
  PetscInt maxiter; // max ksp iterations

  // store pointer to equation systems
  EquationSystems *m_eqnSystems; 

private:
  PetscErrorCode MatMultipleRHSSolve(Mat,Mat,Mat);

  // petsc profile logging tools
  PetscInt  kfLogSetup       ,
            kfLogStatePredict,
            kfLogCovarPredict,
            kfLogStateUpdate ,
            kfLogInvSumCov   ,
            kfLogCovarUpdate ;
            
};

/**
 * Sparse Data Structures in PETSc are used for the Kalman Algorithm Matrices.
 * A @ref ForwardEulerKalman "Forward Euler" time stepping scheme is used for
 * the covariance propagation.  A @ref CrankNicolsonKalman "Crank Nicolson"
 * time stepping scheme is used for the state propagation. 
 *
 * @tparam MathematicalModel The template parameter is assumed to be an
 * interface to the \ref MathematicalModelBaseClass Mathematical Model Base
 * class. The MathematicalModel defines the underlying PDE and assembles the
 * system dynamics matrix and load vector.
 *
 */
class SparseKalmanFilter : public KalmanFilter 
{

public:
  // constructor
  SparseKalmanFilter(EquationSystems *);

  /* Evaluate initial conditions for state and covariance */
  virtual PetscErrorCode InitializeStateCov();

  // invert the covariance sum
  virtual PetscErrorCode InvertSumCov(Mat&) ; 

  // solve for the auxillary measurement vector
  virtual PetscErrorCode MeasurementSolve(char*,Vec) ; 
private:
  PetscInt m_maxSpread ; ///< control sparsity of state cov matrix
};
/**
 * Dense Data Structures in PETSc are used for the Kalman Algorithm Matrices.
 * A @ref ForwardEulerKalman "Forward Euler" time stepping scheme is used for
 * the covariance propagation.  A @ref CrankNicolsonKalman "Crank Nicolson"
 * time stepping scheme is used for the state propagation. 
 */
class DenseKalmanFilter : public KalmanFilter 
{

public:
  // constructor
  DenseKalmanFilter(EquationSystems *);

  /**
   *  Evaluate initial conditions for state and covariance */
  virtual PetscErrorCode InitializeStateCov();

  // solve for the auxillary measurement vector
  // factored matrix is created here and destroyed in InvertSumCov
  virtual PetscErrorCode MeasurementSolve(char*,Vec) ; 

  // invert the covariance sum
  // factored matrix is created in MeasurementSolve and destroyed here 
  virtual PetscErrorCode InvertSumCov(Mat&) ; 

private:
  Mat SumCovFact;  // matrix factorization managed by this class
};

/**
 * Dense linear algebra becomes computationally expensive very quickly. 
 * To achieve real time results, an approximation 
 * of eqn \eqn{ForwardEulerCovariancePrediction}  that ignores the 
 * mechanism of diffusing the uncertainty is used
 * \f[
 * P(t^-_i) = 
 *     P(t^+_{i-1}) 
 *      + 
 *     \left(F P(t^+_{i-1}) + P(t^+_{i-1}) F^T + Q(t_{i-1}) \right) \; \Delta t
 *        \approx 
 *     P(t^+_{i-1}) + Q(t_{i-1}) \; \Delta t
 * \f]
 * Using this approximation and beginning with a diagonal, covariance matrix
 * with diagonal modeling error, Q, the state covariance, \f$ P(t) \f$, can maintain
 * a diagonal uncorrelated structure throughout the time history.
 * \f[ \begin{split}
 * K(t_i)   & = P(t_i^-)[P(t_i^-) + R(t_i)]^{-1} 
 *         = 
 *     \begin{Bmatrix}
 *         \frac{p^i_{11}}{p^i_{11} +\gamma^i_1} & 0 & 0 & ... \\
 *         0 & \frac{p^i_{22}}{p^i_{22} +\gamma^i_2} & 0 & ... \\
 *           &   &  . &   \\
 *           &   &    & . \\
 *         0 & 0 & ... & 0 &\frac{p^i_{nn}}{p^i_{nn} +\gamma^i_{nn}} \\
 *     \end{Bmatrix}
 *         = 
 *          \text{diag}\left( \frac{p^i_{11}}{p^i_{11} +\gamma^i_1},       
 *                \frac{p^i_{22}}{p^i_{22} +\gamma^i_2}, ...  \right)  \\
 *   \\
 * \hat{x}(t_i^+)  & = \hat{x}(t_i^-)  + K(t_i)[z_i - \hat{x}(t_i^-)] 
 *                   = 
 *          \text{diag}\left( 
 *                \frac{\gamma^i_1}{p^i_{11} +\gamma^i_1},       
 *                \frac{\gamma^i_2}{p^i_{22} +\gamma^i_2}, ...  \right)  
 *                     \hat{x}(t_i^-)  + 
 *          \text{diag}\left( 
 *                \frac{p^i_{11}}{p^i_{11} +\gamma^i_1},       
 *                \frac{p^i_{22}}{p^i_{22} +\gamma^i_2}, ...  \right)  
 *                      z_i  \\
 * P(t_i^+) & = P(t_i^-)  - K(t_i)P(t_i^-)                       
 *            = 
 *          \text{diag}\left( 
 *              \frac{\gamma^i_1 p^i_{11}}{p^i_{11} +\gamma^i_1},       
 *              \frac{\gamma^i_2 p^i_{22}}{p^i_{22} +\gamma^i_2}, ...  \right) 
 * \end{split} \f]
 * 
 */
class UncorrelatedKalmanFilter : public DenseKalmanFilter 
{

public:
  // constructor
  UncorrelatedKalmanFilter(EquationSystems *);

  // only update diagonal
  virtual PetscErrorCode CovariancePredict(); 

private:
  // none
};
#endif
