
//constructor
template< typename PDEModelSolver >
onlyMRTIQOI< PDEModelSolver >::
onlyMRTIQOI(AppSolve *user,GetPot &controlfile, PetscInt iii):
qoiBaseClass(user,controlfile,iii) 
{
   m_pdeSolver= dynamic_cast< PDEModelSolver*>( user->pdeSolver() ) ;
   m_sigma = 1.0e9; 
};
// mrti data only
#undef __FUNCT__
#define __FUNCT__ "onlyMRTIQOI::solveqoi"
template< typename PDEModelSolver >
PetscErrorCode onlyMRTIQOI< PDEModelSolver >::
solveqoi(AppSolve *user, GetPot &)
{
  PetscErrorCode info; 
  PetscFunctionBegin; 
  /* load the data structures for the QOI
     PriorLoadQoiData is a function pointer to  */
  info = this->PriorLoadQoiData(user);CHKERRQ(info);

  // plot data
  this->PlotStoredData(user);
    
  PetscFunctionReturn(0);
}

// mrti data only
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::solveqoi"
template< typename PDEModelSolver >
PetscErrorCode noOptQOI< PDEModelSolver >::
solveqoi(AppSolve *user, GetPot &)
{
  PetscErrorCode info; 
  PetscFunctionBegin; 
  /* load the data structures for the QOI
     PriorLoadQoiData is a function pointer to  */
  info = this->PriorLoadQoiData(user);CHKERRQ(info);

  // plot data
  this->PlotInterleaveCompute(user);
    
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "onlyMRTIQOI::SetupAuxSystems"
template< typename PDEModelSolver >
PetscErrorCode onlyMRTIQOI< PDEModelSolver >::
SetupAuxSystems(EquationSystems &EqnSystems)
{
  PetscFunctionBegin; 
  AppSolve *user = AppSolve::getInstance();

  // systems for storing MRTI data and verification problems
  TransientFEMSystem & ideal_system =
   EqnSystems.add_system<TransientFEMSystem>("IdealSystem");
  ideal_system.add_variable ("u0*", FIRST);
  // reserve some space for storage
  ideal_system.ReserveVector(user->get_max_steps()+2);
  ideal_system.ReserveOldVector(1);

  TransientFEMSystem & ideal_uncertainty_system =
   EqnSystems.add_system<TransientFEMSystem>("IdealUncertaintySystem");
  // reserve some space for storage
  ideal_uncertainty_system.add_variable ("du0*", FIRST);
  ideal_uncertainty_system.ReserveVector(user->get_max_steps()+2);
  ideal_uncertainty_system.ReserveOldVector(1);

  // set for compatibility w/ FEMSystem
  ideal_system.time_solver = 
      AutoPtr<TimeSolver>(new SteadySolver(ideal_system));
  ideal_uncertainty_system.time_solver = 
      AutoPtr<TimeSolver>(new SteadySolver(ideal_uncertainty_system));

  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::SetupAuxSystems"
template< typename PDEModelSolver >
PetscErrorCode noOptQOI< PDEModelSolver >::
SetupAuxSystems(EquationSystems &EqnSystems)
{
  PetscErrorCode info; 
  PetscFunctionBegin; 
  AppSolve *user = AppSolve::getInstance();

  // systems for storing MRTI data and verification problems
  TransientFEMSystem & ideal_system =
   EqnSystems.add_system<TransientFEMSystem>("IdealSystem");
  ideal_system.add_variable ("u0*", FIRST);
  // reserve some space for storage
  ideal_system.ReserveVector(user->get_max_steps()+2);
  ideal_system.ReserveOldVector(1);

  // set for compatibility w/ FEMSystem
  ideal_system.time_solver = 
      AutoPtr<TimeSolver>(new SteadySolver(ideal_system));
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- 
   PriorLoadQoiData - compute a fake ideal temperature field
                            used for code verification

   Input Parameters:
   ctx  - user-defined context 
   
   Required Routines: inittempfield must be called before this routine
*/
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::PriorLoadQoiData"
template< typename PDEModelSolver >
PetscErrorCode noOptQOI< PDEModelSolver >::
PriorLoadQoiData(AppSolve *user)
{
  EquationSystems    &EqnSystem = user->get_equation_systems();

  // used for error handling
  PetscFunctionBegin; 

  this->echoOptWindow(std::cout,user->qoiOptimizer()->Nsteplo(),user->qoiOptimizer()->Nstephi());

  if(user->get_num_exact() )
   {
    // get the system for the ideal problem.
    TransientFEMSystem & ideal_system =
     EqnSystem.get_system<TransientFEMSystem>("IdealSystem");
    TransientFEMSystem & state_system = 
     EqnSystem.get_system<TransientFEMSystem> ("StateSystem");

    // number of time steps
    PetscInt nsteps = user->qoiOptimizer()->Nstephi() - user->qoiOptimizer()->Nsteplo() + 1 ;
    // solution buffer for temp storage..
    std::vector< NumericVector<Number>* > solution_buffer(nsteps);

    for(AppSolve::ISTEP = user->qoiOptimizer()->Nsteplo() ; 
        AppSolve::ISTEP <= user->getMaxFEMTimeStep() ; AppSolve::ISTEP++)
     {
       // project the solution 
       // ideal_system.project_solution(pennes_forward_exact,NULL,
       //                                     user->_eqnSystem->parameters);
       state_system.user_initialization();

       // store the solution history
       *ideal_system.vector_solution.at(AppSolve::ISTEP)=
                                   *state_system.current_local_solution;
     }
    // reset IC on state for solve
    AppSolve::ISTEP = user->qoiOptimizer()->Nsteplo() ; 
    state_system.user_initialization();
   }

  PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------- 
   not needed as dirichlet data should take care of this.
   use for verification problems only...
*/
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::DirichletDomains"
template< typename PDEModelSolver >
void spaceTimeFDVerification < PDEModelSolver >::
DirichletDomains(AppSolve *user, Vec yOut, const PetscScalar bcValue )
{
 PetscFunctionBegin;
 PetscErrorCode info; // used to check for functions returning nonzeros 

 std::vector < PetscInt > domainArray; // build id's of all dirichlet domains
 for( this->IterParam=this->Parameters.begin();
this->IterParam!=this->Parameters.end(); this->IterParam++)
  {
   int iii = distance(this->Parameters.begin(),this->IterParam);

   // loop over variables and domains dirichlet domains should be stored
   for( unsigned int i_var = 0 ; i_var < this->m_pdeSolver->n_vars() ; i_var++)
    for ( unsigned int i_domain = 0 ; i_domain < user->get_num_elem_blk() ; i_domain++ ) 
      if( this->m_pdeSolver->dirichletData(i_var,i_domain) )  // TODO add i_var dependence
        domainArray.push_back( iii * user->get_num_elem_blk() + i_domain ) ;
  }

 // set the corresponding entries in the yOut
 for ( PetscInt Ii =0 ; Ii < domainArray.size() ; Ii++)
   info = VecSetValue (yOut, domainArray[Ii] , bcValue, INSERT_VALUES);

 //Communicate Assembly of the yOut
 info = VecAssemblyBegin(yOut); CHKERRV(info);
 info = VecAssemblyEnd(yOut); CHKERRV(info)

 PetscFunctionReturnVoid();
}
#undef __FUNCT__
#define __FUNCT__ "spaceTimeQOI::SolveTAO"
template< typename PDEModelSolver >
PetscErrorCode spaceTimeQOI< PDEModelSolver >::
SolveTAO(AppSolve *user) // run TAO solvers
{
  // used for error handling
  PetscFunctionBegin; 
  PetscErrorCode info;

  /* Optimize the temperature field using FormObjective FormGradient
     Tao call various optimization routines: - TaoSolve_BLMVM
                                             - TaoSolve_BNLS 
                                             - TaoSolve_LMVM 
                                             - TaoSolve_NelderMead
                                             - TaoSolve_FD 
                                             - etc... */
  std::string methodType = this->GetTaoSolverMethod();
  PetscPrintf(PETSC_COMM_WORLD,"\n using TAO Solver:%s\n\n",methodType.c_str());

  /* need to initialize data structures for TaoSolve_FD*/
  PetscTruth  flg,
              testGradient=PETSC_FALSE,// Gradient flag
              testHessian =PETSC_FALSE;// Hessian flag
  if(methodType.find("tao_fd_test")!=std::string::npos)
   {
/* Finite Difference Practical considerations
http://en.wikipedia.org/wiki/Numerical_differentiation

An important consideration in practice when the function is approximated using
floating point arithmetic is how small a value of h to choose. If chosen too
small, the subtraction will yield a large rounding error and in fact all the
finite difference formulae are ill-conditioned[2] and due to cancellation will
produce a value of zero if h is small enough[3]. If too large, the calculation
of the slope of the secant line will be more accurate, but the estimate of the
slope of the tangent by using the secant could be worse.

As discussed in Chapter 5.7 of Numerical Recipes in C
(http://www.nrbook.com/a/bookcpdf/c5-7.pdf), a suitable choice for h is
sqrt(eps) * x where the machine epsilon eps is typically of the order 2.2e-16.
Another important consideration is to make sure that h and x+h are representable
in floating point precision so that the difference between x+h and x is exactly
h. This can be accomplished by placing their values into and out of memory as
follows: h = sqrt(eps) * x, temp = x + h and h = temp - x. It may be necessary
to declare temp as a volatile variable to ensure the steps are not undone by
compiler optimization.
*/
    info=PetscOptionsGetTruth(PETSC_NULL,"-tao_test_gradient",&testGradient,&flg); 
    CHKERRQ(info);
    info=PetscOptionsGetTruth(PETSC_NULL,"-tao_test_hessian",&testHessian,&flg);
    CHKERRQ(info);

    Vec            solution;
    info = TaoAppGetSolutionVec(this->taoapp,&solution); CHKERRQ(info);
    if( testGradient || testHessian )
     { // data structures should have the perturbed solution
       info=TaoAppSetGradientRoutine(this->taoapp,NULL,&user);CHKERRQ(info);
       info=TaoAppSetObjectiveAndGradientRoutine(this->taoapp,FormObjectiveAndGradient
                                                 ,&user);CHKERRQ(info);
       info=TaoAppSetHessianRoutine( this->taoapp,FormHessian,&user);CHKERRQ(info);
       info=TaoAppSetHessianMat(this->taoapp,this->Hessian,this->Hessian);
     }
    if( testHessian )
     { // initialize data structures for hessian computation
       Vec            gradient;
       info = VecDuplicate(solution,&gradient); CHKERRQ(info);
       info = TaoAppComputeGradient(this->taoapp,solution,gradient); CHKERRQ(info);
       PetscPrintf(PETSC_COMM_WORLD, "Adjoint Gradient... \n"); 
       info = VecView(gradient,0); CHKERRQ(info);
       info = VecDestroy(gradient); CHKERRQ(info);
     }
   }

  // print convergence tolerances before solve
  info = TaoView(this->tao); CHKERRQ(info);

  // solve
  info = TaoSolveApplication(this->taoapp,this->tao); CHKERRQ(info);

  /* Get termination information */
  TaoTerminateReason reason;
  info = TaoGetTerminationReason(this->tao,&reason);CHKERRQ(info);
  if (reason <= 0) PetscPrintf(PETSC_COMM_WORLD,
                   "Try a different TAO method, \n \
                    adjust some parameters, \n \
                    or check the function evaluation routines \n ");
  switch(reason){
   case 0:
      PetscPrintf(PETSC_COMM_WORLD,
      "continue iterating...\n"); break;
   case 2:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated: TAO_CONVERGED_ATOL\n"); break;
   case 3:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated: TAO_CONVERGED_RTOL\n"); break;
   case 4:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated: TAO_TRTOL\n"); break;
   case 5:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated: TAO_MINF\n"); break;
   case 6:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated: TAO_CONVERGED_USER\n"); break;
   case -2:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated: TAO_DIVERGED_MAXITS\n"); break;
   case -4:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated: TAO_DIVERGED_NAN\n"); break;
   case -5:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated: TAO_DIVERGED_MAXFCN\n"); break;
   case -6:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated : TAO_DIVERGED_LS_FAILURE\n"); break;
      /* if line search fails, i think that the wrong paramters values are
         stored in the hp3d data structures. The parameters of the last
         fcn evaluation that caused the line search to fail will be 
         in the hp3d data structures. restore the original parameter values
         of to the hp3d the data structures */
   case -7:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated: TAO_DIVERGED_TR_REDUCTION\n");break;
   case -8:
      PetscPrintf(PETSC_COMM_WORLD,
      "Reason TAO terminated: TAO_DIVERGED_USER\n"); break;
  }


  // ensure final solution is in data structures
  info =  VecScatterBegin( this->GLOCSCAT_CNTL, this->QCONTRL, this->CNTL_LOC,
                           INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(info);
  info =  VecScatterEnd(   this->GLOCSCAT_CNTL, this->QCONTRL, this->CNTL_LOC,
                           INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(info);
  info = this->PutLocalCntrlVars(this->CNTL_LOC); 


  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "spaceTimeQOI::DestroyTAO"
template< typename PDEModelSolver >
PetscErrorCode spaceTimeQOI< PDEModelSolver >::
DestroyTAO()
{
  // used for error handling
  PetscFunctionBegin; 
  PetscErrorCode info;

  /* Free TAO data structures */
  info = TaoDestroy(this->tao);       CHKERRQ(info);
  info = TaoAppDestroy(this->taoapp); CHKERRQ(info);
  /* looks like 
   info = MatDestroy(this->Hessian); CHKERRQ(info);
     is called by TaoDestroy
   */

  PetscFunctionReturn(0);
}
//FIXME will need to be modified if # of optimization variables
// change at each optimization step (time varying variables)
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::DestroyPetsc"
template< typename PDEModelSolver >
PetscErrorCode noOptQOI< PDEModelSolver >::
DestroyPetsc()
{
  // used for error handling
  PetscFunctionBegin; 

  /* Free PETSc data structures */
  if(this->QCONTRL)      VecDestroy(       this->QCONTRL);
  if(this->CNTL_LOC)     VecDestroy(       this->CNTL_LOC);   
  if(this->GLOCSCAT_CNTL)VecScatterDestroy(this->GLOCSCAT_CNTL); 

  PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------- */
template< typename PDEModelSolver >
void  noOptQOI< PDEModelSolver >::
GetInitialSolutionVec(AppSolve *user)
{
 PetscFunctionBegin;
 PetscErrorCode info; 
  if( user->Restart() && !user->IDOPT )
    {// read control vars from disk
     this->m_pdeSolver->GetDataFromDisk( user->get_mesh(),
                                          AppSolve::restartFile );

     info = this->GetLocalCntrlVars(this->CNTL_LOC); CHKERRV(info); 

     // scatter to global vec
     info = VecScatterBegin(this->GLOCSCAT_CNTL,this->CNTL_LOC,this->QCONTRL,
                            INSERT_VALUES, SCATTER_FORWARD); CHKERRV(info);
     info = VecScatterEnd(  this->GLOCSCAT_CNTL,this->CNTL_LOC,this->QCONTRL,
                            INSERT_VALUES, SCATTER_FORWARD); CHKERRV(info);
    }
  else
    {// initialize control variables
     this->GetCntrlVars( user->get_mesh() ,this->QCONTRL);

     // ensure field variables are initialized 
     // to the proper initial guess and for initial plot verification 
     info = VecScatterBegin(this->GLOCSCAT_CNTL,this->QCONTRL,this->CNTL_LOC,
                            INSERT_VALUES, SCATTER_REVERSE); CHKERRV(info);
     info = VecScatterEnd(  this->GLOCSCAT_CNTL,this->QCONTRL,this->CNTL_LOC,
                            INSERT_VALUES, SCATTER_REVERSE); CHKERRV(info);
     info = this->PutLocalCntrlVars(this->CNTL_LOC); CHKERRV(info); 
    }
 PetscFunctionReturnVoid();
} 
/* ------ accumulate the load for the sensitivity solve pde deriv ------*/
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::accumulateSensitivityLoad"
template< typename PDEModelSolver >
void spaceTimeQOI< PDEModelSolver >::
accumulateSensitivityLoad(AppSolve *user,
                          const QGauss &qrule,
                          std::vector<PetscScalar> &elemParameter,
                          const unsigned int &field_id, 
                          const std::vector<Real>& JxW,
                          const std::vector<Point>& qpoint,
                          std::vector< DenseSubVector<Number> > &Fi,
                          TransientFEMSystem &state_system)
{
 PetscFunctionBegin; 

 // floor fem time step by ISTEPS_PER_IDEAL  to get power step
 PetscInt idpow= user->get_id_power();

 // loop over params
 for( this->IterParam=this->Parameters.begin();
this->IterParam!=this->Parameters.end(); this->IterParam++)
  {
   int iii = distance(this->Parameters.begin(),this->IterParam);
   optimizationParameter* optParam = *(this->IterParam);
   for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
     // evaluate Solution used in the dpde_dm functer
     this->m_pdeSolver->evaluateSolutionForState(qp,state_system);
     /*   evaluate the sensitivity load at the guass points */
     // loop over shape functions
     for (unsigned int Ii=0; Ii< this->m_pdeSolver->n_vars(); Ii++)
      for (unsigned int Jj=0; Jj< this->m_pdeSolver->n_u_dofs(Ii); Jj++)
       {
         // adjoint variable used in dpde_dm functer and used as storage space
         // to pass in the shape function values FIXME - altering the interface
         // to pass adjoint variable in gradient accumulation and shape
         // functions in sensitivity solve will be a lot of work
         this->m_pdeSolver->evaluateAdjointSolutionAsShapeFunction(Jj,qp);

         // accumulate
         Fi[Ii](Jj) += JxW[qp] *
                     #if defined(PETSC_USE_DEBUG)
                        elemParameter.at(iii) 
                     #else
                        elemParameter[iii] 
                     #endif
                             * (
              - this->m_pdeSolver->dpde_dm(optParam,field_id,qpoint[qp])
                               );
       }
     }
  }

 PetscFunctionReturnVoid(); 
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::accumulateGradient"
template< typename PDEModelSolver >
void spaceTimeQOI< PDEModelSolver >::
accumulateGradient(AppSolve *user, const QGauss &qrule,
                   const unsigned int &field_id, 
                   const std::vector<Real>& JxW,
                   const std::vector<Point>& q_point,
                   DenseVector<Number> &Grad)
{
 PetscFunctionBegin;

 // Get a reference to the NonlinearImplicitSystem we are solving
 EquationSystems    &EqnSystem = user->get_equation_systems();
 TransientFEMSystem& state_system = 
   EqnSystem.get_system<TransientFEMSystem>("StateSystem");

 // Get a reference to the LinearImplicitSystem for adjoint problem
 TransientFEMSystem& adjoint_system = 
   EqnSystem.get_system<TransientFEMSystem>("AdjointSystem");

 for( this->IterParam=this->Parameters.begin();
this->IterParam!=this->Parameters.end(); this->IterParam++)
  {
   int iii = distance(this->Parameters.begin(),this->IterParam);
   optimizationParameter* optParam = *(this->IterParam);
   // build the element contributions
   for (unsigned int qp=0; qp<qrule.n_points(); qp++)
     {

      this->m_pdeSolver->evaluateSolutionForGradient(qp,AppSolve::ISTEP,
                                             state_system,adjoint_system);

      // element gradient contribution
      Grad(iii) += JxW[qp]*(
             CALL_MEMBER_FN(this,optParam->dqoi_dm)(field_id,q_point[qp])
                                          -
             this->m_pdeSolver->dpde_dm(optParam,field_id,q_point[qp])
        	           );

     } // end loop over interior guass points
  } // end loop over parameters

 PetscFunctionReturnVoid();
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::accumulateHessian"
template< typename PDEModelSolver >
void spaceTimeQOI< PDEModelSolver >::
accumulateHessian(AppSolve *user,
                                            const QGauss &qrule,
                                 std::vector<PetscScalar> &elemParameter,
                                            const unsigned int &field_id, 
                                            const std::vector<Real>& JxW,
                                            const std::vector<Point>& q_point,
                                            DenseVector<Number> &HessVecProd)
{
 PetscErrorCode info; /* used to check for functions returning nonzeros */
 PetscFunctionBegin;

 PetscLogEventBegin(AppSolve::logevents[11],0,0,0,0); // elemjac assble

 // Get a reference to the NonlinearImplicitSystem we are solving
 EquationSystems    &EqnSystem = user->get_equation_systems();
 TransientFEMSystem& state_system = 
   EqnSystem.get_system<TransientFEMSystem>("StateSystem");

 // Get a reference to the LinearImplicitSystem for adjoint problem
 TransientFEMSystem& adjoint_system = 
   EqnSystem.get_system<TransientFEMSystem>("AdjointSystem");

 // Get a reference to the LinearImplicitSystem for sensitivity problem
 TransientFEMSystem& sensitivity_system = 
   EqnSystem.get_system<TransientFEMSystem>("SensitivitySystem");

 // Get a reference to the ExplicitSystem for ideal data
 TransientFEMSystem & ideal_system =
   EqnSystem.get_system<TransientFEMSystem>("IdealSystem");

 // get the ideal system uncertainty
 TransientFEMSystem & ideal_uncertainty_system =
   EqnSystem.get_system<TransientFEMSystem>("IdealUncertaintySystem");

 std::vector<optimizationParameter*>::iterator ColumnParam;
 for( this->IterParam=this->Parameters.begin();
this->IterParam!=this->Parameters.end(); this->IterParam++)
  {
   int iii = distance(this->Parameters.begin(),this->IterParam);
   optimizationParameter* optParam = *(this->IterParam);
   // build the element contributions at each gauss point
   for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
     // evaluate Solution FIXME
     this->m_pdeSolver->evaluateIdealSolution(qp,AppSolve::ISTEP,user->qoiOptimizer()->Nstephi(),
                                      state_system,
                                      ideal_system,ideal_uncertainty_system);
     // evaluate state and adnoint at the guass point
     this->m_pdeSolver->evaluateSolutionForGradient(qp,AppSolve::ISTEP,
                                            state_system,adjoint_system);
     // evaluate sensitivities at the guass point
     this->m_pdeSolver->evaluateSolutionForSensitivity(qp,AppSolve::ISTEP,
                                           sensitivity_system );

     // accumulate element hessian sensitivities 
     PetscScalar hessianTmp = 0.0;

     // mixed partials wrt to solution and control params
     hessianTmp += 
      this->m_pdeSolver->d2qoi_du_dm(optParam,field_id,q_point[qp])
                                       -
      this->m_pdeSolver->d2pde_du_dm(optParam,field_id,q_point[qp]);
                                         // * \sum_jjj du_kmhalf_jjj  q_jjj
     // second derivative wrt to control params
     for( ColumnParam=this->Parameters.begin(); 
         ColumnParam!=this->Parameters.end(); ColumnParam++)
      {
       int jjj = distance(this->Parameters.begin(),ColumnParam);
       hessianTmp += 
         #if defined(PETSC_USE_DEBUG)
            elemParameter.at(jjj) 
         #else
            elemParameter[jjj] 
         #endif
            * (
         this->d2qoi_dmi_dmj(field_id,q_point[qp]) //FIXME-matrix of functers
                                          -
         this->m_pdeSolver->d2pde_dmi_dmj(iii,jjj,field_id,q_point[qp])
              );
      }

     // evaluate the adjoint with the solution stored as the current
     // solution to the adjoint problem w/ sensitivity the rhs specific to the
     // matrix-vector product
     // TODO - dangerous to re-use existing p__now,grad_p__now storage space?
     this->m_pdeSolver->evaluateSolutionForAdjoint( qp, adjoint_system );
     // sensitivity adjoint contribution
     hessianTmp -= this->m_pdeSolver->dpde_dm(optParam,field_id,q_point[qp]) ;

     // multiply by quadrature weight
     HessVecProd(iii) += JxW[qp]* hessianTmp;
    }  // end loop over interior guass points
  }
 PetscLogEventEnd(  AppSolve::logevents[11],0,0,0,0); // elemjac assble

 PetscFunctionReturnVoid( );
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::accumulateSensitivityGradient"
template< typename PDEModelSolver >
PetscScalar spaceTimeQOI< PDEModelSolver >::
accumulateSensitivityGradient(AppSolve *user,
               const QGauss &qrule,
               const unsigned int &field_id, 
               const std::vector<Real>& JxW,
               const std::vector<Point>& q_point)
{
 Vec tmpVec; // file handle for petsc viewer
 //int            fileHandle;
 //PetscScalar    *avec;
 //off_t         currentPos;
 PetscFunctionBegin;

 PetscLogEventBegin(AppSolve::logevents[11],0,0,0,0); // elemjac assble

 // Get a reference to the NonlinearImplicitSystem we are solving
 EquationSystems    &EqnSystem = user->get_equation_systems();
 TransientFEMSystem& state_system = 
   EqnSystem.get_system<TransientFEMSystem>("StateSystem");

 // Get a reference to the LinearImplicitSystem for adjoint problem
 TransientFEMSystem& adjoint_system = 
   EqnSystem.get_system<TransientFEMSystem>("AdjointSystem");

 // Get a reference to the LinearImplicitSystem for sensitivity problem
 TransientFEMSystem& sensitivity_system = 
   EqnSystem.get_system<TransientFEMSystem>("SensitivitySystem");

 // Get a reference to the ExplicitSystem for ideal data
 TransientFEMSystem & ideal_system =
   EqnSystem.get_system<TransientFEMSystem>("IdealSystem");

 // get the ideal system uncertainty
 TransientFEMSystem & ideal_uncertainty_system =
   EqnSystem.get_system<TransientFEMSystem>("IdealUncertaintySystem");

 // build the element contributions at each gauss point
 PetscScalar gradientEntry = 0.0;
 for (unsigned int qp=0; qp<qrule.n_points(); qp++)
   {
    // evaluate Solution FIXME
    this->m_pdeSolver->evaluateIdealSolution(qp,AppSolve::ISTEP,user->qoiOptimizer()->Nstephi(),
                                     state_system,
                                     ideal_system,ideal_uncertainty_system);
    // evaluate state and adnoint at the guass point
    this->m_pdeSolver->evaluateSolutionForGradient(qp,AppSolve::ISTEP,
                                           state_system,adjoint_system);
    // evaluate sensitivities at the guass point
    this->m_pdeSolver->evaluateSolutionForSensitivity(qp,AppSolve::ISTEP,
                                              sensitivity_system );

    // verify sensitivity gradient against adjoint gradient
    gradientEntry += JxW[qp]* (
                           this->m_pdeSolver->dqoidu_dudm(field_id)
       // +CALL_MEMBER_FN(this->m_pdeSolver,iiiParam->dpde_dm)(field_id,q_point[qp])
                              );
   }  // end loop over interior guass points
 PetscLogEventEnd(  AppSolve::logevents[11],0,0,0,0); // elemjac assble

 PetscFunctionReturn( gradientEntry );
}
/* -------------------------------------------------------------------- 
   Do Not wait for MRI/MRTI DATA!!!!
     solve for each time step and store in data structures. */
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::ForwardSolve"
template< typename PDEModelSolver >
PetscErrorCode noOptQOI< PDEModelSolver >::
ForwardSolve(AppSolve *user) 
{
  EquationSystems    &EqnSystem = user->get_equation_systems();

  // used for error handling
  PetscFunctionBegin; 

  // Get a reference to the Convection-Diffusion system object.
  TransientFEMSystem & system =
    EqnSystem.get_system<TransientFEMSystem>("StateSystem");

  // reset IC
  *system.current_local_solution =
  *system.vector_solution.at( user->qoiOptimizer()->Nsteplo() ) ;

  //time stepping and function evaluation
  for(AppSolve::ISTEP=user->qoiOptimizer()->Nsteplo()+1 ; AppSolve::ISTEP<=user->qoiOptimizer()->Nstephi() ;
                                                           AppSolve::ISTEP++)
    {
       // update the solution vector from the previous time step
       *system.old_vector_solution[0] = *system.current_local_solution;

       /* !!!NOTE!!! Do not need to store the soln in data structures 
                     after the solve b/c this was already done by the 
                     last function evaluation */
       // Perform NonlinearImplicitSystem::solve () with PetscNonlinearSolver. 
       //      compute_residual  
       //      compute_jacobian 
       // This will put a local copy of solution into current_local_solution.
       system.solve();
   
       // also need to store the solution history for the adjoint
       *system.vector_solution.at(AppSolve::ISTEP)=
                                          *system.current_local_solution;
       
    }

  PetscFunctionReturn(0);
}
// subroutine to extract image data to a NumericVector
#undef __FUNCT__
#define __FUNCT__ "onlyMRTIQOI::ExtractImageDataFromBackground"
template< typename PDEModelSolver >
PetscErrorCode onlyMRTIQOI< PDEModelSolver >::
ExtractImageDataFromBackground(
                  InputImageType::Pointer tempImage, // image to extract from
                  PetscInt idMRTIfile,
                  TransientFEMSystem & local_system ) // system to put soln
{
  PetscFunctionBegin;

  EquationSystems    &EqnSystem = this->m_pdeSolver->get_equation_systems();

  // set pointer
  this->Images->interpolator->SetInputImage( tempImage );

  // background correction tempImage should be phase data
  TransientFEMSystem & background_system =
     EqnSystem.get_system<TransientFEMSystem>("BackgroundSystem");
  // put MRTI thermal image into FEM data structures
  background_system.project_solution(GetITKImageData,NULL,
                                     EqnSystem.parameters);
  // store raw data for later
  *background_system.old_vector_solution[0] = 
                                *background_system.current_local_solution;

  // solve Neumann Problem must remove NULL Space
  PetscDiffSolver* backgroundSolver = 
                  libmesh_cast_ptr<PetscDiffSolver*>(
                  & (*background_system.time_solver->diff_solver()) );
  PetscErrorCode info;
  KSP  snesksp;
  MatNullSpace NullSpace;
  info = SNESGetKSP(backgroundSolver->snes(),&snesksp);CHKERRQ(info);
  info = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,
                            PETSC_NULL,&NullSpace); CHKERRQ(info);
  info = KSPSetNullSpace(snesksp,NullSpace);CHKERRQ(info);

  // assemble and solve 
  background_system.solve();

  // clean up
  info = MatNullSpaceDestroy(NullSpace);CHKERRQ(info);

  // save the total difference to be plotted later...
  if(idMRTIfile)
   {// get previous
    *background_system.vector_solution.at(idMRTIfile) =
    *background_system.vector_solution.at(idMRTIfile-1);
    //background_system.vector_solution.at(idMRTIfile)->zero();
    // store the sum difference
    background_system.vector_solution.at(idMRTIfile)->add(
                  this->Images->GetTmapFactor(),
                  *background_system.current_local_solution );
   // regular phase difference
   *local_system.current_local_solution =
   *local_system.vector_solution.at(idMRTIfile-1);
   local_system.current_local_solution->add(
      // TODO  still need appropriate constant
      this->Images->MagneticPotentialDenominator(),
                       *background_system.old_vector_solution[0] );
   }
  else // initial difference should be zero
   {
   background_system.vector_solution.at(idMRTIfile)->zero();
   // set to initial condition
   local_system.current_local_solution->zero(); Point dum;
   AppSolve *user = AppSolve::getInstance();
   local_system.current_local_solution->add( 
    user->pdeSolver()->InitValues(0,0,dum,EqnSystem.parameters));
   }


  PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "dddasQOI::dddasQOI"
template< typename PDEModelSolver >
dddasQOI < PDEModelSolver >::
dddasQOI(AppSolve *user,GetPot &controlfile,PetscInt iii):
         noOptQOI< PDEModelSolver >::noOptQOI(user,controlfile,iii)
{
  /* for calibration can only have as many optimization steps 
                                      as thermal images will allow */
  // error checking. this enforces that all thermal
  // image requests between intercommunicator groups
  // remain within the bounds [0 , MRTI_ntime]. 
  if(this->IDEAL_NTIME+this->NUMIDEAL * (this->NoptSteps-1) 
                                             > AppSolve::MRTI_ntime)
  {
     if(!this->NUMIDEAL){
        std::cout << "optvars.cpp: calibration input error \n"
             << " IDEAL_NTIME+NUMIDEAL*(NoptSteps-1) > MRTI_ntime\n" 
             << " NUMIDEAL = "<< this->NUMIDEAL <<std::endl; abort();
     }
     //FLOOR property of int division ensures 
     //final time bounded by MRTI_ntime
     this->NoptSteps = ( AppSolve::MRTI_ntime - this->IDEAL_NTIME ) / this->NUMIDEAL + 1;
  }
}
#undef __FUNCT__
#define __FUNCT__ "spaceTimeQOI::SetupAuxSystems"
template< typename PDEModelSolver >
PetscErrorCode spaceTimeQOI< PDEModelSolver >::
SetupAuxSystems(EquationSystems &EqnSystems)
{
  PetscErrorCode info; 
  PetscFunctionBegin; 
  AppSolve *user = AppSolve::getInstance();
  GetPot &controlfile = *EqnSystems.parameters.get<GetPot*>("controlfile") ;

  std::cout << "Setting up Control Vars" <<  std::endl << std::flush;
  info = this->SetupCntrl(user,controlfile); CHKERRQ(info);
  PetscLogDouble memoryUsed;
  info = PetscMemoryGetCurrentUsage(&memoryUsed);
  std::cout << "Memory Used " << memoryUsed <<  std::endl << std::flush;

  // setup the adjoint system variables
  std::cout << "Setting up Adjoint" <<  std::endl << std::flush;
  this->SetupAdjoint(user); // abstract pde solver
  info = PetscMemoryGetCurrentUsage(&memoryUsed);
  std::cout << "Memory Used " << memoryUsed <<  std::endl << std::flush;
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- 
   solveqoi - Use Tao to optimize the parameters 

   Input Parameters:
   user  - user-defined context 
   
   Required Routines: All petsc data structures 
                          must be setup to call this routine
*/
#undef __FUNCT__
#define __FUNCT__ "spaceTimeQOI::solveqoi"
template< typename PDEModelSolver >
PetscErrorCode spaceTimeQOI< PDEModelSolver >::
solveqoi(AppSolve *user, GetPot &)
{
  PetscErrorCode info; 
  PetscFunctionBegin; 

  // setup optimization solver
  info = this->SetupTAO(user); CHKERRQ(info);
                           
  /* load the data structures for the QOI
     PriorLoadQoiData is a function pointer to  */
  info = this->PriorLoadQoiData(user);CHKERRQ(info);
    
  /* Run Optimization Solver */
  info = this->SolveTAO(user);CHKERRQ(info);

  // plot data
  if( this->Nstephi()+1 <= user->get_num_MRTI() * user->IstepsPerIdeal())
      this->PlotInterleaveCompute(user);

  // free data structures
  info = this->DestroyTAO(); CHKERRQ(info);

  PetscFunctionReturn(0);
}
/*----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "spaceTimeQOI::SetupAdjoint"
template< typename PDEModelSolver >
void spaceTimeQOI< PDEModelSolver >::
SetupAdjoint(AppSolve *user)
{
 PetscFunctionBegin;

 // initialize the solution variables
 this->m_pdeSolver->SetupAdjoint( this->Parameters.size() );

 // equation system 
 EquationSystems    &EqnSystem = user->get_equation_systems();

 // allocate memory to hold the solution time history need one more 
 // than Nstep for Initial Condition
 TransientFEMSystem & state_system = 
   EqnSystem.get_system<TransientFEMSystem> ("StateSystem");
 state_system.ReserveVector(user->get_max_steps()+2);

 /* the following system will hold:
          - the exact solution for the verification problems
          - the thermal imaging for the calibration problem
          - the ideal field for the optimal control problem  */
 TransientFEMSystem & ideal_system =
   EqnSystem.add_system<TransientFEMSystem>("IdealSystem");
 
 // allocate memory to hold the adjoint time history
 ideal_system.ReserveVector(user->get_max_steps()+2);
 ideal_system.ReserveOldVector(1);

 /* the following system will hold uncertainty data */
 TransientFEMSystem & ideal_uncertainty_system =
   EqnSystem.add_system<TransientFEMSystem>("IdealUncertaintySystem");
 
 // allocate memory to hold the adjoint time history
 ideal_uncertainty_system.ReserveVector(user->get_max_steps()+2);
 ideal_uncertainty_system.ReserveOldVector(1);

 // Declare the system and its variables.
 // Creates a transient system named "AdjointSystem"
 TransientFEMSystem & adjoint_system = 
   EqnSystem.add_system<TransientFEMSystem> ("AdjointSystem");
 adjoint_system.ReserveVector(user->get_max_steps()+2);
 adjoint_system.ReserveOldVector(1);

 // Declare the system and its variables for the sensitivity problem.
 TransientFEMSystem & sensitivity_system =
   EqnSystem.add_system<TransientFEMSystem>("SensitivitySystem");
 
 //Set system function pointer for the matrix/vector assembly
 adjoint_system.attach_assemble_function (assemble_adjoint);

 for( unsigned int i =0 ; i < state_system.n_vars(); i++)
  {
    std::ostringstream variable_name;
    variable_name << "p" << i;
    if( user->OptimizeMesh() )
      { // Goal oriented error estimation
       adjoint_system.add_variable ( variable_name.str() , SECOND);
      } 
    else
      { 
       adjoint_system.add_variable ( variable_name.str() , FIRST);
      } 
    variable_name.str(""); // reset
    // Add the sensitivity variable "du" to 
    // "SensitivitySystem" as piecewise linear
    variable_name << "du" << i;
    sensitivity_system.add_variable(variable_name.str(),FIRST);

    variable_name.str(""); // reset
    // Add the ideal variable "u*" to "IdealSystem" as piecewise linear 
    variable_name << "u" << i << "*";
    ideal_system.add_variable(variable_name.str(),FIRST);
 
    variable_name.str(""); // reset
    // Add the ideal variable "u*" to "IdealSystem" as piecewise linear 
    variable_name << "du" << i << "*";
    ideal_uncertainty_system.add_variable(variable_name.str(),FIRST);
 
  }

 /* setup verification problem (if necessary) */
 if(user->get_num_exact() )
    {
       /* the following system will hold  the exact solution */
       ExplicitSystem & adjointExact =
          EqnSystem.add_system<ExplicitSystem>("ExactAdjoint");
 
       for( unsigned int i =0 ; i < state_system.n_vars(); i++)
        {
          std::ostringstream variable_name;
          variable_name << "r" << i;
          // add all variables as first order
          adjointExact.add_variable ( variable_name.str() , FIRST);
        }
    }
 
  // allocate memory to hold all sensitivities for a given time step
  //                         AND   at the previous time step
  //sensitivity system has special solve structure that can be exploited 
  sensitivity_system.assemble_before_solve=false;
  sensitivity_system.ReserveVector( user->get_max_steps()+2 );
  sensitivity_system.ReserveOldVector( 2 );
  libMesh::MeshBase &mesh = user->get_mesh();
  libmesh_assert( this->NDOF_CONTROL[1] < mesh.n_nodes() );

  // Give the system a pointer to the initialization function.
  sensitivity_system.attach_init_function(initial_sensitivity_matrix);

  /* reserve and initialize space form element computations of gradient 
    ALWAYS reserve space for error estimate calculation */
  this->GradBeta.assign( mesh.n_local_elem() * this->GetParamSize(),0.0); 

  CHKMEMA; // check for memory corruption use -malloc_debug to enable
  
  // set for compatibility w/ FEMSystem
  ideal_system.time_solver = 
       AutoPtr<TimeSolver>(new SteadySolver(ideal_system));
  sensitivity_system.time_solver       = 
       AutoPtr<TimeSolver>(new SteadySolver(sensitivity_system));
  adjoint_system.time_solver           =
       AutoPtr<TimeSolver>(new SteadySolver(adjoint_system));
  ideal_uncertainty_system.time_solver = 
       AutoPtr<TimeSolver>(new SteadySolver(ideal_uncertainty_system));

  /** set adjoint solver  use below to treat as a linear solver
   * 
   *  -snes_max_it 1                   # Do a maximum of one linear solve
   *  -snes_ls basicnonorms            # Don't do a line search and don't even
   *                                   # compute the residual (3)
   *  -snes_convergence_test skip      # Skip the convergence test before
   *                                   # evaluating the Jacobian.  SNES will
   *                                   # normally bail out early if it starts
   *                                   # with a sufficiently small residual.
   */
  adjoint_system.time_solver->diff_solver() =
        AutoPtr<DiffSolver>(new PetscDiffSolver(adjoint_system));
  sensitivity_system.time_solver->diff_solver() =
        AutoPtr<DiffSolver>(new PetscDiffSolver(sensitivity_system));

  PetscFunctionReturnVoid();
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "spaceTimeQOI::SetupTAO"
template< typename PDEModelSolver >
PetscErrorCode spaceTimeQOI< PDEModelSolver >::
SetupTAO(AppSolve *user)
{
  PetscFunctionBegin;

  PetscErrorCode info;

  /* Create TAO solver with desired solution method */
  info = TaoCreate(PETSC_COMM_WORLD,"tao_blmvm",&this->tao); 
  CHKERRQ(info);
  info = TaoApplicationCreate(PETSC_COMM_WORLD,&this->taoapp);
  CHKERRQ(info);

  /* create Hessian matrix */
  info = MatCreateMPIDense(PETSC_COMM_WORLD,
                           this->NDOF_CONTROL[0],this->NDOF_CONTROL[0],
                           this->NDOF_CONTROL[1],this->NDOF_CONTROL[1],
                           PETSC_NULL,&this->Hessian);CHKERRQ(info);
  /* create shell for hessian vector product */
  info = MatCreateShell(PETSC_COMM_WORLD,
                           this->NDOF_CONTROL[0],this->NDOF_CONTROL[0],
                           this->NDOF_CONTROL[1],this->NDOF_CONTROL[1],
                           &user,&this->HessianShell);CHKERRQ(info);
  info = MatShellSetOperation(this->HessianShell,MATOP_MULT,
                             (void(*)(void))hessianVectorProduct);CHKERRQ(info);
  info = MatShellSetOperation(this->HessianShell,MATOP_GET_DIAGONAL,
                             (void(*)(void))hessianVectorDiagonal);CHKERRQ(info);
  /* 
  FIXME: maybe should use MatCreate() with MatSetType()  instead of creating
  matrix directly.  may also use:
  info = MatCreateMPIAIJ(PETSC_COMM_WORLD,
                 qoiOptimizer->NDOF_CONTROL[0],qoiOptimizer->NDOF_CONTROL[0],
                 qoiOptimizer->NDOF_CONTROL[1],qoiOptimizer->NDOF_CONTROL[1],
                                        n,PETSC_NULL,
                                 N-n,PETSC_NULL,&this->Hessian);CHKERRQ(info);
          - or -
  info = MatCreateSeqBAIJ(PETSC_COMM_SELF,2,user.n,user.n,1,PETSC_NULL,&H);
  CHKERRQ(info); 

  FIXME: MatCreateSeqBAIJ is from rosenbrock example. Is this relates to 
         a block structure to the parameters ? 
  */
  info = MatSetOption(this->Hessian,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(info);
  info = MatSetOption(this->HessianShell,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(info);

  /* Set routines for function, gradient, Hessian, and Variable bounds */
  info = TaoAppSetObjectiveRoutine(this->taoapp,FormObjective,&user);CHKERRQ(info);
  info = TaoAppSetGradientRoutine(this->taoapp,FormGradient,&user);CHKERRQ(info);
  info = TaoAppSetVariableBoundsRoutine(this->taoapp,PhysBnds,&user);CHKERRQ(info);
  info = TaoAppSetHessianRoutine(  this->taoapp,MatrixFreeHessian,&user);CHKERRQ(info);
  info = TaoAppSetHessianMat(this->taoapp,this->HessianShell,this->HessianShell);

  /* Set null preconditioner.  Alternatively, set user-provided 
     preconditioner or explicitly form preconditioning matrix */
  KSP        kspHessian;                 /* PETSc Krylov subspace solver */
  PC          pcHessian;                 /* PETSc Krylov PC */
  info = TaoAppGetKSP(this->taoapp,&kspHessian); CHKERRQ(info);
  if (kspHessian){
    info = KSPGetPC(kspHessian,&pcHessian); CHKERRQ(info);
    info = PCSetType(pcHessian,PCNONE); CHKERRQ(info);
  }

  /* Set Convergence testing routines for the solver */
  info = TaoSetConvergenceTest(this->tao, TaoConverged_CheckPoint,&user); 
  CHKERRQ(info);

  /* Set the pointer to the initial guess*/
  this->GetInitialSolutionVec(user); 
  info = TaoAppSetInitialSolutionVec(this->taoapp,this->QCONTRL); 
  CHKERRQ(info); 

  /* Check for TAO command line options can set different Solvers:
         TaoSolve_NelderMead, TaoSolve_LMVM, TaoSolve_CG, TaoSolve_BLMVM  */
  info = TaoSetOptions(this->taoapp,this->tao); CHKERRQ(info);

  /* Check for linesearch options from command line 
        TaoSetOptions_LineSearch is used to set options for both:
         TaoApply_LineSearch, TaoApply_BoundLineSearch                 */
  if(this->tao->linectx) info = TaoLineSearchSetFromOptions(this->tao); 
  CHKERRQ(info);

  EquationSystems    &EqnSystem = user->get_equation_systems();
  TransientFEMSystem& state_system = 
     EqnSystem.get_system<TransientFEMSystem>("StateSystem");
  if(!user->IDOPT)
     { // set IC
       *state_system.vector_solution.at(user->qoiOptimizer()->Nsteplo()) 
                     = *state_system.current_local_solution;
     }
   else
     { // get the IC
       *state_system.current_local_solution 
             = *state_system.vector_solution.at(user->qoiOptimizer()->Nsteplo());
     }


 PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------
     Setup Parallel Vectors for Optimization solve
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "noOptQOI::SetupCntrl"
template< typename PDEModelSolver >
PetscErrorCode noOptQOI< PDEModelSolver >::
SetupCntrl(AppSolve *user, GetPot &controlfile)
{

  PetscErrorCode info; /* used to check for functions returning nonzeros */
  PetscFunctionBegin;
  std::cerr << "SetupCntrl not ready " << std::endl << std::flush; 
  libmesh_error();
//
//  // setup the variables that will be optimized
//  this->m_pdeSolver->SetupOptimizationVariables(this->Parameters);
//
//  std::vector<PetscTruth>::iterator ConstIter;
//  std::vector<PetscInt> ConstOwnr, FieldOwnr;
//             
//  /* Setup Communication Data Structures for Control variables
//     first count the number of field parameters and constant parameters 
//
//     0 <  subdomain  < n_block 
//     will always completely be a constant portion 
//     of the field associated with a spatial parameter 
//     elements within subdomain == 0 may or may not vary spatially BUT
//     THERE WILL always be a constant portion of the field within subdomain 0
//  */
//  PetscInt nconstparams=0, nfieldparams=0;
//
//  for(this->IterParam=this->Parameters.begin();
//this->IterParam!=this->Parameters.end(); this->IterParam++)
//   {
//    optimizationParameter* optParam = *(this->IterParam);
//    if( optParam->spatial_field ) // spatially varying parameter
//     {
//       // A CONSTANT PORTION ALWAYS EXISTS
//       nconstparams += user->get_num_elem_blk() ; 
//       // there MAY be some spatial varying params
//       nfieldparams += optParam->dofs.size() - user->get_num_elem_blk() ; 
//     }
//    else 
//     {
//       // same for TIME DEPENDENT AND constant
//       nconstparams += optParam->dofs.size();
//     }
//   }
//
//  // initialize local number of owner dofs
//  this->NDOF_CONTROL[0] = 0;
//  this->NDOF_CONTROL[1] = 0;
//  this->NDOF_CONTROL[2] = 0;
//  // all spatially constant parameters belong to rank 0 
//  if(!this->rank) this->NDOF_CONTROL[0]= nconstparams;
//  // add the spatially varying parameters
//  this->NDOF_CONTROL[0] += nfieldparams;
//
//  //the global total is the processor-wise sum
//  info = MPI_Allreduce(&this->NDOF_CONTROL[0],&this->NDOF_CONTROL[1],1,
//                       MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(info);
//
//  //total # of dof OWNED AND SHARED by this processor
//  this->NDOF_CONTROL[2] = nconstparams + nfieldparams;
//
//  /* create the map from the local gradient vector 
//     to the parallel gradient vector */  
//  this->locmap.resize(this->NDOF_CONTROL[2],0);
//
//  /* ---------------- Create parallel control vector   ------------------- */
//  info = VecCreateMPI(PETSC_COMM_WORLD,this->NDOF_CONTROL[0],
//         this->NDOF_CONTROL[1],&this->QCONTRL);CHKERRQ(info);
//  info = VecSetFromOptions(this->QCONTRL);CHKERRQ(info);
//  /* ---------   Create local contribution of parallel vector ---------- */
//  info=VecCreateSeq(PETSC_COMM_SELF,this->NDOF_CONTROL[2],&this->CNTL_LOC); CHKERRQ(info);
//  info = VecSetFromOptions(this->CNTL_LOC);CHKERRQ(info);
//
//  // get ownership range for bounds check
//  PetscInt ilo,ihi;
//  info = VecGetOwnershipRange(this->QCONTRL,&ilo,&ihi); CHKERRQ(info);
//
//
//  // build map from element optimization variable to
//  // global (parallel) optimization variable vector 
//  libMesh::MeshBase &mesh = this->m_pdeSolver->get_mesh();
//
//  CHKMEMA; // check for memory corruption use -malloc_debug to enable
//
//  {
//    // allocate space to store element mapping (scope for the iterator)
//    // looping over ALL elements not just local portion
//    libMesh::MeshBase::const_element_iterator el     = mesh.local_elements_begin();
//    const libMesh::MeshBase::const_element_iterator el_end = mesh.local_elements_end();
//    for ( ; el != el_end; ++el)
//     {
//        // Store a pointer to the element we are currently
//        // working on.  This allows for nicer syntax later.
//        Elem* elem = *el;
//        elem->_param_opt_map.resize( this->GetParamSize() , 0);
//
//        // error check the mesh
//        if( elem->subdomain_id() >   user->get_num_elem_blk() )
//          {
//           std::cout <<": element "<< elem->id() 
//                     <<" has subdomain id " << (int) elem->subdomain_id() 
//                     <<" max subdomain id " << user->get_num_elem_blk();
//           std::cout << std::endl << std::flush ; 
//           libmesh_error();
//          }
//     }
//  }
//
//  CHKMEMA; // check for memory corruption use -malloc_debug to enable
//
//  /* setup dof ordering for optimization global vec should look like
//                        _                       _
//     -----------       | param_0[domain 0]       | (constant field section)
//          |            | param_0[domain 1]       |
//          |            |         .               |
//          |            | param_0[domain n_block] |
//          |            | param_1[domain 0]       |
//          |            | param_1[domain 1]       |
//          |            |         .               |
//          |            | param_1[domain n_block] |
//          |            | param_2[domain 0]       |
//          |            | param_2[domain 1]       |
//          |            |         .               |
//          |            | param_2[domain n_block] |
//     contained on      |         .               |
//      root node        |         .               |
//      (rank = 0)       |         .               |
//          |            | param_N[domain 0]       |
//          |            | param_N[domain 1]       |
//          |            |         .               |
//          |            | param_N[domain n_block] | (end constant field section)
//          |            | ----------------------- |  
//          |            |   param_0[time 0]       | (time varying section)
//          |            |   param_0[time 1]       |
//          |            |           .             |    (only power for now)
//          |            |           .             |
//          |            |   param_0[time n]       | (end time varying section)
//     -----------       | ----------------------- |  
//          |            | param_0[     0      ]   | (varying field section)
//          |            | param_1[     0      ]   |
//          |            | param_2[     0      ]   |
//          |            |         .               |
//          |            | param_N[     0      ]   |
//          |            | param_0[     1      ]   | 
//          |            | param_1[     1      ]   |
//          |            | param_2[     1      ]   |    setup to keep dofs local
//    distributed on     | param_0[     2      ]   |
//       all nodes       | param_1[     2      ]   |
//    (including root)   | param_2[     2      ]   |
//     rank = 0-(N-1)    |         .               |
//          |            | param_N[     2      ]   |
//          |            |         .               |
//          |            |         .               |
//          |            |         .               |
//          |            | param_0[ NfieldElem ]   |
//          |            | param_1[ NfieldElem ]   |
//          |            | param_2[ NfieldElem ]   |
//          |            |         .               |
//          |            | param_N[ NfieldElem ]   | (end varying field section)
//     -----------        -                       -
//  */
//  unsigned int id_spatial_param = 0 , //counter for ith spatiallly varying param
//               rootcount=0, // counter for constant params on root (rank 0)
//               dofcount =0; // counter for this procs total number of dofs
//  for(this->IterParam=this->Parameters.begin();
//this->IterParam!=this->Parameters.end(); this->IterParam++)
//   {
//     optimizationParameter* optParam = *(this->IterParam);
//     int paramID = distance(this->Parameters.begin(),this->IterParam);
//
//     // iterator endpoints for loop over processor local elements
//     libMesh::MeshBase::const_element_iterator       el    =mesh.active_local_elements_begin();
//     const libMesh::MeshBase::const_element_iterator el_end=mesh.active_local_elements_end();
//
//     if( optParam->time_vary ) // time varying parameter
//      {
//        for ( ; el != el_end; ++el)
//         {
//           // Store a pointer to the element we are currently
//           // working on.  This allows for nicer syntax later.
//           Elem* elem = *el;
//           /* -all spatially CONSTANT variables are handled the same
//              -for time varying parameters their is an additional
//               offset added during the hessian/gradient computation  */
//           elem->_param_opt_map.at(paramID) =  paramID;
//         }
//         for( unsigned int iii  = 0 ; iii < optParam->dofs.size(); iii++  )
//         {
//            this->locmap.at( dofcount + iii ) = rootcount;
//            rootcount++;
//         }
//      }
//     else // spatial parameter that MAY ALSO vary spatially
//      {
//        // the first part of the map should always be the constant parameters
//        for( PetscInt iii = 0 ; iii < user->get_num_elem_blk() ; iii++)
//                       this->locmap.at( dofcount + iii ) = rootcount + iii ;
//        // counter for this procs local spatially varying params
//        // begins at user->get_num_elem_blk() b/c constant parameter starts at ZERO
//        // HERE we use the ASSUMPTION that a spatially constant 
//        // portion of the field ALWAYS exists...
//        unsigned int loccount = user->get_num_elem_blk();  
//        for ( ; el != el_end; ++el)
//         {
//           // Store a pointer to the element we are currently
//           // working on.  This allows for nicer syntax later.
//           Elem* elem = *el;
//           if( elem->_mat_data_field >= user->get_num_elem_blk() 
//                                     && optParam->spatial_field) 
//            {/* spatially VARYING parameter field
//                AND this particular element does not belong to the 
//                constant portion of this field 
//
//              NOTE: NfieldElem and _mat_data_field are precomputed variables on 
//              the GLOBAL mesh; can use them w/o values from other processors 
//
//              the first index of the varying field element on the first 
//              processor should start at i = nconstparams  and _mat_data_field
//              begins at user->get_num_elem_blk(), ie the indices are ugly...
//              the numbering is setup to keep dof local
//              */
//            PetscInt dofID =  
//                // begin with # of constant parameters
//                nconstparams           +
//		// _mat_data_field computed over the global mesh acts as a
//		// counter offset by the number of blocks. can use it w/o
//		// needing values from other procs;
//		// subtract one to account for error parameter
//                (elem->_mat_data_field - user->get_num_elem_blk()) * 
//                      (this->GetParamSize() - 1) +
//		// given the element position offset by param id
//                paramID - 1; // subtract one to account for error parameter
//
//            // bounds check used in GetCntrlVars, physicalBounds
//            if( !(dofID >= ilo && dofID < ihi) )
//             {
//              std::cerr << "invalid mapping"        << std::endl 
//                        << "dofID="                 <<   dofID 
//                        << " ilo="                  <<   ilo 
//                        << " ihi="                  <<   ihi 
//                        << " nconstparams="         <<   nconstparams    
//                        << " elem->_mat_data_field="<<   elem->_mat_data_field 
//                        << " user->get_num_elem_blk()="    <<
//user->get_num_elem_blk()
//                        << " this->GetParamSize()=" <<   this->GetParamSize()
//                        << " paramID="              <<   paramID
//                        << std::endl << std::flush; 
//              libmesh_error();
//             }
//
//            // store map
//            elem->_param_opt_map.at(paramID) = dofID;
//            this->locmap.at( dofcount + loccount ) = dofID;
//
//            // update counter
//            loccount++;
//            }
//           else if( elem->_mat_data_field >= user->get_num_elem_blk() )
//	    {/* _mat_data_field >= user->get_num_elem_blk(): belongs to the constant 
//                part of the field. map to original domain */
//            elem->_param_opt_map.at(paramID) = paramID *
//user->get_num_elem_blk() +
//                                                    elem->subdomain_id() - 1;
//            // special case for error estimate always map to zero
//            if( optParam->name().find("error") !=std::string::npos )
//              elem->_param_opt_map.at(paramID) = paramID *
//user->get_num_elem_blk();
//            }
//           else 
//	    {/* _mat_data_field < user->get_num_elem_blk(): belongs to the constant 
//                part of the field of corresponding mesh domain */
//            elem->_param_opt_map.at(paramID) = paramID *
//user->get_num_elem_blk()
//                                                   + elem->_mat_data_field;
//            // special case for error estimate always map to zero
//            if( optParam->name().find("error") !=std::string::npos )
//              elem->_param_opt_map.at(paramID) = paramID *
//user->get_num_elem_blk();
//            }
//         }
//        libmesh_assert( loccount == optParam->dofs.size() );  
//        // keep track of how many spatially varying parameters  we have seen 
//        if(optParam->spatial_field) id_spatial_param++; 
//        //there is always one CONSTANT Part for each mesh subdomain
//        rootcount += user->get_num_elem_blk();
//      }
//     // update counter
//     dofcount = dofcount + optParam->dofs.size();
//   }  // end loop over parameters
//
//  // error check
//  // bounds check used in GetCntrlVars, physicalBounds
//  if( !(nconstparams + user->get_num_field_elem() * id_spatial_param 
//                                              == this->NDOF_CONTROL[1]) )
//   {
//    std::cerr << "unexpect count "        << std::endl 
//              << " nconstparams="         <<nconstparams 
//              << " user->get_num_field_elem()=" <<user->get_num_field_elem() 
//              << " id_spatial_param="     <<id_spatial_param 
//              << " NDOF_CONTROL[1]="      <<this->NDOF_CONTROL[1]
//              << std::endl << std::flush; 
//    libmesh_error();
//   }
//
//  /* ---------   Create local to global scattering context ---------- */
//  IS    locis,globis;
//  // index set of stride one
//  info=ISCreateStride(PETSC_COMM_SELF,this->NDOF_CONTROL[2],0,1,&locis);
//  CHKERRQ(info);
//  // index set of global dofs
//  info=ISCreateGeneral(PETSC_COMM_SELF,this->NDOF_CONTROL[2],
//     	                                    &this->locmap[0],&globis); CHKERRQ(info);
//  // create scattering context 
//  info = VecScatterCreate(this->CNTL_LOC,locis,this->QCONTRL,globis,
//	     	                                 &this->GLOCSCAT_CNTL); CHKERRQ(info);
//  // index sets are no longer needed 
//  info=ISDestroy(locis ); CHKERRQ(info);
//  info=ISDestroy(globis); CHKERRQ(info);
//
//  // hack for computing hessian subset 
//  this->minHessianColumn = controlfile("compexec/minhessiancol", 0  );
//  this->maxHessianColumn = controlfile("compexec/maxhessiancol", 
//                                              this->NDOF_CONTROL[1] );
//  // error check
//  if (this->minHessianColumn < 0 ) this->minHessianColumn = 0 ; 
//  if (this->maxHessianColumn > this->NDOF_CONTROL[1] )
//                   this->maxHessianColumn = this->NDOF_CONTROL[1] ; 
//  // error check
//  if(this->minHessianColumn > this->maxHessianColumn)
//    {
//      std::cerr << "error input hessian columns " << std::endl << std::flush; 
//      libmesh_error();
//    }
//
  PetscFunctionReturn(0);
}
/* --------   read MRTI data for objective function evaluations ---------- */
#undef __FUNCT__
#define __FUNCT__ "onlyMRTIQOI::PriorLoadQoiData"
template< typename PDEModelSolver >
PetscErrorCode onlyMRTIQOI< PDEModelSolver >::
PriorLoadQoiData( AppSolve *user )
{
  PetscFunctionBegin; 

  // determine which files are already in memory and which need to be read in
  for(PetscInt idMRTIfile  = user->IdealNzero()  ; 
               idMRTIfile <= user->IdealNtime()  ; idMRTIfile++)
   {
     this->ExtractThermalImageTimeInstance(user,idMRTIfile);
   }
  PetscFunctionReturn(0);
}
/* --------   read MRTI data for objective function evaluations ---------- */
#undef __FUNCT__
#define __FUNCT__ "onlyMRTIQOI::ExtractThermalImageTimeInstance"
template< typename PDEModelSolver >
void onlyMRTIQOI< PDEModelSolver >::
ExtractThermalImageTimeInstance(AppSolve *user, PetscInt idMRTIfile)
{
 // used for error handling
 PetscFunctionBegin; 
 PetscErrorCode ierr;

 EquationSystems    &EqnSystem = user->get_equation_systems();
 // Get a reference to the Systems for ideal data
 TransientFEMSystem & ideal_system =
    EqnSystem.get_system<TransientFEMSystem>("IdealSystem");
 TransientFEMSystem & ideal_uncertainty_system =
    EqnSystem.get_system<TransientFEMSystem>("IdealUncertaintySystem");
 if(this->MRTI_MEM.at(idMRTIfile)){
     /* if data is already there 
              DO NOTHING                */
 }else{
    /* Data Server will read in the data from disk and broadcast */
    PetscLogEventBegin(AppSolve::logevents[3],0,0,0,0); // read disk
    if(this->Images != NULL)// skip if images data struc not defined
     {
      switch( libMesh::processor_id() )
        {
         case 0:  // read on root then broadcast
          {
           // update imaging data strucs
           this->Images->UpdateImaging(idMRTIfile);

          }// don't break... continue and broadcast
         default:
           InputImageType::PixelContainer* RawData; 
           RawData = this->Images->net_Image->GetPixelContainer();
           MPI_Bcast(RawData->GetBufferPointer(), RawData->Size(),
                     MPIU_SCALAR, 0, PETSC_COMM_WORLD);
           RawData = this->Images->snr_Image->GetPixelContainer();
           MPI_Bcast(RawData->GetBufferPointer(), RawData->Size(),
                     MPIU_SCALAR, 0, PETSC_COMM_WORLD);
        } 
      // put MRTI thermal image into FEM data structures
      if( this->m_BackgroundCorrection ) // use background correction
       {
        EqnSystem.parameters.set<PetscScalar>( "DefaultImageValue") =  0.0;
        this->ExtractImageDataFromBackground(this->Images->net_Image,
                                             idMRTIfile,ideal_system);
       }
      else // conventional CPD
       {
        Point dum;
        EqnSystem.parameters.set<PetscScalar>( "DefaultImageValue") = 
            user->pdeSolver()->InitValues(0,0,dum,EqnSystem.parameters);
        this->ExtractImageData(this->Images->net_Image,ideal_system);
       }
      printf("group %d received mrti time instance %d\n",
                               user->GroupID,idMRTIfile);

      // read the uncertainty data
      EqnSystem.parameters.set<PetscScalar>( "DefaultImageValue") =  this->m_sigma ;
      this->ExtractImageData(this->Images->snr_Image,ideal_uncertainty_system);
      printf("group %d received mrti uncertainty time instance %d\n",
                                                 user->GroupID,idMRTIfile);
      // data in ideal_uncertainty_system should be the 
      // standard deviation. need to square
      ierr = VecPointwiseMult(
               (dynamic_cast< PetscVector<double>* > 
                    (&(*ideal_uncertainty_system.solution)))->vec(),
               (dynamic_cast< PetscVector<double>* > 
                    (&(*ideal_uncertainty_system.solution)))->vec(),
               (dynamic_cast< PetscVector<double>* > 
                    (&(*ideal_uncertainty_system.solution)))->vec());

      // apply dirichlet data if any
      if(this->m_pdeSolver->m_dirichletNodes.size())
        {
          //return probe temp 
          for( unsigned int Ii = 0; 
                  Ii<this->m_pdeSolver->m_dirichletNodes.size(); Ii++) 
            {
               Point dum;
               ideal_system.solution->set(this->m_pdeSolver->m_dirichletNodes[Ii],
            this->m_pdeSolver->InitValues(0,0,dum, EqnSystem.parameters)
                                          );
            ideal_uncertainty_system.solution->set(
               this->m_pdeSolver->m_dirichletNodes[Ii], this->m_sigma );
            }
          ideal_system.solution->localize(
                                        *ideal_system.current_local_solution);
          ideal_uncertainty_system.solution->localize(
                          *ideal_uncertainty_system.current_local_solution);
        }
     }
    else
     { // for regression test get data from command lien
        PetscScalar meascov = 1.0 + 1.0/this->m_pdeSolver->deltat;
        ierr=PetscOptionsGetScalar(PETSC_NULL,"-meascov",&meascov,PETSC_NULL);
        ideal_uncertainty_system.current_local_solution->zero();
        ideal_uncertainty_system.current_local_solution->add(meascov);
        ideal_uncertainty_system.current_local_solution->close();
        ideal_uncertainty_system.solution->zero();
        ideal_uncertainty_system.solution->add(meascov);
        ideal_uncertainty_system.solution->close();
     }

    // store the solution history
    PetscInt nstephi =  idMRTIfile     * user->IstepsPerIdeal() ;
    *ideal_system.vector_solution.at(nstephi)=
                                       *ideal_system.current_local_solution;
    *ideal_uncertainty_system.vector_solution.at(nstephi)=
                      *ideal_uncertainty_system.current_local_solution;
   
    // interpolate the data onto the fem data structures
    if(idMRTIfile > user->IdealNzero() )
    {
       PetscInt nsteplo = (idMRTIfile-1) * user->IstepsPerIdeal();
       for(int jjj = 1; jjj <  user->IstepsPerIdeal()  ; jjj++)
        {
           PetscInt idfemfile = nsteplo + jjj ; 
           // mrti data
           ideal_system.vector_solution.at(idfemfile)->zero();
           ideal_system.vector_solution.at(idfemfile)->add(
                   (double)(user->IstepsPerIdeal()  - jjj )/
                          (double)(user->IstepsPerIdeal() ) ,
                       *ideal_system.vector_solution.at( nsteplo ) );
           ideal_system.vector_solution.at(idfemfile)->add(
                   (double)(                            jjj )/
                          (double)(user->IstepsPerIdeal() ) ,
                       *ideal_system.vector_solution.at( nstephi ) ); 
           // associated uncertainty
           ideal_uncertainty_system.vector_solution.at(idfemfile)->zero();
           ideal_uncertainty_system.vector_solution.at(idfemfile)->add(
                   (double)(user->IstepsPerIdeal()  - jjj )/
                          (double)(user->IstepsPerIdeal() ) ,
                  *ideal_uncertainty_system.vector_solution.at( nsteplo ) );
           ideal_uncertainty_system.vector_solution.at(idfemfile)->add(
                   (double)(                            jjj )/
                          (double)(user->IstepsPerIdeal() ) ,
                  *ideal_uncertainty_system.vector_solution.at( nstephi ) );
        }
    }

    PetscLogEventEnd(AppSolve::logevents[3],0,0,0,0); // read disk
    this->MRTI_MEM.at(idMRTIfile) = 1;
 }
 PetscFunctionReturnVoid();
}

/* -------------  adjoint load for gradient ------------- */
template< typename PDEModelSolver >
void spaceTimeQOI< PDEModelSolver >::
accumulateAdjointQOI(AppSolve *user,
                                     const QGauss &qrule,
                                     std::vector<PetscScalar> &,
                                     const unsigned int &, 
                                     const std::vector<Real>& JxW,
                                     const std::vector<Point>& ,
                                     std::vector< DenseSubVector<Number> > &Fi)
{
 PetscFunctionBegin; 
 EquationSystems    &EqnSystem = user->get_equation_systems();

 // Get a reference to the NonlinearImplicitSystem we are solving
 TransientFEMSystem& state_system = 
   EqnSystem.get_system<TransientFEMSystem>("StateSystem");

 // Get a reference to the ExplicitSystem for ideal data
 TransientFEMSystem & ideal_system =
   EqnSystem.get_system<TransientFEMSystem>("IdealSystem");

 /* the following system will hold uncertainty data */
 TransientFEMSystem & ideal_uncertainty_system =
    EqnSystem.get_system<TransientFEMSystem>("IdealUncertaintySystem");

 // The numeric integration is combined in the same loop
 for (unsigned int qp=0; qp<qrule.n_points(); qp++)
  {
   // evaluate Solution
   this->m_pdeSolver->evaluateIdealSolution(qp,AppSolve::ISTEP, 
                      user->qoiOptimizer()->Nstephi(),
                      state_system, ideal_system, ideal_uncertainty_system);

   // Sub Vector contributions 
   for (unsigned int i=0; i < this->m_pdeSolver->n_p_dofs(0); i++)
     Fi[0](i) += JxW[qp] * this->m_pdeSolver->psi(0,i,qp)  * this->m_pdeSolver->getMass() 
                         * this->m_pdeSolver->getWeightedDifference(0); 
  } // end loop over gauss points
 
 PetscFunctionReturnVoid(); 
}
/* -------------  adjoint load for matrix vector product ------------- */
template< typename PDEModelSolver >
void spaceTimeQOI< PDEModelSolver >::
accumulateAdjointSensitivity(AppSolve *user,
                             const QGauss &qrule,
                             std::vector<PetscScalar> &elemParameter,
                             const unsigned int &field_id, 
                             const std::vector<Real>& JxW,
                             const std::vector<Point>& q_point,
                             std::vector< DenseSubVector<Number> > &Fi)
{
 PetscFunctionBegin; 
 EquationSystems    &EqnSystem = user->get_equation_systems();

 // Get a reference to the NonlinearImplicitSystem we are solving
 TransientFEMSystem& state_system = 
   EqnSystem.get_system<TransientFEMSystem>("StateSystem");

 // Get a reference to the LinearImplicitSystem for adjoint problem
 TransientFEMSystem& adjoint_system = 
   EqnSystem.get_system<TransientFEMSystem>("AdjointSystem");

 // Get a reference to the ExplicitSystem for ideal data
 TransientFEMSystem & ideal_system =
   EqnSystem.get_system<TransientFEMSystem>("IdealSystem");

 /* the following system will hold uncertainty data */
 TransientFEMSystem & ideal_uncertainty_system =
    EqnSystem.get_system<TransientFEMSystem>("IdealUncertaintySystem");

 // Get a reference to the LinearImplicitSystem for sensitivity problem
 TransientFEMSystem& sensitivity_system = 
   EqnSystem.get_system<TransientFEMSystem>("SensitivitySystem");

 // The numeric integration is combined in the same loop
 for (unsigned int qp=0; qp<qrule.n_points(); qp++)
  {
   // evaluate Solution
   this->m_pdeSolver->evaluateIdealSolution(qp,AppSolve::ISTEP,user->qoiOptimizer()->Nstephi(),
                          state_system, ideal_system,ideal_uncertainty_system);

   for (unsigned int i=0; i<this->m_pdeSolver->n_p_dofs(0); i++)
    {
     // evaluate adjoint  Solution
     this->m_pdeSolver->evaluateSolutionForGradient(qp,AppSolve::ISTEP,
                                               state_system,adjoint_system);

     // evaluate sensitivities at the guass point
     this->m_pdeSolver->evaluateSolutionForSensitivity(qp,AppSolve::ISTEP,sensitivity_system);

     // accumulate second derivatives
     PetscScalar secondDerivLoad = this->m_pdeSolver->getMass() *
               this->m_pdeSolver->psi(0,i,qp) *  this->m_pdeSolver->getWeightedSensitivity(0);

     /* 
        REUSE the sensitivity variable to pass in the shape function value
        TODO - altering interface will be a lot of work
     */
     this->m_pdeSolver->evaluateSensitivitySolutionAsShapeFunction(i,qp);

     /* accumulate mixed derivatives at k-1/2 */
     for( this->IterParam=this->Parameters.begin();
this->IterParam!=this->Parameters.end() ; this->IterParam++)
      {
       int jjj = distance(this->Parameters.begin(),this->IterParam);
       optimizationParameter* optParam = *(this->IterParam);
       secondDerivLoad += 
         #if defined(PETSC_USE_DEBUG)
            elemParameter.at(jjj) 
         #else
            elemParameter[jjj] 
         #endif
            * (
          this->m_pdeSolver->d2qoi_du_dm(optParam,field_id,q_point[qp])
                                           -
          0.5 * this->m_pdeSolver->d2pde_du_dm(optParam,field_id,q_point[qp])
              );
      }

     /* accumulate mixed derivatives at k+1/2 */

     // evaluate adjoint  Solution
     this->m_pdeSolver->evaluateSolutionForGradient(qp,AppSolve::ISTEP+1,
                                               state_system,adjoint_system);

     for( this->IterParam=this->Parameters.begin();
this->IterParam!=this->Parameters.end() ; this->IterParam++)
      {
       int jjj = distance(this->Parameters.begin(),this->IterParam);
       optimizationParameter* optParam = *(this->IterParam);
       secondDerivLoad += 
         #if defined(PETSC_USE_DEBUG)
            elemParameter.at(jjj) 
         #else
            elemParameter[jjj] 
         #endif
            * ( 
         - 0.5 * this->m_pdeSolver->d2pde_du_dm(optParam,field_id,q_point[qp])
              );
      }

     // Sub Vector contributions 
       Fi[0](i) += JxW[qp] * secondDerivLoad;
    } // end loop over shape functions
  } // end loop over gauss points
 
 PetscFunctionReturnVoid(); 
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "spaceTimeFDVerification::spaceTimeFDVerification"
template< typename PDEModelSolver >
spaceTimeFDVerification< PDEModelSolver >::
spaceTimeFDVerification(AppSolve *user,GetPot &controlfile, PetscInt iii):
            spaceTimeQOI < PDEModelSolver >::spaceTimeQOI(user,controlfile,iii)
{
}
#undef __FUNCT__
#define __FUNCT__ "VerificationQOI::VerificationQOI"
template< typename PDEModelSolver >
VerificationQOI< PDEModelSolver >::
VerificationQOI(AppSolve *user,GetPot &controlfile, PetscInt iii):
     spaceTimeFDVerification < PDEModelSolver >::
                                 spaceTimeFDVerification(user,controlfile,iii)
{
}
/* ------------------------------------------------------------------- 
   PriorLoadQoiData - compute a fake ideal temperature field
                            used for code verification

   Input Parameters:
   ctx  - user-defined context 
   
   Required Routines: inittempfield must be called before this routine
*/
#undef __FUNCT__
#define __FUNCT__ "VerificationQOI::PriorLoadQoiData"
template< typename PDEModelSolver >
PetscErrorCode VerificationQOI< PDEModelSolver >::
PriorLoadQoiData(AppSolve *user)
{

  // get the system for the ideal problem.
  EquationSystems    &EqnSystem = user->get_equation_systems();
  TransientFEMSystem & ideal_system =
    EqnSystem.get_system<TransientFEMSystem>("IdealSystem");

  PetscErrorCode info;

  // used for error handling
  PetscFunctionBegin; 

  //PetscScalar Tau = AppSolve::getTime(user->qoiOptimizer()->Nstephi());
  //FORTRAN_NAME(setfinaltime)(&Tau);
  
  this->echoOptWindow(std::cout,user->qoiOptimizer()->Nsteplo(),user->qoiOptimizer()->Nstephi());

  /* verification problems making the logic complicated 
     when the Penalty if Fals then we already inside the 
     and DO NOT NEED to compute anything here*/
  if(!this->Penalty) PetscFunctionReturn(0);

  // number of time steps
  PetscInt nsteps = user->qoiOptimizer()->Nstephi() - user->qoiOptimizer()->Nsteplo() + 1 ;
  // solution buffer for temp storage..
  std::vector< NumericVector<Number>* > solution_buffer(nsteps);

  // Set all data structure to 1
  TransientFEMSystem & ideal_uncertainty_system =
     EqnSystem.get_system<TransientFEMSystem>("IdealUncertaintySystem");
  // project one onto FEM data structures
  ideal_uncertainty_system.project_solution(project_one,NULL,
                                        EqnSystem.parameters);
  for(PetscInt iii = 0;iii < nsteps ; iii++)
   {
    PetscInt idstep = user->qoiOptimizer()->Nsteplo() + iii;
    // store onto the data structures
    *ideal_uncertainty_system.vector_solution.at(idstep)=
                          *ideal_uncertainty_system.current_local_solution;
   }

  switch(user->get_num_exact() )
    {
     /* Various Cases for Verification problems*/
     case 9: case 10: case 11: case 12: case 21:
      {
      // get the system for the forward problem.
      TransientFEMSystem & system =
        EqnSystem.get_system<TransientFEMSystem>("StateSystem");

      // get the local control vars to save for later
      info = this->GetLocalCntrlVars(this->CNTL_LOC); 
      CHKERRQ(info);


      for(PetscInt iii = 0;iii < nsteps ; iii++)
       {
        PetscInt idstep = user->qoiOptimizer()->Nsteplo() + iii;
        // save the solution history clone allocates new memory
        solution_buffer.at(iii) = system.vector_solution.at(idstep)->clone().release();
       }

      // put in the control variables for this verfication problem
      if(user->get_num_exact()  != 21){info = this->PutVerifCntrlVars(); CHKERRQ(info);}
      //FORTRAN_NAME(update_iterno)(&qoiOptimizer->NumOptVar,&nzero,&nzero,
      //                        &qoiOptimizer->NDOF_CONTROL[1]);
      this->ForwardSolve(user);

      // return the control variables to the original state   	 	 
      info = this->PutLocalCntrlVars(this->CNTL_LOC);   	 	 
      CHKERRQ(info);  	 	 

      // store solution
      for(PetscInt iii = 0;iii < nsteps ; iii++)
       {
        PetscInt idstep = user->qoiOptimizer()->Nsteplo() + iii;
        // copy the solution into the ideal data structures
        *ideal_system.vector_solution.at(idstep) = *system.vector_solution.at(idstep) ;
        // restore the solution history
        *system.vector_solution.at(idstep) = *solution_buffer.at(iii) ;
        // clean up the memory when done...
        delete solution_buffer.at(iii);
       }

      // reset IC
      *system.current_local_solution = 
      *system.vector_solution.at(user->qoiOptimizer()->Nsteplo()) ;
        
      }
     break;
     default:
       info = this->noOptQOI < PDEModelSolver >::PriorLoadQoiData(user);
       CHKERRQ(info);  	 	 
     break;
    }

  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- 
   verifcomparefdgrad - use a fake ideal temperature field
                        to check adjoint gradient against fd gradient
                        for code verification

   Input Parameters:
   ctx  - user-defined context 
   
*/
#undef __FUNCT__
#define __FUNCT__ "spaceTimeFDVerification::compareFDGrad"
template< typename PDEModelSolver >
PetscErrorCode spaceTimeFDVerification< PDEModelSolver >::
compareFDGrad(TAO_APPLICATION Taoapp, Vec X,Vec G, AppSolve *user)
{
  Vec XCOPY,  *CNTRLFD, GRADFD, DBETA, DIFF, QOI, TMPVEC;
  PetscScalar qoifd,qoi, dum,*gradloc;
  PetscErrorCode info;
  EquationSystems    &EqnSystem = user->get_equation_systems();
  // get the system for the forward problem.
  TransientFEMSystem & system =
    EqnSystem.get_system<TransientFEMSystem>("StateSystem");

  // used for error handling
  PetscFunctionBegin; 

  // rank and size
  //if( !libMesh::processor_id() ) {
  //    for(PetscInt iii= this->IDEAL_NZERO + 1; iii <= this->IDEAL_NTIME; iii++)
  //       printf(" timeerr(%d)= %e\n",iii,AppSolve::TIMEERR.at(iii));
  //    fflush(stdout);
  //}

  switch(user->get_num_exact() ){
    case 9: case 10: case 11: case 14: case 17: 
      PetscPrintf(PETSC_COMM_WORLD,"comparing adjoint grad to fd grad...\n");
      PetscInt nsteps = user->qoiOptimizer()->Nstephi() - user->qoiOptimizer()->Nsteplo() + 1 ;// # of time steps
      // allocate memory
      info=VecDuplicate(X,&XCOPY );CHKERRQ(info); // store initial values 
      info=VecCopy(X,XCOPY );CHKERRQ(info); // store initial values 
      // initialize work vectors
      info=VecDuplicateVecs(X, this->NDOF_CONTROL[1] ,&CNTRLFD);
      CHKERRQ(info); 
      info=VecDuplicate(X,&GRADFD );CHKERRQ(info);
      info=VecDuplicate(X,&DBETA );CHKERRQ(info); 
      info=VecDuplicate(X,&TMPVEC );CHKERRQ(info); 
      info=VecDuplicate(X,&DIFF );CHKERRQ(info); 
      info=VecDuplicate(X,&QOI );CHKERRQ(info); 

      // solution buffer for temp storage..
      std::vector< NumericVector<Number>* > solution_buffer(nsteps);

      for(PetscInt iii = 0;iii < nsteps ; iii++)
       {
        PetscInt idstep = user->qoiOptimizer()->Nsteplo() + iii;
        // save the solution history clone allocates new memory
        solution_buffer.at(iii) = system.vector_solution.at(idstep)->clone().release();
       }


      VecSet(DBETA,0.0);
      // get perturbations
      for (PetscInt iii=0; iii < this->NDOF_CONTROL[1] ; iii++){  
         VecSet(this->CNTL_LOC,0.0);
         this->GetVerifDeltaCntrl(this->CNTL_LOC,iii); 
         CHKERRQ(info);
         info = VecScatterBegin(this->GLOCSCAT_CNTL,
              this->CNTL_LOC,TMPVEC,INSERT_VALUES, SCATTER_FORWARD); 
         CHKERRQ(info);
         info = VecScatterEnd( this->GLOCSCAT_CNTL,
              this->CNTL_LOC,TMPVEC,INSERT_VALUES,SCATTER_FORWARD); 
         CHKERRQ(info);
         /* TMPVEC = [0,0,... ,ith perturbation,...,0,0]*/
         VecPointwiseMult(CNTRLFD[iii],X,TMPVEC);
         VecAXPY(DBETA,1.0,CNTRLFD[iii]); // DBETA is the array of perturbations
         /* CNTRLFD[iii] = 
                  [x_1,x_2,...,x_iii + iiith perturbation,...,x_{n-1},x_n]  */
         VecAXPY(CNTRLFD[iii],1.0,X);
      }


      /*     COMPUTE GRADIENT USING FINITE DIFFERENCES  
             DO NOT compute FD on PENALTY function                  */
      this->Penalty = PETSC_FALSE; 

      // now perturb the input
      for (PetscInt iii=0; iii<this->NDOF_CONTROL[1]; iii++){  

             // compute perturbed qoi
             info=TaoAppComputeObjective(Taoapp,CNTRLFD[iii],&qoifd);
             CHKERRQ(info);

             // Put qoifd in parallel vector
             info = VecSetValue(DIFF,iii,qoifd,INSERT_VALUES); CHKERRQ(info);
      }
      CHKMEMQ; // check for memory corruption use -malloc_debug to enable

      // assemble parallel vector
      info = VecAssemblyBegin(DIFF); CHKERRQ(info);
      info = VecAssemblyEnd(DIFF); CHKERRQ(info)

      /* re-compute initial value of qoi ***NOTE*** this also
          serves to reset the values for the next optimization step  
          AND for the exact gradient of the penalty terms  */
      // reset IC
      *system.current_local_solution = 
      *system.vector_solution.at(user->qoiOptimizer()->Nsteplo()) ;
      info=TaoAppComputeObjective(Taoapp,XCOPY,&qoi);

      for(PetscInt iii=0; iii<this->NDOF_CONTROL[1]; iii++){  
           info = VecSetValue(QOI,iii,qoi,INSERT_VALUES); CHKERRQ(info);
      }
      info = VecAssemblyBegin(QOI); CHKERRQ(info);
      info = VecAssemblyEnd(QOI); CHKERRQ(info)

      // compute FD : (qoifd[iii] - qoi) / deltabeta[iii] 
      VecAXPY(DIFF,-1.0,QOI);
      VecPointwiseDivide(GRADFD,DIFF,DBETA);

      /*---- add penalty term derivatives to the FD gradient---- 
      info = VecSet(this->CNTL_LOC,0.0); CHKERRQ(info);
      info = VecGetArray(this->CNTL_LOC,&gradloc); CHKERRQ(info);
      // get local values of the derivatives of the regularization term 
      this->RegularizationDeriv(gradloc,&GenInfo->scalefact);
      info = VecRestoreArray(this->CNTL_LOC,&gradloc); CHKERRQ(info);
      // accumulate local contributions to global gradient
      info = VecScatterBegin(this->GLOCSCAT_CNTL,
  this->CNTL_LOC,GRADFD,ADD_VALUES, SCATTER_FORWARD);CHKERRQ(info);
      info = VecScatterEnd(this->GLOCSCAT_CNTL,
  this->CNTL_LOC,GRADFD,ADD_VALUES, SCATTER_FORWARD);CHKERRQ(info);
      */

      this->Penalty = PETSC_TRUE; // return back to normal value

      PetscPrintf(PETSC_COMM_WORLD,"Adjoint Gradient\n");
      info = VecView(G,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(info);
      PetscPrintf(PETSC_COMM_WORLD,"FD Gradient\n");
      // divide by zero for error estimate reset first entry to ZERO
      for(PetscInt iii = 0 ; iii < user->get_num_elem_blk() ; iii++)
         info=VecSetValue(GRADFD,iii,0.0,INSERT_VALUES);
      info = VecAssemblyBegin(GRADFD); CHKERRQ(info);
      info = VecAssemblyEnd(GRADFD); CHKERRQ(info);
      info = VecView(GRADFD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(info);
      info=VecWAXPY(TMPVEC,-1.0,G,GRADFD); CHKERRQ(info);
      // normalize by FD
      info=VecPointwiseDivide(G,TMPVEC,GRADFD);CHKERRQ(info);
      // divide by zero for error estimate reset first entry to ZERO
      for(PetscInt iii = 0 ; iii < user->get_num_elem_blk() ; iii++)
        info=VecSetValue(G,iii,0.0,INSERT_VALUES);
      info = VecAssemblyBegin(G); CHKERRQ(info);
      info = VecAssemblyEnd(G); CHKERRQ(info);

      // handle divide by zero errors for static domains
      this->DirichletDomains( user, G , 0.0 ); 

      info = VecNorm(G,NORM_2,&dum); CHKERRQ(info);
      PetscPrintf(PETSC_COMM_WORLD,"qoi = %13.12e, ||grad-grad_fd|| = %e \n",
                                                                       qoi,dum);
      if (dum <= 8.00e-2 ) {
         PetscPrintf(PETSC_COMM_WORLD,"gradient comparision to FD VERIFIED\n");
         PetscPrintf(PETSC_COMM_WORLD,"gradient comparision to FD VERIFIED\n");
         PetscPrintf(PETSC_COMM_WORLD,"gradient comparision to FD VERIFIED\n");
         PetscPrintf(PETSC_COMM_WORLD,"gradient comparision to FD VERIFIED\n");
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"gradient comparision N           \n");
         PetscPrintf(PETSC_COMM_WORLD,"gradient comparision  O          \n");
         PetscPrintf(PETSC_COMM_WORLD,"gradient comparision   T         \n");
         PetscPrintf(PETSC_COMM_WORLD,"gradient comparision     VERIFIED\n");
      }

      // restore solution
      for(PetscInt iii = 0;iii < nsteps ; iii++)
       {
        PetscInt idstep = user->qoiOptimizer()->Nsteplo() + iii;
        // restore the solution history
        *system.vector_solution.at(idstep) = *solution_buffer.at(iii) ;
        // clean up the memory when done...
        delete solution_buffer.at(iii);
       }
      // reset IC
      *system.current_local_solution = 
      *system.vector_solution.at(user->qoiOptimizer()->Nsteplo()) ;
      // free memory 
      info=VecDestroyVecs(CNTRLFD,this->NDOF_CONTROL[1]);
       CHKERRQ(info);
      info=VecDestroy(GRADFD );CHKERRQ(info);
      info=VecDestroy(DBETA );CHKERRQ(info); 
      info=VecDestroy(TMPVEC );CHKERRQ(info); 
      info=VecDestroy(DIFF );CHKERRQ(info); 
      info=VecDestroy(QOI );CHKERRQ(info); 
    break;
  }

  // set the gradient to zero to exit the tao solve
  switch(user->get_num_exact() ){
    case 12: // do not zero out gradient for verif ID 12
      PetscPrintf(PETSC_COMM_WORLD,"Adjoint Gradient\n");
      info = VecView(G,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(info);
      break;  
    default:
      info = VecSet(G,0.0); CHKERRQ(info);
      break;
  }

  PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiSpaceTimeRecompute::getSensitivity"
template< typename PDEModelSolver >
void qoiSpaceTimeRecompute< PDEModelSolver >::
getSensitivity(AppSolve *user,const int &globalDof,
               const std::vector<optimizationParameter*>::iterator idParamIter )
{//default is to recompute sensitivity
 //this->qoiBaseClass::computeSensitivity( user, globalDof, idParamIter );
}
