// dddas.cxx
void build_equation_system_solution_vector(EquationSystems &, std::vector<Number>& );
void build_system_solution_vector(EquationSystems &, System & , std::vector<Number>& );

// setup.cxx
PetscErrorCode SetupVisualization(PetscInt*,PetscInt*,PetscInt*);

// compute_drone.cxx
PetscErrorCode TaoConverged_CheckPoint(TAO_SOLVER tao,void *ctx);
PetscErrorCode SolveOptimizationProb(AppSolve *,GetPot &);
PetscErrorCode WriteControlFile(PetscInt , PetscInt );

// montecarlo.cxx
PetscErrorCode Load_MCdataNopt(PetscTruth);
PetscErrorCode Load_MCdataIdeal(PetscTruth); 

// verif_suite.cxx
PetscErrorCode verifcomputeidealfield(AppSolve *);
PetscErrorCode computeideal_mc(       AppSolve *);
PetscErrorCode verifcomparefdgrad(TAO_APPLICATION ,Vec ,Vec ,void *);

// ???
PetscErrorCode TheRedPill(const char *,const char *);
PetscErrorCode Init_MCdata(PetscTruth); 
