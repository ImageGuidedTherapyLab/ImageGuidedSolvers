PetscScalar    SpaceTime_QOI(AppSolve *ctx);
PetscScalar    Lesbegue2_QOI(AppSolve *ctx);
PetscErrorCode SetControlVariablePointer(CntrlVars *);
PetscErrorCode LoadIdealHSP(   AppSolve *ctx);
PetscErrorCode LoadIdeal_TS(   AppSolve *ctx);
PetscErrorCode LoadIdealArr(   AppSolve *ctx);
PetscErrorCode CurrentIdealArr(AppSolve *ctx);
PetscErrorCode LoadIdealTmp(   AppSolve *ctx);
