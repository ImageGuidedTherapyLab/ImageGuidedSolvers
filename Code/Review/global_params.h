// fortran function signatures
void FORTRAN_NAME(initialize_visualase)(PetscTruth*,PetscTruth*,PetscTruth*,PetscInt*);
void FORTRAN_NAME(reset_iterno)();
PetscScalar FORTRAN_NAME(get_visualase_max_temp)();
void FORTRAN_NAME(cross_product)( PetscScalar*,PetscScalar*,PetscScalar*);
// return PenaltyTerm
PetscScalar  FORTRAN_NAME(getpenalty)();
