#ifndef __SourceTerm_h
#define __SourceTerm_h
/**
 * This class is responsible for assembling the system dynamics matrix and load
 * vector for the prediction of the DCE model.
 * Can template on the  underlying model later when actually have it defined.
 * For now implement as FEMSystem was intended.
 */
class DCESystem : public FEMSystem
{

public:
  // constructor
  DCESystem();

  /** Assemble System Dynamics Matrix and Load Vector */
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext& context);

  // Indices for each variable;
  unsigned int u_var;

private:
  PetscScalar  m_KTrans;
};
#endif
