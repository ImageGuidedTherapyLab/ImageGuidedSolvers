define  es
  set endstall = 1
#  printf "rank = %d groupID = %d \n",rank ,groupID 
  continue
end

define dfpsys
call VecView(((PetscVector<double>)(*system.solution))._vec,0)
end

define id
  shell sh mpdlistjobs
end

define  dfr
  r --gtest_filter=PennesExponentialSourceTerm.NeumannBC
end

define  dfrtp
  r Examples/TreatmentPlanning/treatmentPlanning.py 
end
define  dfr0
  r --gtest_filter=PennesPlanarSourceTerm.CauchyNeumannBC
end

define  dfrbc
  r Examples/ModelAssistedMonitoring/backgroundCorrection.py  -ksp_rtol 1.e-10 -ksp_converged_reason
end
set extension-language .FPP fortran
set breakpoint pending on

# regex breakpoints
#rbreak ^vtkExodusIIReader::

#  use ${DFSRC_DIR}/_obj_1_dtf_${PETSC_ARCH}
#dir dfroutines
#dir dfmodule
#dir ../hp3d/hp_interp
#dir ../hp3d/parallel_version/parallel
#dir ../hp3d/parallel_version/parallel
#dir ../hp3d/parallel_version/parmodule
#dir ../hp3d/parallel_version/parrefine
#dir ../hp3d/parallel_version/parsolve
#dir ../hp3d/parallel_version/parutilities
#dir ../hp3d/meshgen3_2
#dir ../hp3d/meshmodb
#dir ../hp3d/utilities
#dir ../hp3d/datstrb
#dir ../hp3d/constrb
#dir ../hp3d/module
#dir ../hp3d/GMP
#dir ../hp3d/GMP_hpinterp
#dir ../hp3d/elem_util

#  SetupMesh gatherdealloc_hpelemfield setup_fields SetupHp3d adjoint_grad
#  TaoSolve_NelderMead TaoSolve_BLMVM getmrtitemp  verifcomparefdgrad
#  Init_MCdata TaoApply_BoundLineSearch METIS_PartGraphKway SetupVisualization
#  LoadIdealArr output_fem SetupPetscDstruc  Main_Visualization qoieval
#  break "nodgenb.F":213 if (nod .EQ. 101)
#  break "visualization.F":437
#  break "morethuente.c":754
#  break "dddas.cxx":538
#  break "compute_drone.cxx":605
#  break "mesh_communication.C":697
#  break "ucd_io.C":207
#  break "DDDAS_Control.cpp":155
#  break "Compute_Drone.cpp":338
#  break "hp3gendf.F":482
#  break "Setup.c":207
#  break "meshpart.c":44 
#  break "/home/utexas/iv/fuentes/hp3dDF/dfroutines/DDDAS_Control.cpp":342
#  break "/ices/fuentes/DDDAS/3DhpA_work/hp3dDF/dfroutines/Data_Server.cpp":180
#  break "MonteCarlo.cpp":325
#  break "optvars.cpp":1174
#  break "Compute_Drone.cpp":250
#  break "inittempfield.F":368
#  break "neldermead.c":112
#  break "FormFunction.cpp":84
#  break "FormGradient.cpp":204
#  break "evaluations.F":628
#  break "elemjac.FPP":189
#  break "verif_suite.cxx":157
#  break "constit_data.FPP":1059
#break "modelAssistedMonitoring.cxx":209
#break "KalmanFilter.cxx":1314
#break "KalmanFilter.cxx":625
#break "tttkUtilities.cxx":438
#break "thermal_therapy_system.txx":46
#break "pennesSystem.txx":56
#break "pennesSystem.txx":103
#break "pennesSystem.txx":427
#break "ls.c":191
#break "gmres.c":227
#break "pennesSystem.txx":195
#break "PennesRegression.cxx":540
#break "pennesModel.h":378
#break "pennesModel.h":398
#break "PennesRegression.h":265
#break "itkLinearInterpolateImageFunction.txx":86
#break "transient_fem_system.cxx":289
#break "optimizationParameter.cxx":65
