include ${PETSC_DIR}/conf/variables

.PHONY: PennesDMDACudaKernel
PennesDMDACudaKernel: PennesDMDACudaKernel.cu
	${CUDAC} -g -O0 -G0 -gencode arch=compute_20,code=sm_20 -c --compiler-options="${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS}"  $<
	${CLINKER} -o PennesDMDACudaKernel PennesDMDACudaKernel.o ${PETSC_SNES_LIB}
	#${RM} PennesDMDACudaKernel.o
include ${PETSC_DIR}/conf/test
