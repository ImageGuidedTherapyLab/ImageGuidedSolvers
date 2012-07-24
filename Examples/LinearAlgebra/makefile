
CFLAGS 	        = 
FFLAGS	        = 
CPPFLAGS        =
FPPFLAGS        =
LOCDIR          = 
EXAMPLESC       = matmatmult.c 
EXAMPLESF       =
MANSEC          = Mat

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

matmatmult: matmatmult.o  chkopts
	-${CLINKER}  -o matmatmult matmatmult.o ${PETSC_MAT_LIB} 
	./matmatmult -nsize 5000 -log_summary


