#======================================================================
#
#  The makefile was last modified on $Date$ by $Author$
#
#======================================================================
#
#####################
# NOTE
#COMPILER OPTIONS IN THE FILE:
#       ${PETSC_DIR}/bmake/linux-gnu/petscconf
#####################

###############################
#${PETSC_ARCH} in library linking statements to facilitate
#use of multiple compilers 
obj_path     = _obj_$(PETSC_ARCH)
EXEC        := dddas_$(PETSC_ARCH) 
stateEXEC   := treatmentPlanning_$(PETSC_ARCH)
imagingEXEC := image_$(PETSC_ARCH) 

all: 
	make Testing 
	#make pdeopt DEPENDENCY=1
	#make stateSolve  DEPENDENCY=1
	#make imaging DEPENDENCY=1


####### BASE PETSC  INCLUDE FILES are included in TAO base files
include ${TAO_DIR}/bmake/tao_common
include ${DDDAS_SRC}/conf/variables

#If you wish to eliminate the default known suffixes instead of just adding to them, write a rule for .SUFFIXES with no prerequisites. By special dispensation, this eliminates all existing prerequisites of .SUFFIXES. You can then write another rule to add the suffixes you want. For example,
#.SUFFIXES:            # Delete the default suffixes
#.SUFFIXES: .c .o .h   # Define our suffix list
#The `-r' or `--no-builtin-rules' flag causes the default list of suffixes to be empty.
#The variable SUFFIXES is defined to the default list of suffixes before make reads any makefiles. You can change the list of suffixes with a rule for the special target .SUFFIXES, but that does not alter this variable. 
.SUFFIXES:
#    The prerequisites of the special target .PHONY are considered to be phony targets. When it is time to consider such a target, make will run its commands unconditionally, regardless of whether a file with that name exists or what its last-modification time is
.PHONY: build_error banner tar tgz dfclean realclean link done check verify tags dddas_tags mpi_tags petsc_tao_tags libmesh_tags itk_tags Testing $(obj_path)/femLibrary.so
 

###############################
# the order of the source files
# is important for compiling
###############################
SOURCE_FEM    = Code/Common/tttkUtilities.cxx \
                Code/Common/itkVTKImageVariableNameIO.cxx  \
                Code/Common/petsc_fem_system.cxx \
                Code/Common/petsc_fem_context.cxx \
                Code/Common/optimizationParameter.cxx  \
                Code/Common/pdeBaseClass.cxx  
SOURCE_PLAN   = Code/TreatmentPlanning/pennesModel.cxx  \
                Code/TreatmentPlanning/pennesSystem.cxx  
SOURCE_MPAM   = Code/ModelAssistedMonitoring/Imaging.cxx \
                Code/ModelAssistedMonitoring/BackgroundPhase.cxx 
SOURCE_REVIEW = Code/Review/pennesVerification.cxx \
                Code/Review/byte_swapping.cxx \
                Code/Review/string_utils.f90 \
                Code/Review/parse_ini.f90 \
                Code/Review/plot_params.f90 \
                Code/Review/global_params.f90 \
                Code/Review/pennes_model.f90 
SOURCE_TEST   = Testing/PennesRegression.cxx

SOURCE_ALL = $(SOURCE_FEM) $(SOURCE_PLAN) $(SOURCE_MPAM) 

# this is real target if nothing is defined
link: banner exec/$(EXEC) done


# make this as a check to be sure whether we are computing at ACES or TACC
PREC_FLAGS = 

# Setup Preprocessing Flags for different architectures
ifeq ($(PETSC_ARCH),gcc-4.4.3-$(MPI_VERSION)-cxx-dbg)
   FC_FLAGS  += -J$(obj_path) -fbounds-check -x f95-cpp-input 
#   SCRATCH = /workspace   
#   PDE_LIB += -luuid
#   IMG_LIB += -luuid
endif
# RIS
ifeq ($(PETSC_ARCH),gcc-3.4.6-$(MPI_VERSION)-cxx-opt)
   PDE_LIB += -luuid
endif
# intel
ifeq ($(PETSC_ARCH),intel-10.1-$(MPI_VERSION)-cxx-dbg)
   PDE_LIB += -luuid
endif
# shamu
ifeq ($(PETSC_ARCH),intel-10.1-$(MPI_VERSION)-cxx-opt)
   FC_FLAGS  += -module $(obj_path) -CB -traceback -fpp
   PDE_LIB += -luuid
endif
# TACC
ifeq ($(PETSC_ARCH),barcelona-cxx)
   MPI_DIR    = ${MPICH_HOME}
   INCLUDE+=  -I$(TACC_BOOST_INC)/boost-1_37
   # petsc at TACC not compiled w/ X11 libraries (quick fix)
   PDE_LIB += -L/usr/X11R6/lib64  -lX11 -luuid
endif
ifeq ($(PETSC_ARCH),em64t-cxx)
   MPI_DIR    = ${MPICH_HOME}
   # petsc at TACC not compiled w/ X11 libraries (quick fix)
   PDE_LIB += -L/usr/X11R6/lib64  -lX11 -luuid
endif
ifeq ($(PETSC_ARCH),em64t-cxxdebug)
   MPI_DIR    = ${MPICH_HOME}
   # petsc at TACC not compiled w/ X11 libraries (quick fix)
   PDE_LIB += -L/usr/X11R6/lib64  -lX11
endif

#======================================================================

#we'll search for source files in these directories
VPATH = $(sort $(dir $(SOURCE_ALL))):$(sort $(dir $(SOURCE_TEST)))

# dot_o can be changed to .obj for Windows
dot_o = .o
#####################################3
#$(notdir names...)
#    Extracts all but the directory-part of each file name in names. If the file name contains no slash, it is left unchanged. Otherwise, everything through the last slash is removed from it.
#    A file name that ends with a slash becomes an empty string. This is unfortunate, because it means that the result does not always have the same number of whitespace-separated file names as the argument had; but we do not see any other valid alternative.
#    For example,
#
#$(notdir src/foo.c hacks)
#
#    produces the result `foo.c hacks'. 
#####################################3
#$(filter pattern...,text)
#    Returns all whitespace-separated words in text that do match any of the pattern words, removing any words that do not match. The patterns are written using `%', just like the patterns used in the patsubst function above.
#
#    The filter function can be used to separate out different types of strings (such as file names) in a variable. For example:
#sources := foo.c bar.c baz.s ugh.h
#foo: $(sources)
#        cc $(filter %.c %.s,$(sources)) -o foo
#
#    says that `foo' depends of `foo.c', `bar.c', `baz.s' and `ugh.h' 
#    but only `foo.c', `bar.c' and `baz.s' should be specified in
#    the command to the compiler.
#####################################3
#    Substitution references (see section Substitution References) are a simpler way to get the effect of the patsubst function:
#
#$(var:pattern=replacement)
#    is equivalent to
#$(patsubst pattern,replacement,$(var))
#    The second shorthand simplifies one of the most common uses of patsubst: replacing the suffix at the end of file names.
#$(var:suffix=replacement)
#    is equivalent to
#$(patsubst %suffix,%replacement,$(var))
#    For example, you might have a list of object files:
#objects = foo.o bar.o baz.o
#    To get the list of corresponding source files, you could simply write:
#$(objects:.o=.c)
#    instead of using the general form:
#$(patsubst %.o,%.c,$(objects))
#####################################3

#compile all module files first
OBJECT_FEM   := $(addprefix $(obj_path)/, $(SOURCE_ALL:.cxx=$(dot_o)) ) 
#======================================================================
# deprecated
$(obj_path)/femLibrary.so:
	python setup.py build_ext --build-lib=$(obj_path) --build-temp=$(obj_path)
	chmod 755 $(obj_path)
	chmod 755 $(obj_path)/femLibrary.so

# regression using googletest
Testing: $(obj_path)/femLibrary.so
	$(CLINKER) $(LIBMESH_CXXFLAGS) $(INCLUDE) Testing/PennesRegression.cxx -o Testing/PennesRegression -lgtest $(OBJECT_FEM)  $(PDE_LIB)

#TODO dependency info not currently used...
#include ${DDDAS_SRC}/conf/rules
docclean: 
	rm -rf Documentation/TreatmentPlanning/html
	rm -rf Documentation/html
doc: 
	doxygen Documentation/TreatmentPlanning.cfg
	doxygen Documentation/TTTK.cfg
	cd Documentation; pdflatex -jobname latex/references   \
           \\documentclass{article}    \
           \\begin{document}           \
           \\bibliographystyle{plain}  \
           \\bibliography{references}  \
           \\nocite{*}                 \
           \\end{document} ; bibtex latex/references
	cd Documentation/TreatmentPlanning/latex ; make
	sed 's/printindex/input{references.bbl}\n\\printindex/g' Documentation/latex/refman.tex > tmp.tex 
	mv tmp.tex Documentation/latex/refman.tex 
	cd Documentation/latex ; make

grep: 
	@grep $(KEYWORD) $(SOURCE_ALL)

# tags 
tags: 
	@echo "###ASSUMES### that source code for ctags has been hacked"
	@echo "to include keywords PetscScalar/PetscInt/PetscTruth"
	@echo "this may be done by adding to the keywords in fortran.c of "
	@echo "the ctags source "
	ctags -R -h +.blk --langmap=fortran:+.blk.FPP --langmap=c++:+.txx --languages=c,c++,fortran --regex-fortran=/MPI_HEXas/MPI_HEXAS/ --regex-fortran=/MPI_GMVVERTex/MPI_GMVVERTEX/ --regex-fortran=/MPI_CELLINfo/MPI_CELLINFO/ --regex-fortran=/MPI_NODVar/MPI_NODVAR/ $(DDDAS_SRC) $(PETSC_DIR) $(TAO_DIR) $(MPI_SOURCE) $(MPI_DIR) $(ITK_SOURCE) $(LIBMESH_DIR) $(SRC_DIR)
# tags created on local disk for performance. must open tags from local disk
# cd $SCRATCH ; vim ---.c
local_tags: 
	cd $(SCRATCH) ;ctags -R -h +.blk --langmap=fortran:+.blk.FPP --langmap=c++:+.txx --languages=c,c++,fortran --regex-fortran=/MPI_HEXas/MPI_HEXAS/ --regex-fortran=/MPI_GMVVERTex/MPI_GMVVERTEX/ --regex-fortran=/MPI_CELLINfo/MPI_CELLINFO/ --regex-fortran=/MPI_NODVar/MPI_NODVAR/ $(DDDAS_SRC) $(PETSC_DIR) $(TAO_DIR) $(MPI_SOURCE) $(MPI_DIR) $(ITK_SOURCE) $(LIBMESH_DIR) $(SRC_DIR)

sep_tags: dddas_tags  itk_tags  libmesh_tags  mpi_petsc_tao_tags  
	@echo "creating separate tag files "
	@echo "in vim  "
	@echo "       :set tags=<file>"

dddas_tags: 
	ctags -f dddas_tags -R -h +.blk --langmap=fortran:+.blk.FPP --langmap=c++:+.txx --languages=c,c++,fortran --regex-fortran=/MPI_HEXas/MPI_HEXAS/ --regex-fortran=/MPI_GMVVERTex/MPI_GMVVERTEX/ --regex-fortran=/MPI_CELLINfo/MPI_CELLINFO/ --regex-fortran=/MPI_NODVar/MPI_NODVAR/ $(DDDAS_SRC) 

itk_tags: 
	ctags -f itk_tags -R -h +.blk --langmap=fortran:+.blk.FPP --langmap=c++:+.txx --languages=c,c++,fortran  $(ITK_SOURCE) 

libmesh_tags: 
	ctags -f libmesh_tags -R -h +.blk --langmap=fortran:+.blk.FPP --langmap=c++:+.txx --languages=c,c++,fortran $(LIBMESH_DIR)

petsc_tao_tags: 
	ctags -f petsc_tao_tags -R --langmap=c++:+.inl --langmap=c:+.cu --langmap=python:+.pxi --langmap=c++:+.txx --languages=c,c++,python,fortran $(PETSC_DIR) $(TAO_DIR) $(CUDA_DIR)

mpi_tags: 
	ctags -f mpi_tags -R -h +.blk --langmap=fortran:+.blk.FPP --langmap=c++:+.txx --languages=c,c++,fortran  $(MPI_SOURCE) $(MPI_DIR)

check: 
	@echo EXEC=${EXEC} 
	@echo VPATH=${VPATH} 
	@echo MPI_DIR=${MPI_DIR} 
	@echo C=${C} 
	@echo FC_FLAGS=${FC_FLAGS} 
	@echo PETSC_ARCH=${PETSC_ARCH}
	@echo objFPP=$(objFPP)
	@echo MPIEXEC=$(MPIEXEC)
	@echo PETSC_INCLUDE=$(PETSC_INCLUDE)
	@echo TAO_INCLUDE=$(TAO_INCLUDE)
	@echo ITK_INCLUDE=$(ITK_INCLUDE)
	@echo PDE_LIB=${PDE_LIB} 
	@echo IMG_LIB=${IMG_LIB} 
	@echo OBJS=$(OBJS) 
	@echo FLINKER=${FLINKER}
	@echo CC=${CC}
	@echo FC=${FC}
	@echo PCC_FLAGS=${PCC_FLAGS}
	@echo LIBMESH_CXXFLAGS=$(LIBMESH_DIR)
	@echo PREC_FLAGS=${PREC_FLAGS} 
	@echo LIBMESH_LIB=${LIBMESH_LIB} 
	@echo TAO_LIB=${TAO_LIB} 
	@echo CLINKER=${CLINKER}
	@echo PCC=${PCC} 
	@echo SOURCE_ALL=$(SOURCE_ALL) 
	@echo INCLUDE=${INCLUDE} 
	@echo OBJECT_FEM=$(OBJECT_FEM)
