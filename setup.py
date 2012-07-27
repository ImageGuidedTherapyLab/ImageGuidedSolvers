from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import subprocess
import os
libmesh_dir = os.getenv("LIBMESH_DIR")
#tao_dir     = os.getenv("TAO_DIR")
petsc_dir   = os.getenv("PETSC_DIR")
petsc_arch  = os.getenv("PETSC_ARCH")
mpi_dir     = os.getenv("MPI_DIR")
itk_home    = os.getenv("ITK_HOME")
itk_dir     = os.getenv("ITK_DIR")
build_dir   = "_obj_%s" % petsc_arch  

#get and parse libmesh includes,cxxflags,libs
libmesh_include=subprocess.Popen(
    "%s/contrib/bin/libmesh-config --include" % libmesh_dir ,
                  shell=True,stdout=subprocess.PIPE).communicate()[0]
libmesh_include=filter(len,libmesh_include.replace(' ','').replace('\n','').split("-I"))

# append an include to find the PETSc.pxd definition file
# get path from petsc4py path
import petsc4py
libmesh_include.append( "%s/include" % petsc4py.__path__.pop() )

# tao include
#libmesh_include.append( "%s"         % tao_dir)
#libmesh_include.append( "%s/include" % tao_dir)
# local includes
libmesh_include.append( "Code/Common")
libmesh_include.append( "Code/TreatmentPlanning")
libmesh_include.append( "Code/ModelAssistedMonitoring")
# itk includes
libmesh_include.append( "%s/include/InsightToolkit/gdcm/src"                 % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/gdcm"                     % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities"                % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities/vxl/core"       % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities/vxl/vcl"        % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities/vxl/v3p/netlib" % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities/itkExtHdrs"     % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities/nifti/znzlib"   % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities/nifti/niftilib" % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities/expat"          % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities/DICOMParser"    % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities/NrrdIO"         % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Utilities/MetaIO"         % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/SpatialObject"            % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Numerics/NeuralNetworks"  % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Numerics/Statistics"      % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Numerics/FEM"             % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/IO"                       % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Numerics"                 % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Common"                   % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/BasicFilters"             % itk_home)
libmesh_include.append( "%s/include/InsightToolkit/Algorithms"               % itk_home)
libmesh_include.append( "%s/include/InsightToolkit"                          % itk_home)

# cxx flags
libmesh_cxxflags=subprocess.Popen(
    "make --no-print-directory -C %s echo_cxxflags" % libmesh_dir ,
                  shell=True,stdout=subprocess.PIPE).communicate()[0]
libmesh_cxxflags=filter(len,libmesh_cxxflags.replace('\n','').split(' '))
# inline diagnostics 
#libmesh_cxxflags.append( "-Winline" )

# libraries to link w/
libmesh_ldflags=subprocess.Popen(
    "make --no-print-directory -C %s echo_ldflags" % libmesh_dir ,
                  shell=True,stdout=subprocess.PIPE).communicate()[0]
libmesh_ldflags=libmesh_ldflags.replace('\n','')
print "reference libMesh flags", libmesh_ldflags

# Manually set libraries to link with
#TODO- is there a more automagic way to determine the neede libraries???
libmesh_petsc_library_dirs = ["%s/lib/%s" % (libmesh_dir,petsc_arch),
                              "%s/contrib/lib/%s" % (libmesh_dir,petsc_arch),
                              "%s/%s/lib"% (petsc_dir,petsc_arch),
                              "%s/lib"%mpi_dir,itk_dir
                             ]
libmesh_petsc_libraries=["mesh",
               "exodusii",
               "laspack",
               "gmv" ,
#               "gzstream"  ,
               "Hilbert",
               "metis",
               "nemesis",
               "netcdf",
               "parmetis",
               "sfcurves",
               "tetgen",
#               "triangle",
#               "petscsnes",
#               "petscksp",
#               "petscdm",
#               "petscmat",
#               "petscvec",
#               "petsccontrib",
               "petsc",
               "ITKFEM",
               "ITKIO",
               "ITKStatistics",
               "ITKNumerics",
               "ITKBasicFilters",
               "ITKNrrdIO",
               "itkgdcm",
               "itkjpeg12",
               "itkjpeg16",
               "itkopenjpeg",
               "itkpng",
               "itktiff",
               "itkjpeg8",
               "ITKSpatialObject",
               "ITKMetaIO",
               "ITKDICOMParser",
               "ITKEXPAT",
               "ITKniftiio",
               "ITKznz",
               "itkzlib",
               "ITKCommon",
               "itkNetlibSlatec",
               "itkvnl_inst",
               "itkvnl_algo",
               "itkv3p_netlib",
               "itkv3p_lsqr",
               "itkvnl",
               "itkvcl",
               "itksys",
               "mpich",
               "mpichcxx"]
print "configure libmesh with --disable-tecplot"
# source filenames
source_files = ["Code/Common/wrapfemInterface.pyx",  
                "Code/Common/optimizationParameter.cxx",
                "Code/Common/tttkUtilities.cxx",
                "Code/Common/itkVTKImageVariableNameIO.cxx",
                "Code/Common/pdeBaseClass.cxx",
                "Code/Common/petsc_fem_context.cxx",
                "Code/Common/petsc_fem_system.cxx",
                "Code/TreatmentPlanning/pennesModel.cxx",
                "Code/TreatmentPlanning/pennesSystem.cxx",
                "Code/ModelAssistedMonitoring/BackgroundPhase.cxx",
                "Code/ModelAssistedMonitoring/KalmanFilter.cxx",
                "Code/ModelAssistedMonitoring/Imaging.cxx",
               ]
dependency_files = list( source_files) # copy
dependency_files.extend(
               ["Code/Common/wrapfemInterface.pxi",  
                "Code/Common/optimizationParameter.h",
                "Code/Common/tttkUtilities.h",
                "Code/Common/itkVTKImageVariableNameIO.h",
                "Code/Common/pdeBaseClass.h",
                "Code/Common/petsc_fem_context.h",
                "Code/Common/petsc_fem_system.h",
                "Code/Common/thermal_therapy_system.h",
                "Code/TreatmentPlanning/pennesModel.h",
                "Code/TreatmentPlanning/pennesSystem.h",
                "Code/ModelAssistedMonitoring/BackgroundPhase.h",
                "Code/ModelAssistedMonitoring/KalmanFilter.h",
                "Code/ModelAssistedMonitoring/Imaging.h",
               ] )  # header filenames
print source_files
ext_modules_list=[
    Extension("femLibrary", # name of extension
    source_files , # source filenames
    language="c++",              # this causes Cython to create C++ source
    include_dirs= libmesh_include, # include directories
    extra_compile_args = libmesh_cxxflags, # cxx flags
    library_dirs =libmesh_petsc_library_dirs, 
    runtime_library_dirs=libmesh_petsc_library_dirs ,
    libraries=libmesh_petsc_libraries
    ) # Unix-like specific
]


#Cython will generate and compile the signalmodel.cpp file (from the
#signalmodel.pyx), then it will compile signalmodel.cxx (implementation
#of the SignalModel class) and link both objects files together into
#signalmodel.so, which you can then import in Python using import
#signalmodel (if you forget to link the SignalModel.o, you will get
#missing symbols while importing the library in Python).

setup(
  name = "FEMInterface",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules_list
)

