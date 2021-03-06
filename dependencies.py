"""
Setup library dependencies...
"""

import os
import ConfigParser

# read config file
config = ConfigParser.SafeConfigParser({})
config.read("./dependencies.ini")

# dependencies for this version
libmesh_version = config.get('libmesh','version')

def EnvironmentSetup(EnvVariable):
  if(os.getenv(EnvVariable) == None):
   raise ValueError("%s not set?" % EnvVariable)
  else:
   return os.getenv(EnvVariable) 

#default to release w/ debug info
build_type = "Release"
build_type = "Debug"
build_type = "RelWithDebInfo"

# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
parser.add_option( "--petsc", 
                  action="store", dest="build_petsc", default=False,
                  help="build petsc dependency", metavar="bool")
parser.add_option( "--tao", 
                  action="store", dest="build_tao", default=False,
                  help="build tao dependency", metavar="bool")
parser.add_option( "--libmesh", 
                  action="store", dest="build_libmesh", default=False,
                  help="build libmesh dependency", metavar="bool")
parser.add_option( "--itk", 
                  action="store", dest="build_itk", default=False,
                  help="build itk dependency", metavar="bool")
parser.add_option( "--vtk", 
                  action="store", dest="build_vtk", default=False,
                  help="build vtk dependency", metavar="bool")
parser.add_option( "--dakota", 
                  action="store", dest="build_dakota", default=False,
                  help="build dakota dependency", metavar="bool")
(options, args) = parser.parse_args()


# run planning solver w/ default options from ini file
if (options.build_petsc ):
  mpi_dir     = EnvironmentSetup("MPI_DIR")
  petsc_dir   = EnvironmentSetup("PETSC_DIR")
  ## if (os.getenv("MKL_DIR") = None):
  ##    print "MKL FOUND!"
  ##    # read the MKL documentation  for proper setup !!!
  ##    mkl_setup = "--with-blas-lapack-lib=[%s/mkl/lib/intel64/libmkl_rt.so,%s/mkl/lib/intel64/libmkl_intel_thread.so,%s/mkl/lib/intel64/libmkl_core.so,%s/lib/intel64/libiomp5.so] --with-blacs-include=%s/mkl/include --with-blacs-lib=[%s/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so] --download-superlu_dist --download-parmetis --download-mumps --with-scalapack-lib=[%s/mkl/lib/intel64/libmkl_scalapack_lp64.so] --with-scalapack-include=%s/mkl/include " % 
  ## (os.getenv("MKL_DIR"),os.getenv("MKL_DIR"),os.getenv("MKL_DIR"),os.getenv("MKL_DIR"),
  ##  os.getenv("MKL_DIR"),os.getenv("MKL_DIR"),os.getenv("MKL_DIR"),os.getenv("MKL_DIR")
  ## )
  ## 
  ## else: 
  ##    print "no MKL found..."
  ##    mkl_setup = "--download-f-blas-lapack=ifneeded" 
  # setup PETSc
  petsc_root = petsc_dir.split("/")
  petsc_version = petsc_root.pop().split("-")
  petsc_root = "/".join(petsc_root)
  petsc_version.pop(0)
  petsc_version = "-".join(petsc_version)
  print petsc_root ,petsc_version 
  os.system("mkdir -p %s; cd %s; wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-%s.tar.gz;tar -zxf petsc-lite-%s.tar.gz" % (petsc_root,petsc_root,petsc_version,petsc_version) )
  petscConfigOptions="cd %s;config/configure.py  %s --with-mpi-dir=%s --with-shared=1 --with-clanguage=C++ --with-petsc4py=1 --download-petsc4py=yes --with-mpi4py=1 --download-mpi4py=yes;make " % (petsc_dir,mkl_setup,mpi_dir) 
  print petscConfigOptions 
  os.system(petscConfigOptions)
elif (options.build_tao ):
    tao_dir     = EnvironmentSetup("TAO_DIR")

    # setup TAO
    tao_root = tao_dir.split("/")
    tao_root.pop()
    tao_root = "/".join(tao_root)
    os.system("mkdir -p %s; cd %s;wget http://www.mcs.anl.gov/research/projects/tao/download/tao-1.10.1.tar.gz  ;tar -zxf tao-1.10.1.tar.gz " % (tao_root, tao_root) )
    os.system("cd %s;make" % tao_dir)
elif (options.build_libmesh ):
    libMesh_dir = EnvironmentSetup("LIBMESH_DIR")
    # setup libMesh
    libMesh_root = libMesh_dir.split("/")
    libMesh_root.pop()
    libMesh_root = "/".join(libMesh_root)
    os.system("mkdir -p %s;svn co -r 4171  https://libmesh.svn.sourceforge.net/svnroot/libmesh/trunk/libmesh %s" % (libMesh_root,libMesh_dir) )
    #apply cxx patch to libmesh file Makefile.const
    #unwanted_itk_flags=-ansi -pedantic -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
    #libmesh_CXXFLAGS := $(filter-out $(unwanted_itk_flags),$(libmesh_CXXFLAGS))
    #libmesh_CXXFLAGS += -DMPICH_SKIP_MPICXX
    os.system("patch -d %s -p0 < patch_CXX_libMesh_0.7.3.1.diff" % libMesh_dir)
    os.system("cd %s;./configure --enable-everything --disable-tecplot --disable-vtk --disable-gzstreams --disable-triangle;make -j 6" % libMesh_dir)
elif (options.build_itk ):
    # setup ITK
    itk_home    = EnvironmentSetup("ITK_HOME")
    itk_version = EnvironmentSetup("ITK_VERSION")
    itk_root = itk_home.split("/")
    itk_root.pop()
    itk_root = "/".join(itk_root)
    os.system("mkdir -p %s; cd %s;wget http://voxel.dl.sourceforge.net/sourceforge/itk/%s.tar.gz;tar -zxf %s.tar.gz " % (itk_root, itk_root,itk_version ,itk_version ) )
    os.system("mkdir -p %s-build; cd %s-build;cmake -DBUILD_SHARED_LIBS=ON  -DCMAKE_BUILD_TYPE=%s -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=%s ../%s;make -j 4;make install" % (itk_home, itk_home, build_type, itk_home, itk_version ) )
elif (options.build_vtk ):
    # setup VTK
    vtk_home    = EnvironmentSetup("VTK_HOME")
    vtk_version = EnvironmentSetup("VTK_VERSION")
    vtk_root = vtk_home.split("/")
    vtk_root.pop()
    vtk_root = "/".join(vtk_root)
    import sys
    os.system( "mkdir -p %s/lib/python%d.%d/site-packages " % (vtk_home, sys.version_info[0], sys.version_info[1]) )
    os.system( "mkdir -p %s; cd %s; git clone git://vtk.org/VTK.git VTK ;cd VTK; git checkout v%s " % (vtk_root, vtk_root,vtk_version ) )
    os.system( "mkdir -p %s-build; cd %s-build;cmake -DBUILD_SHARED_LIBS=ON  -DCMAKE_BUILD_TYPE=%s -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DCMAKE_VERBOSE_MAKEFILE=ON -DVTK_WRAP_PYTHON=ON -DVTK_USE_PARALLEL=ON -DVTK_USE_N_WAY_ARRAYS=ON -DCMAKE_INSTALL_PREFIX=%s ../VTK;make -j 6;make install" % (vtk_home, vtk_home, build_type ,vtk_home ) )
    
elif (options.build_dakota ):
    # setup Dakota
    dakota_src = EnvironmentSetup("DAKOTA_SRC")
    os.system("patch -d %s -p0 < patch_CXX_Dakota_5_2.diff" % dakota_src )
    os.system("cd %s;./configure --without-graphics --prefix=%s-%s --enable-debug;make -j 6;make install" % (dakota_src,dakota_src,os.getenv(PETSC_ARCH) ) )
    #os.system("cd %s;./configure --without-graphics --prefix=%s-%s --with-lapack=$PETSC_DIR/$PETSC_ARCH/lib/libflapack.a --with-blas=$PETSC_DIR/$PETSC_ARCH/lib/libfblas.a --enable-debug;make -j 6;make install" % (dakota_src,dakota_src,os.getenv(PETSC_ARCH) ) )
else:
  parser.print_help()
  print options



