#include "petscvec.h"

// ITK include files
#include "itkImportImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

// DiffSystem framework files
#include "mesh.h"
#include "equation_systems.h"

// assume 3d
const unsigned int          m_Dimension = 3;

// base class for gradient models
class LinearGradientModel 
{
  public:
  LinearGradientModel(Point&,PetscScalar);
  /** GetPhase should enable functionoid capability */ 
  virtual inline PetscScalar GetPhase(const Point&,PetscScalar);
  protected:
  Point m_linearGradient;
  PetscScalar m_larmorFrequency;
};
// derived class for spiral gradient models
class SpiralGradientModel:public LinearGradientModel
{
  public:
  SpiralGradientModel(Point&,PetscScalar,PetscScalar);
  /** GetPhase should enable functionoid capability */ 
  virtual inline PetscScalar GetPhase(const Point&,PetscScalar);
  private:
  PetscScalar m_angularVelocity;
};
// derived class for epi gradient models
class EPIGradientModel:public LinearGradientModel
{
  public:
  EPIGradientModel(Point&,PetscScalar,int);
  /** GetPhase should enable functionoid capability */ 
  virtual inline PetscScalar GetPhase(const Point&,PetscScalar);
  private:
  int m_kSpaceSkip;
};

enum GradientModelType {LINEARGRADIENT  =0,       
	                SPIRALGRADIENT  =1,      
	                EPIGRADIENT     =2,
                        UNKNOWNGRADIENT  };
//namespace IterativeReconstruction
//{
/**
 * c++ class meant to interface with python
 */
  class SignalModelInterface
  {
  public:

   SignalModelInterface(); ///< constructor
   ~SignalModelInterface(); ///< destructor

   /** print local information */ 
   void printSelf ();

   /** setup libMesh */ 
   PetscErrorCode SetuplibMesh(MPI_Comm);

   /** setup Imaging Data Structures mesh for interpolating imaging data 
     * onto  FEM mesh
     */ 
   PetscErrorCode SetupImaging(int   , int   , int   ,
                                     double, double, double, 
                                     double, double, double );
   
   /** setup FEM mesh for libMesh */ 
   PetscErrorCode SetupStructuredGrid(int   , int   , int   ,
                                      double, double, 
                                      double, double, 
                                      double, double );
   
   /** Project Imaging onto libMesh data structures */ 
   PetscErrorCode ProjectImagingToFEMMesh (char *,Vec);

   /** setup libMesh data structures */ 
   PetscErrorCode AddExplicitSystem (char *,int);

   /** perform matvec product and compute the real and imaginary 
     * part of the signal at a given time instance. assume
     * all data is preload and push the work onto Python
     * TODO for higher compuational speed an ABC for the 
     * gradients can be provided in C++. 
     * This can be coded post development for efficiency.
     * THe current design allows the development to be performed
     * in python rather than C++ :)
     */ 
   void AssembleSignal (std::vector<double>*,
                        std::vector<double>*,
                        std::vector<double>*);

   /** initialize Equation Systems */ 
   PetscErrorCode InitializeEquationSystems();

   /** write libMesh data structures to view in paraview */ 
   PetscErrorCode WriteEquationSystems(char *);

   double m_echoTime; ///< echotime, TE

   PetscErrorCode SetupLinearGradientModel(
                double , double , double , double );
   PetscErrorCode SetupSpiralGradientModel(
                double , double , double , double ,double );
   PetscErrorCode SetupEPIGradientModel(
                double , double , double , double ,int);

   // class typedefs
   typedef PetscScalar                                         PixelType;
   typedef itk::ImportImageFilter< PixelType,m_Dimension >  ImportFilterType;
   typedef ImportFilterType::OutputImageType                   ImageType;
   typedef itk::LinearInterpolateImageFunction<
                                       ImageType, PixelType>  InterpolatorType;

  private:

   Mesh *m_mesh; ///< pointer for mesh object
   EquationSystems *m_eqn_systems; ///< pointer for equation systems object
   LibMeshInit *m_init;///< pointer for equation systems object

   GradientModelType m_InputGradientModel;///< control gradien model
   LinearGradientModel *m_PhaseModel;///< pointer for gradient models

   // ITK data structures
   ImportFilterType::Pointer m_importFilter;
   InterpolatorType::Pointer m_interpolator; 
   // This filter requires the user to specify the size of
   // the image to be produced as output.  The
   // \code{SetRegion()} method is used to this end.  The
   // image size should exactly match the number of pixels
   // available in the locally allocated buffer. 
   
   ImportFilterType::SizeType       m_size; ///< image dimensions
   ImportFilterType::SpacingType m_spacing; ///< image dimensions
   ImageType::PointType           m_origin; ///< image dimensions


  };
//}
