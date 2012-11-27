// libmesh include files
#include "equation_systems.h"
#include "exodusII_io.h"
#include "dof_map.h"
#include "mesh.h"

// local include
#include "pennesModel.h"
#include "petsc_fem_system.h"
#include "tttkUtilities.h"

/** @file pennesModel.cxx
 * 
 * Finite Difference Practical considerations
 * 
 * http://en.wikipedia.org/wiki/Numerical_differentiation
 * 
 * An important consideration in practice when the function is approximated
 * using floating point arithmetic is how small a value of h to choose. If
 * chosen too small, the subtraction will yield a large rounding error and in
 * fact all the finite difference formulae are ill-conditioned[2] and due to
 * cancellation will produce a value of zero if h is small enough[3]. If too
 * large, the calculation of the slope of the secant line will be more accurate,
 * but the estimate of the slope of the tangent by using the secant could be
 * worse.
 * 
 * As discussed in Chapter 5.7 of Numerical Recipes in C
 * (http://www.nrbook.com/a/bookcpdf/c5-7.pdf), a suitable choice for h is
 * sqrt(eps) * x where the machine epsilon eps is typically of the order
 * 2.2e-16.  Another important consideration is to make sure that h and x+h are
 * representable in floating point precision so that the difference between x+h
 * and x is exactly h. This can be accomplished by placing their values into and
 * out of memory as follows: h = sqrt(eps) * x, temp = x + h and h = temp - x.
 * It may be necessary to declare temp as a volatile variable to ensure the
 * steps are not undone by compiler optimization.
 */
const  PetscReal    fdEpsilon = PETSC_SQRT_MACHINE_EPSILON;
/* -------------------------------------------------------------------- 
   Pennes Bioheat Equation is treated as the base class with no source term
   -------------------------------------------------------------------- */ 
void PennesBioheatModel::printSelf(std::ostream& os)
{
  PetscFunctionBegin; 

  os << "PennesModel:               w_1=" <<         w_1       << std::endl;
  os << "PennesModel:               k_1=" <<         k_1       << std::endl;
  os << "PennesModel:               rho=" <<               rho << std::endl;
  os << "PennesModel:     specific_heat=" <<     specific_heat << std::endl;
  os << "PennesModel: bloodspecificheat=" << bloodspecificheat << std::endl;
  os << "PennesModel:          u_artery=" <<          u_artery << std::endl;
  os << "PennesModel:        m_bodyTemp=" <<        m_bodyTemp << std::endl;
  os << "PennesModel:     m_probeTemp="   <<       m_probeTemp << std::endl;
  os << "PennesModel:       w_0_dbeta="   <<       w_0.dbeta() << std::endl;
  os << "PennesModel:       w_0_verif="   <<       w_0.verif() << std::endl;
  os << "PennesModel:       k_0_dbeta="   <<       k_0.dbeta() << std::endl;
  os << "PennesModel:       k_0_verif="   <<       k_0.verif() << std::endl;
  os << "PennesModel:     Power_dbeta="   <<     Power.dbeta() << std::endl;
  os << "PennesModel:     Power_verif="   <<     Power.verif() << std::endl;
  os << "PennesModel:ApplicatorDomain="   <<m_ApplicatorDomain << std::endl;
  for ( PetscInt Ii = 0 ; Ii < this->get_num_elem_blk() ; Ii++)
   {
    os << "PennesModel:          w_0_lb("<<Ii<<")="<<      w_0.lb(Ii)<<std::endl;
    os << "PennesModel:          w_0_ub("<<Ii<<")="<<      w_0.ub(Ii)<<std::endl;
    os << "PennesModel:          k_0_lb("<<Ii<<")="<<      k_0.lb(Ii)<<std::endl;
    os << "PennesModel:          k_0_ub("<<Ii<<")="<<      k_0.ub(Ii)<<std::endl;
    os << "PennesModel:        Power_lb("<<Ii<<")="<<    Power.lb(Ii)<<std::endl;
    os << "PennesModel:        Power_ub("<<Ii<<")="<<    Power.ub(Ii)<<std::endl;
    os << "PennesModel:   BulkFluidFlow("<<Ii<<")="<<m_BulkFluidFlow[Ii]<<std::endl;
    os << "PennesModel:DiffusionDirection("<<Ii<<")="<<m_DiffusionDirection[Ii]<<std::endl;
   }
  //w_0.printStdVector(      os, "PennesModel:    w_0[");
  //k_0.printStdVector(      os, "PennesModel:    k_0[");
  PetscFunctionReturnVoid(); 
}

// constructor
PennesBioheatModel::PennesBioheatModel(
  GetPot &controlfile,EquationSystems &es ):
  _equation_systems(es), // store the pointer
  // thermal conductivity
  k_0( "k_0", es,
       controlfile("thermal_conductivity/k_0_optimize",  false),
       //&PDEModelBaseClass::dpde_dk_0,
       //&PDEModelBaseClass::d2pde_du_dk_0,
       controlfile("thermal_conductivity/k_0_lb"      ,    0.10e0),
       controlfile("thermal_conductivity/k_0_ub"      ,    0.70e0),
       controlfile("thermal_conductivity/k_0_dbeta"   , std::sqrt(fdEpsilon)),
       0.5e0,  /** @todo {determine if this is a bug in adjoint or fd calc} */
       controlfile("thermal_conductivity/k_0_vary"    ,   false) ),
  // perfusion
  w_0( "w_0", es,
       controlfile("perfusion/w_0_optimize",  false),
       //&PDEModelBaseClass::dpde_dw_0,
       //&PDEModelBaseClass::d2pde_du_dw_0,
       controlfile("perfusion/w_0_lb"      ,    0.10e0),
       controlfile("perfusion/w_0_ub"      ,   90.00e0),
       controlfile("perfusion/w_0_dbeta"   , fdEpsilon), 1.0e1,
       controlfile("perfusion/w_0_vary"    ,   false) ),
  // applicator parameters
  Power(     "Power",
             controlfile("probe/power_optimize",  false),
             //&PDEModelBaseClass::dpde_dpow,
             //&PDEModelBaseClass::d2pde_du_dm,
             controlfile("probe/power_lb"      ,   0.00e0),
             controlfile("probe/power_ub"      ,  30.00e0),
             controlfile("probe/power_dbeta"   ,fdEpsilon), 18.8e0)
{
  // density and specific heats
  rho= controlfile("material/rho",1.e3); // [kg/m^3]
  specific_heat= controlfile("material/specific_heat",3840.0); //[J/kg/C]

  // set initial condition parameters for temperature
  m_bodyTemp     = controlfile("initial_condition/u_init",37.0);//celcius
  m_probeTemp    = controlfile("initial_condition/probe_init",21.0); 

  // initial value in healthy tissue
  es.parameters.set<PetscScalar>(   "w_0_healthy") = 
             controlfile( "perfusion/w_0_healthy",0.0e0) ; //[kg/s/m^3]
  es.parameters.set<PetscScalar>(   "k_0_healthy") = 
   controlfile("thermal_conductivity/k_0_healthy",0.57e0) ; //[J/s/m/C]

  // initial value in applicator region
  es.parameters.set<PetscScalar>(   "w_0_probe") = 
             controlfile( "perfusion/w_0_probe",0.0e0) ; //[kg/s/m^3]
  es.parameters.set<PetscScalar>(   "k_0_probe") = 
   controlfile("thermal_conductivity/k_0_probe",0.627e0) ; //[J/s/m/C]

  // initial value in applicator region
  es.parameters.set<PetscScalar>(   "w_0_tumor") = 
             controlfile( "perfusion/w_0_tumor",
  es.parameters.get<PetscScalar>(   "w_0_healthy") ) ; //[kg/s/m^3]
  es.parameters.set<PetscScalar>(   "k_0_tumor") = 
   controlfile("thermal_conductivity/k_0_tumor",
  es.parameters.get<PetscScalar>(   "k_0_healthy") ) ; //[J/s/m/C] 

  // lower bounds
  w_0.lb( 1) = controlfile("perfusion/w_0_tumor_lb"           , w_0.lb(0));
  k_0.lb( 1) = controlfile("thermal_conductivity/k_0_tumor_lb", k_0.lb(0)); 

  // upper bounds
  w_0.ub( 1) = controlfile("perfusion/w_0_tumor_ub"           , w_0.ub(0)); 
  k_0.ub( 1) = controlfile("thermal_conductivity/k_0_tumor_ub", k_0.ub(0)); 

  // nonlinear perfusion parameters
  //  default to linear problem
  w_1= controlfile("perfusion/w_1",
       es.parameters.get<PetscScalar>("w_0_healthy") ) ; //[kg/s/m^3]

  // get perfusion parameters
  u_artery = controlfile("perfusion/u_artery",m_bodyTemp); // [C]
  bloodspecificheat=controlfile("perfusion/c_blood",3840.0);//[J/kg/C]

  // nonlinear thermal conductivity parameters
  k_1= controlfile("thermal_conductivity/k_1",0.0); // [J/s/m/C/C]

  m_ActivationEnergy     = controlfile("arrhenius/activationenergy",3.1e98) ; 
  m_FrequencyFactor      = controlfile("arrhenius/frequencyfactor" ,6.28e5) ; 
  m_GasConstant          = controlfile("arrhenius/gasconstant"     ,8.314472) ; 
  m_BaseTemperature      = controlfile("arrhenius/basetemperature" ,273.0) ; 

  // domains above the dirichletID will have dirichlet BC's applied
  m_ApplicatorDomain = controlfile("probe/applicatordomain",this->get_num_elem_blk()); 
  libmesh_assert(m_ApplicatorDomain > 0 ); 

  // setup power data
  this->GetPowerData( controlfile);
  //// error check
  //if( Power.size() < user->get_max_ideal()
  //                 || 
  //    Power.size() < user->get_num_MRTI()  )
  //  { 
  //    std::cout << "power file of incorrect dimension" << std::endl
  //              << "Power.size() = " << Power.size()  
  //              << " FiniteElementInterface::MAXIDEAL = "   << user->get_max_ideal()
  //              << " FiniteElementInterface::MRTI_ntime = " << user->get_num_MRTI()
  //              << std::endl; abort();
  //  }

  // clear
  m_BulkFluidFlow.clear();
  m_DiffusionDirection.clear();
  for (PetscInt Ii = 0 ; Ii < this->get_num_elem_blk();Ii++)
   {
    Real flowdir [3];
    Real diffdir [3];
   for (PetscInt Jj = 0 ; Jj < 3; Jj++)
    {
     // default to no fluid flow
     OStringStream controlIDKeyBulkFlow;
     controlIDKeyBulkFlow<< "bulkflow/domain_" << Ii << "_"  << Jj; 
     flowdir[Jj] = controlfile(controlIDKeyBulkFlow.str().c_str(),0.0) ;
     // default to normal diffusion in all directions
     OStringStream controlIDKeyDiffusion;
     controlIDKeyDiffusion<< "probe/diffusion_dir_" << Ii << "_"  << Jj; 
     diffdir[Jj] = controlfile(controlIDKeyDiffusion.str().c_str(),1.0) ;
    }
    Gradient FlowDirection(flowdir[0],flowdir[1],flowdir[2]); 
    Gradient DiffDirection(diffdir[0],diffdir[1],diffdir[2]); 
    m_BulkFluidFlow.push_back(FlowDirection);
    m_DiffusionDirection.push_back(DiffDirection);
   }

  // store a pointer to all field parameters for plotting
  _fieldParameters.push_back( &w_0  );  
  _fieldParameters.push_back( &k_0  );  
}
// utility to read the power file
void PennesBioheatModel::GetPowerData(GetPot &powerFile)
{
  PetscFunctionBegin; 
  // power parameters 
  // assume power file has already been preprocessed in a scripting language
  // and able to read in directly
  
  int powerDataSize = powerFile("power/nsize",-1);
  if( powerDataSize == -1)
   { 
     std::cout << "power file format error, size = " << powerDataSize
               << std::endl << std::flush; 
     abort();
   }
  Power.resize(powerDataSize ,0.0); // clear first
  for( PetscInt Ii = 0 ; Ii < Power.size(); Ii++)
  {
     std::ostringstream section_name;
     section_name << "power/power["<< Ii<<"]";
     if ( powerFile.have_variable (section_name.str().c_str())) 
        Power.at(Ii)=powerFile(section_name.str().c_str(),0.0);
     else
      {
        std::cout << "power file format error, "
                  << section_name.str()  << " Not Found !!!"
                  << std::endl << std::flush;  abort();
      }
  }
  Power.printStdVector(    std::cout, "PennesModel:  Power[");
  PetscFunctionReturnVoid(); 
}
/* -------------------------------------------------------------------- 
   Pennes Bioheat Equation with SDA laser is treated as the base class
   for the variety of laser source models
   -------------------------------------------------------------------- */ 
void PennesStandardDiffusionApproximation::printSelf(std::ostream& os)
{
  PetscFunctionBegin; 

  // print base class info
  this->PennesBioheatModel::printSelf(os); 

  os << "PennesSDA:            anfact=" <<        anfact     << std::endl;
  os << "PennesSDA:            mu_a_1=" <<        mu_a_1     << std::endl;
  os << "PennesSDA:            mu_s_1=" <<        mu_s_1     << std::endl;
  os << "PennesSDA:      mu_a_dbeta="   <<      mu_a_0.dbeta() << std::endl;
  os << "PennesSDA:      mu_a_verif="   <<      mu_a_0.verif() << std::endl;
  os << "PennesSDA:      mu_s_dbeta="   <<      mu_s_0.dbeta() << std::endl;
  os << "PennesSDA:      mu_s_verif="   <<      mu_s_0.verif() << std::endl;
  os << "PennesSDA:       x_0_dbeta="   <<       x_0.dbeta() << std::endl;
  os << "PennesSDA:       x_0_verif="   <<       x_0.verif() << std::endl;
  os << "PennesSDA:       y_0_dbeta="   <<       y_0.dbeta() << std::endl;
  os << "PennesSDA:       y_0_verif="   <<       y_0.verif() << std::endl;
  os << "PennesSDA:       z_0_dbeta="   <<       z_0.dbeta() << std::endl;
  os << "PennesSDA:       z_0_verif="   <<       z_0.verif() << std::endl;
  for ( PetscInt Ii = 0 ; Ii < this->get_num_elem_blk() ; Ii++)
   {
    os << "PennesSDA:         mu_a_lb("<<Ii<<")="<<     mu_a_0.lb(Ii)<<std::endl;
    os << "PennesSDA:         mu_a_ub("<<Ii<<")="<<     mu_a_0.ub(Ii)<<std::endl;
    os << "PennesSDA:         mu_s_lb("<<Ii<<")="<<     mu_s_0.lb(Ii)<<std::endl;
    os << "PennesSDA:         mu_s_ub("<<Ii<<")="<<     mu_s_0.ub(Ii)<<std::endl;
    os << "PennesSDA:          x_0_lb("<<Ii<<")="<<      x_0.lb(Ii)<<std::endl;
    os << "PennesSDA:          x_0_ub("<<Ii<<")="<<      x_0.ub(Ii)<<std::endl;
    os << "PennesSDA:          y_0_lb("<<Ii<<")="<<      y_0.lb(Ii)<<std::endl;
    os << "PennesSDA:          y_0_ub("<<Ii<<")="<<      y_0.ub(Ii)<<std::endl;
    os << "PennesSDA:          z_0_lb("<<Ii<<")="<<      z_0.lb(Ii)<<std::endl;
    os << "PennesSDA:          z_0_ub("<<Ii<<")="<<      z_0.ub(Ii)<<std::endl;
   }
  
  //mu_a_0.printStdVector(     os, "          mu_a_0[");
  //mu_s_0.printStdVector(     os, "          mu_s_0[");
  PetscFunctionReturnVoid(); 
}

// constructor
PennesStandardDiffusionApproximation::PennesStandardDiffusionApproximation(
  GetPot &controlfile,EquationSystems &es ) : PennesBioheatModel(controlfile,es),
  // optical laser parameters
  mu_a_0( "mu_a", es,
        controlfile("optical/mu_a_optimize"      ,  false),
        //&PDEModelBaseClass::dpde_dmu_a,
        //&PDEModelBaseClass::d2pde_du_dm,
        controlfile("optical/mu_a_lb"            ,   00.05e2),
        controlfile("optical/mu_a_ub"            ,   400.0e2),
        controlfile("optical/mu_a_dbeta"         , fdEpsilon), 10.0e0,
        controlfile("optical/mu_a_vary"          ,   false) ),
  mu_s_0( "mu_s", es,
        controlfile("optical/mu_s_optimize"      ,  false),
        //&PDEModelBaseClass::dpde_dmu_s,
        //&PDEModelBaseClass::d2pde_du_dm,
        controlfile("optical/mu_s_lb"            , 10.000e2),
        controlfile("optical/mu_s_ub"            , 3000.0e2),
        controlfile("optical/mu_s_dbeta"         ,fdEpsilon), 1400.0e0,
        controlfile("optical/mu_s_vary"          ,  false) ),
  // probe parameters
  x_0("x_0", controlfile("probe/x_0_optimize"      ,  false),
             //&PDEModelBaseClass::dpde_dx_0,
             //&PDEModelBaseClass::d2pde_du_dm,
             controlfile("probe/x_0_lb"            ,-1.00e4),
             controlfile("probe/x_0_ub"            , 1.00e4),
             controlfile("probe/x_0_dbeta"         ,fdEpsilon) , 0.04e0),
  y_0("y_0", controlfile("probe/y_0_optimize"      ,  false),
             //&PDEModelBaseClass::dpde_dy_0,
             //&PDEModelBaseClass::d2pde_du_dm,
             controlfile("probe/y_0_lb"            ,  -1.00e4),
             controlfile("probe/y_0_ub"            ,   1.00e4),
             controlfile("probe/y_0_dbeta"         ,fdEpsilon), 0.06e0),
  z_0("z_0", controlfile("probe/z_0_optimize"      ,  false),
             //&PDEModelBaseClass::dpde_dz_0,
             //&PDEModelBaseClass::d2pde_du_dm,
             controlfile("probe/z_0_lb"            ,  -1.00e4),
             controlfile("probe/z_0_ub"            ,   1.00e4),
             controlfile("probe/z_0_dbeta"         ,fdEpsilon), 0.07e0)
{
  // anisotropy factor
  anfact = controlfile("optical/anfact",0.9); // 

  // initial value in healthy tissue
  es.parameters.set<PetscScalar>("mu_a_healthy") = 
             controlfile("optical/mu_a_healthy",60.00e0) ; // [1/m]
  es.parameters.set<PetscScalar>("mu_s_healthy") = 
             controlfile("optical/mu_s_healthy",453.0e0) ; // [1/m]

  // initial value in applicator region
  es.parameters.set<PetscScalar>("mu_a_probe") = 
             controlfile("optical/mu_a_probe",
  es.parameters.get<PetscScalar>("mu_a_healthy") ) ; // [1/m]
  es.parameters.set<PetscScalar>("mu_s_probe") = 
             controlfile("optical/mu_s_probe",
  es.parameters.get<PetscScalar>("mu_s_healthy") ) ; // [1/m]

  // initial value in tumor tissue
  es.parameters.set<PetscScalar>("mu_a_tumor") = 
             controlfile("optical/mu_a_tumor",
  es.parameters.get<PetscScalar>("mu_a_healthy") ) ; // [1/m]
  es.parameters.set<PetscScalar>("mu_s_tumor") = 
             controlfile("optical/mu_s_tumor",
  es.parameters.get<PetscScalar>("mu_s_healthy") ) ; // [1/m]

  // nonlinear parameters
  mu_a_1 = controlfile("optical/mu_a_1" ,es.parameters.get<PetscScalar>("mu_a_healthy"));
  mu_s_1 = controlfile("optical/mu_s_1" ,es.parameters.get<PetscScalar>("mu_s_healthy"));

  // lower bounds
  mu_a_0.lb(1) = controlfile("optical/mu_a_tumor_lb" ,es.parameters.get<PetscScalar>("mu_a_healthy"));
  mu_s_0.lb(1) = controlfile("optical/mu_s_tumor_lb" ,es.parameters.get<PetscScalar>("mu_s_healthy"));

  // upper bounds
  mu_a_0.ub(1) = controlfile("optical/mu_a_tumor_ub"            ,mu_a_0.ub(0)); 
  mu_s_0.ub(1) = controlfile("optical/mu_s_tumor_ub"            ,mu_s_0.ub(0));

  // setup laser
  this->SetupLaser( controlfile, es, this->get_num_elem_blk() );

  // store a pointer to all field parameters for plotting
  _fieldParameters.push_back( &mu_a_0 );  
  _fieldParameters.push_back( &mu_s_0 );  
}

PetscTruth  PennesStandardDiffusionApproximation::CheckLinearity( const Parameters& parameters )
{
  PetscTruth LinearPDE = PETSC_TRUE; 
  // test for non linear parameters
  if ( mu_a_1 != parameters.get<PetscScalar>("mu_a_healthy") ) LinearPDE = PETSC_FALSE;
  if ( mu_s_1 != parameters.get<PetscScalar>("mu_s_healthy") ) LinearPDE = PETSC_FALSE;
  if ( k_1 != 0.0 )                                            LinearPDE = PETSC_FALSE;
  if ( w_1 != parameters.get<PetscScalar>("w_0_healthy"))      LinearPDE = PETSC_FALSE;
  // option to overwrite
  PetscTruth  nonlinearSolve=PETSC_FALSE;
  PetscOptionsGetTruth(PETSC_NULL,"-fem_linear_solve",&nonlinearSolve,PETSC_NULL);
  if ( nonlinearSolve ) LinearPDE = PETSC_FALSE;

  PetscFunctionReturn( LinearPDE ); 
}
PetscErrorCode PennesStandardDiffusionApproximation::SetupLaser( const GetPot &controlfile, 
                                                                 EquationSystems &es,
                                                                 const PetscInt NumberElementBlocks) 
{
  m_ProbeDomain = controlfile("probe/domain",NumberElementBlocks);

  // Loop over all the elements in the mesh. Store the element centroids and
  libMesh::MeshBase &mesh = es.get_mesh(); // get mesh data
  // volume fractions for each element in the probe domain element.
  libMesh::MeshBase::const_element_iterator el     = mesh.elements_begin();
  const libMesh::MeshBase::const_element_iterator el_end = mesh.elements_end();

  // default to no dirichlet nodes
  m_diffusingradius = 0.0;
  m_diffusinglength = 0.0;
  m_nelemtip        = 0;

  if( m_ProbeDomain >= 0 && m_ProbeDomain < NumberElementBlocks)
   {
     std::cout << "   Setting up WFS Model" << std::endl << std::endl ;
     // reset data structures
     x_0.clear(); y_0.clear(); z_0.clear(); 
     PetscScalar volumeTotal = 0.0;
     for ( ; el != el_end; ++el)
      {
       // Store a pointer to the element we are currently
       // working on.  This allows for nicer syntax later.
       Elem* elem = *el;

       // get centroid and volume
       const unsigned int subdomain_id = elem->subdomain_id() - 1;
       if( subdomain_id == (unsigned int) m_ProbeDomain ) 
         {
           // get centroid
           x_0.push_back(elem->centroid()(0));
           y_0.push_back(elem->centroid()(1));
           z_0.push_back(elem->centroid()(2));

           // get volume
           PetscScalar elemVolume = elem->volume();
           volumeTotal += elemVolume;
           m_volumeFraction.push_back(elemVolume);
           m_elementProbe.push_back(elem);
         }
      }
     // multiply each by the volume fraction
     for( unsigned int Ii = 0 ; Ii < m_volumeFraction.size(); Ii++)
                m_volumeFraction[Ii] = m_volumeFraction[Ii] / volumeTotal;
   }
  else
   { // default is No "m_ProbeDomain"
     m_nelemtip = controlfile("probe/nelemtip",1);// default to 1
     x_0.resize(m_nelemtip, controlfile("probe/x_0",0.0) ); // [m]
     y_0.resize(m_nelemtip, controlfile("probe/y_0",0.0) ); // [m]
     z_0.resize(m_nelemtip, controlfile("probe/z_0",0.0) ); // [m]
     if( m_nelemtip == 1 )
      {
       std::cout << "   Setting up SDA Model" << std::endl << std::endl ;
       m_volumeFraction.resize(1,1.0);
      }
     else if( m_nelemtip > 1 )
      {
       m_diffusinglength = controlfile("probe/length",0.01);// default to 1.0cm
       m_diffusingradius = controlfile("probe/radius",0.00075);// default to diameter 1.5mm
       PetscScalar X_1 = controlfile("probe/x_1", x_0[0]    ), // default to y-direction
                   Y_1 = controlfile("probe/y_1", y_0[0]+1.0), // default to y-direction
                   Z_1 = controlfile("probe/z_1", z_0[0]    ); // default to y-direction
       this->UpdateLaserPosition(x_0[0],y_0[0],z_0[0],X_1,Y_1,Z_1);
      }
     else 
      {
       std::cout << "   error input laser model..." << std::endl << std::endl ;
       abort();
      }
   } 

  PetscFunctionReturn( 0); 
}

// update laser position
PetscErrorCode PennesStandardDiffusionApproximation::UpdateLaserPosition(
                           PetscScalar X0,PetscScalar Y0, PetscScalar Z0, 
                           PetscScalar X1,PetscScalar Y1, PetscScalar Z1) 
{
  PetscFunctionBegin; 
  std::cout << "   Setting up line source WFS " << std::endl << std::endl;
  PetscScalar elemLength = m_diffusinglength / m_nelemtip ;
  // compute direction vector
  // (X_1,Y_1,Z_1) - (X_0,Y_0,Z_0) 
  //  should point towards laser tip as the reference
  PetscScalar differenceNorm = std::sqrt(
                                     std::pow(X1 - X0,2)   +
                                     std::pow(Y1 - Y0,2)   +
                                     std::pow(Z1 - Z0,2)   
                                        );
  m_unitVec[0] = (X1 - X0) / differenceNorm;
  m_unitVec[1] = (Y1 - Y0) / differenceNorm;
  m_unitVec[2] = (Z1 - Z0) / differenceNorm;
  // compute laser starting point 
  // (xhat_dir,yhat_dir,zhat_dir) should point away from laser tip
  m_laserTip[0] = X0 - m_unitVec[0]*0.5*m_diffusinglength; 
  m_laserTip[1] = Y0 - m_unitVec[1]*0.5*m_diffusinglength; 
  m_laserTip[2] = Z0 - m_unitVec[2]*0.5*m_diffusinglength; 

  /**
   *   put points at centroids of elements along the tip
   *    |--x--|--x--|--x--|--x--|--x--|
   */
  for ( PetscInt Ii = 0 ; Ii < m_nelemtip ; Ii++  )
    {
      x_0[Ii]=m_laserTip[0]+m_unitVec[0]*(Ii+0.5) * elemLength;
      y_0[Ii]=m_laserTip[1]+m_unitVec[1]*(Ii+0.5) * elemLength;
      z_0[Ii]=m_laserTip[2]+m_unitVec[2]*(Ii+0.5) * elemLength;
    }
  m_volumeFraction.resize(m_nelemtip,elemLength/m_diffusinglength );
  std::cout << "PennesSDA:      m_laserTip[0] "   <<    m_laserTip[0] << std::endl;
  std::cout << "PennesSDA:      m_laserTip[1] "   <<    m_laserTip[1] << std::endl;
  std::cout << "PennesSDA:      m_laserTip[2] "   <<    m_laserTip[2] << std::endl;
  std::cout << "PennesSDA:      m_unitVec[0]  "   <<    m_unitVec[0]  << std::endl;
  std::cout << "PennesSDA:      m_unitVec[1]  "   <<    m_unitVec[1]  << std::endl;
  std::cout << "PennesSDA:      m_unitVec[2]  "   <<    m_unitVec[2]  << std::endl;
  printStdVector(       std::cout, "volumeFraction[", m_volumeFraction);
  x_0.printStdVector(   std::cout, "           x_0[");
  y_0.printStdVector(   std::cout, "           y_0[");
  z_0.printStdVector(   std::cout, "           z_0[");
  PetscFunctionReturn( 0); 
}
bool PennesStandardDiffusionApproximation::dirichletNodalData(const  Point &NODE )
{
  PetscFunctionBegin; 
  PetscScalar locposRefTip[3]=
                    {NODE(0) - m_laserTip[0],
                     NODE(1) - m_laserTip[1],
                     NODE(2) - m_laserTip[2]};
  PetscScalar axialComponent =
                   locposRefTip[0]* m_unitVec[0] +
                   locposRefTip[1]* m_unitVec[1] +
                   locposRefTip[2]* m_unitVec[2] ;
  PetscScalar locRadius = m_diffusingradius; 
  if( axialComponent > 0.0)
    { // point  is along the fiber tract
      PetscScalar radialVec[3]=
          {locposRefTip[0] - axialComponent * m_unitVec[0],
           locposRefTip[1] - axialComponent * m_unitVec[1],
           locposRefTip[2] - axialComponent * m_unitVec[2]};
      // check if satisfy the radius criteria
      locRadius = std::sqrt( std::pow(radialVec[0],2)   +
                             std::pow(radialVec[1],2)   +
                             std::pow(radialVec[2],2)   
                           );
    }
  PetscFunctionReturn( locRadius < m_diffusingradius ); 
}
void PennesVoltage::printSelf(std::ostream& os)
{
  PetscFunctionBegin; 
  // print base class info
  this->PennesBioheatModel::printSelf(os); 
  //m_ElectricConductivity_0.printStdVector(      os, "RFA.m_ElectricConductivity_0=");
  os << "RFA:      s_0_dbeta="   <<  m_ElectricConductivity_0.dbeta() << std::endl;
  os << "RFA:      s_0_verif="   <<  m_ElectricConductivity_0.verif() << std::endl;
  os << "RFA.m_ElectricConductivity_1="<<m_ElectricConductivity_1 << std::endl;
  PetscFunctionReturnVoid(); 
}


PennesVoltage::PennesVoltage( GetPot &controlfile,EquationSystems &es) : 
                PennesBioheatModel(controlfile,es),
  // electric conductivity
  m_ElectricConductivity_0( "s_0", es,
       controlfile("electric_conductivity/s_0_optimize",  false),
       //&PDEModelBaseClass::dpde_ds_0,
       //&PDEModelBaseClass::d2pde_du_ds_0,
       controlfile("electric_conductivity/s_0_lb"      ,    0.10e0),
       controlfile("electric_conductivity/s_0_ub"      ,    0.70e0),
       controlfile("electric_conductivity/s_0_dbeta"   , std::sqrt(fdEpsilon)),
       0.5e0,  /** @todo {determine if this is a bug in adjoint or fd calc} */
       controlfile("electric_conductivity/s_0_vary"    ,   false) )
{
  PetscFunctionBegin; 

  es.parameters.set<PetscScalar>(   "s_0_healthy") = 
  controlfile("electric_conductivity/s_0_healthy",0.69e0) ; //[S/m]
  es.parameters.set<PetscScalar>(   "s_0_probe") = 
  controlfile("electric_conductivity/s_0_probe",1.0e8) ; //[S/m]
  es.parameters.set<PetscScalar>(   "s_0_tumor") = 
  controlfile("electric_conductivity/s_0_tumor",
  es.parameters.get<PetscScalar>(   "s_0_healthy") ) ; //[S/m] 

  m_ElectricConductivity_1= controlfile("electric_conductivity/s_1",0.0); // [S/m/C]        

  m_InitialVoltage.resize( this->get_num_elem_blk(),
                           controlfile("initial_condition/volt_init",0.0) );

  // store a pointer to all field parameters for plotting
  _fieldParameters.push_back( &m_ElectricConductivity_0);  

  PetscFunctionReturnVoid(); 
}

/* -------------------------------------------------------------------- 
   Delta P1 model
   -------------------------------------------------------------------- */ 
void PennesDeltaP1::printSelf(std::ostream& os)
{
  PennesStandardDiffusionApproximation::printSelf(os);
  os << "PennesDeltaP1: gs=" << m_ScatteringAsymmetry
     <<" f= " << m_ForwardScatterFraction
     <<" m_GuassBeamRadius     " << m_GuassBeamRadius     
     <<" m_SpecularReflectance=" << m_SpecularReflectance 
     <<" m_ReflectionFactor   =" << m_ReflectionFactor  
     << std::endl 
     <<" m_Agar_mu_t_star=    " << m_Agar_mu_t_star
     <<" m_AgarLength    =    " << m_AgarLength    
     << std::endl;
}
PennesDeltaP1::PennesDeltaP1(
   GetPot &controlfile,EquationSystems &es ) : PennesStandardDiffusionApproximation(controlfile,es)
{
  m_ScatteringAsymmetry    = anfact /(anfact + 1.0);
  m_ForwardScatterFraction = anfact * anfact; 
  m_GuassBeamRadius        = controlfile("optical/guass_radius", 0.01e0    );  
  m_SpecularReflectance    = controlfile("optical/specular_reflectance",0.0);  
  PetscScalar RefractionIndex = controlfile("optical/refactive_index", 1.4 );  
  m_ReflectionFactor = -0.13755 * std::pow(RefractionIndex,3)
                       +4.33900 * std::pow(RefractionIndex,2)
                       -4.90366 *          RefractionIndex
                       +1.68960;

  // set up laser position orientation
  m_laserTip[0] =  controlfile("probe/x_0",0.0) ; // [m]
  m_laserTip[1] =  controlfile("probe/y_0",0.0) ; // [m]
  m_laserTip[2] =  controlfile("probe/z_0",0.0) ; // [m]
  m_unitVec[0] = controlfile("probe/x_orientation",0.0) ;
  m_unitVec[1] = controlfile("probe/y_orientation",0.0) ;
  m_unitVec[2] = controlfile("probe/z_orientation",1.0) ;

  // for modeling the separate regions and no attenuation through the agar
  m_Agar_mu_t_star = es.parameters.get<PetscScalar>("mu_a_healthy") + 
                       (1.0 - m_ForwardScatterFraction) * 
                     es.parameters.get<PetscScalar>("mu_s_healthy") ; 
  m_AgarLength     = controlfile("optical/agar_length",0.02) ;
}

/** create attenuation map */
PetscScalar PennesDeltaP1::ComputeRadialLoss( const Point &Location )
{
  // radial loss
  PetscScalar locposRefTip[3]=
                    {Location(0) - m_laserTip[0],
                     Location(1) - m_laserTip[1],
                     Location(2) - m_laserTip[2]};
  PetscScalar axialComponent =
                   locposRefTip[0]* m_unitVec[0] +
                   locposRefTip[1]* m_unitVec[1] +
                   locposRefTip[2]* m_unitVec[2] ;
  PetscScalar radialVec[3]=
      {locposRefTip[0] - axialComponent * m_unitVec[0],
       locposRefTip[1] - axialComponent * m_unitVec[1],
       locposRefTip[2] - axialComponent * m_unitVec[2]};
  return std::sqrt( std::pow(radialVec[0],2)   +
                         std::pow(radialVec[1],2)   +
                         std::pow(radialVec[2],2)   
                   );
}
/** create attenuation map */
PetscScalar PennesDeltaP1::ComputeTotalAttenuation( 
                    const Elem* sourceElem , const Elem* targetElem) 
{
  PetscScalar AttenuationTotal = 0.0;
  // source is the same as the target
  if( sourceElem == targetElem) return AttenuationTotal ;

  // default to the boundary case, ie sourceElem == NULL
  // we are searching to the boundary for a collimated laser
  // input the projection direction from the control file
  Point ProjectionDir( m_unitVec[0], m_unitVec[1], m_unitVec[2]) ;
  // if not looking to a boundary
  // compute the length between the two centroids
  if (sourceElem) {ProjectionDir=targetElem->centroid()-sourceElem->centroid();}

  // ensure that the projection direction is normalized
  // we take the negative because we are search in the reverse
  // direction. ie from the target to the source
  const Point unitProjectionDir = -ProjectionDir.unit(); 

  // keep track of computed length
  PetscScalar approxTotalLength = 0.0;

  // loop until we find the element in question
  // NOTE that the source element may be a NULL
  // indicating that we are search until we reach 
  // the boundary
  Elem *currentElem = const_cast<Elem *>(targetElem);
  libmesh_assert(currentElem != 0 ); 
  while(currentElem != sourceElem ) 
   {
    Point currDistance(ProjectionDir);
    if ( sourceElem )
      { // pointing from source to current element
        currDistance =  currentElem->centroid()- sourceElem->centroid();
      }
    // point from current element to source
    const Point unitcurrDistance = -currDistance.unit();
    PetscScalar maxInnerProduct = 0.0;
    Elem* nextElem = 0;
    Point nextDistance ;
    // loop over all neighbors and find the element most aligned with the
    // direction sought
    for (unsigned int iNeighbor = 0 ; 
                      iNeighbor < currentElem->n_neighbors();iNeighbor++) 
     {
        Elem* neighborElem = currentElem->neighbor (iNeighbor) ;
        Point tmpDistance;
        if( neighborElem ) 
          { // check not a boundary neighbor
           tmpDistance = neighborElem->centroid() - currentElem->centroid();
          }
        else if( neighborElem == NULL ) 
          { // check not a boundary neighbor
           AutoPtr<Elem> boundaryNeighbor = currentElem->build_side(iNeighbor);
           tmpDistance = boundaryNeighbor->centroid() - currentElem->centroid();
          }
        Point unitTmpDistance = tmpDistance.unit() ; 
        // overloaded to dot product
        PetscScalar InnerProduct = unitTmpDistance*unitcurrDistance;
        // take the maximum innter product as the next attenuation
        // direction
        if ( InnerProduct > maxInnerProduct )
          {
            maxInnerProduct = InnerProduct;
            nextDistance(0) = tmpDistance(0) ;
            nextDistance(1) = tmpDistance(1) ;
            nextDistance(2) = tmpDistance(2) ;
            nextElem = neighborElem ; 
          }
     }
    // error check
    if (maxInnerProduct == 0.0 )
      { 
        PetscPrintf(PETSC_COMM_WORLD, "logic error... next elem not set\n");
          PetscPrintf(PETSC_COMM_WORLD, "current id  %d centroid (%22.15e,%22.15e,%22.15e) \n",                 
                   currentElem->id(),
                   currentElem->centroid()(0),
                   currentElem->centroid()(1),
                   currentElem->centroid()(2));
        PetscPrintf(PETSC_COMM_WORLD, "unit current distance     (%22.15e,%22.15e,%22.15e) \n", 
                    unitcurrDistance(0) , unitcurrDistance(1), unitcurrDistance(2));
        PetscPrintf(PETSC_COMM_WORLD, "unit projection direction (%22.15e,%22.15e,%22.15e) \n", 
                    unitProjectionDir(0),unitProjectionDir(1),unitProjectionDir(2));
        if ( sourceElem )
          PetscPrintf(PETSC_COMM_WORLD, "source id   %d centroid (%22.15e,%22.15e,%22.15e) \n",                 
                     sourceElem->id(),
                     sourceElem->centroid()(0),
                     sourceElem->centroid()(1),
                     sourceElem->centroid()(2));
        for (unsigned int iNeighbor = 0 ; 
                          iNeighbor < currentElem->n_neighbors();iNeighbor++) 
         {
            Elem* neighborElem = currentElem->neighbor (iNeighbor) ;
            if( neighborElem ) 
              { // check not a boundary neighbor
               PetscPrintf(PETSC_COMM_WORLD, "neighbor(%d) %d centroid (%22.15e,%22.15e,%22.15e) \n",
                                              iNeighbor,neighborElem->id(),
                                                        neighborElem->centroid()(0),
                                                        neighborElem->centroid()(1),
                                                        neighborElem->centroid()(2));
               Point tmpDistance = neighborElem->centroid() - currentElem->centroid();
               Point unitTmpDistance = tmpDistance.unit() ; 
               // overloaded to dot product
               PetscScalar InnerProduct = unitTmpDistance*unitcurrDistance;
               PetscPrintf(PETSC_COMM_WORLD, "unit neighbor   direction (%22.15e,%22.15e,%22.15e) \n", 
                           unitTmpDistance(0), unitTmpDistance(1), unitTmpDistance(2));
               PetscPrintf(PETSC_COMM_WORLD, "neighbor(%d) IP  %22.15e \n", iNeighbor, InnerProduct );
              }
         }
        libmesh_error();
      }
    //PetscScalar approxElemDiameter = (nextElem->hmax() + nextElem->hmin()) * 0.5;
    // overloaded to dot product
    PetscScalar DistanceProjection = nextDistance*unitProjectionDir; 
    approxTotalLength += DistanceProjection ;

    // get the parameter id
    std::vector<unsigned int> param_dof_indices;
    libMesh::System &template_parameter_system = 
                     this->get_equation_systems().get_system("k_0");
    if(nextElem)
     {//check not at the boundary
      template_parameter_system.get_dof_map().dof_indices (nextElem, param_dof_indices);
     }
    else
     {// if at the boundary use the current element
      template_parameter_system.get_dof_map().dof_indices (currentElem , param_dof_indices);
     }
    // current_solution get the global solution
    // the global solution should already be scatter to a local vector
    // with the same index ordering
    const unsigned int field_id = param_dof_indices[0];
    PetscScalar absorption = mu_a_0.GetGlobalSolution(field_id);
    PetscScalar scattering = mu_s_0.GetGlobalSolution(field_id);
    Real mu_t_star=absorption+scattering*(1.0e0-m_ForwardScatterFraction);

    // update attenuation and move to next element
    AttenuationTotal += mu_t_star * DistanceProjection  ;
    currentElem = nextElem ;
   }
  // error check 
  PetscScalar totalLength = ProjectionDir.size();
  PetscScalar DistanceError = std::abs(totalLength - approxTotalLength );
  if (sourceElem && DistanceError > 1.e-6 ) 
    {
     PetscPrintf(PETSC_COMM_WORLD,"distance between point (%22.15e,%22.15e,%22.15e) & (%22.15e,%22.15e,%22.15e)  %22.15e approximate distance %22.15e error %22.15e \n", 
  targetElem->centroid()(0),targetElem->centroid()(1),targetElem->centroid()(2),
  sourceElem->centroid()(0),sourceElem->centroid()(1),sourceElem->centroid()(2),
          totalLength,approxTotalLength,DistanceError ); 
    }
  return AttenuationTotal ;
}
/** primary fluence irradiance */
PetscScalar PennesDeltaP1::getInterstitialIrradiance(unsigned int,unsigned int,
                                 const Point &qpoint, const Parameters &parameters )
{
  PetscScalar source = 0.0, unitPower = 1.0;
  const Elem* elem = parameters.get<const Elem*>("elem") ;
  // loop over elements in probe domain
  std::vector<Elem*>::iterator probeElemIter;
  for( probeElemIter  = m_elementProbe.begin();
       probeElemIter != m_elementProbe.end()  ; probeElemIter++)
   {
    PetscInt Ii = std::distance(m_elementProbe.begin(), probeElemIter);
    Elem* currentSource = *probeElemIter;
    Point totalDistance = qpoint - currentSource->centroid();
    PetscScalar totalLength = totalDistance.size();  
    //compute total attenuation
    PetscScalar TotalAttenuation =  this->ComputeTotalAttenuation(currentSource,elem); 
    source  += 0.25 * unitPower * m_volumeFraction[Ii] * 
               std::exp(-TotalAttenuation)/(libMesh::pi*totalLength*totalLength);
   } // end loop over probe elements
               
  return source;
}
/** primary fluence irradiance */
PetscScalar PennesDeltaP1::getInterstitialFlux_X(unsigned int,unsigned int,
                                 const Point &qpoint, const Parameters &parameters )
{
  PetscScalar source = 0.0, unitPower = 1.0;
  const Elem* elem = parameters.get<const Elem*>("elem") ;
  // loop over elements in probe domain
  std::vector<Elem*>::iterator probeElemIter;
  for( probeElemIter  = m_elementProbe.begin();
       probeElemIter != m_elementProbe.end()  ; probeElemIter++)
   {
    PetscInt Ii = std::distance(m_elementProbe.begin(), probeElemIter);
    Elem* currentSource = *probeElemIter;
    Point totalDistance = qpoint - currentSource->centroid();
    PetscScalar totalLength = totalDistance.size();  
    //compute total attenuation
    PetscScalar TotalAttenuation =  this->ComputeTotalAttenuation(currentSource,elem); 
    source  += totalDistance(0)/totalLength * 
               0.25 * unitPower * m_volumeFraction[Ii] * 
               std::exp(-TotalAttenuation)/(libMesh::pi*totalLength*totalLength);
   } // end loop over probe elements
               
  return source;
}
/** primary fluence irradiance */
PetscScalar PennesDeltaP1::getInterstitialFlux_Y(unsigned int,unsigned int,
                                 const Point &qpoint, const Parameters &parameters )
{
  PetscScalar source = 0.0, unitPower = 1.0;
  const Elem* elem = parameters.get<const Elem*>("elem") ;
  // loop over elements in probe domain
  std::vector<Elem*>::iterator probeElemIter;
  for( probeElemIter  = m_elementProbe.begin();
       probeElemIter != m_elementProbe.end()  ; probeElemIter++)
   {
    PetscInt Ii = std::distance(m_elementProbe.begin(), probeElemIter);
    Elem* currentSource = *probeElemIter;
    Point totalDistance = qpoint - currentSource->centroid();
    PetscScalar totalLength = totalDistance.size();  
    //compute total attenuation
    PetscScalar TotalAttenuation =  this->ComputeTotalAttenuation(currentSource,elem); 
    source  += totalDistance(1)/totalLength * 
               0.25 * unitPower * m_volumeFraction[Ii] * 
               std::exp(-TotalAttenuation)/(libMesh::pi*totalLength*totalLength);
   } // end loop over probe elements
               
  return source;
}
/** primary fluence irradiance */
PetscScalar PennesDeltaP1::getInterstitialFlux_Z(unsigned int,unsigned int,
                                 const Point &qpoint, const Parameters &parameters )
{
  PetscScalar source = 0.0, unitPower = 1.0;
  const Elem* elem = parameters.get<const Elem*>("elem") ;
  // loop over elements in probe domain
  std::vector<Elem*>::iterator probeElemIter;
  for( probeElemIter  = m_elementProbe.begin();
       probeElemIter != m_elementProbe.end()  ; probeElemIter++)
   {
    PetscInt Ii = std::distance(m_elementProbe.begin(), probeElemIter);
    Elem* currentSource = *probeElemIter;
    Point totalDistance = qpoint - currentSource->centroid();
    PetscScalar totalLength = totalDistance.size();  
    //compute total attenuation
    PetscScalar TotalAttenuation =  this->ComputeTotalAttenuation(currentSource,elem); 
    source  += totalDistance(2)/totalLength * 
               0.25 * unitPower * m_volumeFraction[Ii] * 
               std::exp(-TotalAttenuation)/(libMesh::pi*totalLength*totalLength);
   } // end loop over probe elements
               
  return source;
}

/** primary fluence irradiance */
PetscScalar PennesDeltaP1::getExternalIrradiance(unsigned int,unsigned int,
                                 const Point &qpoint, const Parameters &parameters )
{
  // assume external fiber is collimated and only interested in a
  // certain direction
  Point unitTotalDistance(m_unitVec[0], m_unitVec[1], m_unitVec[2]); 
  const Elem* elem = parameters.get<const Elem*>("elem") ;
  // compute attenuation along collimation direction until reach boundary
  // element should be NULL at boundary
  PetscScalar TotalAttenuation =  this->ComputeTotalAttenuation(NULL,elem); 
  // Compute Loss in radial direction
  PetscScalar locRadius = this->ComputeRadialLoss( qpoint ) ;
  PetscScalar unitPower = 1.0;
  PetscScalar source = 2.0*unitPower/libMesh::pi/m_GuassBeamRadius/m_GuassBeamRadius 
         * (1.0 - m_SpecularReflectance) 
         * std::exp(-TotalAttenuation)
         * std::exp( -2.0 * locRadius * locRadius /
         m_GuassBeamRadius / m_GuassBeamRadius ) ;
  return source;
}
 
/** primary fluence flux */
PetscScalar PennesDeltaP1::getExternalFlux_X(unsigned int,unsigned int,
                                 const Point &qpoint, const Parameters &parameters )
{
  // assume external fiber is collimated and only interested in a
  // certain direction
  Point unitTotalDistance(m_unitVec[0], m_unitVec[1], m_unitVec[2]); 
  const Elem* elem = parameters.get<const Elem*>("elem") ;
  // compute attenuation along collimation direction until reach boundary
  // element should be NULL at boundary
  PetscScalar TotalAttenuation =  this->ComputeTotalAttenuation(NULL,elem); 
  // Compute Loss in radial direction
  PetscScalar locRadius = this->ComputeRadialLoss( qpoint ) ;
  PetscScalar unitPower = 1.0;
  PetscScalar source = 2.0*unitPower/libMesh::pi/m_GuassBeamRadius/m_GuassBeamRadius 
         * (1.0 - m_SpecularReflectance) 
         * std::exp(-TotalAttenuation)
         * std::exp( -2.0 * locRadius * locRadius /
         m_GuassBeamRadius / m_GuassBeamRadius ) ;
  return m_unitVec[0] * source;
}
 
/** primary fluence flux */
PetscScalar PennesDeltaP1::getExternalFlux_Y(unsigned int,unsigned int,
                                 const Point &qpoint, const Parameters &parameters )
{
  // assume external fiber is collimated and only interested in a
  // certain direction
  Point unitTotalDistance(m_unitVec[0], m_unitVec[1], m_unitVec[2]); 
  const Elem* elem = parameters.get<const Elem*>("elem") ;
  // compute attenuation along collimation direction until reach boundary
  // element should be NULL at boundary
  PetscScalar TotalAttenuation =  this->ComputeTotalAttenuation(NULL,elem); 
  // Compute Loss in radial direction
  PetscScalar locRadius = this->ComputeRadialLoss( qpoint ) ;
  PetscScalar unitPower = 1.0;
  PetscScalar source = 2.0*unitPower/libMesh::pi/m_GuassBeamRadius/m_GuassBeamRadius 
         * (1.0 - m_SpecularReflectance) 
         * std::exp(-TotalAttenuation)
         * std::exp( -2.0 * locRadius * locRadius /
         m_GuassBeamRadius / m_GuassBeamRadius ) ;
  return m_unitVec[1] * source;
}

/** primary fluence flux */
PetscScalar PennesDeltaP1::getExternalFlux_Z(unsigned int,unsigned int,
                                 const Point &qpoint, const Parameters &parameters )
{
  // assume external fiber is collimated and only interested in a
  // certain direction
  Point unitTotalDistance(m_unitVec[0], m_unitVec[1], m_unitVec[2]); 
  const Elem* elem = parameters.get<const Elem*>("elem") ;
  // compute attenuation along collimation direction until reach boundary
  // element should be NULL at boundary
  PetscScalar TotalAttenuation =  this->ComputeTotalAttenuation(NULL,elem); 
  // Compute Loss in radial direction
  PetscScalar locRadius = this->ComputeRadialLoss( qpoint ) ;
  PetscScalar unitPower = 1.0;
  PetscScalar source = 2.0*unitPower/libMesh::pi/m_GuassBeamRadius/m_GuassBeamRadius 
         * (1.0 - m_SpecularReflectance) 
         * std::exp(-TotalAttenuation)
         * std::exp( -2.0 * locRadius * locRadius /
         m_GuassBeamRadius / m_GuassBeamRadius ) ;
  return m_unitVec[2] * source;
}
