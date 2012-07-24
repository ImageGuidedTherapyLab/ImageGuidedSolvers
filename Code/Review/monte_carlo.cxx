#include "mri_params.h" // ensure same default values between C and Fortran
#include "parser_defaults.h" // ensure same default values between C and Fortran
#include "fortrandftranslate.h" // header to call fortran routines
#include "CImg.h"
#include "getpot.h" // ini file parser
#include "tao.h" // tao header
#include <vector>  // stl header
#include <numeric> // stl header
// libMesh
#include "libmesh.h"
#include "mesh.h"
#include "equation_systems.h"
#include "exact_solution.h"

using namespace cimg_library;
using namespace std;

#include "baseInfo.h" 
#include "applicationContext.h" 
#include "solvehp3d.h"
#include "pennesVerification.h" 
/* the grid for the monte carlo data should be

        
         -------------> R   radius varies by column (dimy)
         |
         |
         |
         |
       z\|/  
  
     axial direction varies by row (dimx)

  the MC grids are read in from a slight modified form of the CImg .asc file

the first line of the CImg file contains the grid size info

  dx dy dz dv

the modified file format contains the grid spacing in the first line

  dx dy dz dv deltaz deltar

*/ 
static CImg<PetscScalar>     MCdata_ideal     , MCdata     ;
// for grid derivitive data
static CImgList<PetscScalar> MCdata_ideal_grad, MCdata_grad;
static PetscScalar ROTMAT[3][3], // rotation matrix for transformation of grid
                   MC_dR, MC_dZ, // Monte Carlo Grid Spacing 
                   MC_dR_ideal, MC_dZ_ideal, // Monte Carlo Grid Spacing 
                   BEAM_LOC,BEAM_LOC_ideal; /* z-coordinate [meters] wrt to 
                                               mc grid of center of laser tip */
// function prototypes
extern FORTRAN_FUNCTION 
{
PetscScalar FORTRAN_NAME(mcheat_nopt)( PetscScalar *, PetscScalar *);
PetscScalar FORTRAN_NAME(mcheat_ideal)(PetscScalar *, PetscScalar *);
PetscScalar FORTRAN_NAME(dmcheat_noptdvector)(PetscScalar *,PetscScalar *, PetscInt * );
#include "pennes_model.h"
#include "global_params.h"
}

/* --------------------------------------------------------------------
   Initialize monte carlo data. MUST be called after init_control
 */
#undef __FUNCT__
#define __FUNCT__ "Init_MCdata"
PetscErrorCode Init_MCdata(PetscTruth Control_Task){
  PetscInt ndum; // dummy variable
  PetscFunctionBegin; // used for error handling

  GetPot controlfile(dfltINIFile);

  // get grid spacing
  /* Axial and Radial are used as reference points to obtain 
     the axial and radial direction associated w/ the monte carlo grid
  

     default coordinate system is of the form


       Origin
         -------------> Radial  
         |
         |          the remaining coordinate vector is out ot the
         |              plane to make a right handed coordinate system
         |
  Axial \|/  


  */
  vector<PetscScalar> Axial( 3,0.0), 
                      Radial(3,0.0), 
                      Origin(3,0.0); // origin of monte carlo grid
  FORTRAN_NAME(getparam_x_0)(&Origin[0],&ndum)  ;    
  FORTRAN_NAME(getparam_y_0)(&Origin[1],&ndum)  ;    
  FORTRAN_NAME(getparam_z_0)(&Origin[2],&ndum)  ;    
  // define default coordinate system
  Axial.at(0) =controlfile("source_laser/axial_x" ,Origin[0]         );
  Axial.at(1) =controlfile("source_laser/axial_y" ,Origin[1] + MC_dZ );
  Axial.at(2) =controlfile("source_laser/axial_z" ,Origin[2]         );
  Radial.at(0)=controlfile("source_laser/radial_x",Origin[0] + MC_dR );
  Radial.at(1)=controlfile("source_laser/radial_y",Origin[1]         );
  Radial.at(2)=controlfile("source_laser/radial_z",Origin[2]         );
  if(Control_Task)
  {
     cout << "\nInit_MCdata: AXIAL[0]  = " << Axial[0]  ;
     cout << "\nInit_MCdata: AXIAL[1]  = " << Axial[1]  ;
     cout << "\nInit_MCdata: AXIAL[2]  = " << Axial[2]  ;
     cout << "\nInit_MCdata: RADIAL[0] = " << Radial[0] ;
     cout << "\nInit_MCdata: RADIAL[1] = " << Radial[1] ;
     cout << "\nInit_MCdata: RADIAL[2] = " << Radial[2] ;
  }
/* direction cosine rotation matrix

   expect a no
   given a coordinate of a gauss point in the default coordinate system
    x0 = (1,0,0);
    x1 = (0,1,0);
    x2 = (0,0,1);

   we wish to transform this coordinate 
   given wrt to the coordinate system  (x1,x2,x3)
   to the natural coordinate system
   of the monte carlo grid in terms of (y1,y2,y3)

   note xi =    a1i  y1 +    a2i  y2 +    a3i  y3
           = (y1,xi) y1 + (y2,xi) y2 + (y3,xi) y3

      implies b = b1 x1 + b2 x2 + b3 x3
                =   b1 ( a11 y1 + a21 y2 + a31 y3 )
                  + b2 ( a12 y1 + a22 y2 + a32 y3 )
                  + b3 ( a13 y1 + a23 y2 + a33 y3 )
                =   [ a11 , a12 , a13 ]   [b1]
                    [ a21 , a22 , a23 ] * [b2]
                    [ a31 , a32 , a33 ]   [b3]

    where aij = (yi,xj)


*/
  vector<PetscScalar> x0(3,0.0), x1(3,0.0), x2(3,0.0),
                      y0(3,0.0), y1(3,0.0), y2(3,0.0);

  x0.at(0)= 1.0 ; x0.at(1)= 0.0 ; x0.at(2)= 0.0 ;
  x1.at(0)= 0.0 ; x1.at(1)= 1.0 ; x1.at(2)= 0.0 ;
  x2.at(0)= 0.0 ; x2.at(1)= 0.0 ; x2.at(2)= 1.0 ;

  y2.at(0) = Axial.at(0)  - Origin.at(0);
  y2.at(1) = Axial.at(1)  - Origin.at(1);
  y2.at(2) = Axial.at(2)  - Origin.at(2);
  
  y0.at(0) = Radial.at(0) - Origin.at(0);
  y0.at(1) = Radial.at(1) - Origin.at(1);
  y0.at(2) = Radial.at(2) - Origin.at(2);
  
  // y1 = y2 x y0 is perpendicular to plane formed by y0 and y2
  FORTRAN_NAME(cross_product)(&y2[0],&y0[0],&y1[0]);

  // create an orthogonal coordinate system
  FORTRAN_NAME(cross_product)(&y1[0],&y2[0],&y0[0]);

  // normalize
  int iii ;
  PetscScalar sqnorm ;
  sqnorm = sqrt(inner_product(y0.begin(),y0.end(),y0.begin(),0.0)); 
  for (iii=0;iii<3;iii++) y0.at(iii) = y0.at(iii)/sqnorm;
  sqnorm = sqrt(inner_product(y1.begin(),y1.end(),y1.begin(),0.0)); 
  for (iii=0;iii<3;iii++) y1.at(iii) = y1.at(iii)/sqnorm;
  sqnorm = sqrt(inner_product(y2.begin(),y2.end(),y2.begin(),0.0)); 
  for (iii=0;iii<3;iii++) y2.at(iii) = y2.at(iii)/sqnorm;

  // check that have an ON coordinate system
  const PetscScalar tol = 1.e-6;
  if(inner_product(y0.begin(),y0.end(),y1.begin(),0.0) > tol || 
     inner_product(y0.begin(),y0.end(),y2.begin(),0.0) > tol || 
     inner_product(y1.begin(),y1.end(),y2.begin(),0.0) > tol     ){
     printf("(y0,y1) = %e\n",inner_product(y0.begin(),y0.end(),y1.begin(),0.0));
     printf("(y0,y2) = %e\n",inner_product(y0.begin(),y0.end(),y2.begin(),0.0));
     printf("(y1,y2) = %e\n",inner_product(y1.begin(),y1.end(),y2.begin(),0.0));
     printf("non-orthogonal coordinate system") ; abort();
  }
  

  ROTMAT[0][0]=inner_product(y0.begin(),y0.end(),x0.begin(),0.0); 
  ROTMAT[0][1]=inner_product(y0.begin(),y0.end(),x1.begin(),0.0); 
  ROTMAT[0][2]=inner_product(y0.begin(),y0.end(),x2.begin(),0.0);
  ROTMAT[1][0]=inner_product(y1.begin(),y1.end(),x0.begin(),0.0); 
  ROTMAT[1][1]=inner_product(y1.begin(),y1.end(),x1.begin(),0.0); 
  ROTMAT[1][2]=inner_product(y1.begin(),y1.end(),x2.begin(),0.0);
  ROTMAT[2][0]=inner_product(y2.begin(),y2.end(),x0.begin(),0.0); 
  ROTMAT[2][1]=inner_product(y2.begin(),y2.end(),x1.begin(),0.0); 
  ROTMAT[2][2]=inner_product(y2.begin(),y2.end(),x2.begin(),0.0);
  if(Control_Task){
     printf("ROTMAT=[%e %e %e]\n       [%e %e %e]\n       [%e %e %e]\n",
                                 ROTMAT[0][0],ROTMAT[0][1],ROTMAT[0][2],
                                 ROTMAT[1][0],ROTMAT[1][1],ROTMAT[1][2],
                                 ROTMAT[2][0],ROTMAT[2][1],ROTMAT[2][2]);
  }

  PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Load_MCdataIdeal"
PetscErrorCode Load_MCdataIdeal(PetscTruth Control_Task){
  char filename[MAXLEN];
  PetscFunctionBegin; // used for error handling

  GetPot controlfile(dfltINIFile);
 
  sprintf(filename,"%s/%s",
         controlfile("source_laser/mc_filelocation","files"),
         controlfile("source_laser/mc_ideal_file","ideal.asc"));

  // data types converted automatically by CImg
  MCdata_ideal.load(filename);
  // compute fd gradients of image for optimization
  MCdata_ideal_grad = MCdata_ideal.get_gradientXY();

  // get the grid spacing data
  std::FILE *const nfile = cimg::fopen(filename,"rb");
  char line[256] = { 0 };
  std::fscanf(nfile,"%255[^\n]",line);
  unsigned int dx = 0, dy = 1, dz = 1, dv = 1;
  std::sscanf(line,"%u %u %u %u %lf %lf %lf",&dx,&dy,&dz,&dv,
                            &MC_dZ_ideal, &MC_dR_ideal, &BEAM_LOC_ideal);
  cimg::fclose(nfile);
  if(Control_Task){
     printf("loaded file %s \n",filename);
     printf("Load_MCdataIdeal:  MC_dZ_ideal    = %e\n" , MC_dZ_ideal   );
     printf("Load_MCdataIdeal:  MC_dR_ideal    = %e\n" , MC_dR_ideal   );
     printf("Load_MCdataIdeal:  BEAM_LOC_ideal = %e\n" , BEAM_LOC_ideal);
  }

  PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Load_MCdataNopt"
PetscErrorCode Load_MCdataNopt(PetscTruth Control_Task){
  char filename[MAXLEN];
  PetscFunctionBegin; // used for error handling

  GetPot controlfile(dfltINIFile);

  sprintf(filename,"%s/%s",
         controlfile("source_laser/mc_filelocation","files"),
         controlfile("source_laser/mc_file_name","Diff.asc"));

  // data types converted automatically by CImg
  MCdata.load(filename);
  // compute fd gradients of image for optimization
  MCdata_grad = MCdata.get_gradientXY();
  // get the grid spacing data
  std::FILE *const nfile = cimg::fopen(filename,"rb");
  char line[256] = { 0 };
  std::fscanf(nfile,"%255[^\n]",line);
  unsigned int dx = 0, dy = 1, dz = 1, dv = 1; 
  std::sscanf(line,"%u %u %u %u %lf %lf %lf",&dx,&dy,&dz,&dv,
                                          &MC_dZ , &MC_dR , &BEAM_LOC );
  cimg::fclose(nfile);
  if(Control_Task){
     printf("loaded file %s \n",filename);
     printf("Load_MCdataNopt:   MC_dZ    = %e\n"   ,  MC_dZ   );
     printf("Load_MCdataNopt:   MC_dR    = %e\n"   ,  MC_dR   );
     printf("Load_MCdataNopt:   BEAM_LOC = %e\n"   ,  BEAM_LOC);
  }

  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- 
   Register - Register the point with the  natural coordinate system
              of the monte carlo grid

*/
void Register(PetscScalar *point, PetscScalar *origin,
                                  PetscInt *rcoord, PetscInt *zcoord){
  vector<PetscScalar>  point_shift(3,0.0),new_point(3,0.0);
  int  iii,jjj;
  //SHIFT COORDINATES
  for(iii =0 ; iii < 3 ; iii++) point_shift.at(iii) = point[iii] - origin[iii];
  //ROTATE REFERENCE FRAME
  for(iii =0 ; iii < 3 ; iii++) 
      for(jjj =0 ; jjj < 3 ; jjj++) 
          new_point.at(iii) = new_point.at(iii) 
                              + ROTMAT[iii][jjj]*point_shift.at(jjj); 
  //PUT INTO cylindrical
  PetscScalar r = sqrt(pow(new_point.at(0),2)+ pow(new_point.at(1),2));
  PetscScalar z = new_point.at(2) + BEAM_LOC ; 
  *rcoord = int( r / MC_dR);
  *zcoord = int( z / MC_dZ);
  return;
}
/* -------------------------------------------------------------------- 
The functions mcheat_nopt and mcheat_ideal interpolates the Monte Carlo 
source term at the quadrature point (x,y,z). Since the coordinate systems
of the Monte Carlo system and that of the mesh are different, the routine
first shifts the coordinates and then rotates according to the pre-computed
rotation matrix. 

 double *point   -  the quadrature point
 double *origin  -  the location of the laser
 
 the routines are called from:
                   nonlinpennesmonte
                   pennesidealmonte
                   dqmontedpow
*/
PetscScalar FORTRAN_NAME(mcheat_nopt)(PetscScalar *point, PetscScalar *origin){
  PetscScalar heating;
  PetscInt rcoord=0,zcoord=0;

  
  Register(point,origin,&rcoord,&zcoord);

  heating = 0.0;
  // axial dimension varies by row , radial dimension varies by column
  if( rcoord < MCdata.dimy() && zcoord >= 0  && zcoord < MCdata.dimx()) 
                                    heating=MCdata(zcoord,rcoord);

  return heating;
}

PetscScalar FORTRAN_NAME(mcheat_ideal)(PetscScalar *point, PetscScalar *origin){
  PetscScalar heating;
  PetscInt rcoord=0,zcoord=0;

  Register(point,origin,&rcoord,&zcoord);

  heating = 0.0;
  // axial dimension varies by row , radial dimension varies by column
  if( rcoord<MCdata_ideal.dimy() && zcoord>=0  && zcoord<MCdata_ideal.dimx()) 
                            heating=MCdata_ideal(zcoord,rcoord);

  return heating;
}
/* -------------------------------------------------------------------- 
The functions dmcheat_noptdvector computes the directional 
derivative in the direction derivdir of the the Monte Carlo grid
source term at the quadrature point (x,y,z). Since the coordinate systems
of the Monte Carlo system and that of the mesh are different, the routine
first shifts the coordinates and then rotates according to the pre-computed
rotation matrix. 

 double *point   -  the quadrature point
 double *origin  -  the location of the laser
 
 the routines are called from:
                   dqmontedx
                   dqmontedy
                   dqmontedz


the derivatives were done in matlab


 P(x,y,z,o1,o2,o3) = MC(r,zhat)
                   = MC(sqrt(  [a11*(x-o1)+a12*(y-o2)+a13*(z-o3)]^2 +
                               [a21*(x-o1)+a22*(y-o2)+a23*(z-o3)]^2   ) ,
                                a31*(x-o1)+a32*(y-o2)+a33*(z-o3) + BEAM_LOC  )


 dPdo1 = dMCdr drdo1 + dMCdzhat dzhatdo1
 dPdo2 = dMCdr drdo2 + dMCdzhat dzhatdo2
 dPdo3 = dMCdr drdo3 + dMCdzhat dzhatdo3


dPdx = 
  dMCdr (-(a11 (x - o1) + a12 (y - o2) + a13 (z - o3)) a11
         - (a21 (x - o1) + a22 (y - o2) + a23 (z - o3)) a21)/r - dMCdzhat a31

dPdy = 
  dMCdr (-(a11 (x - o1) + a12 (y - o2) + a13 (z - o3)) a12
         - (a21 (x - o1) + a22 (y - o2) + a23 (z - o3)) a22)/r - dMCdzhat a32

dPdz = 
  dMCdr (-(a11 (x - o1) + a12 (y - o2) + a13 (z - o3)) a13
         - (a21 (x - o1) + a22 (y - o2) + a23 (z - o3)) a23)/r - dMCdzhat a33

*/

PetscScalar FORTRAN_NAME(dmcheat_noptdvector)(PetscScalar *point, 
                                              PetscScalar *origin,
                                              PetscInt *icomp){
  PetscScalar DirectionalDeriv;
  PetscInt iii,jjj,rcoord=0,zcoord=0;
  vector<PetscScalar>  point_shift(3,0.0),new_point(3,0.0);

  //SHIFT COORDINATES
  for(iii =0 ; iii < 3 ; iii++) point_shift.at(iii) = point[iii] - origin[iii];
  //ROTATE REFERENCE FRAME
  for(iii =0 ; iii < 3 ; iii++) 
      for(jjj =0 ; jjj < 3 ; jjj++) 
          new_point.at(iii) = new_point.at(iii) 
                              + ROTMAT[iii][jjj]*point_shift.at(jjj); 
  //PUT INTO cylindrical
  PetscScalar r = sqrt(pow(new_point.at(0),2)+ pow(new_point.at(1),2));
  PetscScalar z = new_point.at(2) + BEAM_LOC ; 
  rcoord = int( r / MC_dR);
  zcoord = int( z / MC_dZ);


#if defined(PETSC_USE_DEBUG)
    if(*icomp > 2 || *icomp < 0 ){
        printf("dmcheat_noptdvector: bounds error\n"); abort();
    }
#endif
  DirectionalDeriv = 0.0;
  // axial dimension varies by row , radial dimension varies by column
  if( rcoord < MCdata.dimy() && zcoord >= 0  && zcoord < MCdata.dimx()) 
     DirectionalDeriv =-MCdata_grad[1](zcoord,rcoord)/MC_dR  *
                        (new_point.at(0)*ROTMAT[0][*icomp] + 
                         new_point.at(1)*ROTMAT[1][*icomp] )/r
                       -MCdata_grad[0](zcoord,rcoord)/MC_dZ *ROTMAT[2][*icomp];

  return DirectionalDeriv ;
}


