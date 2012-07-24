#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/iterator/zip_iterator.h>
#include <iostream>

// This example shows how thrust::zip_iterator can be used to create a 
// 'virtual' array of structures.  In this case the structure is a 3d 
// vector type (Float3) whose (x,y,z) components will be stored in 
// three separate float arrays.  The zip_iterator "zips" these arrays
// into a single virtual Float3 array.


// We'll use a 3-tuple to store our 3d vector type
typedef thrust::tuple<PetscScalar,PetscScalar,PetscScalar> PetscScalar3;
// We'll use an 8-tuple to store our element neighbors
typedef thrust::tuple<PetscInt,PetscInt,PetscInt,PetscInt,
                      PetscInt,PetscInt,PetscInt,PetscInt> Neighbor8;


// This functor implements the dot product between 3d vectors
struct DotProduct : public thrust::binary_function<Float3,Float3,float>
{
    __host__ __device__
        float operator()(const Float3& a, const Float3& b) const
        {
            return thrust::get<0>(a) * thrust::get<0>(b) +    // x components
                   thrust::get<1>(a) * thrust::get<1>(b) +    // y components
                   thrust::get<2>(a) * thrust::get<2>(b);     // z components
        }
};

struct tuple_of_tuples_to_ad_hoc_tuple
  : thrust::unary_function<tuple_of_tuples, PetscInt>
{
  __host__ __device__
  PetscInt operator()(const tuple_of_tuples& t)
  {
    return ad_hoc_tuple(thrust::get<0>(thrust::get<0>(t)), thrust::get<0>(thrust::get<1>(t)));
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
        if( neighborElem ) 
          { // check not a boundary neighbor
           Point tmpDistance = neighborElem->centroid() - currentElem->centroid();
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
     }
    // error check
    if (!nextElem )
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
    template_parameter_system.get_dof_map().dof_indices (nextElem, param_dof_indices);

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
  }
};

thrust::host_vector<float> random_vector(const size_t N)
{
    thrust::default_random_engine rng;
    thrust::uniform_real_distribution<float> u01(0.0f, 1.0f);
    thrust::host_vector<float> temp(N);
    for(size_t i = 0; i < N; i++) {
        temp[i] = u01(rng);
    }
    return temp;
}


PetscErrorCode PrimaryFluence(
         thrust::host_vector<PetscScalar> hostCentroidX, //INPUT GLOBAL list
         thrust::host_vector<PetscScalar> hostCentroidY, //INPUT GLOBAL list
         thrust::host_vector<PetscScalar> hostCentroidZ, //INPUT GLOBAL list
         thrust::host_vector<PetscScalar> hostScattering,//INPUT GLOBAL list
         thrust::host_vector<PetscScalar> hostAbsorption,//INPUT GLOBAL list
         thrust::host_vector<PetscInt   > hostNeighbor1, //INPUT GLOBAL list
         thrust::host_vector<PetscInt   > hostNeighbor2, //INPUT GLOBAL list
         thrust::host_vector<PetscInt   > hostNeighbor3, //INPUT GLOBAL list
         thrust::host_vector<PetscInt   > hostNeighbor4, //INPUT GLOBAL list
         thrust::host_vector<PetscInt   > hostNeighbor5, //INPUT GLOBAL list
         thrust::host_vector<PetscInt   > hostNeighbor6, //INPUT GLOBAL list
         thrust::host_vector<PetscInt   > laserList,   //INPUT list of laser probe Id
         thrust::host_vector<PetscScalar> hostNode_x,  //INPUT LOCAL Position
         thrust::host_vector<PetscScalar> hostNode_y,  //INPUT LOCAL Position
         thrust::host_vector<PetscScalar> hostNode_z,  //INPUT LOCAL Position
         thrust::host_vector<PetscScalar> hostFluence, //OUTPUT LOCAL component
         thrust::host_vector<PetscScalar> hostFlux_x,  //OUTPUT LOCAL component
         thrust::host_vector<PetscScalar> hostFlux_y,  //OUTPUT LOCAL component
         thrust::host_vector<PetscScalar> hostFlux_z,  //OUTPUT LOCAL component
                             )
{
    PetscFunctionBegin;
    // number of vectors
    const size_t N = 1000;

    // We'll store the components of the 3d vectors in separate arrays. One set of
    // arrays will store the 'A' vectors and another set will store the 'B' vectors.

    // This 'structure of arrays' (SoA) approach is usually more efficient than the 
    // 'array of structures' (AoS) approach.  The primary reason is that structures,
    // like Float3, don't always obey the memory coalescing rules, so they are not
    // efficiently transferred to and from memory.  Another reason to prefer SoA to
    // AoS is that we don't aways want to process all members of the structure.  For
    // example, if we only need to look at first element of the structure then it 
    // is wasteful to load the entire structure from memory.  With the SoA approach,
    // we can chose which elements of the structure we wish to read.

    thrust::device_vector<PetscScalar> deviceCentroidX=hostCentroidX; // x centroid
    thrust::device_vector<PetscScalar> deviceCentroidY=hostCentroidY; // y centroid
    thrust::device_vector<PetscScalar> deviceCentroidZ=hostCentroidZ; // z centroid

    // element neighbors
    thrust::device_vector<PetscInt> deviceNeighbor1=hostNeighbor1;  
    thrust::device_vector<PetscInt> deviceNeighbor2=hostNeighbor2; 
    thrust::device_vector<PetscInt> deviceNeighbor3=hostNeighbor3;
    thrust::device_vector<PetscInt> deviceNeighbor4=hostNeighbor4;
    thrust::device_vector<PetscInt> deviceNeighbor5=hostNeighbor5;
    thrust::device_vector<PetscInt> deviceNeighbor6=hostNeighbor6;

    // Storage for result of each dot product
    thrust::device_vector<float> result(N);


    // We'll now illustrate two ways to use zip_iterator to compute the dot
    // products.  The first method is verbose but shows how the parts fit together.
    // The second method hides these details and is more concise.
   

    // METHOD #1
    // Defining a zip_iterator type can be a little cumbersome ...
    typedef thrust::device_vector<float>::iterator                     FloatIterator;
    typedef thrust::tuple<FloatIterator, FloatIterator, FloatIterator> FloatIteratorTuple;
    typedef thrust::zip_iterator<FloatIteratorTuple>                   Float3Iterator;

    // Now we'll create some zip_iterators for A and B
    Float3Iterator A_first = thrust::make_zip_iterator(make_tuple(A0.begin(), A1.begin(), A2.begin()));
    Float3Iterator A_last  = thrust::make_zip_iterator(make_tuple(A0.end(),   A1.end(),   A2.end()));
    Float3Iterator B_first = thrust::make_zip_iterator(make_tuple(B0.begin(), B1.begin(), B2.begin()));
                            
    // Finally, we pass the zip_iterators into transform() as if they
    // were 'normal' iterators for a device_vector<Float3>.
    thrust::transform(A_first, A_last, B_first, result.begin(), DotProduct());


    // METHOD #2
    // Alternatively, we can avoid creating variables for X_first, X_last, 
    // and Y_first and invoke transform() directly.
    thrust::transform( thrust::make_zip_iterator(make_tuple(A0.begin(), A1.begin(), A2.begin())),
                       thrust::make_zip_iterator(make_tuple(A0.end(),   A1.end(),   A2.end())),
                       thrust::make_zip_iterator(make_tuple(B0.begin(), B1.begin(), B2.begin())),
                       result.begin(),
                       DotProduct() );
    


    // Finally, we'll print a few results

    // Example output
    // (0.840188,0.45724,0.0860517) * (0.0587587,0.456151,0.322409) = 0.285683
    // (0.394383,0.640368,0.180886) * (0.0138811,0.24875,0.0221609) = 0.168775
    // (0.783099,0.717092,0.426423) * (0.622212,0.0699601,0.234811) = 0.63755
    // (0.79844,0.460067,0.0470658) * (0.0391351,0.742097,0.354747) = 0.389358
    std::cout << std::fixed;
    for(size_t i = 0; i < 4; i++)
    {
        Float3 a = A_first[i];
        Float3 b = B_first[i];
        float dot = result[i];

        std::cout << "(" << thrust::get<0>(a) << "," << thrust::get<1>(a) << "," << thrust::get<2>(a) << ")";
        std::cout << " * ";
        std::cout << "(" << thrust::get<0>(b) << "," << thrust::get<1>(b) << "," << thrust::get<2>(b) << ")";
        std::cout << " = ";
        std::cout << dot << std::endl;
    }   

    PetscFunctionReturn(0);
}
 
// an ad hoc tuple structure
// can have as many elements as necessary
struct ad_hoc_tuple
{
  __host__ __device__
  inline ad_hoc_tuple(float &xx, float &yy)
    : x(xx), y(yy)
  {}

  // allow assignment from float, for expository purposes in the example below
  ad_hoc_tuple &operator=(float f)
  {
    x = f;
    y = f;
    return *this;
  }

  // because these are bare references, you won't be able to
  // dereference the iterator on the host
  float &x;
  float &y;
};

// note the references in the tuples below
typedef thrust::tuple<
  thrust::tuple<float &>,
  thrust::tuple<float &>
> tuple_of_tuples;

struct tuple_of_tuples_to_ad_hoc_tuple
  : thrust::unary_function<tuple_of_tuples, ad_hoc_tuple>
{
  __host__ __device__
  ad_hoc_tuple operator()(tuple_of_tuples t)
  {
    return ad_hoc_tuple(thrust::get<0>(thrust::get<0>(t)), thrust::get<0>(thrust::get<1>(t)));
  }
};

int main(void)
{
  size_t n = 5;
  thrust::device_vector<float> x_vec(n), y_vec(n);

  typedef thrust::device_vector<float>::iterator FloatIterator;

  // first level of zip
  typedef thrust::tuple<FloatIterator> FloatIteratorTuple;
  typedef thrust::zip_iterator<FloatIteratorTuple> ZipIterator1;
  ZipIterator1 x_begin = thrust::make_zip_iterator(thrust::make_tuple(x_vec.begin()));
  ZipIterator1 y_begin = thrust::make_zip_iterator(thrust::make_tuple(y_vec.begin()));

  // second level of zip
  typedef thrust::tuple<ZipIterator1, ZipIterator1> ZipIteratorTuple;
  typedef thrust::zip_iterator<ZipIteratorTuple> ZipIterator2;
  ZipIterator2 zip_begin = thrust::make_zip_iterator(thrust::make_tuple(x_begin, y_begin));

  // transform ZipIterator2 into an iterator returning an ad_hoc_tuple
  typedef thrust::transform_iterator<
    tuple_of_tuples_to_ad_hoc_tuple,
    ZipIterator2,
    ad_hoc_tuple
  > XfrmIterator;

  XfrmIterator iter(zip_begin, tuple_of_tuples_to_ad_hoc_tuple());

  // fill x_vec & y_vec with 13 through the iterator
  thrust::fill(iter, iter + n, 13);

  // print the results on the host
  for(int i = 0; i < n; ++i)
  {
    std::cout << "x[" << i << "]: " << x_vec[i] << std::endl;
    std::cout << "y[" << i << "]: " << y_vec[i] << std::endl;
  }

  return 0;
}

//Interesting, I had not considered nested tuples, I'll look into this. 
//Although, I have got something working but I'm not sure if it is efficient: 
//nested zip_iterators. 
//E.g, 
//  make_zip_iterator(make_tuple(make_zip_iterator(make_tuple(px.begin(),py.beg
//in(),pz.begin())), 
//keys.begin(), indices.begin())); 
//It just seems to work?! 
//For instance, 
struct Mesh 
{ 
   // ... members as first post 
  typedef thrust::tuple< 
     float,float,float, 
     float,float,float, 
     float,float,float 
  > triangle_tuple; 
  typedef typename thrust::device_vector<float>::iterator iter; 
  typedef thrust::tuple< 
     iter,iter,iter, 
     iter,iter,iter, 
     iter,iter,iter 
  > iter_tuple; 
  typedef thrust::zip_iterator<iter_tuple> triangle_iterator; 
  triangle_iterator begin() 
  { 
    return thrust::make_zip_iterator(thrust::make_tuple( 
              px.begin(),py.begin().pz.begin(), 
              qx.begin(),qy.begin(),qz.begin(), 
              rx.begin(),ry.begin(),rz.begin())); 
  } 
  // same as begin() above 
  triangle_iterator end() 
  { 
  } 
  // helper to get points p,q,r from triangle_tuple 
  void get_pqr(triangle_tuple const & t, float3 & p, float3 & q, float3 & r) 
  { 
         p = thrust::make_float3(get<0>(t),get<1>(t),get<2>(t)); 
         q = thrust::make_float3(get<3>(t),get<4>(t),get<5>(t)); 
         r = thrust::make_float3(get<6>(t),get<7>(t),get<8>(t)); 
  } 
}; 

//And usage is as follows: 
thrust::for_each( 
make_zip_iterator(make_tuple(mesh.begin(),indices.begin(),keys.begin()), 
  // nested zip_iterator here 
       make_zip_iterator(make_tuple(mesh.end(),indices.end(),keys.end()), 
       point_in_bbox(nodes,bbox) 
))); 
//And retrieving the triangle_tuple can be down as: 
// point_in_bbox from other post 
 void operator()(Tuple tuple) 
 { 
            Mesh::triangle_tuple tri_tuple = get<0>(t); 
            // Or with the helper: 
            float3 p,q,r; 
            Mesh::get_pqr(get<0>(t),p,q,r); 
            int node_index = get<1>(t); 
            int& key = get<2>(t); 
           // ... do stuff with paramaters ... 
 } 
