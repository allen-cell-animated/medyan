#ifndef MEDYAN_Structure_SurfaceMesh_MTriangle_Hpp
#define MEDYAN_Structure_SurfaceMesh_MTriangle_Hpp

/******************************************************************************

Storing some mechanical properties of the triangle patches

******************************************************************************/

struct MTriangle {

    // Local area elasticity, applicable only in material coordinates
    double kArea = 0.0;
    double eqArea = 0.0;
};


#endif
