#ifndef MEDYAN_Structure_SurfaceMesh_MVertex_hpp
#define MEDYAN_Structure_SurfaceMesh_MVertex_hpp

struct MVoronoiCell {

    MVoronoiCell(short membraneType);

    double kBending; // Local bending modulus
    double eqCurv; // Local spontaneous curvature
};


#endif
