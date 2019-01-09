#ifndef MEDYAN_AdaptiveMeshAttribute_hpp
#define MEDYAN_AdaptiveMeshAttribute_hpp

#include "MathFunctions.h"

// Additional attributes needed for meshwork
struct AdaptiveMeshAttribute {
    struct VertexAttribute {
        double density0;
        double densityAvg;
        double l0Aux;
        mathfunc::Vec3 unitNormal;
    };
    struct HalfEdgeAttribute {
    };
    struct EdgeAttribute {
        double l0;
    };
    struct TriangleAttribute {
        double area0;
    };
};

#endif
