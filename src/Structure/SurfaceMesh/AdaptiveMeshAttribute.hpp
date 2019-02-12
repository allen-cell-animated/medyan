#ifndef MEDYAN_AdaptiveMeshAttribute_hpp
#define MEDYAN_AdaptiveMeshAttribute_hpp

#include "MathFunctions.h"

// Additional attributes needed for meshwork
struct AdaptiveMeshAttribute {
    struct VertexAttribute {
        double size;
        double maxSize;
        double sizeAux; // Used in diffusing
        mathfunc::Vec3 unitNormal;
    };
    struct HalfEdgeAttribute {
    };
    struct EdgeAttribute {
        double eqLength;
    };
    struct TriangleAttribute {
    };
};

#endif
