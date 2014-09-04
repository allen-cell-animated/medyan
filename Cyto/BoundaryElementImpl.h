//
//  BoundaryElementImpl.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryElementImpl__
#define __Cyto__BoundaryElementImpl__

#include <iostream>
#include "BoundaryElement.h"

///SquareBoundaryElement is a square implementation of a boundary element
class SquareBoundaryElement : public BoundaryElement {
    
private:
    std::vector<float> _sides; ///< side lengths of square
    short _orientation; ///< normal vector to the square (can only be unit vectors x(0),y(1),z(2))
    
public:
    SquareBoundaryElement(std::vector<double> coords, std::vector<float> sides, int orientation) :
        _sides(sides), _orientation(orientation), BoundaryElement(coords) {}
    
};

class TriangleBoundaryElement : public BoundaryElement {
    
    ///not yet implemented

};


#endif /* defined(__Cyto__BoundaryElementImpl__) */
