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
    float _size; ///< size of square (from center to side at a right angle)
    
public:
    
    SquareBoundaryElement(std::vector<double>& coords, float size) :
        _size(size), BoundaryElement(coords) {}
    
};

class TriangleBoundaryElement : public BoundaryElement {
    
    ///not yet implemented

};


#endif /* defined(__Cyto__BoundaryElementImpl__) */
