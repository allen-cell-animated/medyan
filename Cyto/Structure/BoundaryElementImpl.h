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
#include <cmath>
#include "BoundaryElement.h"

///PlaneBoundaryElement is a plane implementation of a boundary element
class PlaneBoundaryElement : public BoundaryElement {
    
private:
    ///Parameters of equation (ax + by + cz + d = 0)
    double _a, _b, _c, _d;

public:
    ///Constructor, sets parameters of equation
    PlaneBoundaryElement(std::vector<double> coords, std::vector<double> normal);
    
    ///distance to this plane
    virtual double distance(const std::vector<double>& point);
};

class TriangleBoundaryElement : public BoundaryElement {
    
    ///not yet implemented

};


#endif /* defined(__Cyto__BoundaryElementImpl__) */
