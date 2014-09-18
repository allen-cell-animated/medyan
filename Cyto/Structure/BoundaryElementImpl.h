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

#include "common.h"
#include "BoundaryElement.h"

///PlaneBoundaryElement is a plane implementation of a boundary element
class PlaneBoundaryElement : public BoundaryElement {
    
private:
    ///Parameters of equation (ax + by + cz + d = 0)
    double _a, _b, _c, _d;
    double _k_rep;

public:
    ///Constructor, sets parameters of equation
    PlaneBoundaryElement(std::vector<double> coords, std::vector<double> normal, double repulsConst);
    
    ///Distance to this plane
    virtual double distance(const std::vector<double>& point);
    
    ///Stretched distance to plane
    virtual double stretchedDistance(const std::vector<double>& point, const std::vector<double>& force, double d);
    
    ///get the repulsion constant for this plane
    virtual double getRepulsionConst();
};

class TriangleBoundaryElement : public BoundaryElement {
    
    ///not yet implemented

};


#endif /* defined(__Cyto__BoundaryElementImpl__) */
