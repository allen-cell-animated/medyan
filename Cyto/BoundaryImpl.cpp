//
//  BoundaryImpl.cpp
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryImpl.h"

BoundaryCubic::BoundaryCubic(int nDim, std::vector<std::vector<double>> points, std::vector<int> numDivisions) :
Boundary(nDim, BoundaryType::Cube){
    
    ///Create boundary surfaces, add to vector
    ///Y normal planes
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({points[0], points[1], points[5], points[4]}, numDivisions, 1)));
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({points[3], points[2], points[6], points[7]}, numDivisions, 1)));
    
    ///X normal planes
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({points[0], points[3], points[7], points[4]}, numDivisions, 0)));
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({points[1], points[2], points[6], points[5]}, numDivisions, 0)));
    
    ///Z normal planes
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({points[0], points[1], points[2], points[3]}, numDivisions, 2)));
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({points[4], points[5], points[6], points[7]}, numDivisions, 2)));
    
}