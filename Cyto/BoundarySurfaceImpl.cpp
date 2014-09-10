//
//  BoundarySurfaceImpl.cpp
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundarySurfaceImpl.h"
#include "MathFunctions.h"

using namespace mathfunc;

BasicPlane::BasicPlane(std::vector<std::vector<double>> points, std::vector<int> numDivisions, short orientation) :
_points(points), _orientation(orientation), BoundarySurface(3) {
    
    ///check if valid orientation
    assert(_orientation >= 0 && _orientation <= 2);
    
    ///Calculate side lengths and boundary element sizes
    double L1 = TwoPointDistance(points[0], points[1]);
    double L2 = TwoPointDistance(points[1], points[2]);
    
    float s1 = L1 / numDivisions[0];
    float s2 = L2 / numDivisions[1];
    
    ///Fill this plane with square boundary elements
    
    ///CASE 1: in XY plane
    if(_orientation == 2) {
        
        ///coordinates of boundary element to add
        std::vector<double> BECoordinate;
        
        ///loop over column
        for(int i = 0; i < numDivisions[1]; i++) {
            
            ///reset to next row
            BECoordinate = std::vector<double>{points[0][0] + (s1 / 2), points[0][1] + (i * s2) + (s2 / 2), points[0][2]};
            
            ///loop over row
            for(int j = 0; j < numDivisions[0]; j++) {
                
                _boundaryElements.emplace_back(BoundaryElementDB::Instance(BEDBKey())->CreateSquareBoundaryElement(BECoordinate, {s1, s2}, _orientation));
                BECoordinate[0] += s1;
            }
        }
    }
    
    ///CASE 2: in XZ plane
    if(_orientation == 1) {
        
        ///coordinates of boundary element to add
        std::vector<double> BECoordinate;
        
        ///loop over column
        for(int i = 0; i < numDivisions[1]; i++) {
            
            ///reset to next row
            BECoordinate = std::vector<double>{points[0][0] + (s1 / 2), points[0][1], points[0][2] + (i * s2) + (s2 / 2)};
            
            ///loop over row
            for(int j = 0; j < numDivisions[0]; j++) {
                
                _boundaryElements.emplace_back(BoundaryElementDB::Instance(BEDBKey())->CreateSquareBoundaryElement(BECoordinate, {s1, s2}, _orientation));
                BECoordinate[0] += s1;
            }
        }
    }
    
    ///CASE 3: in YZ plane
    if(_orientation == 0) {
        
        ///coordinates of boundary element to add
        std::vector<double> BECoordinate;
        
        ///loop over column
        for(int i = 0; i < numDivisions[1]; i++) {
            
            ///reset to next row
            BECoordinate = std::vector<double>{points[0][0], points[0][1] + (s1 / 2), points[0][2] + (i * s2) + (s2 / 2)};
            
            ///loop over row
            for(int j = 0; j < numDivisions[0]; j++) {
                
                _boundaryElements.emplace_back(BoundaryElementDB::Instance(BEDBKey())->CreateSquareBoundaryElement(BECoordinate, {s1, s2}, _orientation));
                BECoordinate[1] += s1;
            }
        }
    }
    

}