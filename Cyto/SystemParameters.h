//
//  SystemParameters.h
//  Cyto
//
//  Created by James Komianos on 9/9/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__SystemParameters__
#define __Cyto__SystemParameters__

#include <iostream>
#include "common.h"

///Struct to hold the read mechanics parameters
struct MechanicsParameters {
    
    ///Filament parameters
    double FStretchingK = 0;
    double FBendingK = 0;
    double FBendingTheta = 0;
    double FTwistingK = 0;
    double FTwistingPhi = 0;
    
    ///Linker parameters
    double LStretchingK = 0;
    double LStretchingL = 0;
    double LBendingK = 0;
    double LBendingTheta = 0;
    double LTwistingK = 0;
    double LTwistingPhi = 0;
    
    ///Motor parameters
    double MStretchingK = 0;
    double MStretchingL = 0;
    double MBendingK = 0;
    double MBendingTheta = 0;
    double MTwistingK = 0;
    double MTwistingPhi = 0;
    
    ///Volume parameters
    double VolumeK = 0;
};

///Struct to hold the read geometry parameters
struct GeometryParameters {
    short nDim = 0;
    
    int NX = 0;
    int NY = 0;
    int NZ = 0;
    
    double compartmentSizeX = 0;
    double compartmentSizeY = 0;
    double compartmentSizeZ = 0;
    
    double monomerSize = 0;
    double cylinderSize = 0;
};

///Struct to hold the read boundary parameters
struct BoundaryParameters {
    
    double boundaryCutoff = 0;
    double boundaryK = 0;
};


///This class holds all system-wide parameters, initialized by the Parser
class SystemParameters {
friend class SystemParser;
    
#ifdef TESTING ///Public access if testing only
public:
#endif
    static MechanicsParameters MParams; ///< mechanical parameters
    static GeometryParameters GParams;  ///< geometry parameters
    static BoundaryParameters BParams; ///< boundary parameters
    
public:
    ///Const getters for all parameters
    static const MechanicsParameters& Mechanics() {return MParams;}
    static const GeometryParameters& Geometry() {return GParams;}
    static const BoundaryParameters& Boundaries() {return BParams;}
};






#endif /* defined(__Cyto__SystemParameters__) */
