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
#include <vector>
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
    std::vector<double> LStretchingK = {};
    std::vector<double> LBendingK = {};
    std::vector<double> LBendingTheta = {};
    std::vector<double> LTwistingK = {};
    std::vector<double> LTwistingPhi = {};
    
    ///Motor parameters
    std::vector<double> MStretchingK = {};
    std::vector<double> MBendingK = {};
    std::vector<double> MBendingTheta = {};
    std::vector<double> MTwistingK = {};
    std::vector<double> MTwistingPhi = {};
    
    ///Volume parameters
    double VolumeK = 0;
};

///Struct to hold the read chemistry parameters
struct ChemistryParameters {
    
    ///number of general species
    short numBulkSpecies = 0;
    short numDiffusingSpecies = 0;
    
    ///number of filament related species
    short numFilamentSpecies = 0;
    short numPlusEndSpecies = 0;
    short numMinusEndSpecies = 0;
    short numBoundSpecies = 0;
    short numLinkerSpecies = 0;
    short numMotorSpecies = 0;
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
    double screenLength = 0;
};


///This class holds all system-wide parameters, initialized by the Parser
class SystemParameters {
friend class SystemParser;
    
#ifdef TESTING ///Public access if testing only
public:
#endif
    static MechanicsParameters MParams; ///< mechanical parameters
    static ChemistryParameters CParams; ///< chemistry parameters
    static GeometryParameters GParams;  ///< geometry parameters
    static BoundaryParameters BParams; ///< boundary parameters
    
public:
    ///Const getters for all parameters
    static const MechanicsParameters& Mechanics() {return MParams;}
    static const ChemistryParameters& Chemistry() {return CParams;}
    static const GeometryParameters& Geometry() {return GParams;}
    static const BoundaryParameters& Boundaries() {return BParams;}
};






#endif /* defined(__Cyto__SystemParameters__) */
