
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_SystemParameters_h
#define M3SYM_SystemParameters_h

#include <vector>

#include "common.h"

/// MechanicsParameters is a struct to hold mechanical parameters for the system
struct MechanicsParameters {
    
    //@{
    /// Filament parameter
    double FStretchingK = 0;
    double FBendingK = 0;
    double FBendingTheta = 0;
    double FTwistingK = 0;
    double FTwistingPhi = 0;
    //@}
    
    //@{
    /// Linker parameter
    vector<double> LStretchingK = {};
    vector<double> LBendingK = {};
    vector<double> LBendingTheta = {};
    vector<double> LTwistingK = {};
    vector<double> LTwistingPhi = {};
    //@}
    
    //@{
    /// Motor parameter
    vector<double> MStretchingK = {};
    vector<double> MBendingK = {};
    vector<double> MBendingTheta = {};
    vector<double> MTwistingK = {};
    vector<double> MTwistingPhi = {};
    //@}
    
    //@{
    /// Volume parameter
    double VolumeK = 0;
    double VolumeCutoff = 0;
    //@}
};

/// ChemistryParameters is a struct to hold chemistry parameters for the system
struct ChemistryParameters {
    
    //@{
    /// Number of general species
    short numBulkSpecies = 0;
    short numDiffusingSpecies = 0;
    //@}
    
    //@{
    /// Number of filament related species
    short numFilamentSpecies = 0;
    short numPlusEndSpecies = 0;
    short numMinusEndSpecies = 0;
    short numBoundSpecies = 0;
    short numLinkerSpecies = 0;
    short numMotorSpecies = 0;
    //@}
};

/// GeometryParameters is a struct to hold geometry parameters for the system
struct GeometryParameters {
    
    //@{
    /// Geometry parameter
    short nDim = 0;
    
    int NX = 0;
    int NY = 0;
    int NZ = 0;
    
    double compartmentSizeX = 0;
    double compartmentSizeY = 0;
    double compartmentSizeZ = 0;
    
    double largestCompartmentSide = 0;
    
    double monomerSize = 0;
    double cylinderSize = 0;
    int cylinderIntSize = 0;
    //@}
    
};

/// BoundaryParameters is a struct to hold boundary parameters for the system
struct BoundaryParameters {
    
    double boundaryCutoff = 0;
    double boundaryK = 0;
    double screenLength = 0;
    
    double diameter = 0;
};


///SystemParameters is a static class that holds all simulation parameters, initialized by the [SystemParser] (@ref SystemParser)
class SystemParameters {
friend class SystemParser;
    
#ifdef TESTING ///Public access if testing only
public:
#endif
    static MechanicsParameters MParams; ///< The mechanical parameters
    static ChemistryParameters CParams; ///< The chemistry parameters
    static GeometryParameters GParams;  ///< The geometry parameters
    static BoundaryParameters BParams; ///< The boundary parameters
    
public:
    //@{
    ///Const getter
    static const MechanicsParameters& Mechanics() {return MParams;}
    static const ChemistryParameters& Chemistry() {return CParams;}
    static const GeometryParameters& Geometry() {return GParams;}
    static const BoundaryParameters& Boundaries() {return BParams;}
    //@}
};






#endif /* defined(__Cyto__SystemParameters__) */