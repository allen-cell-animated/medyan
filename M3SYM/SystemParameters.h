
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

/// Struct to hold mechanical parameters for the system
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
    /// MotorGhost parameter
    vector<double> MStretchingK = {};
    vector<double> MBendingK = {};
    vector<double> MBendingTheta = {};
    vector<double> MTwistingK = {};
    vector<double> MTwistingPhi = {};
    //@}
    
    //@{
    /// BranchingPoint parameter
    vector<double> BrStretchingK = {};
    vector<double> BrStretchingL = {};
    vector<double> BrBendingK = {};
    vector<double> BrBendingTheta = {};
    vector<double> BrDihedralK = {};
    vector<double> BrPositionK = {};
    //@}
    
    //@{
    /// Volume parameter
    double VolumeK = 0;
    double VolumeCutoff = 0;
    //@}
};

/// Struct to hold chemistry parameters for the system
struct ChemistryParameters {
    
    //@{
    /// Number of general species
    short numBulkSpecies = 0;
    short numDiffusingSpecies = 0;
    //@}
    
    //@{
    /// Number of filament related species
    short numFilamentSpecies = 0;
    short numPlusEndSpecies  = 0;
    short numMinusEndSpecies = 0;
    short numBoundSpecies    = 0;
    short numLinkerSpecies   = 0;
    short numMotorSpecies    = 0;
    short numBrancherSpecies = 0;
    
    short numBindingSites = 0;
    //@}
};

/// Struct to hold geometry parameters for the system
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
    
    /// Minimum monomer length of a cylinder
    int minCylinderIntSize = 3;
    double minCylinderSize = 0;
    //@}
    
};

/// Struct to hold Boundary parameters for the system
struct BoundaryParameters {
    
    double boundaryCutoff = 0;
    double boundaryK = 0;
    double screenLength = 0;
    
    double diameter = 0;
};

/// Struct to hold dynamic rate changing parameters
struct DynamicRateParameters {
    
    /// Option for dynamic polymerization rate of filaments
    double FDPLength = 0.0;
    /// Option for dynamic unbinding rate of linkers
    vector<double> LDULength = {};
    /// Option for dynamic unbinding rate of motors
    vector<double> MDULength = {};
    /// Option for dynamic walking rate of motors
    vector<double> MDWLength = {};
    
};

/// Static class that holds all simulation parameters,
/// initialized by the SystemParser
class SystemParameters {
friend class SystemParser;
    
#ifdef TESTING ///Public access if testing only
public:
#endif
    static MechanicsParameters MParams;    ///< The mechanical parameters
    static ChemistryParameters CParams;    ///< The chemistry parameters
    static GeometryParameters GParams;     ///< The geometry parameters
    static BoundaryParameters BParams;     ///< The boundary parameters
    static DynamicRateParameters DRParams; ///< The dynamic rate parameters
    
public:
    //@{
    ///Const getter
    static const MechanicsParameters& Mechanics() {return MParams;}
    static const ChemistryParameters& Chemistry() {return CParams;}
    static const GeometryParameters& Geometry() {return GParams;}
    static const BoundaryParameters& Boundaries() {return BParams;}
    static const DynamicRateParameters& DynamicRates() {return DRParams;}
    //@}
};

#endif
