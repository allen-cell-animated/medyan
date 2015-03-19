
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_SysParams_h
#define M3SYM_SysParams_h

#include <vector>
#include <list>

#include "common.h"
#include "Parser.h"

/// Struct to hold mechanical parameters for the system
struct MechParams {
    
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
struct ChemParams {
    
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
    
    
    //@{
    ///Extra motor parameters
    vector<short> motorNumHeadsMin = {};
    vector<short> motorNumHeadsMax = {};
    
    vector<double> motorStepSize = {};
    //@}
};

/// Struct to hold geometry parameters for the system
struct GeoParams {
    
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
    
    int cylinderIntSize = 0;
    double cylinderSize = 0;
    
    /// Minimum monomer length of a cylinder is preset
    int minCylinderIntSize = 3;
    double minCylinderSize = 0;
    //@}
    
};

/// Struct to hold Boundary parameters for the system
struct BoundParams {
    
    double BoundaryK = 0;
    double BScreenLength = 0;
    
    ///This parameter is preset
    double BCeiling = 100.0;
    
    double BoundaryCutoff = 0;
    double diameter = 0;
};

/// Struct to hold dynamic rate changing parameters
struct DyRateParams {
    
    /// Option for dynamic polymerization rate of filaments
    double dFilPolymerizationCharLength = 0.0;
    
    /// Option for dynamic unbinding rate of linkers
    vector<double> dLinkerUnbindingCharLength = {};
    vector<double> dLinkerUnbindingAmplitude = {};
    
    /// Option for dynamic unbinding rate of motors
    vector<double> dMotorUnbindingCharForce = {};
    
    /// Option for dynamic walking rate of motors
    vector<double> dMotorWalkingCharForce = {};
    
};

/// Static class that holds all simulation parameters,
/// initialized by the SystemParser
class SysParams {
friend class SystemParser;
    
#ifdef TESTING ///Public access if testing only
public:
#endif
    static MechParams MParams;    ///< The mechanical parameters
    static ChemParams CParams;    ///< The chemistry parameters
    static GeoParams GParams;     ///< The geometry parameters
    static BoundParams BParams;   ///< The boundary parameters
    static DyRateParams DRParams; ///< The dynamic rate parameters
    
public:
    //@{
    ///Const getter
    static const MechParams& Mechanics() {return MParams;}
    static const ChemParams& Chemistry() {return CParams;}
    static const GeoParams& Geometry() {return GParams;}
    static const BoundParams& Boundaries() {return BParams;}
    static const DyRateParams& DynamicRates() {return DRParams;}
    //@}
    
    //@{
    //Check for consistency of parameters. Done at runtime by the Controller.
    static bool checkChemParameters(ChemistryData& chem);
    static bool checkMechParameters(MechanicsFFType& mech);
    static bool checkDyRateParameters(DynamicRateTypes& dy);
    static bool checkGeoParameters();
    //@}
    
};

#endif
