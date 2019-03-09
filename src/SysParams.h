
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_SysParams_h
#define MEDYAN_SysParams_h

#include <vector>
#include <list>

#include "common.h"
#include "Parser.h"

/// Struct to hold mechanical parameters for the system
struct MechParams {
    
    //@{
    /// Filament parameter
    vector<double> FStretchingK  = {};
    vector<double> FBendingK     = {};
    vector<double> FBendingTheta = {};
    vector<double> FTwistingK    = {};
    vector<double> FTwistingPhi  = {};
    //@}
    
    //@{
    /// Linker parameter
    vector<double> LStretchingK  = {};
    vector<double> LBendingK     = {};
    vector<double> LBendingTheta = {};
    vector<double> LTwistingK    = {};
    vector<double> LTwistingPhi  = {};
    //@}
    
    //@{
    /// MotorGhost parameter
    vector<double> MStretchingK  = {};
    vector<double> MBendingK     = {};
    vector<double> MBendingTheta = {};
    vector<double> MTwistingK    = {};
    vector<double> MTwistingPhi  = {};
    //@}
    
    //@{
    /// BranchingPoint parameter
    vector<double> BrStretchingK  = {};
    vector<double> BrStretchingL  = {};
    vector<double> BrBendingK     = {};
    vector<double> BrBendingTheta = {};
    vector<double> BrDihedralK    = {};
    vector<double> BrPositionK    = {};
    //@}
    
    //@{
    /// Volume parameter
    vector<double> VolumeK = {};
    double VolumeCutoff = 0.0;
    //@}
    
    //@{
    /// Bubble parameter
    vector<double> BubbleK = {};
    vector<double> BubbleRadius = {};
    vector<double> BubbleScreenLength = {};
    
    double BubbleCutoff = 0.0;
    
    ///If using more than one bubble
    short numBubbleTypes = 1;
    //@}
    
    
    //@{
    /// SPECIAL MECHANICAL PROTOCOLS
    
    //Qin
    bool pinLowerBoundaryFilaments = false;
    double pinFraction = 1.0; //test 
    
    ///To pin filaments on boundary via an attractive potential
    bool pinBoundaryFilaments = false;
    double pinDistance = 250; ///< 250nm pinning distance for now
    double pinK = 0.0;       ///< Tethered stiffness
    double pinTime = 0.0;    ///< Time at which to pin the filaments
    //int transfershareaxis=-1;       ///Axis along which activate/deactivate protocols
/// should be executed.
    //@}

    //vectorization
    double *coord;
    vector<vector<bool>> speciesboundvec;
    double* cylsqmagnitudevector;
    vector<int> ncylvec;
    vector<int>bsoffsetvec;

    
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
    /// Vector corresponds to number of species on each filament type
    vector<short> numFilamentSpecies = {};
    vector<short> numPlusEndSpecies  = {};
    vector<short> numMinusEndSpecies = {};
    vector<short> numBoundSpecies    = {};
    vector<short> numLinkerSpecies   = {};
    vector<short> numMotorSpecies    = {};
    vector<short> numBrancherSpecies = {};
    //@}

    /// Number of different filament types
    short numFilaments = 1;
    
    /// Number of binding sites per cylinder
    /// Vector corresponds to each filament type
    vector<short> numBindingSites = {};

    short maxbindingsitespercylinder = 0;
    short minbindingsitespercylinder = 32767;
    
    //@{
    ///Extra motor parameters
    /// Vector corresponds to each filament type
    vector<short> motorNumHeadsMin = {};
    vector<short> motorNumHeadsMax = {};
    
    vector<double> motorStepSize = {};
    //@}
    
    /// Binding sites on filaments
    /// 2D Vector corresponds to each filament type
    vector<vector<short>> bindingSites = {};
    
    //@{
    /// Positions of all bound molecules in species vectors
    /// Vector corresponds to each filament type
    vector<short> brancherBoundIndex = vector<short>(MAX_FILAMENT_TYPES);
    vector<short> linkerBoundIndex   = vector<short>(MAX_FILAMENT_TYPES);
    vector<short> motorBoundIndex    = vector<short>(MAX_FILAMENT_TYPES);
    
    vector<vector<short>> bindingIndices = vector<vector<short>>(MAX_FILAMENT_TYPES);
    //@}
    
    
    //@{
    /// SPECIAL CHEMICAL PROTOCOLS
    
    /// Option to make Filament objects static after a certain time.
    /// This passivates any polymerization or depolymerization
    /// reactions, resulting in constant length filaments for the rest of simulation.
    bool makeFilamentsStatic = false;
    double makeFilamentsStaticTime = 0.0;
    
    /// Option to make Linker objects static after a certain time.
    /// This passivates any unbinding reactions, resulting in permanently
    /// bound linkers for the rest of simulation.
    bool makeLinkersStatic = false;
    double makeLinkersStaticTime = 0.0;
    
    //@}
#ifdef CUDAACCL_NL
    string bindingmanagerlist = "";
    vector<double> bmanagerdistances = {};
#endif
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

    double largestCylinderSize = 0;
    
    vector<double> monomerSize = {};
    
    ///Number of monomers in a cylinder
    vector<int> cylinderNumMon = {};
    
    vector<double> cylinderSize = {};
    
    vector<double> minCylinderSize = {};
    
    /// Minimum monomer length of a cylinder is preset
    int minCylinderNumMon = 3;
    //@}


};

/// Struct to hold Boundary parameters for the system
struct BoundParams {
    
    //@{
    /// Mechanical parameter
    double BoundaryK = 0;
    double BScreenLength = 0;
    //@}
    
    /// Cutoff for force calculation
    double BoundaryCutoff = 0;
    double diameter = 0;
    
    /// Moving speed (if any)
    double moveSpeed = 0;
    
    //@{
    /// Moving times
    double moveStartTime = 0;
    double moveEndTime = 0;
    //@}
    int transfershareaxis=-1;       ///Axis along which activate/deactivate protocols should be executed.
    int planestomove = -1; //tracks if both (2), left/bottom/front (1), or
    // right/top/back(0) planes should be moved.
    vector<vector<double>> fraccompartmentspan = { { 0, 0, 0 },
                                                   { 0.999, 0.999, 0.999 } };
  
};

/// Struct to hold dynamic rate changing parameters
struct DyRateParams {
    
    /// Option for dynamic polymerization rate of filaments
    vector<double> dFilPolymerizationCharLength = {};
    
    /// Option for dynamic unbinding rate of linkers
    vector<double> dLinkerUnbindingCharLength = {};
    vector<double> dLinkerUnbindingAmplitude = {};
    
    /// Option for dynamic unbinding rate of motors
    vector<double> dMotorUnbindingCharForce = {};
    
    /// Option for dynamic walking rate of motors
    vector<double> dMotorWalkingCharForce = {};
    
    //Qin
    /// Option for dynamic branching point unbinding rate
    vector<double> dBranchUnbindingCharLength = {};
};

/// Static class that holds all simulation parameters,
/// initialized by the SystemParser
class SysParams {
friend class Controller;
friend class SystemParser;
friend class ChemManager;
friend class SubSystem;
friend class Cylinder;
    
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
#ifdef NLSTENCILLIST
    static short numcylcylNL;
#endif
    //counter to check excluded volume.
    static int exvolcounter[3]; //positions 0, 1, and 2 correspond to parallel,
    // in-plane and regular cases.
    static long long exvolcounterz[3];
    ///Const getter
    static bool RUNSTATE; //0 refers to restart phase and 1 refers to run phase.
    static bool INITIALIZEDSTATUS; // true refers to sucessful initialization. false
    // corresponds to an on-going initialization state.
    static int numthreads;
    //aravind July11,2016
    static vector<float> MUBBareRate;
    static vector<float> LUBBareRate;
    static vector<float> BUBBareRate;
    //@
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
    static bool checkDyRateParameters(DynamicRateType& dy);
    static bool checkGeoParameters();
    //@}
    
};

#endif
