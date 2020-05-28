
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_SysParams_h
#define MEDYAN_SysParams_h

#include <filesystem>
#include <vector>
#include <list>

#include "common.h"
#include "Parser.h"

//Did not minimize structure
#ifdef TRACKDIDNOTMINIMIZE
struct MinimizationParams{
    vector<floatingpoint> maxF;
    vector<floatingpoint> TotalE;
    vector<vector<floatingpoint>> Energyvec;
    vector<floatingpoint> Lambda;
    vector<floatingpoint> beta;
    vector<bool> safeModeORnot;
    vector<floatingpoint> tempEnergyvec;
    vector<vector<floatingpoint>> gradientvec;
};
#endif
/// Struct to hold mechanical parameters for the system
struct MechParams {

    /// Struct to hold the ForceField types
    struct MechanicsFFType {
        
        //@{
        /// FilamentFF type
        string FStretchingType = "";
        string FBendingType    = "";
        string FTwistingType   = "";
        //@}
        
        //@{
        /// LinkerFF type
        string LStretchingType = "";
        string LBendingType    = "";
        string LTwistingType   = "";
        //@}
        
        //@{
        /// MotorFF type
        string MStretchingType = "";
        string MBendingType    = "";
        string MTwistingType   = "";
        //@}
        
        //@{
        /// BranchingFF type
        string BrStretchingType = "";
        string BrBendingType    = "";
        string BrDihedralType   = "";
        string BrPositionType   = "";
        //@}
        
        /// VolumeFF type
        string VolumeFFType = "";
        
        /// BoundaryFF Type
        string BoundaryFFType = "";
        
        /// Bubble Type
        string BubbleFFType = "";
        
        /// MTOC Type
        string MTOCFFType = "";
        
        /// AFM Type
        string AFMFFType = "";
        
    };

    /// Struct to hold mechanics algorithm information
    struct MechanicsAlgorithm {
        string ConjugateGradient = "";
        
        /// Tolerance and cg parameters
        floatingpoint gradientTolerance = 1.0;
        floatingpoint maxDistance = 1.0;
        floatingpoint lambdaMax = 1.0;
        floatingpoint lambdarunningaverageprobability = 0.0;
        string linesearchalgorithm = "BACKTRACKING";
        
        /// Not yet used
        string MD = "";
    };

    MechanicsFFType       mechanicsFFType;
    MechanicsAlgorithm    mechanicsAlgorithm;

    //@{
    /// Filament parameter
    vector<floatingpoint> FStretchingK  = {};
    vector<floatingpoint> FBendingK     = {};
    vector<floatingpoint> FBendingTheta = {};
    vector<floatingpoint> FTwistingK    = {};
    vector<floatingpoint> FTwistingPhi  = {};
    //@}

    //@{
    /// Linker parameter
    vector<floatingpoint> LStretchingK  = {};
    vector<floatingpoint> LBendingK     = {};
    vector<floatingpoint> LBendingTheta = {};
    vector<floatingpoint> LTwistingK    = {};
    vector<floatingpoint> LTwistingPhi  = {};
    //@}

    //@{
    /// MotorGhost parameter
    vector<floatingpoint> MStretchingK  = {};
    vector<floatingpoint> MBendingK     = {};
    vector<floatingpoint> MBendingTheta = {};
    vector<floatingpoint> MTwistingK    = {};
    vector<floatingpoint> MTwistingPhi  = {};
    //@}

    //@{
    /// BranchingPoint parameter
    vector<floatingpoint> BrStretchingK  = {};
    vector<floatingpoint> BrStretchingL  = {};
    vector<floatingpoint> BrBendingK     = {};
    vector<floatingpoint> BrBendingTheta = {};
    vector<floatingpoint> BrDihedralK    = {};
    vector<floatingpoint> BrPositionK    = {};
    //@}

    //@{
    /// Volume parameter
    vector<floatingpoint> VolumeK = {};
    floatingpoint VolumeCutoff = 0.0;
    //@}

    //@{
    /// Bubble parameter
    vector<floatingpoint> BubbleK = {};
    vector<floatingpoint> BubbleRadius = {};
    vector<floatingpoint> BubbleScreenLength = {};
	vector<floatingpoint> MTOCBendingK = {};
    vector<floatingpoint> AFMBendingK = {};

	floatingpoint BubbleCutoff = 0.0;

    ///If using more than one bubble
    short numBubbleTypes = 1;
    //@}


    //@{
    /// SPECIAL MECHANICAL PROTOCOLS

    bool pinLowerBoundaryFilaments = false;
    floatingpoint pinFraction = 1.0; //test

    ///To pin filaments on boundary via an attractive potential
    bool pinBoundaryFilaments = false;
    floatingpoint pinDistance = 250; ///< 250nm pinning distance for now
    floatingpoint pinK = 0.0;       ///< Tethered stiffness
    floatingpoint pinTime = 0.0;    ///< Time at which to pin the filaments
    //@}

    //vectorization
    floatingpoint *coord;
    vector<vector<bool>> speciesboundvec;
    floatingpoint* cylsqmagnitudevector;
    vector<int> ncylvec;
    vector<int>bsoffsetvec;
    
    // parameters controlling the calculation of the Hessian matrix
    bool hessTracking = false;
    float hessDelta = 0.0001;
    bool denseEstimation = true;
    int hessSkip = 20;

    int sameFilBindSkip = 2;


};

/// Struct to hold chemistry parameters for the system
struct ChemParams {

    /// Struct to hold chemistry algorithm information
    struct ChemistryAlgorithm {
        
        string algorithm = "";
        
        //@{
        /// User can specify total run time of the simulation, as well as
        /// frequency of snapshots, neighbor list updates and minimizations.
        floatingpoint runTime = 0.0;
        
        floatingpoint snapshotTime = 0.0;
        floatingpoint datadumpTime = 1000.0;
        
        floatingpoint minimizationTime = 0.0;
        floatingpoint neighborListTime = 0.0;
        //@}
        
        //@{
        /// Can also specify a frequency in terms of number of chemical reaction steps
        /// Useful for smaller systems and debugging
        int runSteps = 0;
        
        int snapshotSteps = 0;
        
        int minimizationSteps = 0;
        int neighborListSteps = 0;
        //@}
    };

    /// Struct to hold chem setup information
    struct ChemistrySetup {
        
        std::filesystem::path inputFile;
    };

    ChemistryAlgorithm chemistryAlgorithm;
    ChemistrySetup     chemistrySetup;

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

    //Bindingsites are stored as packed 32bit integers. To ensure that there is adequate
    // memory space to store the binding sites, we need to shift based on the maximum
    // number of binding sites per cylinder.
    uint32_t shiftbybits = 0;
    uint32_t maxStableIndex = 0;

    //@{
    ///Extra motor parameters
    /// Vector corresponds to each filament type
    vector<short> motorNumHeadsMin = {};
    vector<short> motorNumHeadsMax = {};

    floatingpoint dutyRatio = 0.1;

    vector<floatingpoint> motorStepSize = {};
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
    floatingpoint makeFilamentsStaticTime = 0.0;

    /// Option to make Linker objects static after a certain time.
    /// This passivates any unbinding reactions, resulting in permanently
    /// bound linkers for the rest of simulation.
    bool makeLinkersStatic = false;

    floatingpoint makeLinkersStaticTime = 0.0;

    bool dissTracking = false;
    bool eventTracking = false;
    int linkerbindingskip = 2;
    
    
    /// Make (de)polymerization depends on stress
    bool makeRateDepend = false;
    double makeRateDependTime = 0.0;
    double makeRateDependForce = 0.0;
    
    /// Make (de)polymerization depends on stress
    bool makeAFM = false;
    double AFMStep1 = 0.0;
    double AFMStep2 = 0.0;
    double IterChange = 0.0;
    double StepTotal = 0.0;
    double StepTime = 0.0;
    float originalPolyPlusRate;
    

    //@}
#ifdef CUDAACCL_NL
    string bindingmanagerlist = "";
    vector<floatingpoint> bmanagerdistances = {};
#endif
};

/// Struct to hold geometry parameters for the system
struct GeoParams {

    //@{
    /// Geometry parameter
    short nDim = 3;

    int NX = 0;
    int NY = 0;
    int NZ = 0;

    floatingpoint compartmentSizeX = 0;
    floatingpoint compartmentSizeY = 0;
    floatingpoint compartmentSizeZ = 0;

    floatingpoint largestCompartmentSide = 0;

    floatingpoint largestCylinderSize = 0;

    vector<floatingpoint> monomerSize = {};

    ///Number of monomers in a cylinder
    vector<int> cylinderNumMon = {};

    vector<floatingpoint> cylinderSize = {};

    vector<floatingpoint> minCylinderSize = {};

    /// Minimum monomer length of a cylinder is preset
    int minCylinderNumMon = 3;
    //@}


};

/// Struct to hold Boundary parameters for the system
struct BoundParams {

    /// Struct to hold the parameters of the Boundary
    struct BoundaryType {
        
        string boundaryShape = "";
        vector<string> boundaryMove = {};
        //string scaleDiffusion = "";
    };

    BoundaryType  boundaryType;

    //@{
    /// Mechanical parameter
    floatingpoint BoundaryK = 0;
    floatingpoint BScreenLength = 0;
    //@}

    /// Cutoff for force calculation
    floatingpoint BoundaryCutoff = 0;
    floatingpoint diameter = 0;

    /// Moving speed (if any)
    floatingpoint moveSpeed = 0;

    //@{
    /// Moving times
    floatingpoint moveStartTime = 0;
    floatingpoint moveEndTime = 0;
    //@}
    int transfershareaxis=-1;       ///Axis along which activate/deactivate protocols should be executed.
    int planestomove = -1; //tracks if both (2), left/bottom/front (1), or
    // right/top/back(0) planes should be moved.
    vector<vector<floatingpoint>> fraccompartmentspan = { { 0, 0, 0 },
                                                   { 0.999, 0.999, 0.999 } };

};

/// Struct to hold dynamic rate changing parameters
struct DyRateParams {

    ///Struct to hold dynamic rate changer type
    struct DynamicRateType {
        
        ///Polymerization rate changing
        vector<string> dFPolymerizationType = {};
        
        ///Linker rate changing
        vector<string> dLUnbindingType = {};
        
        //@{
        ///Motor rate changing
        vector<string> dMUnbindingType = {};
        vector<string> dMWalkingType = {};

        //Qin----
        vector<string> dBUnbindingType = {};
        //@}
    };

    DynamicRateType dynamicRateType;

    /// Option for dynamic polymerization rate of filaments
    vector<floatingpoint> dFilPolymerizationCharLength = {};

    /// Option for dynamic unbinding rate of linkers
    vector<floatingpoint> dLinkerUnbindingCharLength = {};
    vector<floatingpoint> dLinkerUnbindingAmplitude = {};

    /// Option for dynamic unbinding rate of motors
    vector<floatingpoint> dMotorUnbindingCharForce = {};

    /// Option for dynamic walking rate of motors
    vector<floatingpoint> dMotorWalkingCharForce = {};

    //Qin
    /// Option for dynamic branching point unbinding rate
    vector<floatingpoint> dBranchUnbindingCharLength = {};

    /// Option for addinh dynamics branching point unbinding rate based on a
    // characteristic Force.
    vector<floatingpoint> dBranchUnbindingCharForce = {};

    /// Option for manual (de)polymerization rate changer
    /// Start time
    double manualCharStartTime = 100000.0;
    /// Plusend Polymerization Rate Ratio
    float manualPlusPolyRate = 1.0;
    /// Plusend Depolymerization Rate Ratio
    float manualPlusDepolyRate = 1.0;
    /// Minusend Polymerization Rate Ratio
    float manualMinusPolyRate = 1.0;
    /// Minusend Depolymerization Rate Ratio
    float manualMinusDepolyRate = 1.0;
};

struct SpecialParams {

    /// Struct to hold any special setup types
    struct SpecialSetupType {
        
        ///MTOC configuration
        bool mtoc = false;

        //@{
        ///MTOC Parameters
        short mtocFilamentType = 0;
        int mtocNumFilaments   = 0;
        int mtocFilamentLength = 0;
        short mtocBubbleType   = 0;
        //@}
        
        ///AFM configuration
        bool afm = false;
        
        //@{
        ///MTOC Parameters
        short afmFilamentType = 0;
        int afmNumFilaments   = 0;
        int afmFilamentLength = 0;
        short afmBubbleType   = 0;
        //@}
        vector<float> mtocInputCoordXYZ = {};
    };

    SpecialSetupType specialSetupType;

    /// Parameters for initializing MTOC attached filaments
    float mtocTheta1 = 0.0;
    float mtocTheta2 = 1.0;
    float mtocPhi1 = 0.0;
    float mtocPhi2 = 1.0;
};

using BubbleData = vector<tuple<short, vector<floatingpoint>>>;
using FilamentData = tuple<
    vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>>,
    vector<tuple<string, short, vector<vector<floatingpoint>>>>,
    vector<tuple<string, short, vector<floatingpoint>>>,
    vector<vector<floatingpoint>> >;

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
    #ifdef TRACKDIDNOTMINIMIZE
    static MinimizationParams MinParams;
	#endif

    static SpecialParams SParams; ///< Other parameters
    
public:
    //@{
#ifdef NLSTENCILLIST
    static short numcylcylNL;//Tracks total number of neighborlists in the system
#endif
    //counter to check excluded volume.
    static int exvolcounter[3]; //positions 0, 1, and 2 correspond to parallel,
    // in-plane and regular cases.
    static long long exvolcounterz[3];
    ///Const getter
    static bool RUNSTATE; //0 refers to restart phase and 1 refers to run phase.
    static bool USECHEMCOPYNUM; // if set to 0, restart file copy numbers are used. If
    // not, chemistry file copy numbers are used.
    static bool INITIALIZEDSTATUS; // true refers to sucessful initialization. false
    static bool DURINGCHEMISTRY; //true if MEDYAN is running chemistry, false otherwise.
    // corresponds to an on-going initialization state.

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
    static const SpecialParams& SpecialInputs() {return SParams;}
    //@}

    //@{
    //Check for consistency of parameters. Done at runtime by the Controller.
    static bool checkChemParameters(ChemistryData& chem);
    static bool checkMechParameters(MechParams::MechanicsFFType& mech);
    static bool checkDyRateParameters(DyRateParams::DynamicRateType& dy);
    static bool checkGeoParameters();
    //@}

    static void addChemParameters(ChemistryData& chem);
    #ifdef TRACKDIDNOTMINIMIZE
    static MinimizationParams& Mininimization() { return MinParams;}
	#endif

};

#endif
