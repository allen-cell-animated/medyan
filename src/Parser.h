
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

#ifndef MEDYAN_Parser_h
#define MEDYAN_Parser_h

#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <ios>

#include "common.h"
#include "utility.h"

/// Struct to hold output types;
struct OutputTypes {
    
    //@{
    /// Possible output type
    
    /// A basic snapshot
    bool basicSnapshot = false;
    /// Birth times
    bool birthTimes = false;
    /// Forces
    bool forces = false;
    /// Tensions
    bool tensions = false;
    /// Chemistry
    bool chemistry = false;
    //@}
};

/// Struct to hold mechanics algorithm information
struct MechanicsAlgorithm {
    string ConjugateGradient = "";
    
    /// Tolerance and cg parameters
    double gradientTolerance = 1.0;
    double maxDistance = 1.0;
    double lambdaMax = 1.0;
    
    /// Not yet used
    string MD = "";
};

/// Struct to hold chemistry algorithm information
struct ChemistryAlgorithm {
    
    string algorithm = "";
    
    //@{
    /// User can specify total run time of the simulation, as well as
    /// frequency of snapshots, neighbor list updates and minimizations.
    double runTime = 0.0;
    
    double snapshotTime = 0.0;
    
    double minimizationTime = 0.0;
    double neighborListTime = 0.0;
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

/// Struct to hold Species and Reaction information
/// @note - all filament-related reactions and species are 2D vectors corresponding
///         to the filament type specified in the input file.
struct ChemistryData {
    
    /// Reaction happening between SpeciesBulk and SpeciesDiffusing ONLY
    vector<tuple<vector<string>, vector<string>, double, double>> genReactions = {};
    
    /// Reaction happening between SpeciesBulk ONLY
    vector<tuple<vector<string>, vector<string>, double>> bulkReactions = {};
    
    /// Filament nucleation reaction
    vector<vector<tuple<vector<string>, vector<string>, double>>> nucleationReactions;
    
    //@{
    /// Filament reactions
    /*!
     *  All Filament reactions are held using a vector containing a tuple with the 
     *  string of reactants, string of products, and the reaction rate.
     */
    /// Polymerization reactions
    vector<vector<tuple<vector<string>, vector<string>, double, double>>> polymerizationReactions;
    /// Depolymerization reactions
    vector<vector<tuple<vector<string>, vector<string>, double, double>>> depolymerizationReactions;
    /// Aging reactions
    vector<vector<tuple<vector<string>, vector<string>, double, double>>> agingReactions;
    /// Destruction reactions
    vector<vector<tuple<vector<string>, vector<string>, double>>> destructionReactions;
    
    /// Branching reactions
    /// This reaction also contains the off rate, and a string
    /// specifying the nucleation zone and relevant distance parameter
    vector<vector<tuple<vector<string>, vector<string>, double, double, string, double>>> branchingReactions;
    
    /// Severing reactions
    vector<vector<tuple<string, double>>> severingReactions;
    //@}
    
    //@{
    /// Cross Filament binding and unbinding reactions
    /*!
     *  All cross Filament reactions are held using a vector containing a tuple with 
     *  the string of reactants, string of products, the reaction rate, and binding 
     *  range.
     */
    /// Linker reactions
    vector<vector<tuple<vector<string>, vector<string>, double, double, double, double, double>>> linkerReactions;
    /// MotorGhost reactions
    vector<vector<tuple<vector<string>, vector<string>, double, double, double, double, double>>> motorReactions;
    //@}
    
    /// MotorGhost walking reactions
    vector<vector<tuple<vector<string>, vector<string>, double, double>>> motorWalkingReactions;
    
    /// SpeciesBulk parsed, in the form of a tuple which contains the name and
    /// initial copy number, release time, removal time, and CONST/REG qualifier
    vector<tuple<string, int, double, double, string>> speciesBulk = {};
    
    /// SpeicesDiffusing parsed, in the form of a tuple which contains name,
    /// initial copy number per compartment, the rate of diffusion, release time,
    /// AVG/REG qualifier, and number of events to average if applicable.
    vector<tuple<string, int, double, double, double, string, int>> speciesDiffusing = {};
    
    //@{
    /// Filament species parsed
    vector<vector<string>> speciesFilament;
    vector<vector<string>> speciesPlusEnd;
    vector<vector<string>> speciesMinusEnd;
    
    vector<vector<string>> speciesBound;
    vector<vector<string>> speciesLinker;
    vector<vector<string>> speciesMotor;
    vector<vector<string>> speciesBrancher;
    //@}
    
    //@{
    /// Binding sites parsed
    vector<string> B_BINDING_INDEX;
    vector<string> L_BINDING_INDEX;
    vector<string> M_BINDING_INDEX;
    //@}
    
    ///Constructor initializes memory of vector members
    ChemistryData()
    
    : nucleationReactions(MAX_FILAMENT_TYPES),
      polymerizationReactions(MAX_FILAMENT_TYPES),
      depolymerizationReactions(MAX_FILAMENT_TYPES),
      agingReactions(MAX_FILAMENT_TYPES),
      destructionReactions(MAX_FILAMENT_TYPES),
    
      branchingReactions(MAX_FILAMENT_TYPES),
      severingReactions(MAX_FILAMENT_TYPES),
      linkerReactions(MAX_FILAMENT_TYPES),
      motorReactions(MAX_FILAMENT_TYPES),
      motorWalkingReactions(MAX_FILAMENT_TYPES),
    
      speciesFilament(MAX_FILAMENT_TYPES),
      speciesPlusEnd(MAX_FILAMENT_TYPES),
      speciesMinusEnd(MAX_FILAMENT_TYPES),
      speciesBound(MAX_FILAMENT_TYPES),
      speciesLinker(MAX_FILAMENT_TYPES),
      speciesMotor(MAX_FILAMENT_TYPES),
      speciesBrancher(MAX_FILAMENT_TYPES),
    
      B_BINDING_INDEX(MAX_FILAMENT_TYPES),
      L_BINDING_INDEX(MAX_FILAMENT_TYPES),
      M_BINDING_INDEX(MAX_FILAMENT_TYPES) {}
    
};

/// Struct to hold the parameters of the Boundary
struct BoundaryType {
    
    string boundaryShape = "";
    string boundaryMove = "";
};

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
    
};

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
    //@}
};


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
};

/// Struct to hold chem setup information
struct ChemistrySetup {
    
    string inputFile = "";
};

/// Struct to hold Filament setup information
struct FilamentSetup {
    
    string inputFile = "";
    
    ///If want a random distribution, used if inputFile is left blank
    int numFilaments = 0;
    ///Filament length, in number of cylinders
    int filamentLength = 1;
    ///Filament type to create
    short filamentType = 0;
    ///Filament projection type.
    string projectionType="STRAIGHT";
    
    ///For resetting pin positions in restart phase
    string pinRestartFile = "";
};

/// Struct to hold Bubble setup information
struct BubbleSetup {
    
    string inputFile = "";
    
    ///If want a random distribution, used if inputFile is left blank
    int numBubbles = 0;
    ///Bubble type to create
    short bubbleType = 0;
};

/// A general parser
/*!
 *  A parser object, when initialized, opens an input file. Upon destruction, it 
 *  closes the file.
 */
class Parser {
protected:
    fstream _inputFile; ///< input file being used
    
public:
    Parser(string inputFileName) {
        _inputFile.open(inputFileName);
        if(!_inputFile.is_open()) {
            cout << "There was an error parsing file " << inputFileName
                 << ". Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        cout << "Loading file " << inputFileName << endl;
    }
    ~Parser() {_inputFile.close();}
};


/// To parse a system input file, initialized by the Controller.
class SystemParser : public Parser{
public:
    SystemParser(string inputFileName) : Parser(inputFileName) {}
    ~SystemParser() {}
    
    //@{
    /// Parameter parser. Reads input directly into system parameters
    /// @note - does not check for correctness and consistency here.
    void readMechParams();
    void readChemParams();
    void readGeoParams();
    void readBoundParams();
    void readDyRateParams();
    //@}
    
    //@{
    /// Algorithm parser
    MechanicsAlgorithm readMechanicsAlgorithm();
    ChemistryAlgorithm readChemistryAlgorithm();
    //@}
    
    //@{
    /// Type parser
    MechanicsFFType readMechanicsFFType();
    DynamicRateType readDynamicRateType();
    BoundaryType readBoundaryType();
    SpecialSetupType readSpecialSetupType();
    //@}
    
    /// Read Filament information
    FilamentSetup readFilamentSetup();
    
    /// Read Bubble information
    BubbleSetup readBubbleSetup();
    
    /// Chemistry information
    ChemistrySetup readChemistrySetup();
};

/// Used to parse initial Filament information, initialized by the Controller.
class FilamentParser : public Parser {
    
public:
    FilamentParser(string inputFileName) : Parser(inputFileName) {}
    ~FilamentParser() {}
    
    /// Reads filament input file. Returns a vector of tuples containing
    /// filament type and positions (start and end points).
    /// @note - Does not check for coordinate correctness.
     tuple< vector<tuple<short, vector<double>, vector<double>>> , vector<tuple<string, short, vector<vector<double>>>> , vector<tuple<string, short, vector<double>>> , vector<vector<double>> >  readFilaments();
};

/// Used to parse initial Bubble information, initialized by the Controller.
class BubbleParser : public Parser {
    
public:
    BubbleParser(string inputFileName) : Parser(inputFileName) {}
    ~BubbleParser() {}
    
    /// Reads bubble input file. Returns a vector of tuples containing
    /// bubble type and position.
    /// @note - Does not check for coordinate correctness.
    vector<tuple<short, vector<double>>> readBubbles();
};


/// Used to parse all chemical information, initialized by the Controller.
class ChemistryParser: public Parser {
    
public:
    ChemistryParser(string inputFileName) : Parser(inputFileName) {}
    ~ChemistryParser() {}
    
    /// Reads chemical reactions and species from input file. Returns a
    /// ChemistryData struct containing this data
    /// @note - this does not check for consistency and correctness, the only
    ///         sanity check here is that there are no duplicate species names.
    ChemistryData readChemistryInput();
};


/// Used to parse pin positions if needed upon restart
class PinRestartParser: public Parser {
    
public:
    PinRestartParser(string inputFileName) : Parser(inputFileName) {}
    ~PinRestartParser() {}
    
    /// Reads pin positions from file, and sets filaments
    void resetPins();
};



#endif
