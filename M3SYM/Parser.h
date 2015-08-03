
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

#ifndef M3SYM_Parser_h
#define M3SYM_Parser_h

#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <ios>

#include "common.h"

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
    /// Stresses
    bool stresses = false;
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
    
    /// Not yet used
    string MD = "";
};

/// Struct to hold chemistry algorithm information
struct ChemistryAlgorithm {
    
    string algorithm = "";
    
    //@{
    /// User can specify either the total number of chemical steps to perform,
    /// or the total run time of the simulation.
    int numTotalSteps = 0;
    double runTime = 0.0;
    
    int numStepsPerSnapshot = 0;
    double snapshotTime = 0;
    //@}
    
    int numChemSteps = 0; ///< Specifying number of chemical steps at a time
    int numStepsPerNeighbor = 0; ///< Number of chemical steps per a neighbor
                                 ///< list update. Will affect efficiency.
};

/// Struct to hold Species and Reaction information
struct ChemistryData {
    
    /// Reaction happening between SpeciesBulk and SpeciesDiffusing ONLY
    vector<tuple<vector<string>, vector<string>, double>> genReactions = {};
    
    /// Reaction happening between SpeciesBulk ONLY
    vector<tuple<vector<string>, vector<string>, double>> bulkReactions = {};
    
    /// Filament nucleation reaction
    vector<tuple<vector<string>, vector<string>, double>> nucleationReactions = {};
    
    //@{
    /// Filament reactions
    /*!
     *  All Filament reactions are held using a vector containing a tuple with the 
     *  string of reactants, string of products, and the reaction rate.
     */
    /// Polymerization reactions
    vector<tuple<vector<string>, vector<string>, double>> polymerizationReactions = {};
    
    /// Depolymerization reactions
    vector<tuple<vector<string>, vector<string>, double>> depolymerizationReactions = {};
    
    /// Aging reactions
    vector<tuple<vector<string>, vector<string>, double>> agingReactions = {};
    
    /// Destruction reactions
    vector<tuple<vector<string>, vector<string>, double>> destructionReactions = {};
    
    /// Branching reactions
    /// This reaction also contains the off rate, and a string
    /// specifying the nucleation zone and relevant distance parameter
    vector<tuple<vector<string>, vector<string>, double, double, string, double>> branchingReactions = {};
    
    /// Severing reactions
    vector<tuple<string, double>> severingReactions = {};
    //@}
    
    //@{
    /// Cross Filament binding and unbinding reactions
    /*!
     *  All cross Filament reactions are held using a vector containing a tuple with 
     *  the string of reactants, string of products, the reaction rate, and binding 
     *  range.
     */
    /// Linker reactions
    vector<tuple<vector<string>, vector<string>, double, double, double, double>> linkerReactions = {};
    /// MotorGhost reactions
    vector<tuple<vector<string>, vector<string>, double, double, double, double>> motorReactions = {};
    //@}
    
    /// MotorGhost walking reactions
    vector<tuple<vector<string>, vector<string>, double>> motorWalkingReactions = {};
    
    /// SpeciesBulk parsed, in the form of a tuple which contains the name and
    /// initial copy number, release time, and CONST/REG qualifier
    vector<tuple<string, int, double, string>> speciesBulk = {};
    
    /// SpeicesDiffusing parsed, in the form of a tuple which contains name,
    /// initial copy number per compartment, the rate of diffusion, release time,
    /// AVG/REG qualifier, and number of events to average if applicable.
    vector<tuple<string, int, double, double, string, int>> speciesDiffusing = {};
    
    //@{
    /// Filament species parsed
    vector<string> speciesFilament = {};
    vector<string> speciesPlusEnd  = {};
    vector<string> speciesMinusEnd = {};
    
    vector<string> speciesBound  = {};
    vector<string> speciesLinker = {};
    vector<string> speciesMotor  = {};
    vector<string> speciesBrancher = {};
    //@}
    
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
    
};

///Struct to hold dynamic rate changer types
struct DynamicRateTypes {
    
    ///Polymerization rate changing
    string dFPolymerizationType = "";
    
    ///Linker rate changing
    vector<string> dLUnbindingType = {};
    
    //@{
    ///Motor rate changing
    vector<string> dMUnbindingType = {};
    vector<string> dMWalkingType = {};
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
    
    /// Output parser
    OutputTypes readOutputTypes();
    
    //@{
    /// Parameter parser. Reads input directly into system parameters
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
    DynamicRateTypes readDynamicRateTypes();
    BoundaryType readBoundaryType();
    //@}
    
    /// Read Filament information
    FilamentSetup readFilamentSetup();
    
    /// Chemistry information
    ChemistrySetup readChemistrySetup();
};

/// Used to parse initial Filament information, initialized by the Controller.
class FilamentParser : public Parser {
    
public:
    FilamentParser(string inputFileName) : Parser(inputFileName) {}
    ~FilamentParser() {}
    
    /// Reads filament input file. Returns a vector of Filament positions,
    /// all containing starting and ending points.
    vector<vector<vector<double>>> readFilaments();
};

/// Used to parse all chemical information, initialized by the Controller.
class ChemistryParser: public Parser {
    
public:
    ChemistryParser(string inputFileName) : Parser(inputFileName) {}
    ~ChemistryParser() {}
    
    /// Reads chemical reactions and species from input file. Returns a
    /// ChemistryData struct containing this data
    ChemistryData readChemistryInput();
};


#endif
