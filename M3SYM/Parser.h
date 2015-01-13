
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

#ifndef M3SYM_Parser_h
#define M3SYM_Parser_h

#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <ios>

#include "common.h"

/// Struct to hold output types;
struct OutputTypes {
    
    //@{
    /// Possible output type
    bool basicSnapshot = false;
    bool birthTimes = false;
    bool forces = false;
    bool stresses = false;
    //@}
};

/// Struct to hold mechanics algorithm information
struct MechanicsAlgorithm {
    string ConjugateGradient = "";
    double gradientTolerance = 0.001;
    
    //not yet used
    string MD = "";
};

/// Struct to hold chemistry algorithm information
struct ChemistryAlgorithm {
    
    string algorithm = "";
    int numSteps = 0;
    int numStepsPerMech = 0;
    int numStepsPerSnapshot = 0;
    int numStepsPerNeighbor = 0;
};

/// Struct to hold Species and Reaction information
struct ChemistryData {
    
    /// Reaction happening between SpeciesBulk and SpeciesDiffusing ONLY
    vector<tuple<vector<string>,
                 vector<string>,
                        double>> genReactions = {};
    /// Reaction happening between SpeciesBulk ONLY
    vector<tuple<vector<string>,
                 vector<string>,
                        double>> bulkReactions = {};
    
    /// Filament nucleation reaction
    vector<tuple<vector<string>,
                 vector<string>,
                        double>> nucleationReactions = {};
    
    //@{
    /// Filament reactions
    /*!
     *  All Filament reactions are held using a vector containing a tuple with the 
     *  string of reactants, string of products, and the reaction rate.
     */
    /// Polymerization reactions
    vector<tuple<vector<string>,
                 vector<string>,
                        double>> polymerizationReactions = {};
    /// Depolymerization reactions
    vector<tuple<vector<string>,
                 vector<string>,
                        double>> depolymerizationReactions = {};
    /// Aging reactions
    vector<tuple<vector<string>,
                 vector<string>,
                        double>> agingReactions = {};
    
    /// Destruction reactions
    vector<tuple<vector<string>,
                 vector<string>,
                        double>> destructionReactions = {};
    
    /// Branching reactions
    vector<tuple<vector<string>,
                 vector<string>,
                        double, double>> branchingReactions = {};
    
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
    vector<tuple<vector<string>,
                 vector<string>,
                 double, double, double, double>> linkerReactions = {};
    /// MotorGhost reactions
    vector<tuple<vector<string>,
                 vector<string>,
                 double, double, double, double>> motorReactions = {};
    //@}
    
    /// MotorGhost walking reactions
    vector<tuple<vector<string>,
                 vector<string>,
                        double>> motorWalkingReactions = {};
    
    /// SpeciesBulk parsed, in the form of a tuple which contains the name and
    /// initial copy number, and whether this is a constant species
    vector<tuple<string, int, string>> speciesBulk = {};
    
    /// SpeicesDiffusing parsed, in the form of a tuple which contains name,
    /// initial copy number per compartment, and the rate of diffusion.
    vector<tuple<string, int, double>> speciesDiffusing = {};
    
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
};

/// Struct to hold the ForceField types
struct MechanicsFFType {
    
    //@{
    /// FilamentFF type
    string FStretchingType = "";
    string FBendingType = "";
    string FTwistingType = "";
    //@}
    
    //@{
    /// LinkerFF type
    string LStretchingType = "";
    string LBendingType = "";
    string LTwistingType = "";
    //@}
    
    //@{
    /// MotorFF type
    string MStretchingType = "";
    string MBendingType = "";
    string MTwistingType = "";
    //@}
    
    //@{
    /// BranchingFF type
    string BrStretchingType = "";
    string BrBendingType = "";
    string BrDihedralType = "";
    //@}
    
    /// VolumeFF type
    string VolumeFFType = "";
    
    /// BoundaryFF Type
    string BoundaryFFType = "";
    
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


/// Split a string by whitespace into generic type
template<typename T>
vector<T> split(const string& line) {
    istringstream is(line);
    return vector<T>(istream_iterator<T>(is), istream_iterator<T>());
}


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
            cout << "There was an error parsing file " << inputFileName <<
                ". Exiting" << endl;
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
    void readMechanicsParameters();
    void readChemistryParameters();
    void readGeometryParameters();
    void readBoundaryParameters();
    //@}
    
    //@{
    /// Algorithm parser
    MechanicsAlgorithm readMechanicsAlgorithm();
    ChemistryAlgorithm readChemistryAlgorithm();
    //@}
    
    //@{
    /// Type parser
    MechanicsFFType readMechanicsFFType();
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

/// Used to parse all chemical information
class ChemistryParser: public Parser {
    
public:
    ChemistryParser(string inputFileName) : Parser(inputFileName) {}
    ~ChemistryParser() {}
    
    /// Reads chemical reactions and species from input file. Returns a
    /// ChemistryData struct containing this data
    ChemistryData readChemistryInput();
};


#endif
