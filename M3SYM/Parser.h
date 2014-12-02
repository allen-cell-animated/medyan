
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

/// Struct to hold mechanics algorithm information
struct MechanicsAlgorithm {
    string ConjugateGradient = "";
    string MD = "";
};

/// Struct to hold chemistry algorithm information
struct ChemistryAlgorithm {
    
    string algorithm = "";
    int numSteps = 0;
    int numStepsPerMech = 0;
};

/// Struct to hold chemistry species and reaction information
struct ChemistryData {
    
    /// Reactions happening between bulk and diffusing species ONLY
    vector<tuple<vector<string>, vector<string>, double>> genReactions = {};
    /// Reactions happening between bulk species ONLY
    vector<tuple<vector<string>, vector<string>, double>> bulkReactions = {};
    
    //@{
    /// Filament reactions
    /*!
     *  All filament reactions are held using a vector containing a tuple with the string
     *  of reactants, string of products, and the reaction rate.
     */
    /// Polymerization reactions
    vector<tuple<vector<string>, vector<string>, double>> polymerizationReactions = {};
    /// Depolymerization reactions
    vector<tuple<vector<string>, vector<string>, double>> depolymerizationReactions = {};
    /// Binding reactions
    vector<tuple<vector<string>, vector<string>, double>> bindingReactions = {};
    /// Unbinding reactions
    vector<tuple<vector<string>, vector<string>, double>> unbindingReactions = {};
    //@}
    
    //@{
    /// Cross filament binding reactions
    /*!
     *  All cross filament reactions are held using a vector containing a tuple with the string
     *  of reactants, string of products, the reaction rate, and binding range.
     */
    /// Linker binding reactions
    vector<tuple<vector<string>, vector<string>, double, double, double>> linkerBindingReactions = {};
    /// Motor binding reactions
    vector<tuple<vector<string>, vector<string>, double, double, double>> motorBindingReactions = {};
    //@}
    
    /// Motor walking reactions
    vector<tuple<vector<string>, vector<string>, double>> motorWalkingReactions = {};
    
    ///[SpeciesBulk] (@ref SpeciesBulk) parsed, in the form of a tuple which contains the name and
    /// initial copy number.
    vector<tuple<string, int>> speciesBulk = {};
    
    /// [SpeciesDiffusing] (@ref SpeciesDiffusing) parsed, in the form of a tuple which contains name,
    /// initial copy number per compartment, and the rate of diffusion.
    vector<tuple<string, int, double>> speciesDiffusing = {};
    
    //@{
    /// Filament species parsed
    vector<string> speciesFilament = {};
    vector<string> speciesBound = {};
    vector<string> speciesLinker = {};
    vector<string> speciesMotor = {};
    vector<string> speciesPlusEnd = {};
    vector<string> speciesMinusEnd = {};
    //@}
    
};

/// Struct to hold the parameters of the boundary
struct BoundaryType {
    
    string boundaryShape = "";
};

/// Struct to hold the force field types
struct MechanicsFFType {
    
    //@{
    /// Filament FF types
    string FStretchingType = "";
    string FBendingType = "";
    string FTwistingType = "";
    //@}
    
    //@{
    /// Linker FF types
    string LStretchingType = "";
    string LBendingType = "";
    string LTwistingType = "";
    //@}
    
    //@{
    /// Motor FF type
    string MStretchingType = "";
    string MBendingType = "";
    string MTwistingType = "";
    //@}
    
    /// Volume FF type
    string VolumeFFType = "";
    
    /// Boundary FF Type
    string BoundaryFFType = "";
    
};

/// Struct to hold chem setup information
struct ChemistrySetup {
    
    string inputFile = "";
};

/// Struct to hold filament setup information
struct FilamentSetup {
    
    string inputFile = "";
    
    ///If want a random distribution, used if inputFile is left blank
    int numFilaments = 0;
};


/// Split a string by whitespace into generic type
template<typename T>
vector<T> split(const string& line) {
    istringstream is(line);
    return vector<T>(istream_iterator<T>(is), istream_iterator<T>());
}


/// A general parser
/*!
 *  A parser object, when initialized, opens an input file. Upon destruction, it closes the file.
 */
class Parser {
protected:
    fstream _inputFile; ///< input file being used
    
public:
    Parser(string inputFileName) {
        _inputFile.open(inputFileName);
        if(!_inputFile.is_open()) {
            cout << "There was an error parsing file " << inputFileName << ". Exiting" << endl;
            exit(EXIT_FAILURE);
        }
        cout << "Loading file " << inputFileName << endl;
    }
    ~Parser() {_inputFile.close();}
};


/// To parse a system input file, initialized by the [Controller] (@ref Controller).
class SystemParser : public Parser{
public:
    SystemParser(string inputFileName) : Parser(inputFileName) {}
    ~SystemParser() {}
    
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

/// FilamentParser class is used to parse initial filament information, initialized by the [Controller] (@ref Controller).
class FilamentParser : public Parser {
    
public:
    FilamentParser(string inputFileName) : Parser(inputFileName) {}
    ~FilamentParser() {}
    
    /// Reads filament input file. Returns a vector of filament positions,
    /// all containing starting and ending points.
    vector<vector<vector<double>>> readFilaments();
};

/// ChemistryParser class is used to parse all chemical information
class ChemistryParser: public Parser {
    
public:
    ChemistryParser(string inputFileName) : Parser(inputFileName) {}
    ~ChemistryParser() {}
    
    /// Reads chemical reactions and species from input file. Returns a
    /// chemistry setup struct containing this data
    ChemistryData readChemistryInput();
};




#endif
