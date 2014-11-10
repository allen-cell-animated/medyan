//
//  Parser.h
//  Cyto
//
//  Created by James Komianos on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Parser__
#define __Cyto__Parser__

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <ios>
#include "common.h"

///Struct to hold mechanics algorithm information
struct MechanicsAlgorithm {
    string ConjugateGradient = "";
    string MD = "";
};

///Struct to hold chemistry algorithm (simple for now)
struct ChemistryAlgorithm {
    
    string algorithm = "";
    int numSteps = 0;
    int numStepsPerMech = 0;
};

///Struct to hold chemical info in system
struct ChemistryData {
    
    ///Reactions happening between bulk and diffusing species ONLY
    vector<tuple<vector<string>, vector<string>, double>> genReactions = {};
    ///Reactions happening between bulk species ONLY
    vector<tuple<vector<string>, vector<string>, double>> bulkReactions = {};
    
    ///Filament reactions
    ///Polymerization reactions
    vector<tuple<vector<string>, vector<string>, double>> polymerizationReactions = {};
    ///Depolymerization reactions
    vector<tuple<vector<string>, vector<string>, double>> depolymerizationReactions = {};
    ///Binding reactions
    vector<tuple<vector<string>, vector<string>, double>> bindingReactions = {};
    ///unbinding reactions
    vector<tuple<vector<string>, vector<string>, double>> unbindingReactions = {};
    
    ///Cross filament binding reactions
    ///Linker binding
    vector<tuple<vector<string>, vector<string>, double, double, double>> linkerBindingReactions = {};
    ///Motor binding
    vector<tuple<vector<string>, vector<string>, double, double, double>> motorBindingReactions = {};

    ///Motor walking reactions
    vector<tuple<vector<string>, vector<string>, double>> motorWalkingReactions = {};
    
    ///SpeciesBulk parsed, in the form of a tuple which contains the name and initial copy number
    vector<tuple<string, int>> speciesBulk = {};
    
    ///SpeciesDiffusing parsed, in the form of a tuple which contains name, initial copy number
    /// per compartment, and the rate of diffusion.
    vector<tuple<string, int, double>> speciesDiffusing = {};
    
    ///All filament species parsed
    vector<string> speciesFilament = {};
    vector<string> speciesBound = {};
    vector<string> speciesLinker = {};
    vector<string> speciesMotor = {};
    vector<string> speciesPlusEnd = {};
    vector<string> speciesMinusEnd = {};
    
};

///Struct to hold boundary parameters (simple for now)
struct BoundaryType {
    
    string boundaryShape = "";
};

///Struct to hold mechanics ff type
struct MechanicsFFType {
    ///Filament FF types
    string FStretchingType = "";
    string FBendingType = "";
    string FTwistingType = "";
    
    ///Linker FF types
    string LStretchingType = "";
    string LBendingType = "";
    string LTwistingType = "";
    
    ///Motor FF type
    string MStretchingType = "";
    string MBendingType = "";
    string MTwistingType = "";
    
    ///Volume FF type
    string VolumeFFType = "";
    
    ///Boundary FF Type
    string BoundaryFFType = "";
    
};

struct ChemistrySetup {
    
    ///If Reading in
    string inputFile = "";
};


struct FilamentSetup {
    
    ///If reading in
    string inputFile = "";
    
    ///If want a random distribution
    int numFilaments = 0;
};


///Function to split a string by whitespace into generic type
template<typename T>
vector<T> split(const string& line) {
    istringstream is(line);
    return vector<T>(istream_iterator<T>(is), istream_iterator<T>());
}


///General parser class for io
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


///SystemParser class, to parse a system input file
class SystemParser : public Parser{
public:
    SystemParser(string inputFileName) : Parser(inputFileName) {}
    ~SystemParser() {}
    
    ///Checks if mechanics is activated
    bool mechanics();
    
    ///Checks if chemistry is activated
    bool chemistry();
    
    ///Parameter parsers. These read input directly into system parameters
    void readMechanicsParameters();
    void readChemistryParameters();
    void readGeometryParameters();
    void readBoundaryParameters();
    
    ///Algorithm parsers
    MechanicsAlgorithm readMechanicsAlgorithm();
    ChemistryAlgorithm readChemistryAlgorithm();
    
    ///Type parsers
    MechanicsFFType readMechanicsFFType();
    BoundaryType readBoundaryType();
    
    //Filament information
    FilamentSetup readFilamentSetup();
    
    //Chemistry information
    ChemistrySetup readChemistrySetup();

};

class FilamentParser : public Parser {
    
public:
    FilamentParser(string inputFileName) : Parser(inputFileName) {}
    ~FilamentParser() {}
    
    ///Reads filament input file. Returns a vector of filament positions,
    ///all containing starting and ending points.
    vector<vector<vector<double>>> readFilaments();
};


class ChemistryParser: public Parser {
    
public:
    ChemistryParser(string inputFileName) : Parser(inputFileName) {}
    ~ChemistryParser() {}
    
    ///Reads chemical reactions and species from input file. Returns a
    ///chemistry setup struct containing this data
    ChemistryData readChemistryInput();
};




#endif /* defined(__Cyto__Parser__) */
