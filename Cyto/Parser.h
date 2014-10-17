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
    std::string ConjugateGradient = "";
    std::string MD = "";
};

///Struct to hold chemistry algorithm (simple for now)
struct ChemistryAlgorithm {
    
    std::string algorithm = "";
    int numSteps = 0;
    int numStepsPerMech = 0;
};

///Struct to hold chemical species and reactions in system
struct ChemistrySpeciesAndReactions {
    
    ///Reactions parsed, in the form of a tuple which contains reactants, products and rate
    ///For cross-filament reactions, also contains rMin, rMax denoting the reaction range
    
    ///Reactions happening between bulk and diffusing species ONLY
    std::vector<std::tuple<std::vector<std::string>, std::vector<std::string>, double>> generalReactions = {};
    ///Polymerization reactions
    std::vector<std::tuple<std::vector<std::string>, std::vector<std::string>, double>> polymerizationReactions = {};
    ///Depolymerization reactions
    std::vector<std::tuple<std::vector<std::string>, std::vector<std::string>, double>> depolymerizationReactions = {};
    ///Binding reactions
    std::vector<std::tuple<std::vector<std::string>, std::vector<std::string>, double>> bindingReactions = {};
    ///unbinding reactions
    std::vector<std::tuple<std::vector<std::string>, std::vector<std::string>, double>> unbindingReactions = {};
    
    ///Cross filament binding reactions
    std::vector<std::tuple<std::vector<std::string>, std::vector<std::string>, double, double, double>> crossFilamentBindingReactions = {};
    
    ///SpeciesBulk parsed, in the form of a tuple which contains the name and initial copy number
    std::vector<std::tuple<std::string, int>> speciesBulk = {};
    
    ///SpeciesDiffusing parsed, in the form of a tuple which contains name, initial copy number
    /// per compartment, and the rate of diffusion.
    std::vector<std::tuple<std::string, int, double>> speciesDiffusing = {};
    
    ///All filament species parsed
    std::vector<std::string> speciesFilament = {};
    std::vector<std::string> speciesBound = {};
    std::vector<std::string> speciesLinker = {};
    std::vector<std::string> speciesMotor = {};
    std::vector<std::string> speciesPlusEnd = {};
    std::vector<std::string> speciesMinusEnd = {};
    
};

///Struct to hold boundary parameters (simple for now)
struct BoundaryType {
    
    std::string boundaryShape = "";
};

///Struct to hold mechanics ff type
struct MechanicsFFType {
    ///Filament FF types
    std::string FStretchingType = "";
    std::string FBendingType = "";
    std::string FTwistingType = "";
    
    ///Linker FF types
    std::string LStretchingType = "";
    std::string LBendingType = "";
    std::string LTwistingType = "";
    
    ///Motor FF type
    std::string MStretchingType = "";
    std::string MBendingType = "";
    std::string MTwistingType = "";
    
    ///Volume FF type
    std::string VolumeFFType = "";
    
    ///Boundary FF Type
    std::string BoundaryFFType = "";
    
};

struct ChemistrySetup {
    
    ///If Reading in
    std::string inputFile = "";
};


struct FilamentSetup {
    
    ///If reading in
    std::string inputFile = "";
    
    ///If want a random distribution
    int numFilaments = 0;
};


///Function to split a string by whitespace into generic type
template<typename T>
std::vector<T> split(const std::string& line) {
    std::istringstream is(line);
    return std::vector<T>(std::istream_iterator<T>(is), std::istream_iterator<T>());
}


///General parser class for io
class Parser {
protected:
    std::fstream _inputFile; ///< input file being used
    
public:
    Parser(std::string inputFileName) {
        _inputFile.open(inputFileName);
        if(!_inputFile.is_open()) {
            std::cout << "There was an error parsing file " << inputFileName << ". Exiting" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Loading file " << inputFileName << std::endl;
    }
    ~Parser() {_inputFile.close();}
};


///SystemParser class, to parse a system input file
class SystemParser : public Parser{
public:
    SystemParser(std::string inputFileName) : Parser(inputFileName) {}
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
    FilamentParser(std::string inputFileName) : Parser(inputFileName) {}
    ~FilamentParser() {}
    
    ///Reads filament input file. Returns a vector of filament positions,
    ///all containing starting and ending points.
    std::vector<std::vector<std::vector<double>>> readFilaments();
};


class ChemistryParser: public Parser {
    
public:
    ChemistryParser(std::string inputFileName) : Parser(inputFileName) {}
    ~ChemistryParser() {}
    
    ///Reads chemical reactions and species from input file. Returns a
    ///chemistry setup struct containing this data
    ChemistrySpeciesAndReactions readChemistryInput();
};




#endif /* defined(__Cyto__Parser__) */
