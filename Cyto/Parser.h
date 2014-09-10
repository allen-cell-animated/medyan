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
#include "Mcommon.h"

///Struct to hold mechanics algorithm information
struct MechanicsAlgorithm {
    std::string algorithm = "";
};

///Struct to hold chemistry algorithm (simple for now)
struct ChemistryAlgorithm {
    
    std::string algorithm = "";
    std::string setup = "";
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
        std::cout << "Parsing file " << inputFileName << std::endl;
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

};

class FilamentParser : public Parser {
    
public:
    FilamentParser(std::string inputFileName) : Parser(inputFileName) {}
    ~FilamentParser() {}
    
    ///Reads filament input file. Returns a vector of filament positions,
    ///all containing starting and ending points.
    std::vector<std::vector<std::vector<double>>> readFilaments();
};



#endif /* defined(__Cyto__Parser__) */
