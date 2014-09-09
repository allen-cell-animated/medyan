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


///Function to split a string by whitespace into generic type
template<typename T>
std::vector<T> split(const std::string& line) {
    std::istringstream is(line);
    return std::vector<T>(std::istream_iterator<T>(is), std::istream_iterator<T>());
}

///Parser class, to parse an input file
class Parser {
    
private:
    std::fstream _inputFile; ///< input file being used
    
public:

    Parser(std::string inputFileName) {
        _inputFile.open(inputFileName);
        if(!_inputFile.is_open()) {
            std::cout << "There was an error parsing input file " << inputFileName << ". Exiting" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Parsing input file " << inputFileName << std::endl;
    }
    
    ~Parser() {_inputFile.close();}
    
    ///Checks if mechanics is activated
    bool mechanics();
    
    ///Checks if chemistry is activated
    bool chemistry();
    
    ///Chemistry parameters parser
    ChemistryAlgorithm readChemistryAlgorithm();
    
    ///Mechanics parsers
    MechanicsAlgorithm readMechanicsAlgorithm();
    MechanicsFFType readMechanicsFFType();
    void readMechanicsParameters();
    
    ///Boundary parser
    BoundaryType readBoundaryType();
    void readBoundaryParameters();

    ///Geometry parser
    void readGeometryParameters();
    
};


#endif /* defined(__Cyto__Parser__) */
