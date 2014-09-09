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

///Struct to hold chemistry parameters (simple for now)
struct ChemistryParameters {
    
    std::string algorithm = "";
    std::string setup = "";
};


///Struct to hold boundary parameters (simple for now)
struct BoundaryParameters {
    
    std::string boundaryType = "";
};

///Struct to hold the read geometry parameters
struct GeometryParameters {
    
    short nDim = 0;
    
    int NX = 0;
    int NY = 0;
    int NZ = 0;
    
    double compartmentSizeX = 0;
    double compartmentSizeY = 0;
    double compartmentSizeZ = 0;
    
    double monomerSize = 0;
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
    std::string VolumeType = "";
};

///Struct to hold the read mechanics parameters
struct MechanicsParameters {
    
    ///Filament parameters
    double FStretchingK = 0;
    double FStretchingL = 0;
    double FBendingK = 0;
    double FBendingTheta = 0;
    double FTwistingK = 0;
    double FTwistingPhi = 0;
    
    ///Linker parameters
    double LStretchingK = 0;
    double LStretchingL = 0;
    double LBendingK = 0;
    double LBendingTheta = 0;
    double LTwistingK = 0;
    double LTwistingPhi = 0;
    
    ///Motor parameters
    double MStretchingK = 0;
    double MStretchingL = 0;
    double MBendingK = 0;
    double MBendingTheta = 0;
    double MTwistingK = 0;
    double MTwistingPhi = 0;
 
    ///Volume parameters
    double VolumeK = 0;
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
    ChemistryParameters readChemistryParameters();
    
    ///Boundary Parameters parser
    BoundaryParameters readBoundaryParameters();

    ///Geometry parameters parser
    GeometryParameters readGeometryParameters();
    
    ///Mechanics FF Types parser
    MechanicsFFType readMechanicsFFType();
    
    ///Mechanics Parameters parser
    MechanicsParameters readMechanicsParameters();
    
    ///Function to check consistency between all desired forcefields, boundaries, and constants
    bool checkInput(ChemistryParameters &CParams, BoundaryParameters &BParams, GeometryParameters &GParams,
                    MechanicsFFType &MTypes, MechanicsParameters &MParams);
    
};


#endif /* defined(__Cyto__Parser__) */
