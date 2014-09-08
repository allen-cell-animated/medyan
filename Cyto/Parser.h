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
    }
    
    ~Parser() {_inputFile.close();}
    
    
    ///GEOMETRY CONSTANT PARSER
    std::vector<std::vector<double>> readGeometryParameters() {
        
        std::vector<std::vector<double>> returnVector;
        
        std::vector<double> gridTemp;
        std::vector<double> compartmentTemp;
        
        //find grid size lines
        std::string line;
        while(getline(_inputFile, line)) {
            if (line.find("NX") != std::string::npos
                || line.find("NY") != std::string::npos
                || line.find("NZ") != std::string::npos) {
                
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() != 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                gridTemp.push_back(double(std::atoi(lineVector[1].c_str())));
            }

            else if (line.find("COMPARTMENTSIZEX") != std::string::npos
                || line.find("COMPARTMENTSIZEY") != std::string::npos
                || line.find("COMPARTMENTSIZEZ") != std::string::npos) {
                
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() != 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                compartmentTemp.push_back(double(std::atoi(lineVector[1].c_str())));
            }
            
            else {}
        }

        return std::vector<std::vector<double>>(returnVector.begin(), returnVector.end());
    }
    
    ///MECHANICS CONSTANT PARSER
    std::vector<std::vector<double>> readMechanicsParameters() {
        
        std::vector<std::vector<double>> returnVector;
        
        ///Filament parameters
        std::vector<double> fStretchingParams;
        std::vector<double> fBendingParams;
        std::vector<double> fTwistingParams;
        
        ///Linker parameters
        std::vector<double> lStretchingParams;
        std::vector<double> lBendingParams;
        std::vector<double> lTwistingParams;
        
        ///Motor parameters
        std::vector<double> mStretchingParams;
        std::vector<double> mBendingParams;
        std::vector<double> mTwistingParams;
        
        ///Volume parameters
        std::vector<double> volumeParams;
        
        
        std::string line;
        while(getline(_inputFile, line)) {
            if (line.find("FSTRETCHINGK") != std::string::npos
                || line.find("FSTRETCHINGL") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    fStretchingParams.push_back(double(std::atoi(lineVector[1].c_str())));
                }
            }
            else if (line.find("FBENDINGK") != std::string::npos
                || line.find("FBENDINGTHETA") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    fBendingParams.push_back(double(std::atoi(lineVector[1].c_str())));
                }
            }
            else if (line.find("FTWISTINGK") != std::string::npos
                     || line.find("FTWISTINGPHI") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    fTwistingParams.push_back(double(std::atoi(lineVector[1].c_str())));
                }
            }
            else if (line.find("LSTRETCHINGK") != std::string::npos
                || line.find("LSTRETCHINGL") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    lStretchingParams.push_back(double(std::atoi(lineVector[1].c_str())));
                }
            }
            else if (line.find("LBENDINGK") != std::string::npos
                     || line.find("LBENDINGTHETA") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    lBendingParams.push_back(double(std::atoi(lineVector[1].c_str())));
                }
            }
            else if (line.find("LTWISTINGK") != std::string::npos
                     || line.find("LTWISTINGPHI") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    lTwistingParams.push_back(double(std::atoi(lineVector[1].c_str())));
                }
            }
            else if (line.find("MSTRETCHINGK") != std::string::npos
                     || line.find("MSTRETCHINGL") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    mStretchingParams.push_back(double(std::atoi(lineVector[1].c_str())));
                }
            }
            else if (line.find("MBENDINGK") != std::string::npos
                     || line.find("MBENDINGTHETA") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    mBendingParams.push_back(double(std::atoi(lineVector[1].c_str())));
                }
            }
            else if (line.find("MTWISTINGK") != std::string::npos
                     || line.find("MTWISTINGPHI") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    mTwistingParams.push_back(double(std::atoi(lineVector[1].c_str())));
                }
            }
            else if (line.find("VOLUMEK") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    volumeParams.push_back(double(std::atoi(lineVector[1].c_str())));
                }
            }
            else {}
        }
        return {fStretchingParams, fBendingParams, fTwistingParams,
                lStretchingParams, lBendingParams, lTwistingParams,
                mStretchingParams, mBendingParams, mTwistingParams, volumeParams};
    }

    ///Mechanics objects to be included
    std::vector<std::string> readMechanicsObjects() {
        
        std::string filament, linker, motor, volume;
        
        std::string line;
        while(getline(_inputFile, line)) {
            
            if (line.find("FILAMENTS") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() != 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                filament = lineVector[1];
            }
            else if (line.find("LINKERS") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() != 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                linker = lineVector[1];
            }
            else if (line.find("MOTORS") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() != 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                motor = lineVector[1];
            }
            else if (line.find("VOLUME") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() != 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                volume = lineVector[1];
            }
            else {}
        }
        return {filament, linker, motor, volume};
    }
    
    ///Mechanics force field types
    std::vector<std::string> readMechanicsFFType() {
        
        ///Filament type
        std::string fStretchingType;
        std::string fBendingType;
        std::string fTwistingType;
        
        ///Linker type
        std::string lStretchingType;
        std::string lBendingType;
        std::string lTwistingType;
        
        ///Motor type
        std::string mStretchingType;
        std::string mBendingType;
        std::string mTwistingType;
        
        ///Volume type
        std::string volumeType;

        std::string line;
        while(getline(_inputFile, line)) {
            
            if (line.find("FSTRETCHINGTYPE") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    fStretchingType = lineVector[1];
                }
            }
            else if (line.find("FBENDINGTYPE") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    fBendingType = lineVector[1];
                }
            }
            else if (line.find("FTWISTINGTYPE") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    fTwistingType = lineVector[1];
                }
            }
            else if (line.find("LSTRETCHINGTYPE") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    lStretchingType = lineVector[1];
                }
            }
            else if (line.find("LBENDINGTYPE") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    lBendingType = lineVector[1];
                }
            }
            else if (line.find("LTWISTINGTYPE") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    lTwistingType = lineVector[1];
                }
            }
            else if (line.find("MSTRETCHINGTYPE") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    mStretchingType = lineVector[1];
                }
            }
            else if (line.find("MBENDINGTYPE") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    mBendingType = lineVector[1];
                }
            }
            else if (line.find("MTWISTINGTYPE") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    mTwistingType = lineVector[1];
                }
            }
            else if (line.find("VOLUMETYPE") != std::string::npos) {
                
                std::vector<std::string> lineVector = split<std::string>(line);
                if(lineVector.size() >= 2) {
                    std::cout << "There was an error parsing input file. Exiting" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    volumeType = lineVector[1];
                }
            }
            else {}
        }
        return {fStretchingType, fBendingType, fTwistingType,
            lStretchingType, lBendingType, lTwistingType,
            mStretchingType, mBendingType, mTwistingType, volumeType};
        
    }

    
};


#endif /* defined(__Cyto__Parser__) */
