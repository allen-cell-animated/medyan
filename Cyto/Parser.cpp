//
//  Parser.cpp
//  Cyto
//
//  Created by James Komianos on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "Parser.h"

bool SystemParser::mechanics() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("MECHANICS:") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "Need to specify Mechanics ON/OFF. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                if(lineVector[1] == "ON") return true;
            }
        }
    }
    ///default is false
    return false;
}

bool SystemParser::chemistry() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("CHEMISTRY:") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "Need to specify Chemistry ON/OFF. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                if(lineVector[1] == "ON") return true;
            }
        }
    }
    ///default is false
    return false;
}

///CHEMISTRY PARSER
ChemistryAlgorithm SystemParser::readChemistryAlgorithm() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    ChemistryAlgorithm CAlgorithm;
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("CALGORITHM") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Chemistry parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.algorithm = lineVector[1];
            }
        }
        
        if (line.find("SETUP") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Chemistry parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.setup = lineVector[1];
            }
        }
    }
    return CAlgorithm;
}

///Mechanics force field types
MechanicsFFType SystemParser::readMechanicsFFType() {
    
    MechanicsFFType MTypes;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("FSTRETCHINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament stretching type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FStretchingType = lineVector[1];
            }
        }
        else if (line.find("FBENDINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament bending type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FBendingType = lineVector[1];
            }
        }
        else if (line.find("FTWISTINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament twisting type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FTwistingType = lineVector[1];
            }
        }
        else if (line.find("LSTRETCHINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker stretching type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LStretchingType = lineVector[1];
            }
        }
        else if (line.find("LBENDINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker bending type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LBendingType = lineVector[1];
            }
        }
        else if (line.find("LTWISTINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker twisting type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LTwistingType = lineVector[1];
            }
        }
        else if (line.find("MSTRETCHINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor stretching type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MStretchingType = lineVector[1];
            }
        }
        else if (line.find("MBENDINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor bending type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MBendingType = lineVector[1];
            }
        }
        else if (line.find("MTWISTINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor twisting type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MTwistingType = lineVector[1];
            }
        }
        else if (line.find("BOUNDARYTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Boundary type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.BoundaryFFType = lineVector[1];
            }
        }
        else if (line.find("VOLUMETYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Volume type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.VolumeFFType = lineVector[1];
            }
        }
        
        else {}
    }
    return MTypes;
}


///MECHANICS CONSTANT PARSER
void SystemParser::readMechanicsParameters() {
    
    MechanicsParameters MParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        ///Filament stretching
        if (line.find("FSTRETCHINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FStretchingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        if (line.find("FSTRETCHINGL") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FStretchingL = double(std::atoi(lineVector[1].c_str()));
            }
        }
        
        ///Filament bending
        else if (line.find("FBENDINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FBendingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("FBENDINGTHETA") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FBendingTheta = double(std::atoi(lineVector[1].c_str()));
            }
        }
        
        ///Filament twisting
        else if (line.find("FTWISTINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FTwistingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("FTWISTINGPHI") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FTwistingPhi = double(std::atoi(lineVector[1].c_str()));
            }
        }
        
        ///Linker stretching
        if (line.find("LSTRETCHINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LStretchingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        if (line.find("LSTRETCHINGL") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LStretchingL = double(std::atoi(lineVector[1].c_str()));
            }
        }
        
        ///Linker bending
        else if (line.find("LBENDINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LBendingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("LBENDINGTHETA") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LBendingTheta = double(std::atoi(lineVector[1].c_str()));
            }
        }
        
        ///Linker twisting
        else if (line.find("LTWISTINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LTwistingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("LTWISTINGPHI") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LTwistingPhi = double(std::atoi(lineVector[1].c_str()));
            }
        }
        
        ///Motor stretching
        if (line.find("MSTRETCHINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MStretchingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        if (line.find("MSTRETCHINGL") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MStretchingL = double(std::atoi(lineVector[1].c_str()));
            }
        }
        
        ///Motor bending
        else if (line.find("MBENDINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MBendingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("FBENDINGTHETA") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MBendingTheta = double(std::atoi(lineVector[1].c_str()));
            }
        }
        
        ///Motor twisting
        else if (line.find("MTWISTINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MTwistingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("MTWISTINGPHI") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MTwistingPhi = double(std::atoi(lineVector[1].c_str()));
            }
        }
        
        ///Volume parameter
        else if (line.find("VOLUMEK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Volume parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.VolumeK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else {}
    }
    ///Set system parameters
    SystemParameters::MParams = MParams;
}

MechanicsAlgorithm SystemParser::readMechanicsAlgorithm() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    MechanicsAlgorithm MAlgorithm;
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("CONJUGATEGRADIENT") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "A conjugate gradient method must be specified. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MAlgorithm.ConjugateGradient = lineVector[1];
            }
        }
        else if (line.find("MD") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "A Mechanics algorithm must be specified. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MAlgorithm.MD = lineVector[1];
            }
        }
    }
    return MAlgorithm;
}

    
///BOUNDARY PARSERS
void SystemParser::readBoundaryParameters() {
    
    BoundaryParameters BParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("BINTERACTIONCUTOFF") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "There was an error parsing input file at Boundary parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.interactionCutoff = double(std::atoi(lineVector[1].c_str()));
            }
        }
    }
    ///Set system parameters
    SystemParameters::BParams = BParams;
}

BoundaryType SystemParser::readBoundaryType() {
        
    _inputFile.clear();
    _inputFile.seekg(0);
    
    BoundaryType BType;
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("BOUNDARYSHAPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "A boundary shape needs to be specified. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BType.boundaryShape = lineVector[1];
            }
        }
    }
    return BType;
}
    
    
///GEOMETRY PARSERS
void SystemParser::readGeometryParameters() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    GeometryParameters GParams;
    
    std::vector<double> gridTemp;
    std::vector<double> compartmentTemp;
    double monomerSize = 0;
    double cylinderSize = 0;
    short nDim = 0;
    
    //find grid size lines
    std::string line;
    while(getline(_inputFile, line)) {
        if (line.find("NX") != std::string::npos
            || line.find("NY") != std::string::npos
            || line.find("NZ") != std::string::npos) {
            
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at grid dimensions. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2)
                gridTemp.push_back(double(std::atoi(lineVector[1].c_str())));
            else {}
        }
        
        else if (line.find("COMPARTMENTSIZEX") != std::string::npos
                 || line.find("COMPARTMENTSIZEY") != std::string::npos
                 || line.find("COMPARTMENTSIZEZ") != std::string::npos) {
            
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at compartment size. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                compartmentTemp.push_back(double(std::atoi(lineVector[1].c_str())));
            else {}
        }
        
        else if(line.find("MONOMERSIZE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "A monomer size needs to be specified. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                monomerSize = double(std::atoi(lineVector[1].c_str()));
            else {}
        }
        
        else if(line.find("CYLINDERSIZE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "A cylinder size needs to be specified. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                cylinderSize = double(std::atoi(lineVector[1].c_str()));
            else {}
        }
        
        else if(line.find("NDIM") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() !=  2) {
                std::cout << "Number of dimensions needs to be specified. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                nDim = short(std::atoi(lineVector[1].c_str()));
            }
            else {}
        }
        else {}
    }
    ///set geometry parameters and return
    GParams.nDim = nDim;
    GParams.cylinderSize = cylinderSize;
    GParams.monomerSize = monomerSize;
    if(gridTemp.size() >= 1) GParams.NX = gridTemp[0];
    if(gridTemp.size() >= 2) GParams.NY = gridTemp[1];
    if(gridTemp.size() >= 3) GParams.NZ = gridTemp[2];
    if(compartmentTemp.size() >= 1) GParams.compartmentSizeX = compartmentTemp[0];
    if(compartmentTemp.size() >= 2) GParams.compartmentSizeY = compartmentTemp[1];
    if(compartmentTemp.size() >= 3) GParams.compartmentSizeZ = compartmentTemp[2];
    
    SystemParameters::GParams = GParams;
}

///FILAMENT SETUP PARSER
FilamentSetup SystemParser::readFilamentSetup() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    FilamentSetup FSetup;
    
    std::string line;
    while(getline(_inputFile, line)) {
        if(line.find("FILAMENTFILE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "Error reading filament input file. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.inputFile = lineVector[1];
            else {}
        }
        else if(line.find("NUMFILAMENTS") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "Error reading number of filaments. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.numFilaments = std::atoi(lineVector[1].c_str());
            else {}
        }
    }
    return FSetup;
}

///FILAMENT DATA PARSER
std::vector<std::vector<std::vector<double>>> FilamentParser::readFilaments() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::vector<std::vector<std::vector<double>>> returnVector;
    
    std::string line;
    int filamentIndex = 0;
    
    while(getline(_inputFile, line)) {
        std::vector<double> lineVector = split<double>(line);
        if(lineVector.size() == 7) {
            
            std::vector<double> coord1;
            std::vector<double> coord2;
            for(auto it = lineVector.begin() + 1; it != lineVector.begin() + 4; it++)
                coord1.push_back(*it);
            for(auto it = lineVector.begin() + 5; it != lineVector.end(); it++)
                coord2.push_back(*it);
            
            returnVector[filamentIndex] = {coord1, coord2};
        }
        filamentIndex++;
    }
    return returnVector;
}



