//
//  Parser.cpp
//  Cyto
//
//  Created by James Komianos on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "SystemParameters.h"
#include "Parser.h"
//
//bool SystemParser::mechanics() {
//    
//    _inputFile.clear();
//    _inputFile.seekg(0);
//    
//    std::string line;
//    while(getline(_inputFile, line)) {
//        
//        if (line.find("MECHANICS:") != std::string::npos) {
//            
//            std::vector<std::string> lineVector = split<std::string>(line);
//            if(lineVector.size() != 2) {
//                std::cout << "Need to specify Mechanics ON/OFF. Exiting" << std::endl;
//                exit(EXIT_FAILURE);
//            }
//            else if (lineVector.size() == 2) {
//                if(lineVector[1] == "ON") return true;
//            }
//        }
//    }
//    ///default is false
//    return false;
//}
//
//bool SystemParser::chemistry() {
//    
//    _inputFile.clear();
//    _inputFile.seekg(0);
//    
//    std::string line;
//    while(getline(_inputFile, line)) {
//        
//        if (line.find("CHEMISTRY:") != std::string::npos) {
//            
//            std::vector<std::string> lineVector = split<std::string>(line);
//            if(lineVector.size() != 2) {
//                std::cout << "Need to specify Chemistry ON/OFF. Exiting" << std::endl;
//                exit(EXIT_FAILURE);
//            }
//            else if (lineVector.size() == 2) {
//                if(lineVector[1] == "ON") return true;
//            }
//        }
//    }
//    ///default is false
//    return false;
//}

///CHEMISTRY PARSER
ChemistryAlgorithm SystemParser::readChemistryAlgorithm() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    ChemistryAlgorithm CAlgorithm;
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != std::string::npos) { continue; }
        
        
        if (line.find("CALGORITHM") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "There was an error parsing input file at Chemistry parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.algorithm = lineVector[1];
            }
        }
        if (line.find("NUMSTEPS:") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "There was an error parsing input file at Chemistry parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.numSteps = std::atoi(lineVector[1].c_str());
            }
        }
        if (line.find("NUMSTEPSPERM:") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "There was an error parsing input file at Chemistry parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.numStepsPerMech = std::atoi(lineVector[1].c_str());
            }
        }
        
        
    }
    return CAlgorithm;
}

///CHEMISTRY SETUP PARSER
ChemistrySetup SystemParser::readChemistrySetup() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    ChemistrySetup CSetup;
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != std::string::npos) { continue; }
        
        if(line.find("CHEMISTRYFILE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "Error reading chemistry input file. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                CSetup.inputFile = lineVector[1];
            else {}
        }
    }
    return CSetup;
}


///Mechanics force field types
MechanicsFFType SystemParser::readMechanicsFFType() {
    
    MechanicsFFType MTypes;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != std::string::npos) { continue; }
        
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
        
        if(line.find("#") != std::string::npos) { continue; }
        
        ///Filament stretching
        if (line.find("FSTRETCHINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FStretchingK = std::atof(lineVector[1].c_str());
            }
        }
//        if (line.find("FSTRETCHINGL") != std::string::npos) {
//            
//            std::vector<std::string> lineVector = split<std::string>(line);
//            if(lineVector.size() > 2) {
//                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
//                exit(EXIT_FAILURE);
//            }
//            else if (lineVector.size() == 2) {
//                MParams.FStretchingL = std::atof(lineVector[1].c_str());
//            }
//        }
        
        ///Filament bending
        else if (line.find("FBENDINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FBendingK = std::atof((lineVector[1].c_str()));
            }
        }
        else if (line.find("FBENDINGTHETA") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FBendingTheta = std::atof((lineVector[1].c_str()));
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
                MParams.FTwistingK = std::atof((lineVector[1].c_str()));
            }
        }
        else if (line.find("FTWISTINGPHI") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FTwistingPhi = std::atof((lineVector[1].c_str()));
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
                MParams.LStretchingK = std::atof((lineVector[1].c_str()));
            }
        }
        if (line.find("LSTRETCHINGL") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LStretchingL = std::atof((lineVector[1].c_str()));
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
                MParams.LBendingK = std::atof((lineVector[1].c_str()));
            }
        }
        else if (line.find("LBENDINGTHETA") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LBendingTheta = std::atof((lineVector[1].c_str()));
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
                MParams.LTwistingK = std::atof((lineVector[1].c_str()));
            }
        }
        else if (line.find("LTWISTINGPHI") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LTwistingPhi = std::atof((lineVector[1].c_str()));
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
                MParams.MStretchingK = std::atof((lineVector[1].c_str()));
            }
        }
        if (line.find("MSTRETCHINGL") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MStretchingL = std::atof((lineVector[1].c_str()));
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
                MParams.MBendingK = std::atof((lineVector[1].c_str()));
            }
        }
        else if (line.find("MBENDINGTHETA") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MBendingTheta = std::atof((lineVector[1].c_str()));
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
                MParams.MTwistingK = std::atof((lineVector[1].c_str()));
            }
        }
        else if (line.find("MTWISTINGPHI") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MTwistingPhi = std::atof((lineVector[1].c_str()));
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
                MParams.VolumeK = std::atof((lineVector[1].c_str()));
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
        
        if(line.find("#") != std::string::npos) { continue; }
        
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
        
        if(line.find("#") != std::string::npos) { continue; }
        
        if (line.find("BCUTOFF") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "There was an error parsing input file at Boundary parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.boundaryCutoff = std::atof((lineVector[1].c_str()));
            }
            ///Default value to be half compartment size
            else {
                BParams.boundaryCutoff = SystemParameters::Geometry().compartmentSizeX / 2;
            }
        }
        else if (line.find("BINTERACTIONK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "There was an error parsing input file at Boundary parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.boundaryK = std::atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("BSCREENLENGTH") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "There was an error parsing input file at Boundary parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.screenLength = std::atof((lineVector[1].c_str()));
            }
            else {}
        }
        else {}
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
        
        if(line.find("#") != std::string::npos) { continue; }
        
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
        
        if(line.find("#") != std::string::npos) { continue; }
        
        if (line.find("NX") != std::string::npos
            || line.find("NY") != std::string::npos
            || line.find("NZ") != std::string::npos) {
            
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at grid dimensions. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2)
                gridTemp.push_back(std::atof((lineVector[1].c_str())));
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
                compartmentTemp.push_back(std::atof((lineVector[1].c_str())));
            else {}
        }
        
        else if(line.find("MONOMERSIZE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "A monomer size needs to be specified. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                monomerSize = std::atof(lineVector[1].c_str());
            else {}
        }
        
        else if(line.find("CYLINDERSIZE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() != 2) {
                std::cout << "A cylinder size needs to be specified. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                cylinderSize = std::atof(lineVector[1].c_str());
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
        
        if(line.find("#") != std::string::npos) { continue; }
        
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
    
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != std::string::npos) { continue; }
        
        std::vector<std::string> lineVector = split<std::string>(line);
        if(lineVector.size() == 7) {
            std::vector<double> coord1;
            std::vector<double> coord2;
            for(auto it = lineVector.begin() + 1; it != lineVector.begin() + 4; it++) {
                coord1.push_back(std::atof(((*it).c_str())));
            }
            for(auto it = lineVector.begin() + 4; it != lineVector.end(); it++) {
                coord2.push_back(std::atof(((*it).c_str())));
            }
            
            returnVector.push_back({coord1, coord2});
        }
    }
    return returnVector;
}

///CHEMISTRY INPUT PARSER
ChemistrySpeciesAndReactions ChemistryParser::readChemistryInput() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    ChemistrySpeciesAndReactions chemSR;
    std::string line;

    while(getline(_inputFile, line)) {
        
        if(line.find("#") != std::string::npos) { continue; }
        
        if(line.find("SPECIESBULK") != std::string::npos) {
        
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() !=  3) {
                std::cout << "Error reading a bulk species. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3)
                chemSR.speciesBulk.push_back(std::tuple<std::string, int>(lineVector[1], std::atoi(lineVector[2].c_str())));
            else {}
        }
        if(line.find("SPECIESDIFFUSING") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() !=  4) {
                std::cout << "Error reading a diffusing species. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 4)
                chemSR.speciesDiffusing.push_back(std::tuple<std::string, int, double>
                        (lineVector[1], std::atoi(lineVector[2].c_str()), std::atof(lineVector[3].c_str())));
            else {}
        }
        
        if(line.find("SPECIESFILAMENT") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() !=  2) {
                std::cout << "Error reading a filament species. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chemSR.speciesFilament.push_back(lineVector[1]);
            else {}
        }
        if(line.find("SPECIESBOUND") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() !=  2) {
                std::cout << "Error reading a filament bound species. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chemSR.speciesBound.push_back(lineVector[1]);
            else {}
        }
        if(line.find("SPECIESPLUSEND") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() !=  2) {
                std::cout << "Error reading a filament plus end species. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chemSR.speciesPlusEnd.push_back(lineVector[1]);
            else {}
        }
        if(line.find("SPECIESMINUSEND") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() !=  2) {
                std::cout << "Error reading a filament minus end species. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chemSR.speciesMinusEnd.push_back(lineVector[1]);
            else {}
        }
        
        ///loop through a reaction
        if(line.find("REACTION") != std::string::npos) {
            
            std::vector<std::string> reactants;
            std::vector<std::string> products;
            
            std::vector<std::string> lineVector = split<std::string>(line);
            
            auto arrowIt = std::find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
            
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chemSR.reactions.push_back(std::tuple<std::vector<std::string>, std::vector<std::string>, double>
                                        (reactants, products, std::atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                std::cout << "Error reading reaction. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
    }
    
    return chemSR;
}



