//
//  Parser.cpp
//  Cyto
//
//  Created by James Komianos on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "Parser.h"

bool mechanics() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("MECHANICS") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Chemistry parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                if(lineVector[1] == "ON") return true;
            }
        }
    ///default is false
    return false;
}

bool chemistry() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("CHEMISTRY") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Chemistry parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                if(lineVector[1] == "ON") return true;
            }
        }
    ///default is false
    return false;
}

///CHEMISTRY PARSER
ChemistryParameters Parser::readChemistryParameters() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    ChemistryParameters CParams;
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("SETUP") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Chemistry parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.algorithm = lineVector[1];
            }
        }
        
        if (line.find("SETUP") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Chemistry parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.setup = lineVector[1];
            }
        }
    }
    return CParams;
}

BoundaryParameters Parser::readBoundaryParameters() {
    
    BoundaryParameters BParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("BOUNDARYTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Boundary parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.boundaryType = lineVector[1];
            }
        }
    }
    return BParams;
}

///GEOMETRY CONSTANT PARSER
GeometryParameters Parser::readGeometryParameters() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    GeometryParameters GParams;
    
    std::vector<double> gridTemp;
    std::vector<double> compartmentTemp;
    double monomerSize = 0;
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
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at monomer size. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                monomerSize = double(std::atoi(lineVector[1].c_str()));
            else {}
        }
        
        else if(line.find("NDIM") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() > 2) {
                std::cout << "There was an error parsing input file at compartment size. Exiting" << std::endl;
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
    GParams.monomerSize = monomerSize;
    if(gridTemp.size() >= 1) GParams.NX = gridTemp[0];
    if(gridTemp.size() >= 2) GParams.NY = gridTemp[1];
    if(gridTemp.size() >= 3) GParams.NZ = gridTemp[2];
    if(compartmentTemp.size() >= 1) GParams.compartmentSizeX = compartmentTemp[0];
    if(compartmentTemp.size() >= 2) GParams.compartmentSizeY = compartmentTemp[1];
    if(compartmentTemp.size() >= 3) GParams.compartmentSizeZ = compartmentTemp[2];
    
    return GParams;
    
}


///Mechanics force field types
MechanicsFFType Parser::readMechanicsFFType() {
    
    MechanicsFFType MTypes;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        if (line.find("FSTRETCHINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Filament stretching type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FStretchingType = lineVector[1];
            }
        }
        else if (line.find("FBENDINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Filament bending type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FBendingType = lineVector[1];
            }
        }
        else if (line.find("FTWISTINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Filament twisting type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FTwistingType = lineVector[1];
            }
        }
        else if (line.find("LSTRETCHINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Linker stretching type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LStretchingType = lineVector[1];
            }
        }
        else if (line.find("LBENDINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Linker bending type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LBendingType = lineVector[1];
            }
        }
        else if (line.find("LTWISTINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Linker twisting type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LTwistingType = lineVector[1];
            }
        }
        else if (line.find("MSTRETCHINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Motor stretching type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MStretchingType = lineVector[1];
            }
        }
        else if (line.find("MBENDINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Motor bending type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MBendingType = lineVector[1];
            }
        }
        else if (line.find("MTWISTINGTYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Motor twisting type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MTwistingType = lineVector[1];
            }
        }
        else if (line.find("VOLUMETYPE") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Volume type. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.VolumeType = lineVector[1];
            }
        }
        else {}
    }
    return MTypes;
}


///MECHANICS CONSTANT PARSER
MechanicsParameters Parser::readMechanicsParameters() {
    
    MechanicsParameters MParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    std::string line;
    while(getline(_inputFile, line)) {
        
        ///Filament stretching
        if (line.find("FSTRETCHINGK") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FStretchingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        if (line.find("FSTRETCHINGL") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
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
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FBendingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("FBENDINGTHETA") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
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
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Filament parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FTwistingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("FTWISTINGPHI") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
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
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LStretchingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        if (line.find("LSTRETCHINGL") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
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
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LBendingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("LBENDINGTHETA") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
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
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Linker parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.LTwistingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("LTWISTINGPHI") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
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
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MStretchingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        if (line.find("MSTRETCHINGL") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
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
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MBendingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("FBENDINGTHETA") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
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
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Motor parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MTwistingK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else if (line.find("MTWISTINGPHI") != std::string::npos) {
            
            std::vector<std::string> lineVector = split<std::string>(line);
            if(lineVector.size() >= 2) {
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
            if(lineVector.size() >= 2) {
                std::cout << "There was an error parsing input file at Volume parameters. Exiting" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.VolumeK = double(std::atoi(lineVector[1].c_str()));
            }
        }
        else {}
    }
    return MParams;
}


///Function to check consistency between all desired forcefields, boundaries, and constants
bool Parser::checkInput(ChemistryParameters &CParams, BoundaryParameters &BParams, GeometryParameters &GParams,
                MechanicsFFType &MTypes, MechanicsParameters &MParams) {
    
    std::cout << "Checking parsed data..." << std::endl;
    
    ///Check chemistry input
    if(CParams.algorithm == "" || CParams.setup == "") {
        std::cout << "Did not specify all chemistry parameters." << std::endl;
        return false;
    }
    
      ///Check mechanics input consistency
    //Filaments
    if(MTypes.FStretchingType != "") {
        if(MParams.FStretchingK == 0 || MParams.FStretchingL == 0) {
            std::cout << "Did not specify non-zero FStretching constants." << std::endl;
            return false;
        }
    }
    else {
        if(MParams.FStretchingK != 0 || MParams.FStretchingL != 0) {
            std::cout << "Specified FStretching constants but did not choose FF type." << std::endl;
            return false;
        }
    }
    
    if(MTypes.FBendingType != "") {
        if(MParams.FBendingK == 0 || MParams.FBendingTheta == 0) {
            std::cout << "Did not specify non-zero FBending constants." << std::endl;
            return false;
        }
    }
    else {
        if(MParams.FBendingK != 0 || MParams.FBendingTheta != 0) {
            std::cout << "Specified FBending constants but did not choose FF type." << std::endl;
            return false;
        }
    }
    
    if(MTypes.FTwistingType != "") {
        if(MParams.FTwistingK == 0 || MParams.FTwistingPhi == 0) {
            std::cout << "Did not specify non-zero FTwisting constants." << std::endl;
            return false;
        }
    }
    else {
        if(MParams.FTwistingK != 0 || MParams.FTwistingPhi != 0) {
            std::cout << "Specified FTwisting constants but did not choose FF type." << std::endl;
            return false;
        }
    }
    
    //Linkers
    if(MTypes.LStretchingType != "") {
        if(MParams.LStretchingK == 0 || MParams.LStretchingL == 0) {
            std::cout << "Did not specify non-zero LStretching constants." << std::endl;
            return false;
        }
    }
    else {
        if(MParams.LStretchingK != 0 || MParams.LStretchingL != 0) {
            std::cout << "Specified LStretching constants but did not choose FF type." << std::endl;
            return false;
        }
    }
    
    if(MTypes.LBendingType != "") {
        if(MParams.LBendingK == 0 || MParams.LBendingTheta == 0) {
            std::cout << "Did not specify non-zero LBending constants." << std::endl;
            return false;
        }
    }
    else {
        if(MParams.LBendingK != 0 || MParams.LBendingTheta != 0) {
            std::cout << "Specified FBending constants but did not choose FF type." << std::endl;
            return false;
        }
    }
    
    if(MTypes.LTwistingType != "") {
        if(MParams.LTwistingK == 0 || MParams.LTwistingPhi == 0) {
            std::cout << "Did not specify non-zero LTwisting constants." << std::endl;
            return false;
        }
    }
    else {
        if(MParams.LTwistingK != 0 || MParams.LTwistingPhi != 0) {
            std::cout << "Specified LTwisting constants but did not choose FF type." << std::endl;
            return false;
        }
    }
    
    //Motors
    if(MTypes.MStretchingType != "") {
        if(MParams.MStretchingK == 0 || MParams.MStretchingL == 0) {
            std::cout << "Did not specify non-zero MStretching constants." << std::endl;
            return false;
        }
    }
    else {
        if(MParams.MStretchingK != 0 || MParams.MStretchingL != 0) {
            std::cout << "Specified MStretching constants but did not choose FF type." << std::endl;
            return false;
        }
    }
    
    if(MTypes.MBendingType != "") {
        if(MParams.MBendingK == 0 || MParams.MBendingTheta == 0) {
            std::cout << "Did not specify non-zero MBending constants." << std::endl;
            return false;
        }
    }
    else {
        if(MParams.MBendingK != 0 || MParams.MBendingTheta != 0) {
            std::cout << "Specified MBending constants but did not choose FF type." << std::endl;
            return false;
        }
    }
    
    if(MTypes.MTwistingType != "") {
        if(MParams.MTwistingK == 0 || MParams.MTwistingPhi == 0) {
            std::cout << "Did not specify non-zero MTwisting constants." << std::endl;
            return false;
        }
    }
    else {
        if(MParams.MTwistingK != 0 || MParams.MTwistingPhi != 0) {
            std::cout << "Specified MTwisting constants but did not choose FF type." << std::endl;
            return false;
        }
    }
    
    ///Volume
    if(MTypes.VolumeType != "") {
        if(MParams.VolumeK == 0) {
            std::cout << "Did not specify non-zero volume constant." << std::endl;
            return false;
        }
    }
    else {
        if(MParams.VolumeK != 0 ) {
            std::cout << "Specified MStretching constants but did not choose FF type." << std::endl;
            return false;
        }
    }
    
    
    ///Check Geometry consistency
    if(GParams.nDim > 3 || GParams.nDim <= 0)
        std::cout << "Specified nDim must be between 1 and 3!" <<std::endl;
    
    if(GParams.nDim == 1) {
        if((GParams.NY != 0 || GParams.NZ != 0) || (GParams.NX == 0)) {
            std::cout << "Specified grid dimensions do not match NDIM." << std::endl;
            return false;
        }
        if((GParams.compartmentSizeY != 0 || GParams.compartmentSizeZ != 0) || (GParams.compartmentSizeX == 0)) {
            std::cout << "Specified compartment dimensions do not match NDIM." << std::endl;
            return false;
        }
        
        
    }
    else if(GParams.nDim == 2) {
        if(GParams.NZ != 0 || (GParams.NX == 0 || GParams.NY == 0)) {
            std::cout << "Specified grid dimensions do not match NDIM." << std::endl;
            return false;
        }
        if(GParams.compartmentSizeZ != 0 || (GParams.compartmentSizeX == 0 || GParams.compartmentSizeY == 0)) {
            std::cout << "Specified compartment dimensions do not match NDIM." << std::endl;
            return false;
        }
        
    }
    else {
        if(GParams.NX == 0 || GParams.NY == 0 || GParams.NZ == 0) {
            std::cout << "Specified grid dimensions do not match NDIM." << std::endl;
            return false;
        }
        if(GParams.compartmentSizeX == 0 || GParams.compartmentSizeY == 0 || GParams.compartmentSizeZ == 0) {
            std::cout << "Specified compartment dimensions do not match NDIM." << std::endl;
            return false;
        }
    }
    
    if(GParams.monomerSize == 0) {
        std::cout << "Must specify a non-zero monomer size" << std::endl;
        return false;
    }
    
    
    std::cout << "Parsing checks all passed!" <<std::endl;
    return true;
}
