
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

#include "Parser.h"

#include "SystemParameters.h"

void SystemParser::readChemistryParameters() {
    
    ChemistryParameters CParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("NUMBULKSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numBulkSpecies = atof(lineVector[1].c_str());
            }
        }
        
        if (line.find("NUMDIFFUSINGSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numDiffusingSpecies = atof(lineVector[1].c_str());
            }
        }
    
        if (line.find("NUMFILAMENTSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numFilamentSpecies = atof(lineVector[1].c_str());
            }
        }
    
        if (line.find("NUMPLUSENDSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numPlusEndSpecies = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMMINUSENDSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numMinusEndSpecies = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMBOUNDSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numBoundSpecies = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMLINKERSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numLinkerSpecies = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMMOTORSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numMotorSpecies = atof(lineVector[1].c_str());
            }
        }
    }
    //set system parameters
    SystemParameters::CParams = CParams;
}

ChemistryAlgorithm SystemParser::readChemistryAlgorithm() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    ChemistryAlgorithm CAlgorithm;
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        
        if (line.find("CALGORITHM") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.algorithm = lineVector[1];
            }
        }
        if (line.find("NUMSTEPS:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.numSteps = atoi(lineVector[1].c_str());
            }
        }
        if (line.find("NUMSTEPSPERM:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "There was an error parsing input file at Chemistry parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.numStepsPerMech = atoi(lineVector[1].c_str());
            }
        }
        
        
    }
    return CAlgorithm;
}

ChemistrySetup SystemParser::readChemistrySetup() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    ChemistrySetup CSetup;
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if(line.find("CHEMISTRYFILE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading chemistry input file. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                CSetup.inputFile = lineVector[1];
            else {}
        }
    }

    return CSetup;
}

MechanicsFFType SystemParser::readMechanicsFFType() {
    
    MechanicsFFType MTypes;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("FSTRETCHINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Filament stretching type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FStretchingType = lineVector[1];
            }
        }
        else if (line.find("FBENDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Filament bending type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FBendingType = lineVector[1];
            }
        }
        else if (line.find("FTWISTINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Filament twisting type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FTwistingType = lineVector[1];
            }
        }
        else if (line.find("LSTRETCHINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Linker stretching type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LStretchingType = lineVector[1];
            }
        }
        else if (line.find("LBENDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Linker bending type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LBendingType = lineVector[1];
            }
        }
        else if (line.find("LTWISTINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Linker twisting type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LTwistingType = lineVector[1];
            }
        }
        else if (line.find("MSTRETCHINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Motor stretching type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MStretchingType = lineVector[1];
            }
        }
        else if (line.find("MBENDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Motor bending type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MBendingType = lineVector[1];
            }
        }
        else if (line.find("MTWISTINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Motor twisting type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MTwistingType = lineVector[1];
            }
        }
        if (line.find("BRSTRETCHINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Filament stretching type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.BrStretchingType = lineVector[1];
            }
        }
        else if (line.find("BRBENDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Filament bending type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.BrBendingType = lineVector[1];
            }
        }
        else if (line.find("BRTWISTINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Filament twisting type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.BrTwistingType = lineVector[1];
            }
        }
        else if (line.find("BOUNDARYTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Boundary type. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.BoundaryFFType = lineVector[1];
            }
        }
        else if (line.find("VOLUMETYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Volume type. Exiting" << endl;
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

void SystemParser::readMechanicsParameters() {
    
    MechanicsParameters MParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        //Filament stretching
        if (line.find("FSTRETCHINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Mechanics parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FStretchingK = atof(lineVector[1].c_str());
            }
        }

        //Filament bending
        else if (line.find("FBENDINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Mechanics parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FBendingK = atof((lineVector[1].c_str()));
            }
        }
        else if (line.find("FBENDINGTHETA") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Mechanics parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FBendingTheta = atof((lineVector[1].c_str()));
            }
        }
        
        //Filament twisting
        else if (line.find("FTWISTINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Mechanics parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FTwistingK = atof((lineVector[1].c_str()));
            }
        }
        else if (line.find("FTWISTINGPHI") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Mechanics parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FTwistingPhi = atof((lineVector[1].c_str()));
            }
        }
        
        //Linker stretching
        if (line.find("LSTRETCHINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.LStretchingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        
        //Linker bending
        else if (line.find("LBENDINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.LBendingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("LBENDINGTHETA") != string::npos) {
            
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.LBendingTheta.push_back(atof((lineVector[i].c_str())));
            }
            
        }
        
        //Linker twisting
        else if (line.find("LTWISTINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.LTwistingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("LTWISTINGPHI") != string::npos) {
            
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.LTwistingPhi.push_back(atof((lineVector[i].c_str())));
            }
        }
        
        //Motor stretching
        if (line.find("MSTRETCHINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.MStretchingK.push_back(atof((lineVector[i].c_str())));
            }
        }

        //Motor bending
        else if (line.find("MBENDINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.MBendingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("MBENDINGTHETA") != string::npos) {
            
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.MBendingTheta.push_back(atof((lineVector[i].c_str())));
            }
        }
        
        //Motor twisting
        else if (line.find("MTWISTINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.LTwistingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("MTWISTINGPHI") != string::npos) {
            
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.MTwistingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        
        //Branch stretching
        if (line.find("BRSTRETCHINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BrStretchingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("BRSTRETCHINGL") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BrStretchingL.push_back(atof((lineVector[i].c_str())));
            }
        }
        
        
        //Branch bending
        else if (line.find("BRBENDINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BrBendingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("BRBENDINGTHETA") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BrBendingTheta.push_back(atof((lineVector[i].c_str())));
            }
            
        }
        
        //Branch twisting
        else if (line.find("BRTWISTINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BrTwistingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("BRTWISTINGPHI") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BrTwistingPhi.push_back(atof((lineVector[i].c_str())));
            }
        }
        
        //Volume parameter
        else if (line.find("VOLUMEK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Mechanics parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.VolumeK = atof((lineVector[1].c_str()));
            }
        }
        
        else if (line.find("VOLUMECUTOFF") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at Mechanics parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.VolumeCutoff = atof((lineVector[1].c_str()));
            }
        }
        else {}
    }
    //Set system parameters
    SystemParameters::MParams = MParams;
}

MechanicsAlgorithm SystemParser::readMechanicsAlgorithm() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    MechanicsAlgorithm MAlgorithm;
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("CONJUGATEGRADIENT") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "A conjugate gradient method must be specified. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MAlgorithm.ConjugateGradient = lineVector[1];
            }
        }
        else if (line.find("MD") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "A Mechanics algorithm must be specified. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MAlgorithm.MD = lineVector[1];
            }
        }
    }
    return MAlgorithm;
}

void SystemParser::readBoundaryParameters() {
    
    BoundaryParameters BParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("BCUTOFF") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "There was an error parsing input file at Boundary parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.boundaryCutoff = atof((lineVector[1].c_str()));
            }
            //Default value to be half compartment size
            else {
                BParams.boundaryCutoff = SystemParameters::Geometry().compartmentSizeX / 2;
            }
        }
        else if (line.find("BINTERACTIONK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "There was an error parsing input file at Boundary parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.boundaryK = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("BSCREENLENGTH") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "There was an error parsing input file at Boundary parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.screenLength = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("BDIAMETER") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "There was an error parsing input file at Boundary parameters. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.diameter = atof((lineVector[1].c_str()));
            }
            else {}
        }
        
        else {}
    }
    //Set system parameters
    SystemParameters::BParams = BParams;
}

BoundaryType SystemParser::readBoundaryType() {
        
    _inputFile.clear();
    _inputFile.seekg(0);
    
    BoundaryType BType;
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("BOUNDARYSHAPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "A boundary shape needs to be specified. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BType.boundaryShape = lineVector[1];
            }
        }
    }
    return BType;
}

void SystemParser::readGeometryParameters() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    GeometryParameters GParams;
    
    vector<double> gridTemp;
    vector<double> compartmentTemp;
    double monomerSize = 0;
    double cylinderSize = 0;
    short nDim = 0;
    
    //find grid size lines
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("NX") != string::npos
            || line.find("NY") != string::npos
            || line.find("NZ") != string::npos) {
            
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at grid dimensions. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2)
                gridTemp.push_back(atof((lineVector[1].c_str())));
            else {}
        }
        
        else if (line.find("COMPARTMENTSIZEX") != string::npos
                 || line.find("COMPARTMENTSIZEY") != string::npos
                 || line.find("COMPARTMENTSIZEZ") != string::npos) {
            
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at compartment size. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                compartmentTemp.push_back(atof((lineVector[1].c_str())));
            else {}
        }
        
        else if(line.find("MONOMERSIZE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "A monomer size needs to be specified. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                monomerSize = atof(lineVector[1].c_str());
            else {}
        }
        
        else if(line.find("CYLINDERSIZE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "A cylinder size needs to be specified. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                cylinderSize = atof(lineVector[1].c_str());
            else {}
        }
        
        else if(line.find("NDIM") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Number of dimensions needs to be specified. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                nDim = short(atoi(lineVector[1].c_str()));
            }
            else {}
        }
        else {}
    }
    //set geometry parameters and return
    GParams.nDim = nDim;
    GParams.cylinderSize = cylinderSize;
    GParams.monomerSize = monomerSize;
    
#ifdef CHEMISTRY
    if(cylinderSize / monomerSize < 5) {
        cout << "With chemistry, cylinder size specified needs to be at least 5 monomers long. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    GParams.cylinderIntSize = int(cylinderSize / monomerSize);
#endif
    
    if(gridTemp.size() >= 1) GParams.NX = gridTemp[0];
    if(gridTemp.size() >= 2) GParams.NY = gridTemp[1];
    if(gridTemp.size() >= 3) GParams.NZ = gridTemp[2];
    if(compartmentTemp.size() >= 1) GParams.compartmentSizeX = compartmentTemp[0];
    if(compartmentTemp.size() >= 2) GParams.compartmentSizeY = compartmentTemp[1];
    if(compartmentTemp.size() >= 3) GParams.compartmentSizeZ = compartmentTemp[2];
    
    //find max compartment side
    if(GParams.compartmentSizeX > GParams.largestCompartmentSide) GParams.largestCompartmentSide = GParams.compartmentSizeX;
    if(GParams.compartmentSizeY > GParams.largestCompartmentSide) GParams.largestCompartmentSide = GParams.compartmentSizeY;
    if(GParams.compartmentSizeZ > GParams.largestCompartmentSide) GParams.largestCompartmentSide = GParams.compartmentSizeZ;
    
    SystemParameters::GParams = GParams;
}

FilamentSetup SystemParser::readFilamentSetup() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    FilamentSetup FSetup;
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if(line.find("FILAMENTFILE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading filament input file. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.inputFile = lineVector[1];
            else {}
        }
        else if(line.find("NUMFILAMENTS") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading number of filaments. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.numFilaments = atoi(lineVector[1].c_str());
            else {}
        }
    }
    return FSetup;
}

vector<vector<vector<double>>> FilamentParser::readFilaments() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    vector<vector<vector<double>>> returnVector;
    string line;
    
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        vector<string> lineVector = split<string>(line);
        if(lineVector.size() == 7) {
            vector<double> coord1;
            vector<double> coord2;
            for(auto it = lineVector.begin() + 1; it != lineVector.begin() + 4; it++) {
                coord1.push_back(atof(((*it).c_str())));
            }
            for(auto it = lineVector.begin() + 4; it != lineVector.end(); it++) {
                coord2.push_back(atof(((*it).c_str())));
            }
            
            returnVector.push_back({coord1, coord2});
        }
    }
    return returnVector;
}

ChemistryData ChemistryParser::readChemistryInput() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    ChemistryData chem;
    string line;

    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        else if(line.find("SPECIESBULK") != string::npos) {
        
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a bulk species. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                chem.speciesBulk.push_back(tuple<string, int>(lineVector[1], atoi(lineVector[2].c_str())));
            }
            else {}
        }
        else if(line.find("SPECIESDIFFUSING") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  4) {
                cout << "Error reading a diffusing species. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 4) {
                chem.speciesDiffusing.push_back(tuple<string, int, double>
                        (lineVector[1], atoi(lineVector[2].c_str()), atof(lineVector[3].c_str())));
            }
            else {}
        }
        
        else if(line.find("SPECIESFILAMENT") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament species. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chem.speciesFilament.push_back(lineVector[1]);
            else {}
        }
        else if(line.find("SPECIESBOUND") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament bound species. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chem.speciesBound.push_back(lineVector[1]);
            else {}
        }
        
        else if(line.find("SPECIESLINKER") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament linker species. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chem.speciesLinker.push_back(lineVector[1]);
            else {}
        }
        else if(line.find("SPECIESMOTOR") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament linker species. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chem.speciesMotor.push_back(lineVector[1]);
            else {}
        }
        
        else if(line.find("SPECIESPLUSEND") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament plus end species. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) 
                chem.speciesPlusEnd.push_back(lineVector[1]);
            else {}
        }
        else if(line.find("SPECIESMINUSEND") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament minus end species. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chem.speciesMinusEnd.push_back(lineVector[1]);
            else {}
        }
        
        //loop through a reaction
       else if(line.find("GENREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
   
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.genReactions.push_back(tuple<vector<string>, vector<string>, double>
                                                  (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a general reaction. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("BULKREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.bulkReactions.push_back(tuple<vector<string>, vector<string>, double>
                                                  (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a bulk reaction. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("DEPOLYMERIZATIONREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.depolymerizationReactions.push_back(tuple<vector<string>, vector<string>, double>
                                                           (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a depolymerization reaction. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("POLYMERIZATIONREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.polymerizationReactions.push_back(tuple<vector<string>, vector<string>, double>
                                                  (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a polymerization reaction. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(line.find("LINKERBINDINGREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 3; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.linkerBindingReactions.push_back(tuple<vector<string>, vector<string>, double, double, double>
                     (reactants, products, atof(lineVector[lineVector.size() - 3].c_str()),
                     atof(lineVector[lineVector.size() - 2].c_str()), atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a linker binding reaction. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("MOTORBINDINGREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 3; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.motorBindingReactions.push_back(tuple<vector<string>, vector<string>, double, double, double>
                     (reactants, products, atof(lineVector[lineVector.size() - 3].c_str()),
                     atof(lineVector[lineVector.size() - 2].c_str()), atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a motor binding reaction. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        
        else if(line.find("UNBINDINGREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.unbindingReactions.push_back(tuple<vector<string>, vector<string>, double>
                     (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading an unbinding reaction. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("BINDINGREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.bindingReactions.push_back(tuple<vector<string>, vector<string>, double>
                     (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a binding reaction. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("MOTORWALKINGREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.motorWalkingReactions.push_back(tuple<vector<string>, vector<string>, double>
                     (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a motor walking reaction. Exiting" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    return chem;
}



