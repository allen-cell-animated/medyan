
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "Parser.h"

#include "SysParams.h"

OutputTypes SystemParser::readOutputTypes() {
    
    OutputTypes oTypes;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("OUTPUTTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 3) {
                cout <<
                "There was an error parsing input file at output types. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector[1] == "SNAPSHOT")   oTypes.basicSnapshot = true;
            
            else if (lineVector[1] == "BIRTHTIMES") oTypes.birthTimes = true;
            
            else if (lineVector[1] == "FORCES")     oTypes.forces = true;
            
            else if (lineVector[1] == "STRESSES")   oTypes.stresses = true;
        }
    }
    return oTypes;
}

void SystemParser::readChemParams() {
    
    ChemParams CParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("NUMBULKSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numBulkSpecies = atof(lineVector[1].c_str());
            }
        }
        
        if (line.find("NUMDIFFUSINGSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numDiffusingSpecies = atof(lineVector[1].c_str());
            }
        }
    
        if (line.find("NUMFILAMENTSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numFilamentSpecies = atof(lineVector[1].c_str());
            }
        }
    
        if (line.find("NUMPLUSENDSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numPlusEndSpecies = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMMINUSENDSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numMinusEndSpecies = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMBOUNDSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numBoundSpecies = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMLINKERSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numLinkerSpecies = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMMOTORSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numMotorSpecies = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMBRANCHERSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numBrancherSpecies = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMBINDINGSITES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numBindingSites = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMMOTORHEADSMIN") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.motorNumHeadsMin.push_back(atoi(lineVector[i].c_str()));
            }
        }
        if (line.find("NUMMOTORHEADSMAX") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.motorNumHeadsMax.push_back(atoi(lineVector[i].c_str()));
            }
        }
        if (line.find("MOTORSTEPSIZE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.motorStepSize.push_back(atof(lineVector[i].c_str()));
            }
        }
        
    }
    //Figure out the binding sites
    int deltaBinding = SysParams::Geometry().cylinderIntSize /
                       CParams.numBindingSites;
    
    int firstBindingSite = deltaBinding / 2 + 1;
    int bindingCount = firstBindingSite;
    
    //add all other binding sites
    while(bindingCount < SysParams::Geometry().cylinderIntSize) {
        
        CParams.bindingSites.push_back(bindingCount);
        bindingCount += deltaBinding;
    }
    
    //set system parameters
    SysParams::CParams = CParams;
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
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.algorithm = lineVector[1];
            }
        }
        if (line.find("NUMTOTALSTEPS:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.numTotalSteps = atoi(lineVector[1].c_str());
            }
        }
        
        if (line.find("RUNTIME:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.runTime = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NUMSTEPSPERS:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.numStepsPerSnapshot = atoi(lineVector[1].c_str());
            }
        }
        if (line.find("SNAPSHOTTIME:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.snapshotTime = atof(lineVector[1].c_str());
            }
        }
        
        if (line.find("NUMCHEMSTEPS:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.numChemSteps = atoi(lineVector[1].c_str());
            }
        }
        if (line.find("NUMSTEPSPERN:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.numStepsPerNeighbor = atoi(lineVector[1].c_str());
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
                cout << "Error reading chemistry input file. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            
            else if (lineVector.size() == 2)
                CSetup.inputFile = lineVector[1];
            
            else if(lineVector.size() < 2) {
                cout << "Must specify a chemistry input file. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
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
                cout <<
                "There was an error parsing input file at Filament stretching type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FStretchingType = lineVector[1];
            }
        }
        else if (line.find("FBENDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Filament bending type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FBendingType = lineVector[1];
            }
        }
        else if (line.find("FTWISTINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Filament twisting type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.FTwistingType = lineVector[1];
            }
        }
        else if (line.find("LSTRETCHINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Linker stretching type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LStretchingType = lineVector[1];
            }
        }
        else if (line.find("LBENDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Linker bending type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LBendingType = lineVector[1];
            }
        }
        else if (line.find("LTWISTINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Linker twisting type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.LTwistingType = lineVector[1];
            }
        }
        else if (line.find("MSTRETCHINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Motor stretching type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MStretchingType = lineVector[1];
            }
        }
        else if (line.find("MBENDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Motor bending type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MBendingType = lineVector[1];
            }
        }
        else if (line.find("MTWISTINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Motor twisting type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.MTwistingType = lineVector[1];
            }
        }
        if (line.find("BRSTRETCHINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Branch stretching type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.BrStretchingType = lineVector[1];
            }
        }
        else if (line.find("BRBENDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Branch bending type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.BrBendingType = lineVector[1];
            }
        }
        else if (line.find("BRDIHEDRALTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Branch dihedral type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.BrDihedralType = lineVector[1];
            }
        }
        else if (line.find("BRPOSITIONTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Branch dihedral type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.BrPositionType = lineVector[1];
            }
        }
        else if (line.find("BOUNDARYTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Boundary type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MTypes.BoundaryFFType = lineVector[1];
            }
        }
        else if (line.find("VOLUMETYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Volume type. Exiting."
                << endl;
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

void SystemParser::readMechParams() {
    
    MechParams MParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        //Filament stretching
        if (line.find("FSTRETCHINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Mechanics parameters. Exiting."
                << endl;
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
                cout <<
                "There was an error parsing input file at Mechanics parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FBendingK = atof((lineVector[1].c_str()));
            }
        }
        else if (line.find("FBENDINGTHETA") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Mechanics parameters. Exiting."
                << endl;
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
                cout <<
                "There was an error parsing input file at Mechanics parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.FTwistingK = atof((lineVector[1].c_str()));
            }
        }
        else if (line.find("FTWISTINGPHI") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Mechanics parameters. Exiting."
                << endl;
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
        
        //Branch dihedral
        else if (line.find("BRDIHEDRALK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BrDihedralK.push_back(atof((lineVector[i].c_str())));
            }
        }
        //Branch position
        else if (line.find("BRPOSITIONK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BrPositionK.push_back(atof((lineVector[i].c_str())));
            }
        }
        
        //Volume parameter
        else if (line.find("VOLUMEK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Mechanics parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.VolumeK = atof((lineVector[1].c_str()));
            }
        }
        
        else if (line.find("VOLUMECUTOFF") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Mechanics parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.VolumeCutoff = atof((lineVector[1].c_str()));
            }
        }
        else {}
    }
    //Set system parameters
    SysParams::MParams = MParams;
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
                cout <<
                "A conjugate gradient method must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MAlgorithm.ConjugateGradient = lineVector[1];
            }
        }
        else if (line.find("MD") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "A Mechanics algorithm must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MAlgorithm.MD = lineVector[1];
            }
        }
        else if (line.find("GRADIENTTOLERANCE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                MAlgorithm.gradientTolerance = atof(lineVector[1].c_str());
            }
        }
        else if (line.find("ENERGYTOLERANCE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                MAlgorithm.energyTolerance = atof(lineVector[1].c_str());
            }
        }
        else if (line.find("MAXDISTANCE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                MAlgorithm.maxDistance = atof(lineVector[1].c_str());
            }
        }
    }
    return MAlgorithm;
}

void SystemParser::readBoundParams() {
    
    BoundParams BParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("BCUTOFF") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Boundary parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.BoundaryCutoff = atof((lineVector[1].c_str()));
            }
            //Default value to be half compartment size
            else {
                BParams.BoundaryCutoff =
                SysParams::Geometry().compartmentSizeX / 2;
            }
        }
        else if (line.find("BINTERACTIONK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Boundary parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.BoundaryK = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("BSCREENLENGTH") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Boundary parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BParams.BScreenLength = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("BDIAMETER") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                BParams.diameter = atof((lineVector[1].c_str()));
            }
            else {}
        }
        
        else {}
    }
    //Set system parameters
    SysParams::BParams = BParams;
}

void SystemParser::readDyRateParams() {
    
    DyRateParams DRParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }

        if (line.find("DFPOLYMERIZATIONLEN") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at dynamic rate parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                DRParams.dFilPolymerizationCharLength =
                          atof((lineVector[1].c_str()));
            }
            else {}
        }

        else if (line.find("DMUNBINDINGFORCE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRParams.dMotorUnbindingCharForce.push_back(
                                  atof((lineVector[i].c_str())));
            }
            else {}
        }
        
        else if (line.find("DMWALKINGFORCE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRParams.dMotorWalkingCharForce.push_back(
                                atof((lineVector[i].c_str())));
            }
            else {}
        }
        
        else if (line.find("DLUNBINDINGLEN") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRParams.dLinkerUnbindingCharLength.push_back(
                                    atof((lineVector[i].c_str())));
            }
            else {}
        }
        
        else if (line.find("DLUNBINDINGAMP") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRParams.dLinkerUnbindingAmplitude.push_back(
                                   atof((lineVector[i].c_str())));
            }
            else {}
        }
    }
    
    //set system parameters
    SysParams::DRParams = DRParams;
}

DynamicRateTypes SystemParser::readDynamicRateTypes() {
    
    DynamicRateTypes DRTypes;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("DFPOLYMERIZATIONTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at dynamic polymerization type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                DRTypes.dFPolymerizationType = lineVector[1];
            }
        }
        
        else if (line.find("DMUNBINDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRTypes.dMUnbindingType.push_back(lineVector[i]);
            }
        }
        
        else if (line.find("DMWALKINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRTypes.dMWalkingType.push_back(lineVector[i]);
            }
        }
        
        else if (line.find("DLUNBINDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRTypes.dLUnbindingType.push_back(lineVector[i]);
            }
        }
    }
    return DRTypes;
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
                cout << "A boundary shape needs to be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BType.boundaryShape = lineVector[1];
            }
        }
    }
    return BType;
}

void SystemParser::readGeoParams() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    GeoParams GParams;
    
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
                cout << "There was an error parsing input file at grid dimensions. Exiting." << endl;
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
                cout << "There was an error parsing input file at compartment size. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                compartmentTemp.push_back(atof((lineVector[1].c_str())));
            else {}
        }
        
        else if(line.find("MONOMERSIZE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "A monomer size needs to be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                monomerSize = atof(lineVector[1].c_str());
            else {}
        }
        
        else if(line.find("CYLINDERSIZE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "A cylinder size needs to be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                cylinderSize = atof(lineVector[1].c_str());
            else {}
        }
        
        else if(line.find("NDIM") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Number of dimensions needs to be specified. Exiting." << endl;
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
    if(cylinderSize / monomerSize < SysParams::Geometry().minCylinderIntSize) {
        cout <<
        "With chemistry, cylinder size specified is too short. Exiting."
        << endl;
        exit(EXIT_FAILURE);
    }
#endif
    GParams.cylinderIntSize = int(cylinderSize / monomerSize);
    
    GParams.minCylinderSize =
        SysParams::Geometry().minCylinderIntSize * GParams.monomerSize;
    
    if(gridTemp.size() >= 1) GParams.NX = gridTemp[0];
    if(gridTemp.size() >= 2) GParams.NY = gridTemp[1];
    if(gridTemp.size() >= 3) GParams.NZ = gridTemp[2];
    
    if(compartmentTemp.size() >= 1) GParams.compartmentSizeX = compartmentTemp[0];
    if(compartmentTemp.size() >= 2) GParams.compartmentSizeY = compartmentTemp[1];
    if(compartmentTemp.size() >= 3) GParams.compartmentSizeZ = compartmentTemp[2];
    
    //find max compartment side
    if(GParams.compartmentSizeX > GParams.largestCompartmentSide)
        GParams.largestCompartmentSide = GParams.compartmentSizeX;
    
    if(GParams.compartmentSizeY > GParams.largestCompartmentSide)
        GParams.largestCompartmentSide = GParams.compartmentSizeY;
    
    if(GParams.compartmentSizeZ > GParams.largestCompartmentSide)
        GParams.largestCompartmentSide = GParams.compartmentSizeZ;
    
    SysParams::GParams = GParams;
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
                cout << "Error reading filament input file. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.inputFile = lineVector[1];
            else {}
        }
        else if(line.find("NUMFILAMENTS") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading number of filaments. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.numFilaments = atoi(lineVector[1].c_str());
            else {}
        }
        else if(line.find("FILAMENTLENGTH") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading filament length. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.filamentLength = atoi(lineVector[1].c_str());
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
            if(lineVector.size() !=  5) {
                cout << "Error reading a bulk species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 5) {
                
                if(lineVector[4] != "CONST" && lineVector[4] != "REG") {
                    
                    cout << "Option for bulk species not valid. Exiting." << endl;
                    cout << lineVector[4] << endl;
                    exit(EXIT_FAILURE);
                }
                
                chem.speciesBulk.push_back(tuple<string, int, double, string>
                    (lineVector[1], atoi(lineVector[2].c_str()),
                     atof(lineVector[3].c_str()), lineVector[4]));
            }
            else {}
        }
        else if(line.find("SPECIESDIFFUSING") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() >  7 || lineVector.size() < 6) {
                cout << "Error reading a diffusing species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 7) {
                
                if(lineVector[5] != "AVG") {
                    
                    cout << "Too many arguments for a non AVG-qualified diffusing species. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                
                chem.speciesDiffusing.push_back(tuple<string, int, double, double, string, int>
                    (lineVector[1], atoi(lineVector[2].c_str()),
                     atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                     lineVector[5], atoi(lineVector[6].c_str())));
            }
            else if (lineVector.size() == 6) {
                
                if(lineVector[5] != "REG") {
                    
                    cout << "Not enough arguments for a non REG-qualified diffusing species. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                
                chem.speciesDiffusing.push_back(tuple<string, int, double, double, string, int>
                     (lineVector[1], atoi(lineVector[2].c_str()),
                     atof(lineVector[3].c_str()), atof(lineVector[4].c_str()), lineVector[5], 0));
            }
            else {}
        }
        
        else if(line.find("SPECIESFILAMENT") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chem.speciesFilament.push_back(lineVector[1]);
            else {}
        }
        else if(line.find("SPECIESBOUND") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament bound species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chem.speciesBound.push_back(lineVector[1]);
            else {}
        }
        
        else if(line.find("SPECIESLINKER") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament linker species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chem.speciesLinker.push_back(lineVector[1]);
            else {}
        }
        else if(line.find("SPECIESMOTOR") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament motor species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chem.speciesMotor.push_back(lineVector[1]);
            else {}
        }
        else if(line.find("SPECIESBRANCHER") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament brancher species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                chem.speciesBrancher.push_back(lineVector[1]);
            else {}
        }
        else if(line.find("SPECIESPLUSEND") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament plus end species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) 
                chem.speciesPlusEnd.push_back(lineVector[1]);
            else {}
        }
        else if(line.find("SPECIESMINUSEND") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  2) {
                cout << "Error reading a filament minus end species. Exiting." << endl;
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
                
                chem.genReactions.push_back(
                tuple<vector<string>, vector<string>, double>
                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a general reaction. Exiting." << endl;
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
                
                chem.bulkReactions.push_back(
                tuple<vector<string>, vector<string>, double>
                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a bulk reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("NUCLEATIONREACTION") != string::npos) {
            
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
                
                chem.nucleationReactions.push_back(
                tuple<vector<string>, vector<string>, double>
                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
            }
            else {
                cout << "Error reading a nucleation reaction. Exiting." << endl;
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
                
                chem.depolymerizationReactions.push_back(
                tuple<vector<string>, vector<string>, double>
                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a depolymerization reaction. Exiting." << endl;
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
                    
                    if(*it != "+") products.push_back((*it));
                }
                
                chem.polymerizationReactions.push_back(
                tuple<vector<string>, vector<string>, double>
                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a polymerization reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(line.find("LINKERREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "<->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 4; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.linkerReactions.push_back(
                tuple<vector<string>, vector<string>, double, double, double, double>
                (reactants, products, atof(lineVector[lineVector.size() - 4].c_str()),
                                      atof(lineVector[lineVector.size() - 3].c_str()),
                                      atof(lineVector[lineVector.size() - 2].c_str()),
                                      atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a linker reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(line.find("MOTORREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "<->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 4; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.motorReactions.push_back(
                tuple<vector<string>, vector<string>, double, double, double, double>
                (reactants, products, atof(lineVector[lineVector.size() - 4].c_str()),
                                    atof(lineVector[lineVector.size() - 3].c_str()),
                                    atof(lineVector[lineVector.size() - 2].c_str()),
                                    atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a motor reaction. Exiting." << endl;
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
                
                chem.motorWalkingReactions.push_back(
                tuple<vector<string>, vector<string>, double>
                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a motor walking reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("AGINGREACTION") != string::npos) {
            
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
                
                chem.agingReactions.push_back(
                tuple<vector<string>, vector<string>, double>
                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading an aging reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("DESTRUCTIONREACTION") != string::npos) {
            
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
                
                chem.destructionReactions.push_back(
                tuple<vector<string>, vector<string>, double>
                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
            }
            else {
                cout << "Error reading a destruction reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("BRANCHINGREACTION") != string::npos) {
            
            vector<string> reactants;
            vector<string> products;
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "<->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 4; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.branchingReactions.push_back(
                tuple<vector<string>, vector<string>, double, double, string, double>
                (reactants, products, atof(lineVector[lineVector.size() - 4].c_str()),
                 atof(lineVector[lineVector.size() - 3].c_str()),
                      lineVector[lineVector.size() - 2].c_str(),
                 atof(lineVector[lineVector.size() - 1].c_str())));
            }
            else {
                cout << "Error reading a branching reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(line.find("SEVERINGREACTION") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "AT");
            if(arrowIt != lineVector.end()) {
                
                auto it = arrowIt + 1;
                
                chem.severingReactions.push_back(tuple<string, double>
                ((*it), atof(lineVector[lineVector.size() - 1].c_str())));
            }
            else {
                cout << "Error reading a severing reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
    }
    return chem;
}
