
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "Parser.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

#include "SysParams.h"

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
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.numFilamentSpecies.push_back(atoi(lineVector[i].c_str()));
            }
        }
    
        if (line.find("NUMPLUSENDSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.numPlusEndSpecies.push_back(atoi(lineVector[i].c_str()));
            }
        }
        if (line.find("NUMMINUSENDSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.numMinusEndSpecies.push_back(atoi(lineVector[i].c_str()));
            }
        }
        if (line.find("NUMBOUNDSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.numBoundSpecies.push_back(atoi(lineVector[i].c_str()));
            }
        }
        if (line.find("NUMLINKERSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.numLinkerSpecies.push_back(atoi(lineVector[i].c_str()));
            }
        }
        if (line.find("NUMMOTORSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.numMotorSpecies.push_back(atoi(lineVector[i].c_str()));
            }
        }
        if (line.find("NUMBRANCHERSPECIES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.numBrancherSpecies.push_back(atoi(lineVector[i].c_str()));
            }
        }
        if (line.find("NUMFILAMENTTYPES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.numFilaments = atoi(lineVector[1].c_str());
            }
        }
        
        if (line.find("NUMBINDINGSITES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    CParams.numBindingSites.push_back(atoi(lineVector[i].c_str()));
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
        
        if (line.find("SPECIALPROTOCOL") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if(lineVector.size() > 4) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                
                if(lineVector[1] == "MAKELINKERSSTATIC") {
                    CParams.makeLinkersStatic = true;
                    CParams.makeLinkersStaticTime = atof(lineVector[2].c_str());
                }
                if(lineVector[1] == "MAKEFILAMENTSSTATIC") {
                    CParams.makeFilamentsStatic = true;
                    CParams.makeFilamentsStaticTime = atof(lineVector[2].c_str());
                }
            }
        }
    }
    
    //Figure out the binding sites
    for(int i = 0; i < CParams.numBindingSites.size(); i++) {
    
        vector<short> tempBindingSites;
        
        int deltaBinding = SysParams::Geometry().cylinderNumMon[i] /
                           CParams.numBindingSites[i];
    
        int firstBindingSite = deltaBinding / 2 + 1;
        int bindingCount = firstBindingSite;
        
        //add all other binding sites
        while(bindingCount < SysParams::Geometry().cylinderNumMon[i]) {
            tempBindingSites.push_back(bindingCount);
            bindingCount += deltaBinding;
        }
        
        
        //push to CParams
        CParams.bindingSites.push_back(tempBindingSites);
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
        if (line.find("RUNSTEPS:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.runSteps = atoi(lineVector[1].c_str());
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
        if (line.find("SNAPSHOTSTEPS:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.snapshotSteps = atoi(lineVector[1].c_str());
            }
        }
        if (line.find("MINIMIZATIONTIME:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.minimizationTime = atof(lineVector[1].c_str());
            }
        }
        if (line.find("MINIMIZATIONSTEPS:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.minimizationSteps = atoi(lineVector[1].c_str());
            }
        }
        if (line.find("NEIGHBORLISTTIME:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.neighborListTime = atof(lineVector[1].c_str());
            }
        }
        if (line.find("NEIGHBORLISTSTEPS:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.neighborListSteps = atoi(lineVector[1].c_str());
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
    
    MechanicsFFType MType;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("FSTRETCHINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Filament stretching FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.FStretchingType = lineVector[1];
            }
        }
        else if (line.find("FBENDINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Filament bending FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.FBendingType = lineVector[1];
            }
        }
        else if (line.find("FTWISTINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Filament twisting FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.FTwistingType = lineVector[1];
            }
        }
        else if (line.find("LSTRETCHINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Linker stretching FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.LStretchingType = lineVector[1];
            }
        }
        else if (line.find("LBENDINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Linker bending FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.LBendingType = lineVector[1];
            }
        }
        else if (line.find("LTWISTINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Linker twisting FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.LTwistingType = lineVector[1];
            }
        }
        else if (line.find("MSTRETCHINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Motor stretching FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.MStretchingType = lineVector[1];
            }
        }
        else if (line.find("MBENDINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Motor bending FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.MBendingType = lineVector[1];
            }
        }
        else if (line.find("MTWISTINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Motor twisting FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.MTwistingType = lineVector[1];
            }
        }
        if (line.find("BRSTRETCHINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Branch stretching FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.BrStretchingType = lineVector[1];
            }
        }
        else if (line.find("BRBENDINGFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Branch bending FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.BrBendingType = lineVector[1];
            }
        }
        else if (line.find("BRDIHEDRALFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Branch dihedral FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.BrDihedralType = lineVector[1];
            }
        }
        else if (line.find("BRPOSITIONFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Branch position FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.BrPositionType = lineVector[1];
            }
        }
        else if (line.find("BOUNDARYFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Boundary FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.BoundaryFFType = lineVector[1];
            }
        }
        else if (line.find("VOLUMEFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Volume FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.VolumeFFType = lineVector[1];
            }
        }
        else if (line.find("BUBBLEFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at Bubble FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.BubbleFFType = lineVector[1];
            }
        }
        else if (line.find("MTOCFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at MTOC FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.MTOCFFType = lineVector[1];
            }
        }
        else if(line.find("MEM_STRETCHING_FF_TYPE") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at membrane stretching FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2) {
                MType.MemStretchingFFType = lineVector[1];
            }
        }
        else if(line.find("MEM_BENDING_FF_TYPE") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at membrane bending FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2) {
                MType.MemBendingFFType = lineVector[1];
            }
        }
        else if(line.find("MEM_CYLINDER_VOLUME_FF_TYPE") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at membrane cylinder volume FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2) {
                MType.MemCylinderVolumeFFType = lineVector[1];
            }
        }
        
        else {}
    }
    return MType;
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
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.FStretchingK.push_back(atof((lineVector[i].c_str())));
            }
        }

        //Filament bending
        else if (line.find("FBENDINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.FBendingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("FBENDINGTHETA") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.FBendingTheta.push_back(atof((lineVector[i].c_str())));
            }
        }
        
        //Filament twisting
        else if (line.find("FTWISTINGK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.FTwistingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("FTWISTINGPHI") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.FTwistingPhi.push_back(atof((lineVector[i].c_str())));
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
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.VolumeK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("VOLUMECUTOFF") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "Error reading Volume cutoff. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.VolumeCutoff = atoi(lineVector[1].c_str());
            }
        }
        else if (line.find("MEM_CYLINDER_VOLUME_K") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.MemCylinderVolumeK.push_back(atof(lineVector[i].c_str()));
            }
        }
        else if (line.find("MEM_CYLINDER_VOLUME_CUTOFF") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "Error reading Volume cutoff. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.MemCylinderVolumeCutoff = atoi(lineVector[1].c_str());
            }
        }
        
        //Bubble parameter
        else if (line.find("BUBBLEINTERACTIONK") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BubbleK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("BUBBLESCREENLENGTH") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BubbleScreenLength.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("BUBBLECUTOFF") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "Error reading Bubble cutoff. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.BubbleCutoff = atoi(lineVector[1].c_str());
            }
        }
        else if (line.find("BUBBLERADIUS") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.BubbleRadius.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("NUMBUBBLETYPES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "Error reading number of Bubble types. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.numBubbleTypes = atoi(lineVector[1].c_str());
            }
        }

        // Membrane parameter
        else if (line.find("MEM_ELASTIC_K") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() >= 2) {
                for(size_t i = 1; i < lineVector.size(); ++i)
                    MParams.MemElasticK.push_back(atof(lineVector[i].c_str()));
            }
        }
        else if (line.find("MEM_BENDING_K") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() >= 2) {
                for(size_t i = 1; i < lineVector.size(); ++i)
                    MParams.MemBendingK.push_back(atof(lineVector[i].c_str()));
            }
        }
        else if (line.find("MEM_EQ_CURV") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() >= 2) {
                for(size_t i = 1; i < lineVector.size(); ++i)
                    MParams.MemEqCurv.push_back(atof(lineVector[i].c_str()));
            }
        }
        
        if (line.find("SPECIALPROTOCOL") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if(lineVector.size() > 4) {
                cout <<
                "There was an error parsing input file at Chemistry parameters. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 4) {
                
                if(lineVector[1] == "PINBOUNDARYFILAMENTS") {
                    
                    MParams.pinBoundaryFilaments = true;
                    MParams.pinK = atof(lineVector[2].c_str());
                    MParams.pinTime = atof(lineVector[3].c_str());
                    
                }
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
        else if (line.find("MAXDISTANCE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                MAlgorithm.maxDistance = atof(lineVector[1].c_str());
            }
        }
        else if (line.find("LAMBDAMAX") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                MAlgorithm.lambdaMax = atof(lineVector[1].c_str());
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
        
        if (line.find("BOUNDARYCUTOFF") != string::npos) {
            
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
            //Default value to be compartment size
            else {
                BParams.BoundaryCutoff =
                SysParams::Geometry().compartmentSizeX;
            }
        }
        else if (line.find("BOUNDARYINTERACTIONK") != string::npos) {
            
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
        else if (line.find("BOUNDARYSCREENLENGTH") != string::npos) {
            
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
        else if (line.find("BOUNDARYDIAMETER") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                BParams.diameter = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("BMOVESPEED") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                BParams.moveSpeed = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("BMOVESTARTTIME") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                BParams.moveStartTime = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("BMOVEENDTIME") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                BParams.moveEndTime = atof((lineVector[1].c_str()));
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
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRParams.dFilPolymerizationCharLength.push_back(
                                  atof((lineVector[i].c_str())));
            }
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

DynamicRateType SystemParser::readDynamicRateType() {
    
    DynamicRateType DRType;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("DFPOLYMERIZATIONTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRType.dFPolymerizationType.push_back(lineVector[i]);
            }
        }
        
        else if (line.find("DMUNBINDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRType.dMUnbindingType.push_back(lineVector[i]);
            }
        }
        
        else if (line.find("DMWALKINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRType.dMWalkingType.push_back(lineVector[i]);
            }
        }
        
        else if (line.find("DLUNBINDINGTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRType.dLUnbindingType.push_back(lineVector[i]);
            }
        }
    }
    return DRType;
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
        else if (line.find("BOUNDARYMOVE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "A boundary move type needs to be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                BType.boundaryMove = lineVector[1];
            }
        }
    }
    return BType;
}

SpecialSetupType SystemParser::readSpecialSetupType() {
    
    SpecialSetupType SType;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if (line.find("SPECIALSETUP") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 3) {
                cout <<
                "There was an error parsing input file at special setup types. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector[1] == "MTOC") SType.mtoc = true;
        }
        else if (line.find("MTOCFILAMENTTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "A filament type to connect to an MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                SType.mtocFilamentType = atoi(lineVector[1].c_str());
            }
        }
        else if (line.find("MTOCNUMFILAMENTS") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "A number of filaments to connect to an MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                SType.mtocNumFilaments = atoi(lineVector[1].c_str());
            }
        }
        else if (line.find("MTOCFILAMENTLENGTH") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "A filament length for filaments connected to an MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                SType.mtocFilamentLength = atoi(lineVector[1].c_str());
            }
        }
        
        else if (line.find("MTOCBUBBLETYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "A bubble type to connect to an MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                SType.mtocBubbleType = atoi(lineVector[1].c_str());
            }
        }
    }
    return SType;
}


void SystemParser::readGeoParams() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    GeoParams GParams;
    
    vector<double> gridTemp;
    vector<double> compartmentTemp;
    vector<double> monomerSize = {};
    vector<double> cylinderSize = {};
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
            
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    monomerSize.push_back(atof((lineVector[i].c_str())));
            }
            else {}
        }
        
        else if(line.find("CYLINDERSIZE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    cylinderSize.push_back(atof((lineVector[i].c_str())));
            }
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
    
    if(GParams.cylinderSize.size() != GParams.monomerSize.size()) {
        
        cout << "Must specify an equivalent number of cylinder and monomer sizes. Exiting."
             << endl;
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0; i < GParams.cylinderSize.size(); i++) {
        
#ifdef CHEMISTRY
        if(cylinderSize[i] / monomerSize[i] < SysParams::Geometry().minCylinderNumMon) {
            cout <<
            "With chemistry, cylinder size specified is too short. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
#endif
        GParams.cylinderNumMon.push_back(int(cylinderSize[i] / monomerSize[i]));
        
        GParams.minCylinderSize.push_back(
        SysParams::Geometry().minCylinderNumMon * GParams.monomerSize[i]);
        
    }
        
    if(gridTemp.size() >= 1) GParams.NX = gridTemp[0];
    if(gridTemp.size() >= 2) GParams.NY = gridTemp[1];
    if(gridTemp.size() >= 3) GParams.NZ = gridTemp[2];
    
    if(compartmentTemp.size() >= 1) GParams.compartmentSizeX = compartmentTemp[0];
    if(compartmentTemp.size() >= 2) GParams.compartmentSizeY = compartmentTemp[1];
    if(compartmentTemp.size() >= 3) GParams.compartmentSizeZ = compartmentTemp[2];
    
    //find max compartment side
    GParams.largestCompartmentSide = max(GParams.compartmentSizeX,
                                     max(GParams.compartmentSizeY, GParams.compartmentSizeZ));
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
        else if(line.find("NUMFILAMENTS") != string::npos &&
                line.find("NUMFILAMENTSPECIES") == string::npos &&
                line.find("MTOCNUMFILAMENTS") == string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading number of filaments. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.numFilaments = atoi(lineVector[1].c_str());
            else {}
        }
        else if(line.find("FILAMENTLENGTH") != string::npos &&
                line.find("MTOCFILAMENTLENGTH") == string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading filament length. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.filamentLength = atoi(lineVector[1].c_str());
            else {}
        }
        else if(line.find("FILAMENTTYPE") != string::npos &&
                line.find("NUMFILAMENTTYPES") == string::npos &&
                line.find("MTOCFILAMENTTYPE") == string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading filament type. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.filamentType = atoi(lineVector[1].c_str());
            else {}
        }
        else if (line.find("RESTARTPHASE") != string::npos){SysParams::RUNSTATE=false;}
        else if(line.find("PROJECTIONTYPE")!=string::npos){
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading filament projection type. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.projectionType = lineVector[1];
            else {}
        }
        else if(line.find("PINRESTARTFILE")!=string::npos){
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading filament projection type. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                FSetup.pinRestartFile = lineVector[1];
            else {}
        }
    }
    return FSetup;
}

MembraneSetup SystemParser::readMembraneSetup() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    MembraneSetup MemSetup;
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if(line.find("MEMBRANE_FILE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "Error reading membrane input file. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                MemSetup.inputFile = lineVector[1];
        }
        // number of membranes will be automatically inferred from the input file
    }
    return MemSetup;
}


BubbleSetup SystemParser::readBubbleSetup() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    BubbleSetup BSetup;
    
    string line;
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        if(line.find("BUBBLEFILE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading bubble input file. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                BSetup.inputFile = lineVector[1];
            else {}
        }
        else if(line.find("NUMBUBBLES") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading number of bubbles. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                BSetup.numBubbles = atoi(lineVector[1].c_str());
            else {}
        }
        else if(line.find("BUBBLETYPE") != string::npos &&
                line.find("NUMBUBBLETYPES") == string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading bubble type. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                BSetup.bubbleType = atoi(lineVector[1].c_str());
            else {}
        }
    }
    return BSetup;
}
    tuple< vector<tuple<short, vector<double>, vector<double>>> , vector<tuple<string, short, vector<vector<double>>>> , vector<tuple<string, short, vector<double>>> , vector<vector<double>> > FilamentParser::readFilaments() {
    _inputFile.clear();
    _inputFile.seekg(0);
     vector<tuple<short, vector<double>, vector<double>>> filamentVector;
     vector<vector<vector<double>>> linkerVector;
     vector<vector<vector<double>>> motorVector;
     vector<vector<double>> staticVector;
     vector<tuple<string, short, vector<vector<double>>>> boundVector;
     vector<tuple<string, short, vector<double>>> branchVector;
     string line;
    
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        vector<string> lineVector = split<string>(line);
        if(lineVector.size() >= 8) {
            vector<double> coord1;
            vector<double> coord2;
            vector<vector<double>> coord3;
            short type;
            //aravind parse linkers, motors. June 30,2016.
            if(lineVector[0]=="FILAMENT"){
            type = atoi((*(lineVector.begin() + 1)).c_str());
            for(auto it = lineVector.begin() + 2; it != lineVector.begin() + 5; it++) {
                coord1.push_back(atof(((*it).c_str())));
            }
            for(auto it = lineVector.begin() + 5; it != lineVector.end(); it++) {
                coord2.push_back(atof(((*it).c_str())));
                
            }
                filamentVector.emplace_back(type, coord1, coord2);}
            else
            {
                type = atoi((*(lineVector.begin() + 1)).c_str());
                string boundType=lineVector[0];
                for(auto it = lineVector.begin() + 2; it != lineVector.begin() + 5; it++) {
                    coord1.push_back(atof(((*it).c_str())));
                }
                for(auto it = lineVector.begin() + 5; it != lineVector.end(); it++) {
                    coord2.push_back(atof(((*it).c_str())));
                }
                coord3.push_back(coord1);
                coord3.push_back(coord2);
                boundVector.emplace_back(boundType, type, coord3);
            }
        }
        //aravind Feb 19, 2016. Parase Linkers, Motors.
        else if(lineVector.size()==5) {
            vector<double> coord1;
            vector<vector<double>> coord3;
            if(lineVector[0]=="STATIC"){
                for(auto it = lineVector.begin() + 1; it != lineVector.begin() + 5; it++) {
                    coord1.push_back(atof(((*it).c_str()))); //FORMAT FILAMENTTYPE COORDx COORDy COORDz.
                }
                staticVector.push_back({coord1});}
            else{ // BRANCHER
                short type = atoi((*(lineVector.begin() + 1)).c_str()); //FILAMENT TYPE THAT IT BINDS TO.
                string boundType=lineVector[0];//BRANCHER BOUND NAME
                for(auto it = lineVector.begin() + 2; it != lineVector.begin() + 5; it++) {
                    coord1.push_back(atof(((*it).c_str())));
                }
                branchVector.emplace_back(boundType,type,coord1);
            }
        }
    }
      tuple< vector<tuple<short, vector<double>, vector<double>>> , vector<tuple<string, short, vector<vector<double>>>> , vector<tuple<string, short, vector<double>>> , vector<vector<double>> > returnVector=make_tuple(filamentVector,boundVector,branchVector, staticVector);
    return returnVector;
}

vector<MembraneParser::membraneInfo> MembraneParser::readMembranes() {

    vector<membraneInfo> res;
    
    bool wasEmpty = true;
    
    while(getline(_inputFile, line)) {
        
        bool isEmpty = (line.empty() || line.find("#") != string::npos); // empty line or line with '#'

        if(wasEmpty && !isEmpty) { // New membrane
            res.emplace_back(membraneInfo());
        }

        wasEmpty = isEmpty;
        if(isEmpty) continue;

        if(res.empty()) {
            cout << "This line should never be executed. Something is wrong when reading membrane information. "
                 << "Exiting." <<endl;
            exit(EXIT_FAILURE);
        }

        auto& activeMem = res.back();

        vector<string> lineVector = split<string>(line);
        size_t lineVectorSize = lineVector.size();

        if(lineVectorSize < 3) {
            cout << "Error occured when reading membrane files. "
                 << "Each vertex should have 3 doubles for position. "
                 << "Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        else {
            activeMem.emplace_back(vertexInfo());
            auto& activeVertex = activeMem.back();
            auto& activePosition = get<0>(activeVertex);
            auto& activeNeighbor = get<1>(activeVertex);

            // Parse coordinate information
            for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                activePosition[coordIdx] = lineVector[coordIdx];
            }

            // Parse neighbor indices
            activeNeighbor.reserve(lineVectorSize - 3);
            for(size_t nIdx = 0; nIdx < lineVectorSize - 3; ++nIdx) {
                activeNeighbor.push_back((size_t)atoi(lineVector[nIdx + 3].c_str()));
            }
        }
    }

    return res;
}

vector<tuple<short, vector<double>>> BubbleParser::readBubbles() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    vector<tuple<short, vector<double>>> returnVector;
    string line;
    
    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        vector<string> lineVector = split<string>(line);
        if(lineVector.size() == 5) {
            vector<double> coord;
            
            short type = atoi((*(lineVector.begin() + 1)).c_str());
            
            for(auto it = lineVector.begin() + 2; it != lineVector.end(); it++) {
                coord.push_back(atof(((*it).c_str())));
            }
            returnVector.emplace_back(type, coord);
        }
    }
    return returnVector;
}

ChemistryData ChemistryParser::readChemistryInput() {
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    ///To keep track of duplicate names
    vector<string> allSpeciesNames;
    
    ChemistryData chem; string line;

    while(getline(_inputFile, line)) {
        
        if(line.find("#") != string::npos) { continue; }
        
        else if(line.find("SPECIESBULK") != string::npos) {
        
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  6) {
                cout << "Error reading a bulk species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 6) {
                
                if(lineVector[5] != "CONST" && lineVector[5] != "REG") {
                    
                    cout << "Option for bulk species not valid. Exiting." << endl;
                    cout << lineVector[5] << endl;
                    exit(EXIT_FAILURE);
                }

                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                
                chem.speciesBulk.push_back(tuple<string, int, double, double, string>
                    (lineVector[1], atoi(lineVector[2].c_str()),
                     atof(lineVector[3].c_str()), atof(lineVector[4].c_str()), lineVector[5]));
            }
            else {}
        }
        else if(line.find("SPECIESDIFFUSING") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() >  8 || lineVector.size() < 7) {
                cout << "Error reading a diffusing species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 8) {
                
                if(lineVector[6] != "AVG") {
                    
                    cout << "Too many arguments for a non AVG-qualified diffusing species. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                
                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                
                chem.speciesDiffusing.push_back(tuple<string, int, double, double, double, string, int>
                    (lineVector[1], atoi(lineVector[2].c_str()),
                     atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                     atof(lineVector[5].c_str()), lineVector[6], atoi(lineVector[7].c_str())));
            }
            else if (lineVector.size() == 7) {
                
                if(lineVector[6] != "REG") {
                    
                    cout << "Not enough arguments for a non REG-qualified diffusing species. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                
                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                
                chem.speciesDiffusing.push_back(tuple<string, int, double, double, double, string, int>
                     (lineVector[1], atoi(lineVector[2].c_str()),
                     atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                     atof(lineVector[5].c_str()), lineVector[6], 0));
            }
            else {}
        }
        
        else if(line.find("SPECIESFILAMENT") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                
                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                
                chem.speciesFilament[atoi(lineVector[2].c_str())].push_back(lineVector[1]);
                
            }
            else {}
        }
        else if(line.find("SPECIESBOUND") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament bound species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                
                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                
                chem.speciesBound[atoi(lineVector[2].c_str())].push_back(lineVector[1]);
            }
            else {}
        }
        
        else if(line.find("SPECIESLINKER") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament linker species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                
                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                
                chem.speciesLinker[atoi(lineVector[2].c_str())].push_back(lineVector[1]);
            }
            else {}
        }
        else if(line.find("SPECIESMOTOR") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament motor species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                
                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                
                chem.speciesMotor[atoi(lineVector[2].c_str())].push_back(lineVector[1]);
            }
            else {}
        }
        else if(line.find("SPECIESBRANCHER") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament brancher species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
            
                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                
                chem.speciesBrancher[atoi(lineVector[2].c_str())].push_back(lineVector[1]);
            }
            else {}
        }
        else if(line.find("SPECIESPLUSEND") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament plus end species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                
                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                
                chem.speciesPlusEnd[atoi(lineVector[2].c_str())].push_back(lineVector[1]);
            }
            else {}
        }
        else if(line.find("SPECIESMINUSEND") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament minus end species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                
                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                
                chem.speciesMinusEnd[atoi(lineVector[2].c_str())].push_back(lineVector[1]);
            }
            else {}
        }
        else if(line.find("BRANCHERBINDINGSITE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a brancher binding site. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3)
                chem.B_BINDING_INDEX[atoi(lineVector[2].c_str())] = lineVector[1];
            else {}
        }
        else if(line.find("LINKERBINDINGSITE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a linker binding site. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3)
                chem.L_BINDING_INDEX[atoi(lineVector[2].c_str())] = lineVector[1];
            else {}
        }
        
        else if(line.find("MOTORBINDINGSITE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() !=  3) {
                cout << "Error reading a motor binding site. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3)
                chem.M_BINDING_INDEX[atoi(lineVector[2].c_str())] = lineVector[1];
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
            
            int filType = atoi(lineVector[1].c_str());
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.nucleationReactions[filType].push_back(
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
            
            int filType = atoi(lineVector[1].c_str());
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.depolymerizationReactions[filType].push_back(
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
            
            int filType = atoi(lineVector[1].c_str());
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    
                    if(*it != "+") products.push_back((*it));
                }
                
                chem.polymerizationReactions[filType].push_back(
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

            int filType = atoi(lineVector[1].c_str());
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "<->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 4; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.linkerReactions[filType].push_back(
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
            
            int filType = atoi(lineVector[1].c_str());
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "<->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 4; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.motorReactions[filType].push_back(
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
            
            int filType = atoi(lineVector[1].c_str());
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.motorWalkingReactions[filType].push_back(
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
            
            int filType = atoi(lineVector[1].c_str());
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.agingReactions[filType].push_back(
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
            
            int filType = atoi(lineVector[1].c_str());
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.destructionReactions[filType].push_back(
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
            
            int filType = atoi(lineVector[1].c_str());
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "<->");
            if(arrowIt != lineVector.end()) {
                
                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                
                for(auto it = arrowIt + 1; it != lineVector.end() - 4; it++) {
                    if(*it != "+")  products.push_back((*it));
                }
                
                chem.branchingReactions[filType].push_back(
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
            
            int filType = atoi(lineVector[1].c_str());
            
            auto arrowIt = find(lineVector.begin(), lineVector.end(), "AT");
            if(arrowIt != lineVector.end()) {
                
                auto it = arrowIt + 1;
                
                chem.severingReactions[filType].push_back(tuple<string, double>
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

void PinRestartParser::resetPins() {
    
    //loop through filaments
    for(auto &f: Filament::getFilaments()) {
    
        _inputFile.clear();
        _inputFile.seekg(0);
        
        // Get minus end bead
        auto b1 = f->getMinusEndCylinder()->getFirstBead();
        auto b2 = f->getPlusEndCylinder()->getSecondBead();
        
        int filID = f->getID();
        string searchID = "FILAMENT " + std::to_string(filID) + ":";
        
        string line;
        
        while(getline(_inputFile, line)) {
            
            if(line.find("#") != string::npos) { continue; }
            
            else if(line.find(searchID) != string::npos) {
                
                vector<string> lineVector = split<string>(line);
                if(lineVector.size() !=  8) {
                    cout << "Error reading a restart pin position. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 8) {
                    
                    b1->pinnedPosition = vector<double>{atof(lineVector[2].c_str()), atof(lineVector[3].c_str()), atof(lineVector[4].c_str())};
                    b2->pinnedPosition = vector<double>{atof(lineVector[5].c_str()), atof(lineVector[6].c_str()), atof(lineVector[7].c_str())};
                    
                    if(!areEqual(b1->pinnedPosition[0],0.0) && !areEqual(b1->pinnedPosition[1],0.0) && !areEqual(b1->pinnedPosition[2],0.0)) {
                        b1->addAsPinned();
                    
//                        cout << "Pinned filament! coordinates = " << b1->coordinate[0] << " " << b1->coordinate[1] << " " << b1->coordinate[2] << endl;
//                        cout << "Pin position = " << b1->pinnedPosition[0] << " " << b1->pinnedPosition[1] << " " << b1->pinnedPosition[2] << endl;                        
                    }
                    
                    if(!areEqual(b2->pinnedPosition[0],0.0) && !areEqual(b2->pinnedPosition[1],0.0) && !areEqual(b2->pinnedPosition[2],0.0)) {
                        b2->addAsPinned();
                       
//                        cout << "Pinned filament! coordinates = " << b2->coordinate[0] << " " << b2->coordinate[1] << " " << b2->coordinate[2] << endl;
//                        cout << "Pin position = " << b2->pinnedPosition[0] << " " << b2->pinnedPosition[1] << " " << b2->pinnedPosition[2] << endl;
                    }
                }
            }
        }
    }
}

