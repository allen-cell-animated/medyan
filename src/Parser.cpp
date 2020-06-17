
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include <stdexcept> // runtime_error
#include <utility> // move

#include "Parser.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

#include "SysParams.h"
#include "Util/Io/Log.hpp"

void SystemParser::readChemParams() {

    ChemParams CParams;

    _inputFile.clear();
    _inputFile.seekg(0);

    string line;
    while(getline(_inputFile, line)) {

        if(line.find("#") != string::npos) { continue; }

    //        if (line.find("NUMBULKSPECIES") != string::npos) {
    //
    //            vector<string> lineVector = split<string>(line);
    //            if(lineVector.size() > 2) {
    //                cout <<
    //                     "There was an error parsing input file at Chemistry parameters. Exiting."
    //                     << endl;
    //                exit(EXIT_FAILURE);
    //            }
    //            else if (lineVector.size() == 2) {
    //                CParams.numBulkSpecies = atof(lineVector[1].c_str());
    //            }
    //        }
    //
    //        if (line.find("NUMDIFFUSINGSPECIES") != string::npos) {
    //
    //            vector<string> lineVector = split<string>(line);
    //            if(lineVector.size() > 2) {
    //                cout <<
    //                     "There was an error parsing input file at Chemistry parameters. Exiting."
    //                     << endl;
    //                exit(EXIT_FAILURE);
    //            }
    //            else if (lineVector.size() == 2) {
    //                CParams.numDiffusingSpecies = atof(lineVector[1].c_str());
    //            }
    //        }

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
            //the vector size can be 5 for PINLOWERBOUNDARYFILAMENTS
            if(lineVector.size() > 7) {
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
            else if (lineVector.size() == 4) {
                if(lineVector[1]  == "RATEDEPEND") {
                    CParams.makeRateDepend = true;
                    CParams.makeRateDependTime = atof(lineVector[2].c_str());
                    CParams.makeRateDependForce = atof(lineVector[3].c_str());
                }
            }
            else if (lineVector.size() == 7) {
                if(lineVector[1]  == "AFM") {
                    CParams.makeAFM = true;
                    //displacement of each pull
                    CParams.AFMStep1 = atof(lineVector[2].c_str());
                    CParams.AFMStep2 = atof(lineVector[3].c_str());
                    //change dispalcement from 1 to 2
                    CParams.IterChange = atof(lineVector[4].c_str());
                    //total step of each AFM pull
                    CParams.StepTotal = atof(lineVector[5].c_str());
                    //time between each pull
                    CParams.StepTime = atof(lineVector[6].c_str());
                }
            }
        }

        if (line.find("DISSIPATIONTRACKING:") != string::npos) {

            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {

                const char * testStr1 = "ON";
                const char * testStr2 = lineVector[1].c_str();
                if(strcmp(testStr1, testStr2) == 0){
                    CParams.dissTracking = true;

                }

            }
        }

        if (line.find("EVENTTRACKING:") != string::npos) {

            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {

                const char * testStr1 = "ON";
                const char * testStr2 = lineVector[1].c_str();
                if(strcmp(testStr1, testStr2) == 0){
                    CParams.eventTracking = true;

                }

            }
        }

        if (line.find("LINKERBINDINGSKIP:") != string::npos) {

            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CParams.linkerbindingskip = atoi(lineVector[1].c_str());

            }
        }

    }

    //Figure out the binding sites
    for(int i = 0; i < CParams.numBindingSites.size(); i++) {

        CParams.maxbindingsitespercylinder = max(CParams.maxbindingsitespercylinder,
                                                 CParams.numBindingSites[i]);

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
    //Find the maximum allowed Cindex and shift operator
    auto np2 = mathfunc::nextPowerOf2(uint32_t(CParams
            .maxbindingsitespercylinder));

    if(np2 == CParams.maxbindingsitespercylinder)
        np2 *= 2;

	CParams.shiftbybits = log2(np2);
    CParams.maxStableIndex = numeric_limits<uint32_t>::max()/CParams.shiftbybits -1;
//	cout<<"shiftbybits "<<CParams.shiftbybits<<" maxbindingsitespercylinder "<<CParams
//	.maxbindingsitespercylinder<<endl;
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
        if (line.find("DATADUMPTIME:") != string::npos) {

            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                CAlgorithm.datadumpTime = atof(lineVector[1].c_str());
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
        else if(line.find("MEM_TENSION_TYPE") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at membrane stretching accumulation type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2) {
                MType.memTensionFFType = lineVector[1];
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
        else if(line.find("MEM_BEAD_VOLUME_FF_TYPE") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at membrane cylinder volume FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2) {
                MType.MemBeadVolumeFFType = lineVector[1];
            }
        }
        else if(line.find("VOLUME_CONSERVATION_FF_TYPE") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at membrane cylinder volume FF type. Exiting."
                    << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2) {
                MType.volumeConservationFFType = lineVector[1];
            }
        }
        

        else if (line.find("AFMFFTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout <<
                "There was an error parsing input file at AFM FF type. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MType.AFMFFType = lineVector[1];
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
        else if (line.find("MEM_BEAD_VOLUME_K") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.memBeadVolumeK = stod(lineVector[i]);
            }
        }
        else if (line.find("MEM_BEAD_VOLUME_CUTOFF") != string::npos) {
            
            if (line.find("MEM_BEAD_VOLUME_CUTOFF_MECH") != string::npos) {
            
                vector<string> lineVector = split<string>(line);
                if(lineVector.size() != 2) {
                    LOG(ERROR) << "Error reading membrane-bead volume cutoff for mech.";
                    throw std::runtime_error("Error reading volume cutoff");
                }
                else {
                    MParams.MemBeadVolumeCutoffMech = std::stod(lineVector[1]);
                }
            }
            else {
                vector<string> lineVector = split<string>(line);
                if(lineVector.size() != 2) {
                    LOG(ERROR) << "Error reading membrane-bead volume cutoff.";
                    throw std::runtime_error("Error reading volume cutoff");
                }
                else {
                    MParams.MemBeadVolumeCutoff = std::stod(lineVector[1]);
                }
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

        else if (line.find("MTOCBENDINGK") != string::npos){
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                MParams.MTOCBendingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        else if (line.find("AFMBENDINGK") != string::npos){
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    MParams.AFMBendingK.push_back(atof((lineVector[i].c_str())));
            }
        }
        
        if (line.find("HESSIANTRACKING:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 3) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                MParams.hessTracking = true;
                //MParams.hessDelta = atof(lineVector[1].c_str());
                MParams.hessSkip = atof(lineVector[1].c_str());
                int dense = atoi(lineVector[2].c_str());
                if(dense == 0){
                    MParams.denseEstimation = true;
                }else{
                    MParams.denseEstimation = false;
                }
                
            }
        }
        
        if (line.find("SAMEFILBINDINGSKIP:") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                MParams.sameFilBindSkip = atoi(lineVector[1].c_str());
                
            }
        }
        

        // Membrane parameter
        else if (line.find("MEM_AREA_K") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() >= 2) {
                for(size_t i = 1; i < lineVector.size(); ++i)
                    MParams.memAreaK.push_back(atof(lineVector[i].c_str()));
            }
        }
        else if (line.find("MEM_EQ_AREA_FACTOR") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() >= 2) {
                for(size_t i = 1; i < lineVector.size(); ++i)
                    MParams.memEqAreaFactor.push_back(atof(lineVector[i].c_str()));
            }
        }
        else if (line.find("MEM_TENSION") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() >= 2) {
                for(size_t i = 1; i < lineVector.size(); ++i)
                    MParams.MemTension.push_back(atof(lineVector[i].c_str()));
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

        // Water compressibility
        else if(line.find("BULK_MODULUS") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() >= 2) {
                MParams.bulkModulus = atof(lineVector[1].c_str());
            }
        }
        
        if (line.find("SPECIALPROTOCOL") != string::npos) {

            vector<string> lineVector = split<string>(line);

            if(lineVector[1] == "PINBOUNDARYFILAMENTS") {
                if(lineVector.size() > 5) {
                    cout <<
                         "There was an error parsing input file at Chemistry parameters. Exiting."
                         << endl;
                    exit(EXIT_FAILURE);
                }

                else{
                    MParams.pinBoundaryFilaments = true;
                    MParams.pinK = atof(lineVector[2].c_str());
                    MParams.pinDistance = atof(lineVector[3].c_str());
                    MParams.pinTime = atof(lineVector[4].c_str());
                }

            }
            else if(lineVector[1] == "PIN_BUBBLES") {
                if(lineVector.size() != 4) {
                    LOG(ERROR) << "Error parameter specification for pinning bubbles";
                    LOG(INFO) << "PIN_BUBBLES <pin-k> <pin-time>";
                    throw std::runtime_error("Error pinning bubbles");
                }

                else {
                    MParams.pinBubbles = true;
                    MParams.pinK = std::stod(lineVector[2]);
                    MParams.pinTime = std::stod(lineVector[3]);
                }
            }
            else if(lineVector[1] == "PIN_MEMBRANE_BORDER_VERTICES") {
                if(lineVector.size() != 3) {
                    LOG(ERROR) << "Error parameter specification for pinning membrane border vertices";
                    LOG(INFO) << "PIN_MEMBRANE_BORDER_VERTICES <pin-k>";
                    throw std::runtime_error("Error pinning membrane border vertices");
                }

                else {
                    MParams.pinMembraneBorderVertices = true;
                    MParams.pinK = std::stod(lineVector[2]);
                }
            }
            else if(lineVector[1] == "PIN_INITIAL_FILAMENT_BELOW_Z") {
                if(lineVector.size() != 4) {
                    LOG(ERROR) << "Error parameter specification for pinning filament beads below z";
                    LOG(INFO) << "PIN_INITIAL_FILAMENT_BELOW_Z <pin-k> <z-max>";
                    throw std::runtime_error("Error pinning initial filament beads");
                }

                else {
                    MParams.pinInitialFilamentBelowZ      = true;
                    MParams.pinK                          = std::stod(lineVector[2]);
                    MParams.pinInitialFilamentBelowZValue = std::stod(lineVector[3]);
                }
            }
            else if (lineVector.size() == 5) {

                //Qin
                if(lineVector[1] == "PINLOWERBOUNDARYFILAMENTS") {

                    MParams.pinLowerBoundaryFilaments = true;
                    MParams.pinK = atof(lineVector[2].c_str());
                    MParams.pinTime = atof(lineVector[3].c_str());
                    MParams.pinFraction = atof(lineVector[4].c_str());
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
        else if (line.find("LAMBDARUNNINGAVERAGEPROBABILITY") != string::npos) {

            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                MAlgorithm.lambdarunningaverageprobability = atof(lineVector[1].c_str());
            }
        }
        else if (line.find("LINESEARCHALGORITHM")!= string::npos){
            vector<string> lineVector = split<string>(line);
            if (lineVector.size() == 2) {
                MAlgorithm.linesearchalgorithm = (lineVector[1].c_str());
            }
        }
    }
    return MAlgorithm;
}

void SystemParser::readBoundParams() {

    BoundParams BParams;

    _inputFile.clear();
    _inputFile.seekg(0);
    vector<int> leftfrontbottom = {0,0,0};
    vector<int> rightbacktop = {0,0,0};
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
        else if(line.find("BOUNDARYMOVE") != string::npos){
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout << "A boundary move type needs to be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                //set planes to move and transfershare axis
                if(lineVector[1] == "LEFT")
                    leftfrontbottom[0] = 1;
                else if(lineVector[1] == "BOTTOM")
                    leftfrontbottom[1] = 1;
                else if(lineVector[1] == "FRONT")
                    leftfrontbottom[2] = 1;
                else if(lineVector[1] == "RIGHT")
                    rightbacktop[0] = 1;
                else if(lineVector[1] == "TOP")
                    rightbacktop[1] = 1;
                else if(lineVector[1] == "BACK")
                    rightbacktop[2] = 1;
            }
        }
        else if(line.find("TRANSFERSHAREAXIS") != string::npos){
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 3) {
                cout <<
                     "There was an error parsing input file at Chemistry parameters. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }

            else{
                cout<<"TRANSFERSHARE AXIS "<<lineVector[2]<<endl;
                if(lineVector[2]=="X")
                    BParams.transfershareaxis=0;
                else if(lineVector[2]=="Y")
                    BParams.transfershareaxis=1;
                else if(lineVector[2]=="Z")
                    BParams.transfershareaxis=2;
                else if(lineVector[2]=="RADIAL") {
                    BParams.transfershareaxis = 3;
                    cout<<"RADIAL transfer not implemented. Change paramters. Exiting"
                            "."<<endl;
                    exit(EXIT_FAILURE);
                }
                else{
                    cout << "There was an error parsing input file at Chemistry parameters. Exiting."
                         << endl;
                    exit(EXIT_FAILURE);}
            }
        }
        else if(line.find("FILCREATIONBOUNDS") != string::npos){
            vector<string> lineVector = split<string>(line);
            if(lineVector.size()!=7) {
                cout << "FILCREATIONBOUNDS should have 6 elements. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else{
                vector<floatingpoint> tempvec;
                vector<vector<floatingpoint>> tempbounds;
                for(int i = 1;i<4;i++)
                    tempvec.push_back(atof((lineVector[i].c_str())));
                tempbounds.push_back(tempvec);
                tempvec.clear();
                for(int i = 4;i<7;i++)
                    tempvec.push_back(atof((lineVector[i].c_str())));
                tempbounds.push_back(tempvec);
                tempvec.clear();
                BParams.fraccompartmentspan = tempbounds;
            }
        }

        else {}
    }

    for(int i = 0; i < 3; i++){
        int addthemup = leftfrontbottom[i] + rightbacktop[i];
        if(addthemup > 0)
            BParams.transfershareaxis = i;
        if(addthemup == 2)
            BParams.planestomove = 2;
        else if(leftfrontbottom[i] == 1)
            BParams.planestomove = 1;
        else if(rightbacktop[i] == 1)
            BParams.planestomove = 0;
    }
//    std::cout<<BParams.transfershareaxis<<" "<<BParams.planestomove<<endl;
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
        else if (line.find("DBUNBINDINGLEN") != string::npos) {
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRParams.dBranchUnbindingCharLength.push_back(
                            atof((lineVector[i].c_str())));
            }
            else {}
        }
        else if (line.find("DBUNBINDINGF") != string::npos) {
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRParams.dBranchUnbindingCharForce.push_back(
                            atof((lineVector[i].c_str())));
            }
            else {}
        }
        /// Manual Rate Changer
        // It takes 5 inputs as start_time, plusend_poly, plusend_depoly, minusend_poly, minusend_depoly
        // Currently it applies type 0 to all filament types
        else if (line.find("MANUALSTARTTIME") != string::npos) {
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                DRParams.manualCharStartTime = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("MANUALPLUSPOLYRATIO") != string::npos) {
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                DRParams.manualPlusPolyRate = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("MANUALPLUSDEPOLYRATIO") != string::npos) {
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                DRParams.manualPlusDepolyRate = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("MANUALMINUSPOLYRATIO") != string::npos) {
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                DRParams.manualMinusPolyRate = atof((lineVector[1].c_str()));
            }
            else {}
        }
        else if (line.find("MANUALMINUSDEPOLYRATIO") != string::npos) {
            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                DRParams.manualMinusDepolyRate = atof((lineVector[1].c_str()));
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

            //adding branching dy type
        else if (line.find("DBUNBINDINGTYPE") != string::npos) {

            vector<string> lineVector = split<string>(line);

            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    DRType.dBUnbindingType.push_back(lineVector[i]);
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
                BType.boundaryMove.push_back(lineVector[1]);
            }
        }
        //add Compartment Scaling
        //else if (line.find("DIFFUSIONSCALE") != string::npos) {

        //  vector<string> lineVector = split<string>(line);
        //  if(lineVector.size() != 2) {
        //      cout << "Diffusion scaling needs to be specified. Exiting." << endl;
        //      exit(EXIT_FAILURE);
        //  }
        //  else if (lineVector.size() == 2) {
        //      BType.scaleDiffusion = lineVector[1];
        //  }
        //}
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
            else if (lineVector[1] == "AFM") SType.afm = true;
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
        
        else if (line.find("AFMFILAMENTTYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "A filament type to connect to an AFM must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                SType.afmFilamentType = atoi(lineVector[1].c_str());
            }
        }
        else if (line.find("AFMNUMFILAMENTS") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "A number of filaments to connect to an AFM must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                SType.afmNumFilaments = atoi(lineVector[1].c_str());
            }
        }
        else if (line.find("AFMFILAMENTLENGTH") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "A filament length for filaments connected to an AFM must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                SType.afmFilamentLength = atoi(lineVector[1].c_str());
            }
        }
        
        else if (line.find("AFMBUBBLETYPE") != string::npos) {
            
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                cout <<
                "A bubble type to connect to an MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                SType.afmBubbleType = atoi(lineVector[1].c_str());
            }
        }
        else if (line.find("MTOCXYZCOORD") != string::npos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 4) {
                cout <<
                "Coordinates of MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 4) {
                SType.mtocInputCoordXYZ.push_back(stof(lineVector[1].c_str()));
                SType.mtocInputCoordXYZ.push_back(stof(lineVector[2].c_str()));
                SType.mtocInputCoordXYZ.push_back(stof(lineVector[3].c_str()));
            }
        }
    }
    return SType;
}

void SystemParser::readSpecialParams() {
    
    SpecialParams SParams;
    
    _inputFile.clear();
    _inputFile.seekg(0);
    
    string line;
    
    if (line.find("MTOCFILAMENTCOORD") != string::npos) {
        vector<string> lineVector = split<string>(line);
        if(lineVector.size() != 5) {
            cout << "4 coordinates of MTOC filaments must be specified. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        else if (lineVector.size() == 5) {
            SParams.mtocTheta1 = stof(lineVector[1].c_str());
            SParams.mtocTheta2 = stof(lineVector[2].c_str());
            SParams.mtocPhi1 = stof(lineVector[3].c_str());
            SParams.mtocPhi2 = stof(lineVector[4].c_str());
        }
    }

    SysParams::SParams = SParams;
    
}

void SystemParser::readSimulParams() {
    auto& sp = SysParams::simulParamsMut();

    _inputFile.clear();
    _inputFile.seekg(0);

    std::string line;
    while(getline(_inputFile, line)) {
        const auto hashPos = line.find('#');

        if(line.find("TRACK_FORCES") < hashPos) {
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() != 2) {
                LOG(ERROR) << "1 value (on/off) is required as the argument for this option";
                throw std::runtime_error("Invalid number of arguments for TRACK_FORCES");
            }
            else {
                if(lineVector[1] == "on") sp.trackForces = true;
                else if(lineVector[1] == "off") sp.trackForces = false;
                else {
                    LOG(ERROR) << "Unrecognized argument " << lineVector[1] << " for TRACK_FORCES";
                    throw std::runtime_error("Invalid argument for TRACK_FORCES");
                }
            }
        }
    }
}

void SystemParser::readGeoParams() {

    _inputFile.clear();
    _inputFile.seekg(0);

    GeoParams GParams;

    vector<floatingpoint> gridTemp;
    vector<floatingpoint> compartmentTemp;
    vector<floatingpoint> monomerSize = {};
    vector<floatingpoint> cylinderSize = {};
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
    //find max Cylinder size
    GParams.largestCylinderSize = 0;
    for(auto l:GParams.cylinderSize)
        GParams.largestCylinderSize = max(GParams.largestCylinderSize, l);
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

                line.find("MTOCNUMFILAMENTS") == string::npos &&
                line.find("AFMNUMFILAMENTS") == string::npos) {

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

                line.find("MTOCFILAMENTLENGTH") == string::npos &&
                line.find("AFMFILAMENTLENGTH") == string::npos) {
            

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
                line.find("MTOCFILAMENTTYPE") == string::npos &&
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
        else if (line.find("RESTARTPHASE") != string::npos){SysParams::RUNSTATE=false;
            vector<string> lineVector = split<string>(line);
            if(lineVector.size() > 2) {
                cout << "Error reading restart params. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2){
                if(lineVector[1].find("USECHEMCOPYNUM"))
                SysParams::USECHEMCOPYNUM = true;
            }
        }
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
    
    MembraneSetup memSetup;
    
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
                memSetup.inputFile = lineVector[1];
        }
        else if(line.find("MEMBRANE_INIT") != string::npos) {
            std::vector< std::string > lineVector = split< std::string >(line);
            if(lineVector.size() < 2) {
                LOG(ERROR) << "Error reading membrane initialization parameters.";
                throw std::runtime_error("Error reading membrane initialization");
            } else {
                lineVector.erase(lineVector.begin()); // remove first element
                memSetup.meshParam.push_back(std::move(lineVector));
            }
        }
    }
    return memSetup;
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
    tuple< vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>> , vector<tuple<string, short, vector<vector<floatingpoint>>>> , vector<tuple<string, short, vector<floatingpoint>>> , vector<vector<floatingpoint>> > FilamentParser::readFilaments() {
    _inputFile.clear();
    _inputFile.seekg(0);
     vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>> filamentVector;
     vector<vector<vector<floatingpoint>>> linkerVector;
     vector<vector<vector<floatingpoint>>> motorVector;
     vector<vector<floatingpoint>> staticVector;
     vector<tuple<string, short, vector<vector<floatingpoint>>>> boundVector;
     vector<tuple<string, short, vector<floatingpoint>>> branchVector;
     string line;
    
    while(getline(_inputFile, line)) {

        if(line.find("#") != string::npos) { continue; }

        vector<string> lineVector = split<string>(line);
        if(lineVector.size() >= 8) {
            vector<floatingpoint> coord1;
            vector<floatingpoint> coord2;
            vector<vector<floatingpoint>> coord3;
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
                /*Linker Motor*/
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
            vector<floatingpoint> coord1;
            vector<vector<floatingpoint>> coord3;
            //USED ONLY TO RESTART PINNED TRAJECTORIES.
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
      tuple< vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>> , vector<tuple<string, short, vector<vector<floatingpoint>>>> , vector<tuple<string, short, vector<floatingpoint>>> , vector<vector<floatingpoint>> > returnVector=make_tuple(filamentVector,boundVector,branchVector, staticVector);
    return returnVector;
}

vector<MembraneParser::MembraneInfo> MembraneParser::readMembranes() {

    vector<MembraneInfo> res;
    
    bool wasEmpty = true;
    int stage; // 0: init, number of vertices; 1: vertex coordinate; 2: triangle vertex index
    size_t numVertices;
    
    string line;
    while(getline(_inputFile, line)) {
        
        bool isEmpty = (line.empty() || line.find("#") != string::npos); // empty line or line with '#'

        if(wasEmpty && !isEmpty) { // New membrane
            res.emplace_back();
            stage = 0;
        }

        wasEmpty = isEmpty;
        if(isEmpty) continue;

        auto& activeMem = res.back();

        vector<string> lineVector = split<string>(line);
        size_t lineVectorSize = lineVector.size();

        switch(stage) {
        case 0: // init, number of vertices
            if(lineVectorSize != 1) cout << "First line of membrane should be number of vertices." << endl;
            numVertices = atoi(lineVector[0].c_str());
            activeMem.vertexCoordinateList.reserve(numVertices);
            stage = 1;
            break;
        case 1: // vertex coordinate
            {
                if(lineVectorSize != 3) {
                    cout << "Each vertex should have 3 coordinates" << endl;
                }

                activeMem.vertexCoordinateList.emplace_back();
                auto& activeCoord = activeMem.vertexCoordinateList.back();

                // Parse coordinate information
                for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                    activeCoord[coordIdx] = atof(lineVector[coordIdx].c_str());
                }
            }

            if(activeMem.vertexCoordinateList.size() == numVertices) stage = 2;
            break;
        case 2: // triangle vertex indices
            {
                if(lineVectorSize != 3) {
                    cout << "Each triangle should have 3 indices" << endl;
                }

                activeMem.triangleVertexIndexList.emplace_back();
                auto& activeTriangle = activeMem.triangleVertexIndexList.back();

                for(size_t i = 0; i < 3; ++i)
                    activeTriangle[i] = atoi(lineVector[i].c_str());
            }
            break;
        }
    }

    return res;
}

vector<tuple<short, vector<floatingpoint>>> BubbleParser::readBubbles() {

    _inputFile.clear();
    _inputFile.seekg(0);
    
    vector<tuple<short, vector<floatingpoint>>> returnVector;
    string line;

    while(getline(_inputFile, line)) {

        if(line.find("#") != string::npos) { continue; }

        vector<string> lineVector = split<string>(line);
        if(lineVector.size() == 5) {
            vector<floatingpoint> coord;
            
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
            if(lineVector.size() !=  6 && lineVector.size() !=  8) {
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

                chem.speciesBulk.push_back(tuple<string, int, floatingpoint, floatingpoint,
                        string, string, floatingpoint>(lineVector[1], atoi(lineVector[2].c_str()),
                                atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                                lineVector[5], "NONE", 0.0));
            }
            else if (lineVector.size() == 8) {

                if(lineVector[5] != "CONST" && lineVector[5] != "REG") {

                    cout << "Option for bulk species not valid. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }

                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }

                chem.speciesBulk.push_back(tuple<string, int, floatingpoint, floatingpoint,
                        string, string, floatingpoint>(lineVector[1], atoi(lineVector[2].c_str()),
                                atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                                lineVector[5],lineVector[6], atof(lineVector[7].c_str())));
            }
            else {}
        }
        else if(line.find("SPECIESDIFFUSING") != string::npos) {

            vector<string> lineVector = split<string>(line);


            if(lineVector.size() >  9 || lineVector.size() < 7) {
                cout << "Error reading a diffusing species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 8) {

                if(lineVector[6] != "AVG") {

                    cout << "Too many arguments for a non AVG-qualified diffusing "
                            "species. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }

                if(find(allSpeciesNames.begin(), allSpeciesNames.end(), lineVector[1]) != allSpeciesNames.end()) {
                    cout << "Duplicate species names are not allowed. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    allSpeciesNames.push_back(lineVector[1]);
                }
                if(lineVector[6] == "AVG")
                    chem.speciesDiffusing.push_back(tuple<string, int, floatingpoint, floatingpoint,
                            floatingpoint, string, int, string, floatingpoint>
                    (lineVector[1], atoi(lineVector[2].c_str()),
                     atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                     atof(lineVector[5].c_str()), lineVector[6], atoi(lineVector[7].c_str
                             ()),"NONE", 0.0));
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
                
                chem.speciesDiffusing.push_back(tuple<string, int, floatingpoint, floatingpoint,
                        floatingpoint, string, int, string, floatingpoint>
                     (lineVector[1], atoi(lineVector[2].c_str()),
                     atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                     atof(lineVector[5].c_str()), lineVector[6], 0, "NONE", 0.0));
            }
            else if (lineVector.size() == 9) {

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

                chem.speciesDiffusing.push_back(tuple<string, int, floatingpoint, floatingpoint,
                        floatingpoint, string, int, string, floatingpoint>
                                                        (lineVector[1], atoi(lineVector[2].c_str()),
                                                         atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                                                         atof(lineVector[5].c_str()),
                                                         lineVector[6], 0, lineVector[7],
                                                         atof(lineVector[8].c_str())));
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

            // check parameters related to dissipation tracking if it is enabled
            float gnum = 0.0;
            int dissOffSet = 0;
            string HRCDID = "NA";
            if(SysParams::Chemistry().dissTracking){
            string dissString = lineVector[1].c_str();
            istringstream iss(dissString);
            string token;
            vector<string> dissTokens;
            dissOffSet = 1;

            while (std::getline(iss, token, ':')) {
                if (!token.empty())
                    dissTokens.push_back(token);
            }

            gnum = atof(dissTokens[0].c_str());

            if(dissTokens.size()!=1){
                HRCDID = dissTokens[1];
            } else {
                HRCDID = "NA";
            }
            }


            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                for(auto it  = lineVector.begin() + 1 + dissOffSet; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }

                chem.genReactions.push_back(
                        tuple<vector<string>, vector<string>, floatingpoint, floatingpoint,
                        string>(reactants, products, atof(lineVector[lineVector.size() - 1].c_str()),gnum,HRCDID));
                
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
                tuple<vector<string>, vector<string>, floatingpoint>
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
                tuple<vector<string>, vector<string>, floatingpoint>
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

            int filType = atoi(lineVector[2].c_str());
            // check parameters related to dissipation tracking if it is enabled
            float gnum = 0.0;
            int dissOffSet = 0;
            string HRCDID = "NA";
            if(SysParams::Chemistry().dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffSet = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {

                for(auto it  = lineVector.begin() + 2 + dissOffSet; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }

                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }

                chem.depolymerizationReactions[filType].push_back(
                        tuple<vector<string>, vector<string>, floatingpoint,floatingpoint, string>
                                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str()),gnum,HRCDID));
                
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

            int filType = atoi(lineVector[2].c_str());
            // check parameters related to dissipation tracking if it is enabled
            floatingpoint gnum = 0.0;
            int dissOffSet = 0;
            string HRCDID = "NA";
            if(SysParams::Chemistry().dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffSet = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {

                for(auto it  = lineVector.begin() + 2 + dissOffSet; it != arrowIt; it++) {

                    if(*it != "+") reactants.push_back((*it));
                }

                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {

                    if(*it != "+") products.push_back((*it));
                }


                chem.polymerizationReactions[filType].push_back(
                        tuple<vector<string>, vector<string>, floatingpoint,floatingpoint,string>
                                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str()),gnum,HRCDID));
                
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

            int filType = atoi(lineVector[2].c_str());
            // check parameters related to dissipation tracking if it is enabled
            floatingpoint gnum = 0.0;
            int dissOffSet = 0;
            string HRCDID = "NA";
            if(SysParams::Chemistry().dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffSet = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            auto arrowIt = find(lineVector.begin(), lineVector.end(), "<->");
            if(arrowIt != lineVector.end()) {

                for(auto it  = lineVector.begin() + 2 + dissOffSet; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }

                for(auto it = arrowIt + 1; it != lineVector.end() - 4; it++) {
                    if(*it != "+")  products.push_back((*it));
                }

                chem.linkerReactions[filType].push_back(
                        tuple<vector<string>, vector<string>, floatingpoint, floatingpoint, floatingpoint, floatingpoint, floatingpoint,string>
                                (reactants, products, atof(lineVector[lineVector.size() - 4].c_str()),
                                 atof(lineVector[lineVector.size() - 3].c_str()),
                                 atof(lineVector[lineVector.size() - 2].c_str()),
                                 atof(lineVector[lineVector.size() - 1].c_str()),gnum, HRCDID));
                
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

            int filType = atoi(lineVector[2].c_str());
            // check parameters related to dissipation tracking if it is enabled
            floatingpoint gnum = 0.0;
            int dissOffSet = 0;
            string HRCDID = "NA";
            if(SysParams::Chemistry().dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffSet = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            auto arrowIt = find(lineVector.begin(), lineVector.end(), "<->");
            if(arrowIt != lineVector.end()) {

                for(auto it  = lineVector.begin() + 2 + dissOffSet; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }

                for(auto it = arrowIt + 1; it != lineVector.end() - 4; it++) {
                    if(*it != "+")  products.push_back((*it));
                }

                chem.motorReactions[filType].push_back(
                        tuple<vector<string>, vector<string>, floatingpoint, floatingpoint, floatingpoint, floatingpoint, floatingpoint,string>
                                (reactants, products, atof(lineVector[lineVector.size() - 4].c_str()),
                                 atof(lineVector[lineVector.size() - 3].c_str()),
                                 atof(lineVector[lineVector.size() - 2].c_str()),
                                 atof(lineVector[lineVector.size() - 1].c_str()),gnum,HRCDID));
                
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

            int filType = atoi(lineVector[2].c_str());
            // check parameters related to dissipation tracking if it is enabled
            floatingpoint gnum = 0.0;
            int dissOffSet = 0;
            string HRCDID = "NA";
            if(SysParams::Chemistry().dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffSet = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {

                for(auto it  = lineVector.begin() + 2 + dissOffSet; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }

                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }

                chem.motorWalkingReactions[filType].push_back(
                        tuple<vector<string>, vector<string>, floatingpoint, floatingpoint,string>
                                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str()),gnum,HRCDID));
                
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

            int filType = atoi(lineVector[2].c_str());
            // check parameters related to dissipation tracking if it is enabled 
            floatingpoint gnum = 0.0;
            int dissOffSet = 0;
            string HRCDID = "NA";
            if(SysParams::Chemistry().dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffSet = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {

                for(auto it  = lineVector.begin() + 2 + dissOffSet; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }

                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }

                chem.agingReactions[filType].push_back(
                        tuple<vector<string>, vector<string>, floatingpoint, floatingpoint,string>
                                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str()),gnum,HRCDID));
                
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
                tuple<vector<string>, vector<string>, floatingpoint>
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
                tuple<vector<string>, vector<string>, floatingpoint, floatingpoint, string, floatingpoint>
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
                
                chem.severingReactions[filType].push_back(tuple<string, floatingpoint>
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

        int filID = f->getId();
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
                    
                    b1->pinnedPosition = vector<floatingpoint>{stof(lineVector[2].c_str()
                                                               ), stof(lineVector[3].c_str()), stof(lineVector[4].c_str())};
                    b2->pinnedPosition = vector<floatingpoint>{stof(lineVector[5].c_str()
                                                               ), stof(lineVector[6].c_str()), stof(lineVector[7].c_str())};
                    
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
